"""新写验证工具：多工况空间／相位离散矩阵，不改变计算内核或验收阈值。

分别固定另一离散轴，保存实际输入、原始数组与平衡残差。所谓通过仅指
所列各相邻级的指标变化满足 abs + rel*|细级值|，不是连续解误差上界、
标准认证或实物验证。quick 只检查管道，正式收敛状态始终为 null。
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p / "pyproject.toml").is_file() and (p / "gear_contact").is_dir())
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from gear_contact import Gear, GearPair, SurfaceModification, analyze_config, export_result
from gear_contact.io import write_json

_SPEC = importlib.util.spec_from_file_location("convergence_equilibrium_checks", PROJECT / "examples/validate_cases.py")
_CHECKS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CHECKS)
check_result = _CHECKS.check_result

CASE_LABELS = {
    "combined_spur": "外直齿组合修形，80 N·m",
    "internal_helical": "内斜齿、鼓形与3 μm失准，80 N·m",
    "low_load": "外直齿组合修形低载，0.8 N·m",
    "edge_misalignment": "外直齿鼓形与100 μm强失准边缘承载，80 N·m",
    "iso_love_18_77": "18/77 ISO/Love延伸接触及中心／端点压力",
}

# 各指标绝对门槛采用其SI单位，近零时不使用任意极小分母放大相对误差。
DEFAULT_ABSOLUTE_TOLERANCES = {
    "te_peak_to_peak_m": 1e-10,
    "mean_mesh_stiffness_n_per_m": 100.,
    "max_center_pressure_pa": 1e4,
    "max_endpoint_pressure_pa": 1e4,
    "max_line_load_n_per_m": 1.,
    "mean_active_contact_length_m": 1e-8,
}


def case_config(case_id):
    """返回独立可复算配置；18/77沿用原名义扭矩与LTCA乘子，不重复放大。"""
    if case_id not in CASE_LABELS:
        raise ValueError("未知矩阵案例：" + str(case_id))
    if case_id == "iso_love_18_77":
        case = json.loads((PROJECT / "examples/helical_18_77_case.json").read_text(encoding="utf-8"))
        config = deepcopy(case["analysis"])
        nominal = case["loading"]["nominal_pinion_torque_n_m"]
        config["torque_n_m"] = nominal * case["loading"]["ltca_torque_multiplier"]
        config["iso_settings"]["nominal_pinion_torque_n_m"] = nominal
        return config
    internal = case_id == "internal_helical"
    beta = .18 if internal else 0.
    pair = GearPair(Gear(28, .0025, .025, helix_angle_rad=beta),
                    Gear(84 if internal else 56, .0025, .025,
                         helix_angle_rad=beta if internal else -beta, internal=internal))
    surface = SurfaceModification(lead_crown_m=8e-6)
    if case_id in ("combined_spur", "low_load"):
        surface = SurfaceModification(lead_crown_m=8e-6, tip_relief_m=6e-6,
                                      tip_start=.72, root_relief_m=3e-6, root_end=.18)
    return {"pair": asdict(pair), "torque_n_m": .8 if case_id == "low_load" else 80.,
            "modification_pinion": asdict(surface), "elastic_mode": "coupled", "local_contact": True,
            "misalignment_m": 3e-6 if internal else 100e-6 if case_id == "edge_misalignment" else 0.}


def _count(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 2:
        raise ValueError(name + "须为不小于2的整数")
    return value


def _levels(values, name):
    if not isinstance(values, (tuple, list)) or len(values) < 3:
        raise ValueError(name + "须至少三级")
    levels = [_count(v, name) for v in values]
    if any(b <= a for a, b in zip(levels, levels[1:])):
        raise ValueError(name + "须严格递增，不自动排序或去重")
    return levels


def _tolerances(relative, absolute):
    if isinstance(relative, bool) or not isinstance(relative, (int, float)) or not math.isfinite(relative) or not 0 < relative < 1:
        raise ValueError("相对门槛须为(0,1)内的有限数")
    values = dict(DEFAULT_ABSOLUTE_TOLERANCES)
    if absolute is not None:
        if not isinstance(absolute, dict) or set(absolute) - set(values):
            raise ValueError("绝对门槛须按已定义指标名提供")
        values.update(absolute)
    for name, value in values.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(name + "的绝对门槛须为有限正数")
    return float(relative), values


def metrics_from_arrays(arrays, *, require_endpoints=False):
    """直接从无pickle原数组复算指标；缺少端点表示未知，不补为零。

    active_mask沿cell_offsets切分，每相位受载长度为活跃条带宽度之和，
    最后作相位均值；它不是三维接触面积。
    """
    required = ("phase", "normal_approach_m", "mesh_stiffness_n_per_m", "cell_offsets",
                "pressure_pa", "line_load_n_per_m", "widths_m", "active_mask")
    if any(name not in arrays for name in required):
        raise ValueError("原数组缺少指标所需字段")
    phase = np.asarray(arrays["phase"])
    offsets = np.asarray(arrays["cell_offsets"])
    if (phase.ndim != 1 or not phase.size or offsets.dtype.kind not in "iu" or
            offsets.shape != (len(phase) + 1,) or offsets[0] != 0 or np.any(np.diff(offsets) <= 0)):
        raise ValueError("相位或cell_offsets布局无效")
    count = int(offsets[-1])
    for name in required:
        value = np.asarray(arrays[name])
        expected = len(phase) + 1 if name == "cell_offsets" else len(phase) if name in (
            "phase", "normal_approach_m", "mesh_stiffness_n_per_m") else count
        if value.shape != (expected,) or not np.all(np.isfinite(value)):
            raise ValueError(name + "布局不符或含非有限值")
    if np.asarray(arrays["active_mask"]).dtype.kind != "b":
        raise ValueError("active_mask必须保留布尔数组")
    for name in ("mesh_stiffness_n_per_m", "pressure_pa", "line_load_n_per_m", "widths_m"):
        if np.any(np.asarray(arrays[name]) < 0):
            raise ValueError(name + "不能为负")
    active = np.asarray(arrays["active_mask"])
    widths = np.asarray(arrays["widths_m"])
    lengths = [float(np.sum(widths[a:b][active[a:b]])) for a, b in zip(offsets, offsets[1:])]
    endpoint = None
    if "endpoint_pressure_pa" in arrays:
        pressure = np.asarray(arrays["endpoint_pressure_pa"])
        ends = np.asarray(arrays.get("rotation_endpoint_offsets", []))
        if (pressure.ndim != 1 or not pressure.size or not np.all(np.isfinite(pressure)) or np.any(pressure < 0) or
                ends.dtype.kind not in "iu" or ends.shape != offsets.shape or ends[0] != 0 or
                ends[-1] != len(pressure) or np.any(np.diff(ends) <= 0)):
            raise ValueError("端点压力或逐相位端点布局无效")
        endpoint = float(np.max(pressure))
    elif require_endpoints or "rotation_endpoint_offsets" in arrays:
        raise ValueError("该工况必须保留完整端点压力，缺失值不能补零")
    return {"te_peak_to_peak_m": float(np.ptp(arrays["normal_approach_m"])),
            "mean_mesh_stiffness_n_per_m": float(np.mean(arrays["mesh_stiffness_n_per_m"])),
            "max_center_pressure_pa": float(np.max(arrays["pressure_pa"])),
            "max_endpoint_pressure_pa": endpoint,
            "max_line_load_n_per_m": float(np.max(arrays["line_load_n_per_m"])),
            "mean_active_contact_length_m": float(np.mean(lengths))}


def compare_metrics(coarse, fine, *, relative_tolerance=.02, absolute_tolerances=None):
    """逐指标采用绝对加相对门槛；两级均无定义的端点量保持未评估。"""
    relative, absolute = _tolerances(relative_tolerance, absolute_tolerances)
    rows = {}
    for name, atol in absolute.items():
        a, b = coarse[name], fine[name]
        if a is None or b is None:
            rows[name] = {"coarse": a, "fine": b, "evaluated": False,
                          "passed": None if a is None and b is None else False,
                          "reason": "两级均未定义端点量" if a is None and b is None else "两级指标定义不一致"}
            continue
        if not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in (a, b)):
            raise ValueError("比较指标须为有限数或明确未定义的端点量")
        difference = abs(b-a)
        limit = atol + relative * abs(b)
        rows[name] = {"coarse": a, "fine": b, "evaluated": True, "absolute_difference": difference,
                      "relative_difference": None if b == 0 else difference / abs(b),
                      "absolute_tolerance": atol, "relative_tolerance": relative,
                      "allowed_difference": limit, "passed": difference <= limit}
    return {"metrics": rows, "passed": all(row["passed"] is not False for row in rows.values())}


def run(output, *, cases=None, space_levels=(11, 21, 41), phase_levels=(16, 32, 64),
        fixed_phases=32, fixed_strips=21, relative_tolerance=.02, absolute_tolerances=None, quick=False):
    """执行矩阵并持续写报告；单级求解／校核失败保留输入并继续其他级及案例。

    quick固定使用3/5/7条带和2/4/8相位，以4相位／5条带固定另一轴；
    请求的正式设置仍记录。正式通过要求全部级原数组校核通过且所有相邻级
    的已定义指标变化达标。物理正确性和全参数域保证始终不在本验收范围内。
    """
    space_levels = _levels(space_levels, "空间级数")
    phase_levels = _levels(phase_levels, "相位级数")
    fixed_phases, fixed_strips = _count(fixed_phases, "固定相位数"), _count(fixed_strips, "固定条带数")
    relative_tolerance, absolute_tolerances = _tolerances(relative_tolerance, absolute_tolerances)
    if not isinstance(quick, bool):
        raise ValueError("quick必须为布尔值")
    selected = list(CASE_LABELS) if cases is None else list(cases)
    if not selected or len(set(selected)) != len(selected) or any(c not in CASE_LABELS for c in selected):
        raise ValueError("案例选择须为非空、不重复的已定义案例列表")
    requested = {"space_levels": space_levels, "phase_levels": phase_levels,
                 "fixed_phases": fixed_phases, "fixed_strips": fixed_strips}
    if quick:
        space_levels, phase_levels, fixed_phases, fixed_strips = [3, 5, 7], [2, 4, 8], 4, 5
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("输出目录必须为空，避免旧成功结果混入本次验收")
    started = time.perf_counter()
    report = {"schema": "cylindrical_convergence_matrix_v1", "implementation": "新写验证工具",
              "generated_utc": datetime.now(timezone.utc).isoformat(), "quick_mode": quick,
              "requested_resolution": requested,
              "actual_resolution": {"space_levels": space_levels, "phase_levels": phase_levels,
                                    "fixed_phases": fixed_phases, "fixed_strips": fixed_strips},
              "relative_tolerance": relative_tolerance, "absolute_tolerances": absolute_tolerances,
              "acceptance_rule": "所有级原数组校核通过；全部相邻级已定义指标满足abs+rel*|细级值|",
              "interpretation": "所列合成工况的离散敏感性验收，不是连续解误差上界或全参数域保证",
              "physical_or_experimental_validation": False, "cases": [], "completed": False,
              "formal_convergence_passed": None, "numerical_checks_passed": None,
              "python": sys.version.split()[0], "numpy": np.__version__,
              "kernel_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted((PROJECT / "gear_contact").glob("*.py"))},
              "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

    def save():
        report["elapsed_seconds"] = time.perf_counter() - started
        write_json(output / "convergence_matrix.json", report)

    save()
    for case_id in selected:
        base = case_config(case_id)
        cache = {}
        case = {"case_id": case_id, "label": CASE_LABELS[case_id], "base_input": base,
                "endpoint_pressure_required": "extended_contact_settings" in base, "axes": {}}
        report["cases"].append(case)
        for axis, levels in (("space", space_levels), ("phase", phase_levels)):
            rows = []
            for count in levels:
                strips, phases = (count, fixed_phases) if axis == "space" else (fixed_strips, count)
                key = (strips, phases)
                if key not in cache:
                    folder = output / case_id / f"strips_{strips}_phases_{phases}"
                    folder.mkdir(parents=True)
                    config = dict(deepcopy(base), strips=strips, phases=phases)
                    # 在求解之前保存完整请求，求解器抛错也能复算并定位失败工况。
                    write_json(folder / "requested_input.json", config)
                    row = {"strips": strips, "phases": phases, "path": folder.relative_to(output).as_posix(),
                           "numerical_checks_passed": False, "metrics": None, "status": "未完成"}
                    try:
                        result = analyze_config(config)
                        export_result(result, folder)
                        validation = check_result(result)
                        write_json(folder / "validation.json", validation)
                        with np.load(folder / "arrays.npz", allow_pickle=False) as arrays:
                            metrics = metrics_from_arrays(arrays, require_endpoints=case["endpoint_pressure_required"])
                        row.update(metrics=metrics, validation=validation,
                                   numerical_checks_passed=bool(validation["passed"]),
                                   status="已求解并通过原数组校核" if validation["passed"] else "原数组校核失败")
                        if case_id == "edge_misalignment":
                            edge_phases = 0
                            for rotation in result.rotations:
                                loaded_edge = False
                                for tooth in np.unique(rotation.grid.tooth_pair_index):
                                    ids = np.flatnonzero(rotation.grid.tooth_pair_index == tooth)
                                    z = rotation.grid.face_coordinate_m[ids]
                                    ends = ids[[int(np.argmin(z)), int(np.argmax(z))]]
                                    loaded_edge |= bool(np.any(rotation.solution.forces_n[ends] > 0))
                                edge_phases += int(loaded_edge)
                            row["loaded_edge_phase_count"] = edge_phases
                            row["edge_loading_observed"] = edge_phases > 0
                            if not row["edge_loading_observed"]:
                                row.update(numerical_checks_passed=False, status="未观察到该案例要求的边缘承载")
                    except Exception as error:
                        row.update(status="求解或原数组归档校核失败",
                                   failure={"exception_type": type(error).__name__, "message": str(error),
                                            "traceback": traceback.format_exc()})
                        write_json(folder / "failure.json", row["failure"])
                    cache[key] = row
                    write_json(folder / "level_result.json", row)
                    print(f"{case_id}：条带{strips}／相位{phases}，{row['status']}", flush=True)
                rows.append(dict(cache[key], level=count))
            comparisons = []
            for coarse, fine in zip(rows, rows[1:]):
                entry = {"coarse_level": coarse["level"], "fine_level": fine["level"], "passed": False}
                if coarse["numerical_checks_passed"] and fine["numerical_checks_passed"]:
                    entry.update(compare_metrics(coarse["metrics"], fine["metrics"],
                        relative_tolerance=relative_tolerance, absolute_tolerances=absolute_tolerances))
                else:
                    entry["reason"] = "相邻级至少一级未通过求解或原数组校核，不计算伪变化量"
                comparisons.append(entry)
            axis_passed = all(row["numerical_checks_passed"] for row in rows) and all(c["passed"] for c in comparisons)
            case["axes"][axis] = {"levels": rows, "comparisons": comparisons,
                                  "sampled_changes_passed": axis_passed,
                                  "formal_convergence_passed": None if quick else axis_passed}
            save()
        case["numerical_checks_passed"] = all(row["numerical_checks_passed"] for row in cache.values())
        case["sampled_changes_passed"] = all(axis["sampled_changes_passed"] for axis in case["axes"].values())
        case["formal_convergence_passed"] = None if quick else case["sampled_changes_passed"]
        save()
    report["completed"] = True
    report["numerical_checks_passed"] = all(case["numerical_checks_passed"] for case in report["cases"])
    report["sampled_changes_passed"] = all(case["sampled_changes_passed"] for case in report["cases"])
    report["formal_convergence_passed"] = None if quick else report["sampled_changes_passed"]
    save()
    return report


def exit_code(report, *, require_convergence=False):
    """正式门禁必须明确为True；快速或未完成报告都不能通过。"""
    if report.get("completed") is not True or report.get("numerical_checks_passed") is not True:
        return 1
    if require_convergence and (report.get("quick_mode") is not False or report.get("formal_convergence_passed") is not True):
        return 2
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="新的空输出目录，保存逐级输入、数组与总报告")
    parser.add_argument("--cases", nargs="+", choices=list(CASE_LABELS), help="选择工况；默认全部五案")
    parser.add_argument("--space-levels", type=int, nargs="+", default=[11, 21, 41], help="至少三级递增条带数")
    parser.add_argument("--phase-levels", type=int, nargs="+", default=[16, 32, 64], help="至少三级递增相位数")
    parser.add_argument("--fixed-phases", type=int, default=32, help="空间细化时固定的相位数")
    parser.add_argument("--fixed-strips", type=int, default=21, help="相位细化时固定的条带数")
    parser.add_argument("--relative-tolerance", type=float, default=.02, help="各相邻级变化的相对门槛")
    parser.add_argument("--absolute-tolerance", nargs=2, action="append", metavar=("METRIC", "SI_VALUE"),
                        help="按指标名覆盖SI绝对门槛，可重复；指标见DEFAULT_ABSOLUTE_TOLERANCES")
    parser.add_argument("--quick", action="store_true", help="固定小网格检查完整管道，正式收敛状态保持null")
    parser.add_argument("--require-convergence", action="store_true", help="要求正式收敛明确通过，否则非零退出")
    args = parser.parse_args(argv)
    try:
        absolute = {}
        for name, value in args.absolute_tolerance or []:
            if name in absolute:
                raise ValueError("同一指标不能重复指定绝对门槛：" + name)
            absolute[name] = float(value)
        report = run(args.output, cases=args.cases, space_levels=args.space_levels, phase_levels=args.phase_levels,
                     fixed_phases=args.fixed_phases, fixed_strips=args.fixed_strips,
                     relative_tolerance=args.relative_tolerance, absolute_tolerances=absolute, quick=args.quick)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps({"完成": report["completed"], "原数组校核通过": report["numerical_checks_passed"],
                      "正式收敛通过": report["formal_convergence_passed"],
                      "报告": str(args.output / "convergence_matrix.json")}, ensure_ascii=False))
    return exit_code(report, require_convergence=args.require_convergence)


if __name__ == "__main__":
    raise SystemExit(main())
