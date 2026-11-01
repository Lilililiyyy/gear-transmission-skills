"""新写验证工具：重算工程案例、相似律及指定工况的离散敏感性。

只调用本项目内核。案例为确定性合成输入；方程残差、
相似律和离散敏感性分别登记，不能据此宣称实物或三维接触精度。
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p / "pyproject.toml").is_file() and (p / "gear_contact").is_dir())
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from gear_contact import Gear, GearPair, Material, SurfaceModification, analyze_pair, export_result
from gear_contact.io import write_json


def standard_pair(*, scale=1., beta=0., internal=False):
    """28/56 外啮合或 28/84 内啮合；长度缩放不改变无量纲几何。"""
    return GearPair(Gear(28, .0025 * scale, .025 * scale, helix_angle_rad=beta),
                    Gear(84 if internal else 56, .0025 * scale, .025 * scale,
                         helix_angle_rad=beta if internal else -beta, internal=internal))


def check_result(result):
    """直接复核原数组的平衡、互补、虚功和 Hertz 压力积分恒等式。

    使用各次求解声明的间隙容差，另给出绝对残差；零载荷时切线刚度
    不作正值要求，因为局部线接触的零载右切线可为零。
    """
    force = result.summary["normal_force_n"]
    force_residual = max(abs(float(np.sum(r.solution.forces_n)) - force) for r in result.rotations)
    penetration = max(max(0., -float(np.min(r.solution.residual_gaps_m))) for r in result.rotations)
    active_error = max((float(np.max(np.abs(r.solution.residual_gaps_m[r.solution.active_mask])))
                        if np.any(r.solution.active_mask) else 0.) for r in result.rotations)
    complementarity = max(float(np.max(np.abs(r.solution.forces_n * r.solution.residual_gaps_m)))
                          for r in result.rotations)
    gap_tolerance = max(3e-12, *(1.1 * r.solution.diagnostics.get("gap_tolerance_m", 0.)
                               for r in result.rotations))
    work_a = force * result.normal_approach_m
    work_b = result.inputs["torque_n_m"] * result.angular_te_rad
    work_residual = float(np.max(np.abs(work_a - work_b)))
    hertz_error = max(float(np.max(np.abs(math.pi * r.pressure_pa * r.half_width_m / 2
                                         - r.line_load_n_per_m))) for r in result.rotations)
    q_scale = max(1., result.summary["max_line_load_n_per_m"])
    checks = {
        "all_converged": all(r.solution.diagnostics["converged"] for r in result.rotations),
        "nonnegative_forces": all(np.all(r.solution.forces_n >= 0) for r in result.rotations),
        "finite_response": all(np.all(np.isfinite(a)) for a in (
            result.normal_approach_m, result.angular_te_rad, result.mesh_stiffness_n_per_m))
            and all(np.all(np.isfinite(r.pressure_pa)) for r in result.rotations),
        "force_balance": force_residual <= max(1e-8, force * 1e-10),
        "unilateral_gaps": penetration <= gap_tolerance,
        "active_compatibility": active_error <= gap_tolerance,
        "complementarity": complementarity <= max(force, 1.) * gap_tolerance,
        "virtual_work": work_residual <= max(1e-13, float(np.max(np.abs(work_a))) * 1e-12),
        "hertz_integral": hertz_error <= q_scale * 1e-12,
        "nonnegative_stiffness": bool(np.all(result.mesh_stiffness_n_per_m >= 0)),
    }
    return {"passed": all(checks.values()), "checks": {k: bool(v) for k, v in checks.items()},
            "force_residual_n": force_residual, "maximum_penetration_m": penetration,
            "active_gap_error_m": active_error, "complementarity_n_m": complementarity,
            "gap_tolerance_m": gap_tolerance, "virtual_work_residual_j": work_residual,
            "hertz_integral_residual_n_per_m": hertz_error}


def relative_array_error(actual, expected):
    """用被比数组的最大幅值归一化，避免近零局部值放大相对误差。"""
    a, b = np.asarray(actual), np.asarray(expected)
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), np.finfo(float).tiny))


def run(output, *, quick=False, resolution_tolerance=.02):
    """保存每案可重算输入与原数组；正式模式对末两级敏感性作显式判定。

    quick 仅用于接口验收，不计算空间/相位细化，formal_acceptance_passed
    保持 null。返回的 all_checks_passed 和正式验收状态不得混用。
    """
    if not math.isfinite(resolution_tolerance) or resolution_tolerance <= 0:
        raise ValueError("离散指标相对变化阈值必须为有限正数")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    phases, strips = (4, 5) if quick else (32, 21)
    results, rows, checks = {}, [], []
    crown = SurfaceModification(lead_crown_m=8e-6)
    combined = SurfaceModification(lead_crown_m=8e-6, tip_relief_m=6e-6,
                                   tip_start=.72, root_relief_m=3e-6, root_end=.18)

    def calculate(case_id, label, pair=None, **options):
        params = {"torque_n_m": 80., "phases": phases, "strips": strips, **options}
        result = analyze_pair(standard_pair() if pair is None else pair, **params)
        folder = output / case_id
        export_result(result, folder)
        validation = check_result(result)
        write_json(folder / "validation.json", validation)
        row = {"case": case_id, "label": label, "passed": validation["passed"],
               "phases": params["phases"], "strips": params["strips"],
               "normal_force_n": result.summary["normal_force_n"],
               "te_peak_to_peak_m": result.summary["te_peak_to_peak_m"],
               "max_pressure_pa": result.summary["max_pressure_pa"],
               "mean_mesh_stiffness_n_per_m": result.summary["mean_mesh_stiffness_n_per_m"],
               "force_residual_n": validation["force_residual_n"],
               "gap_residual_m": max(validation["maximum_penetration_m"], validation["active_gap_error_m"])}
        rows.append(row)
        results[case_id] = result
        return result

    calculate("spur", "直齿基准")
    base = calculate("combined", "组合修形", modification_pinion=combined)
    calculate("helical", "20° 外啮合斜齿", standard_pair(beta=math.radians(20)),
              modification_pinion=crown, misalignment_m=3e-6)
    calculate("internal", "28/84 内啮合直齿", standard_pair(internal=True))
    calculate("internal_helical", "28/84 内啮合斜齿", standard_pair(internal=True, beta=.18),
              modification_pinion=crown, misalignment_m=3e-6)
    other_pair = GearPair(Gear(25, .0032, .031, helix_angle_rad=.27, profile_shift=.2),
                          Gear(67, .0032, .037, helix_angle_rad=-.27, profile_shift=.1),
                          axial_offset_m=.003)
    calculate("nondefault", "非默认几何、异材料与载荷", other_pair, torque_n_m=137.,
              material_pinion=Material(180e9, .28), material_wheel=Material(105e9, .32),
              modification_pinion=combined, misalignment_m=-7e-6)
    positive = calculate("misalignment_positive", "鼓形叠加正失准", modification_pinion=crown, misalignment_m=12e-6)
    negative = calculate("misalignment_negative", "鼓形叠加负失准", modification_pinion=crown, misalignment_m=-12e-6)
    zero = calculate("zero_load", "零载荷接触开闭", torque_n_m=0., modification_pinion=combined)
    checks.append({"check": "zero_load", "passed": zero.summary["max_pressure_pa"] == 0.
                   and all(np.all(r.solution.forces_n == 0.) for r in zero.rotations),
                   "scope": "检查零载力和压力；不将零载刚度强制为正"})
    mirror_errors = {"approach": relative_array_error(negative.normal_approach_m, positive.normal_approach_m),
                     "stiffness": relative_array_error(negative.mesh_stiffness_n_per_m, positive.mesh_stiffness_n_per_m),
                     "forces": max(relative_array_error(b.solution.forces_n.reshape(-1, strips),
                                                        a.solution.forces_n.reshape(-1, strips)[:, ::-1])
                                   for a, b in zip(positive.rotations, negative.rotations))}
    checks.append({"check": "misalignment_mirror", "passed": max(mirror_errors.values()) <= 2e-6,
                   "relative_tolerance": 2e-6, "relative_errors": mirror_errors,
                   "scope": "对称直齿、对称鼓形；翻转齿向条带，齿对顺序不变"})

    # E 与 T 同乘 c 时应变不变：位移相同，力/压力/切线刚度均乘 c。
    material_factor = 1.4
    material_scaled = calculate("material_scaling", "模量与载荷同比缩放", torque_n_m=80 * material_factor,
                                material_pinion=Material(210e9 * material_factor, .3),
                                material_wheel=Material(210e9 * material_factor, .3), modification_pinion=combined)
    material_errors = {
        "approach": relative_array_error(material_scaled.normal_approach_m, base.normal_approach_m),
        "stiffness": relative_array_error(material_scaled.mesh_stiffness_n_per_m,
                                           base.mesh_stiffness_n_per_m * material_factor),
        "pressure": max(relative_array_error(b.pressure_pa, a.pressure_pa * material_factor)
                        for a, b in zip(base.rotations, material_scaled.rotations))}
    checks.append({"check": "material_load_similarity", "passed": max(material_errors.values()) <= 2e-6,
                   "relative_tolerance": 2e-6, "relative_errors": material_errors,
                   "modulus_and_torque_factor": material_factor})
    scale = 1.7
    scaled_surface = SurfaceModification(lead_crown_m=8e-6 * scale, tip_relief_m=6e-6 * scale,
                                         tip_start=.72, root_relief_m=3e-6 * scale, root_end=.18)
    larger = calculate("geometry_scaling", "长度与扭矩相似缩放", standard_pair(scale=scale),
                       torque_n_m=80 * scale**3, modification_pinion=scaled_surface)
    scale_errors = {"approach": relative_array_error(larger.normal_approach_m, base.normal_approach_m * scale),
                    "angular_te": relative_array_error(larger.angular_te_rad, base.angular_te_rad),
                    "stiffness": relative_array_error(larger.mesh_stiffness_n_per_m, base.mesh_stiffness_n_per_m * scale),
                    "pressure": max(relative_array_error(b.pressure_pa, a.pressure_pa)
                                    for a, b in zip(base.rotations, larger.rotations))}
    checks.append({"check": "geometry_similarity", "passed": max(scale_errors.values()) <= 2e-6,
                   "relative_tolerance": 2e-6, "relative_errors": scale_errors,
                   "length_factor": scale, "torque_factor": scale**3,
                   "scope": "长度及修形乘s、扭矩乘s³；压力相同、位移及刚度乘s、角度TE相同"})

    resolution = {"performed": not quick, "relative_tolerance": resolution_tolerance,
                  "scope": "仅28/56直齿组合修形、80 N·m；末两级变化不是模型误差或全参数域误差界"}
    if not quick:
        fields = ("te_peak_to_peak_m", "max_pressure_pa", "mean_mesh_stiffness_n_per_m")
        for name, levels, fixed in (("space", (11, 21, 41), {"phases": 32}),
                                    ("phase", (16, 32, 64), {"strips": 21})):
            level_results = []
            for count in levels:
                variable = {"strips": count} if name == "space" else {"phases": count}
                if count == 21 and name == "space" or count == 32 and name == "phase":
                    result = base
                else:
                    result = calculate(f"resolution_{name}_{count}", f"{name}细化{count}",
                                       modification_pinion=combined, **fixed, **variable)
                level_results.append(result)
            changes = {field: abs(level_results[-1].summary[field] - level_results[-2].summary[field])
                             / max(abs(level_results[-1].summary[field]), np.finfo(float).tiny)
                       for field in fields}
            resolution[name] = {"levels": list(levels), "fixed": fixed,
                                "metrics": [{field: r.summary[field] for field in fields} for r in level_results],
                                "last_two_relative_changes": changes,
                                "passed": all(v <= resolution_tolerance for v in changes.values())}

    passed = all(row["passed"] for row in rows) and all(check["passed"] for check in checks)
    if not quick:
        passed = passed and resolution["space"]["passed"] and resolution["phase"]["passed"]
    record = {"generated_utc": datetime.now(timezone.utc).isoformat(), "implementation": "新写验证工具",
              "quick_mode": quick, "all_checks_passed": bool(passed),
              "formal_acceptance_passed": None if quick else bool(passed),
              "input_origin": "确定性参数化合成案例", "case_count": len(rows), "cases": rows,
              "cross_case_checks": checks, "resolution": resolution,
              "physical_or_experimental_validation": False,
              "python": sys.version.split()[0], "numpy": np.__version__,
              "kernel_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted((PROJECT / "gear_contact").glob("*.py"))},
              "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "elapsed_seconds": time.perf_counter() - started}
    write_json(output / "case_validation.json", record)
    with (output / "case_validation.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="输入、原数组与验收报告的输出目录")
    parser.add_argument("--quick", action="store_true", help="小网格接口复核，不执行正式离散验收")
    parser.add_argument("--resolution-tolerance", type=float, default=.02, help="末两级指标变化上限，默认0.02")
    args = parser.parse_args()
    record = run(args.output, quick=args.quick, resolution_tolerance=args.resolution_tolerance)
    print(json.dumps({"案例数": record["case_count"], "检查通过": record["all_checks_passed"],
                      "正式验收通过": record["formal_acceptance_passed"],
                      "报告": str(args.output / "case_validation.json")}, ensure_ascii=False))
    if not record["all_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
