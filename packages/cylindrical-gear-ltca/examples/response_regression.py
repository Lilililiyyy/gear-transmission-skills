"""新写数值回归工具：固定输入重算响应，与本项目旧版本自有基线比较。

基线只用于比较，不参与求解。通过表示本版本响应未超出声明的回归容差，
不改变 MASTA 数学问题等价性、标准评级或物理精度的验证状态。
普通检查不写回基线；更新基线必须另行审查源码身份、输入和原始结果。
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import sys
import time

import numpy as np
import scipy

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p / "pyproject.toml").is_file() and (p / "gear_contact").is_dir())
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from gear_contact import analyze_config, export_result
from gear_contact.io import write_json


# 所有容差按 SI 单位逐项声明。离散活动集合和索引必须完全相同。
FIELD_UNITS = {
    "mean_normal_approach_m": "m", "te_peak_to_peak_m": "m",
    "mean_mesh_stiffness_n_per_m": "N/m", "max_pressure_pa": "Pa",
    "max_endpoint_pressure_pa": "Pa", "phase": "base_pitch_period",
    "normal_approach_m": "m", "angular_te_rad": "rad",
    "mesh_stiffness_n_per_m": "N/m", "phase_max_pressure_pa": "Pa",
    "phase_max_endpoint_pressure_pa": "Pa",
}
EXACT_FIELDS = {"cell_offsets", "tooth_pair_index", "active_mask"}
SCALAR_FIELDS = {"mean_normal_approach_m", "te_peak_to_peak_m",
                 "mean_mesh_stiffness_n_per_m", "max_pressure_pa", "max_endpoint_pressure_pa"}
ENDPOINT_FIELDS = {"max_endpoint_pressure_pa", "phase_max_endpoint_pressure_pa"}
SHA256 = re.compile(r"[a-f0-9]{64}")


def canonical_sha256(value):
    """对实际输入的规范 JSON 绑定身份；拒绝 NaN 和无穷，不依赖键顺序。"""
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True,
                                    separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def environment():
    """记录计算环境；版本差异本身不等于数值回归失败。"""
    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "platform": platform.platform()}


def response_values(result):
    """读取真实求解结果；变长单元按相位和原单元索引展开，不排序或插值。"""
    scalar = SCALAR_FIELDS - ENDPOINT_FIELDS
    values = {name: result.summary[name] for name in scalar}
    values.update({name: np.asarray(getattr(result, name)).tolist() for name in
                   ("phase", "normal_approach_m", "angular_te_rad", "mesh_stiffness_n_per_m")})
    values["phase_max_pressure_pa"] = [float(np.max(r.pressure_pa)) for r in result.rotations]
    values["cell_offsets"] = np.r_[0, np.cumsum([
        len(r.solution.forces_n) for r in result.rotations])].tolist()
    values["active_mask"] = np.concatenate([r.solution.active_mask for r in result.rotations]).tolist()
    values["tooth_pair_index"] = np.concatenate([r.grid.tooth_pair_index for r in result.rotations]).tolist()
    if any(r.endpoint_pressure is not None for r in result.rotations):
        if not all(r.endpoint_pressure is not None for r in result.rotations):
            raise ValueError("同一周期的端点压力不能部分缺失")
        values["max_endpoint_pressure_pa"] = result.summary["max_endpoint_pressure_pa"]
        values["phase_max_endpoint_pressure_pa"] = [
            float(np.max(r.endpoint_pressure["pressure_pa"])) for r in result.rotations]
    return values


def _finite_array(value, name):
    """拒绝空值、字符串、布尔数及非有限响应，不能由 JSON null 转换掩盖失败。"""
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or not raw.size:
        raise ValueError(f"{name} 必须是非空有限数值，不能含布尔数、字符串或空值")
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} 含非有限值")
    return array


def _required_fields(config):
    return (set(FIELD_UNITS) | EXACT_FIELDS) - (
        set() if config.get("extended_contact_settings") is not None else ENDPOINT_FIELDS)


def _validate_values(values, config):
    """校验必需字段、单位索引对应的数组长度和布尔活动集合。"""
    if not isinstance(values, dict) or set(values) != _required_fields(config):
        raise ValueError("响应字段缺失或多余，必须与该工况的中心/端点定义一致")
    phases = config["phases"]
    for name in set(values) - EXACT_FIELDS:
        array = _finite_array(values[name], name)
        if array.shape != (() if name in SCALAR_FIELDS else (phases,)):
            raise ValueError(f"{name} 的形状与相位数或标量定义不一致")
    offsets = np.asarray(values["cell_offsets"])
    if (offsets.dtype.kind not in "iu" or offsets.shape != (phases + 1,)
            or offsets[0] != 0 or np.any(np.diff(offsets) <= 0)):
        raise ValueError("cell_offsets 必须为从零开始、逐相位严格递增的整数索引")
    size = int(offsets[-1])
    labels = np.asarray(values["tooth_pair_index"])
    active = np.asarray(values["active_mask"])
    if labels.dtype.kind not in "iu" or labels.shape != (size,):
        raise ValueError("tooth_pair_index 必须与原单元布局同长且为整数")
    if active.dtype.kind != "b" or active.shape != (size,):
        raise ValueError("active_mask 必须与原单元布局同长且为布尔数组")


def validate_baseline(baseline):
    """检查基线契约；缺字段或非法容差不允许被当作空比较通过。"""
    if not isinstance(baseline, dict) or baseline.get("schema_version") != 1:
        raise ValueError("响应基线 schema_version 必须为 1")
    source = baseline.get("source", {})
    if (not isinstance(source, dict) or not re.fullmatch(r"[a-f0-9]{40}", str(source.get("commit", "")))
            or not isinstance(source.get("kernel_sha256"), dict) or not source["kernel_sha256"]
            or any(not SHA256.fullmatch(str(v)) for v in source["kernel_sha256"].values())):
        raise ValueError("基线缺少有效源提交和内核 SHA-256 身份")
    env = baseline.get("environment", {})
    if not isinstance(env, dict) or any(not isinstance(env.get(k), str) or not env[k]
                                       for k in ("python", "numpy", "scipy", "platform")):
        raise ValueError("基线必须记录 Python、NumPy、SciPy 和平台版本")
    tolerances = baseline.get("tolerances", {})
    if not isinstance(tolerances, dict) or set(tolerances) != set(FIELD_UNITS):
        raise ValueError("基线必须为每个数值响应字段明确提供容差")
    for name, unit in FIELD_UNITS.items():
        spec = tolerances[name]
        if not isinstance(spec, dict) or set(spec) != {"unit", "absolute", "relative"} or spec["unit"] != unit:
            raise ValueError(f"{name} 容差单位或字段不符合契约")
        for key in ("absolute", "relative"):
            number = _finite_array(spec[key], f"{name}.{key}")
            if number.shape != () or float(number) < 0:
                raise ValueError(f"{name} 容差必须为非负有限标量")
        if spec["absolute"] + spec["relative"] <= 0:
            raise ValueError(f"{name} 容差至少有一项为正")
    cases = baseline.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("基线必须包含非空工况列表")
    names = set()
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("每个基线工况必须为对象")
        name = case.get("case_id")
        if (not isinstance(name, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) is None
                or name in names):
            raise ValueError("工况 ID 必须唯一，且为不含路径的 ASCII 文件名")
        names.add(name)
        config = case.get("config")
        if not isinstance(config, dict) or any(k not in config for k in (
                "pair", "torque_n_m", "phases", "strips", "elastic_mode", "local_contact")):
            raise ValueError(f"{name} 缺少显式输入、离散或模型定义")
        for key in ("phases", "strips"):
            if isinstance(config[key], bool) or not isinstance(config[key], int) or config[key] < 2:
                raise ValueError(f"{name}.{key} 必须为至少 2 的整数")
        if case.get("config_sha256") != canonical_sha256(config):
            raise ValueError(f"{name} 输入身份与基线哈希不一致")
        _validate_values(case.get("response"), config)


def compare_response(actual, expected, tolerances, config):
    """逐元素比较 abs(a-b) <= atol + rtol*abs(b)，近零由绝对容差控制。"""
    _validate_values(actual, config)
    _validate_values(expected, config)
    rows = []
    for name in sorted(expected):
        a, b = np.asarray(actual[name]), np.asarray(expected[name])
        if a.shape != b.shape:
            rows.append({"field": name, "passed": False, "reason": "数组形状改变",
                         "actual_shape": list(a.shape), "baseline_shape": list(b.shape)})
            continue
        if name in EXACT_FIELDS:
            mismatch = a != b
            rows.append({"field": name, "passed": not bool(np.any(mismatch)), "comparison": "exact",
                         "shape": list(a.shape), "mismatch_indices": np.argwhere(mismatch).tolist()})
        else:
            spec = tolerances[name]
            with np.errstate(over="ignore", invalid="ignore"):
                difference = a.astype(float) - b.astype(float)
                limit = spec["absolute"] + spec["relative"] * np.abs(b)
            if not np.all(np.isfinite(difference)) or not np.all(np.isfinite(limit)):
                raise ValueError(f"{name} 的差值或容差边界产生非有限值")
            passed = np.abs(difference) <= limit
            rows.append({"field": name, "unit": spec["unit"], "passed": bool(np.all(passed)),
                         "shape": list(a.shape), "absolute_tolerance": spec["absolute"],
                         "relative_tolerance": spec["relative"], "signed_difference": difference.tolist(),
                         "allowed_absolute_difference": np.asarray(limit).tolist(),
                         "maximum_absolute_difference": float(np.max(np.abs(difference))),
                         "failed_indices": np.argwhere(~passed).tolist()})
    return rows


def run(output, *, baseline_path=None, case_ids=None):
    """在新建或空目录重算选定工况，保留请求输入、原数组与逐项差异。

    拒绝复用非空目录，防止失败重跑遗留旧成功数组；不刷新或覆写基线。
    """
    started = time.perf_counter()
    baseline_path = Path(baseline_path or PROJECT / "examples/response_baseline.json").resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    validate_baseline(baseline)
    output = Path(output).resolve()
    if output == baseline_path or output in baseline_path.parents or baseline_path in output.parents:
        raise ValueError("回归输出不能覆盖基线文件及其父路径")
    by_id = {case["case_id"]: case for case in baseline["cases"]}
    selected = list(case_ids) if case_ids is not None else list(by_id)
    if not selected or len(set(selected)) != len(selected) or set(selected) - set(by_id):
        raise ValueError("需要非空、无重复且存在于基线的工况 ID")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("响应回归输出必须为新建或空目录，不能覆盖已有成功或失败记录")
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for name in selected:
        case = by_id[name]
        target = (output / name).resolve()
        if output not in target.parents:
            raise ValueError("工况输出路径必须保持在回归输出目录内")
        target.mkdir(parents=True, exist_ok=True)
        write_json(target / "requested_input.json", case["config"])
        try:
            # 只把独立输入交给内核；响应参照和容差不会进入求解参数。
            result = analyze_config(copy.deepcopy(case["config"]))
            actual = response_values(result)
            _validate_values(actual, case["config"])
            export_result(result, target)
            write_json(target / "response.json", actual)
            rows = compare_response(actual, case["response"], baseline["tolerances"], case["config"])
            converged = bool(result.summary["all_converged"])
            record = {"case_id": name, "config_sha256": case["config_sha256"],
                      "solver_converged": converged,
                      "response_regression_passed": converged and all(row["passed"] for row in rows),
                      "differences": rows}
        except (ValueError, RuntimeError, KeyError, TypeError, np.linalg.LinAlgError) as error:
            record = {"case_id": name, "config_sha256": case["config_sha256"], "response_regression_passed": False,
                      "error": f"{type(error).__name__}: {error}"}
        write_json(target / "response_comparison.json", record)
        records.append(record)
    report = {"schema_version": 1, "checked_at_utc": datetime.now(timezone.utc).isoformat(),
              "response_regression_passed": all(r["response_regression_passed"] for r in records),
              "case_count": len(records), "selected_case_ids": selected,
              "complete_baseline_checked": set(selected) == set(by_id),
              "baseline_sha256": hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
              "baseline_source": baseline["source"], "baseline_environment": baseline["environment"],
              "current_environment": environment(),
              "elapsed_seconds": time.perf_counter() - started,
              "current_kernel_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in sorted((PROJECT / "gear_contact").glob("*.py"))},
              "validation_tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "reference_problem_equivalence_assessed": False,
              "physical_or_experimental_validation": False,
              "scope": "固定输入与离散下的本项目版本间数值回归；不等于原组件等价、离散收敛或物理精度验证",
              "cases": records}
    write_json(output / "response_regression.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description="重算圆柱 LTCA 自有响应基线；不会自动更新基线")
    parser.add_argument("--output", type=Path, required=True, help="原数组和回归差异输出目录")
    parser.add_argument("--baseline", type=Path, help="经审查的基线 JSON，缺省使用仓库固定基线")
    parser.add_argument("--case", action="append", dest="case_ids", help="只检查具名工况，可重复；缺省检查全部")
    parser.add_argument("--require-match", action="store_true", help="任一响应回归失败时以状态码 2 退出")
    args = parser.parse_args(argv)
    try:
        report = run(args.output, baseline_path=args.baseline, case_ids=args.case_ids)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(2, f"响应回归输入无效：{error}\n")
    print(json.dumps({k: report[k] for k in (
        "response_regression_passed", "case_count", "complete_baseline_checked",
        "physical_or_experimental_validation")}, ensure_ascii=False, indent=2))
    if args.require_match and not report["response_regression_passed"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
