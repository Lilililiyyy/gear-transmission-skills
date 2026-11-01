"""新写验证工具：当前圆柱内核求解及MASTA分层对标。

本项目响应只由config调用analyze_config求得；MASTA仅用于对标。
共同定义字段采用显式容差，模型响应差异不自动变成精度验收通过。
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import time

import numpy as np

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p / "pyproject.toml").is_file() and (p / "gear_contact").is_dir())
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from gear_contact import Gear, GearPair, Material, analyze_config, modification_from_config, export_result
from gear_contact.io import json_value, write_json


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalized_config(config):
    """展开数据类缺省值，比较实际数学输入；不读取任何参照响应。"""
    pair = dict(config["pair"])
    pair["pinion"], pair["wheel"] = Gear(**pair["pinion"]), Gear(**pair["wheel"])
    values = {"pair": asdict(GearPair(**pair)), "torque_n_m": config["torque_n_m"]}
    for name in ("material_pinion", "material_wheel"):
        values[name] = asdict(Material(**config.get(name, {})))
    for name in ("modification_pinion", "modification_wheel"):
        values[name] = json_value(modification_from_config(config.get(name, {})))
    for name, default in (("phases", 32), ("strips", 21), ("misalignment_m", 0.),
                          ("phase_offset", 0.),
                          ("elastic_mode", "coupled"), ("local_contact", True),
                          ("local_reference_depths_m", None)):
        values[name] = config.get(name, default)
    if config.get("iso_settings") is not None:
        from gear_contact.iso_contact import ISOContactSettings
        values["iso_settings"] = asdict(ISOContactSettings(**config["iso_settings"]))
    if config.get("extended_contact_settings") is not None:
        from gear_contact.extended_contact import ExtendedContactSettings
        values["extended_contact_settings"] = asdict(ExtendedContactSettings(**config["extended_contact_settings"]))
    return json_value(values)


def common_values(result):
    """仅映射可由当前几何/载荷接口实际求得的量，保持 SI 单位。"""
    pair = result.pair
    values = {key: float(getattr(pair, key)) for key in
              ("center_distance_m", "base_pitch_m", "working_pressure_angle_rad",
               "base_helix_angle_rad", "face_overlap_m", "transverse_contact_ratio", "overlap_ratio")}
    values["normal_force_n"] = result.summary["normal_force_n"]
    for role in ("pinion", "wheel"):
        for key in ("reference_radius_m", "base_radius_m", "transverse_pressure_angle_rad",
                    "tip_radius_m", "root_radius_m", "root_form_radius_m", "tip_form_radius_m"):
            values[f"{role}.{key}"] = float(getattr(getattr(pair, role), key))
    return values


def field_unit(field):
    if field.endswith("_n_per_m"):
        return "N_per_m"
    if field.endswith("_pa"):
        return "Pa"
    if field.endswith("_rad"):
        return "rad"
    if field.endswith("_m"):
        return "m"
    if field.endswith("_n"):
        return "N"
    return "1"


def comparison_row(case, field, project, reference, group):
    """零参照的相对差为空；没有判据的模型差异不标为通过。"""
    if not np.isfinite(project) or not np.isfinite(reference):
        raise ValueError(f"{case['case_id']} 的 {field} 不是有限数，不能参与对标")
    unit = field_unit(field)
    absolute = abs(project - reference)
    relative = absolute / abs(reference) if reference != 0 else None
    # 容差用于公式/数值共同量；不是机械模型的工程精度界限。
    atol = {"m": 1e-11, "rad": 1e-10, "N": 1e-8, "1": 1e-9}.get(unit, 1e-10) if group == "common" else None
    rtol = 1e-8 if group == "common" else None
    matched = absolute <= atol + rtol * abs(reference) if group == "common" else None
    return {"case_id": case["case_id"], "label": case["label"], "group": group,
            "field": field, "unit": unit, "project": float(project), "reference": float(reference),
            "signed_difference": float(project-reference), "absolute_error": float(absolute),
            "relative_error": relative, "absolute_tolerance": atol, "relative_tolerance": rtol,
            "matched": matched}


def _inside(path, directory):
    return path == directory or directory in path.parents


def run(output, *, reference_dir=None, cases_path=None, case_ids=None, phases=None, strips=None, plots=False):
    """实际求解并导出每案原数组，响应不同仍完整保留结果。

    comparison_passed 只有同问题且所有字段有判据并通过才为真；本批
    MASTA 存在模型差异，程序正常执行与完整工程精度通过是两种状态。
    """
    reference_dir = Path(reference_dir or PROJECT / "examples/masta_reference").resolve()
    output = Path(output).resolve()
    if _inside(output, reference_dir) or _inside(reference_dir, output):
        raise ValueError("输出目录不能覆盖MASTA对标数据目录及其父子目录")
    reference_path = reference_dir / "cases.json"
    manifest_path = reference_dir / "manifest.json"
    baseline = read_json(reference_path)
    manifest = read_json(manifest_path) if manifest_path.is_file() else {}
    identities = manifest.get("files", {})
    if isinstance(identities, list):
        identities = {item["file"]: item["sha256"] for item in identities}
    binding = identities.get("cases.json") == sha256(reference_path)
    supplied = read_json(cases_path) if cases_path else baseline
    originals = {c["case_id"]: c for c in baseline["cases"]}
    supplied_cases = {c["case_id"]: c for c in supplied["cases"]}
    if len(originals) != len(baseline["cases"]) or len(supplied_cases) != len(supplied["cases"]):
        raise ValueError("案例ID重复，无法唯一绑定参照")
    selected = list(case_ids) if case_ids else list(originals)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("必须选择至少一个且不重复的案例")
    if set(selected)-originals.keys() or set(selected)-supplied_cases.keys():
        raise ValueError("所选案例缺少MASTA对标数据或本项目计算输入")
    # 自定义参照也必须受相同输出边界约束，不能让案例ID变成路径。
    for case_id in selected:
        if not isinstance(case_id, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", case_id) is None:
            raise ValueError("案例ID须为单一ASCII文件名，仅含字母、数字、点、下划线或连字符")
        target = (output / case_id).resolve()
        if (target == output or output not in target.parents
                or _inside(target, reference_dir) or _inside(reference_dir, target)):
            raise ValueError("案例输出必须严格位于结果目录内，且不能覆盖MASTA对标数据")
    output.mkdir(parents=True, exist_ok=True)
    rows, records, same, equivalent = [], [], {}, {}
    started = time.monotonic()
    for case_id in selected:
        case = originals[case_id]
        config = dict(supplied_cases[case_id]["config"])
        if phases is not None:
            config["phases"] = phases
        if strips is not None:
            config["strips"] = strips
        result = analyze_config(config)
        export_result(result, output / case_id)
        same[case_id] = normalized_config(config) == normalized_config(case["config"])
        comparability = case["comparability"]
        equivalent[case_id] = same[case_id] and comparability.get("input_mapping_complete", False) and comparability.get("ltca_model_equivalent", False)
        if (not result.summary["all_converged"] or result.summary["max_force_residual_n"] > 1e-7
                or result.summary["max_gap_residual_m"] > 3e-12):
            raise RuntimeError(f"{case_id} 的当前内核求解未达到声明的残差要求")
        values = common_values(result)
        for field, reference in case["reference_common_values"].items():
            rows.append(comparison_row(case, field, values[field], reference, "common"))
        for field, reference in case["reference_ltca_metrics"].items():
            rows.append(comparison_row(case, field, result.summary[field], reference, "ltca"))
        records.append({"case_id": case_id, "config": config, "project_summary": result.summary,
                        "common_values": values, "comparability": comparability})
    common = [row for row in rows if row["group"] == "common"]
    response = [row for row in rows if row["group"] == "ltca"]
    summary = {"checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "case_count": len(records), "project_rotations": sum(r["config"].get("phases", 32) for r in records),
        "own_values_recomputed": True, "reference_own_values_used": False,
        "masta_runtime_required": False, "reference_binding_verified": binding,
        "same_inputs_as_reference": same,
        "input_check_scope": "同输入指案例config；完整数学问题等效性另由full_reference_problem_equivalent报告",
        "full_reference_problem_equivalent": equivalent,
        "common_fields": len(common), "common_fields_matched": sum(r["matched"] is True for r in common),
        "common_comparison_passed": bool(common) and binding and all(same.values()) and all(r["matched"] for r in common),
        "response_difference_fields": len(response), "response_accuracy_verified": False,
        "maximum_response_relative_errors": {field: max(r["relative_error"] for r in response
            if r["field"] == field and r["relative_error"] is not None)
            for field in sorted({r["field"] for r in response if r["relative_error"] is not None})},
        "comparison_values_within_tolerance": all(r["matched"] is True for r in rows),
        "comparison_passed": binding and all(equivalent.values()) and all(r["matched"] is True for r in rows),
        "physical_or_experimental_validation": False, "elapsed_seconds": time.monotonic()-started,
        "reference_sha256": sha256(reference_path),
        "manifest_sha256": sha256(manifest_path) if manifest_path.is_file() else None,
        "kernel_sha256": {p.name: sha256(p) for p in sorted((PROJECT/"gear_contact").glob("*.py"))},
        "validation_tool_sha256": sha256(Path(__file__)),
        "scope": "共同几何和法向载荷采用显式容差；按各案例config的几何相位独立求解，汇总TE、端点峰压和刚度差异。响应未设工程精度判据，comparison_passed不代表本次计算是否收敛或数值是否接近。"}
    write_json(output / "case_inputs_and_results.json", records)
    write_json(output / "comparison.json", rows)
    with (output / "comparison.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if plots:
        # 非绘图路径只依赖计算内核及其已有数值依赖。
        from benchmark_plotting import plot_comparison
        summary["figures"] = plot_comparison(rows, output)
    write_json(output / "benchmark_summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="独立结果目录，不得覆盖MASTA对标数据")
    parser.add_argument("--reference-dir", type=Path, help="包含cases.json和manifest.json的参照目录")
    parser.add_argument("--cases", type=Path, help="同结构JSON，可修改config检验参数响应")
    parser.add_argument("--case-id", action="append", help="只运行指定案例，可多次给出")
    parser.add_argument("--phases", type=int, help="覆盖自身相位数，报告会标记映射输入变化")
    parser.add_argument("--strips", type=int, help="覆盖自身条带数，报告会标记映射输入变化")
    parser.add_argument("--plots", action="store_true", help="另外生成PNG/SVG/PDF，需要matplotlib")
    parser.add_argument("--require-common-match", action="store_true", help="共同字段或身份核验失败时以非零退出")
    args = parser.parse_args()
    result = run(args.output, reference_dir=args.reference_dir, cases_path=args.cases,
        case_ids=args.case_id, phases=args.phases, strips=args.strips, plots=args.plots)
    print(json.dumps({k: result[k] for k in ("case_count", "common_fields", "common_fields_matched",
        "common_comparison_passed", "response_difference_fields", "comparison_passed")}, ensure_ascii=False, indent=2))
    if args.require_common_match and not result["common_comparison_passed"]:
        raise SystemExit("共同量验收未通过，请检查输入、参照身份及逐字段差异")


if __name__ == "__main__":
    main()
