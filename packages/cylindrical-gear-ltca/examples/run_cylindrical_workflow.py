"""新写单案例验证入口：几何、显式因子强度、LTCA及同相位接触斑点。

所有本项目数值由当前参数化内核求得；参照文件只用于输入身份核验和
结果比较。原强度因子为外部输入，接触域及柔度与MASTA仍存在差异。
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p/"pyproject.toml").is_file() and (p/"gear_contact").is_dir())
for folder in (PROJECT, PROJECT/"examples"):
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))

from gear_contact import Material, analyze_config, export_result
from gear_contact.io import write_json, json_value
from gear_contact.strength import MeshRatingFactors, GearRatingFactors, StressLimits, evaluate_iso6336_explicit_factors
from run_benchmarks import common_values, normalized_config, comparison_row


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_case(case):
    """统一名义载荷记账，分别构造强度名义输入与LTCA评级后目标。"""
    config = deepcopy(case["analysis"])
    if "torque_n_m" in config:
        raise ValueError("单案例扭矩须由loading统一提供，不能在analysis中重复给出")
    loading = case["loading"]
    nominal = loading["nominal_pinion_torque_n_m"]
    multiplier = loading["ltca_torque_multiplier"]
    if (isinstance(nominal, bool) or not math.isfinite(nominal) or nominal < 0
            or isinstance(multiplier, bool) or not math.isfinite(multiplier) or multiplier <= 0):
        raise ValueError("名义扭矩须非负有限，LTCA目标乘子须正有限")
    config["torque_n_m"] = nominal*multiplier
    if config.get("elastic_mode") in ("iso", "iso_love"):
        settings = dict(config.get("iso_settings", {}))
        if "nominal_pinion_torque_n_m" in settings:
            raise ValueError("单案例ISO名义扭矩须由loading统一提供，不在iso_settings中重复给出")
        settings["nominal_pinion_torque_n_m"] = nominal
        config["iso_settings"] = settings
    args = deepcopy(case["strength"])
    args["nominal_pinion_torque_n_m"] = nominal
    args["mesh_factors"] = MeshRatingFactors(**args["mesh_factors"])
    for role in ("pinion", "wheel"):
        args[role+"_factors"] = GearRatingFactors(**args[role+"_factors"])
        if args.get(role+"_limits") is not None:
            args[role+"_limits"] = StressLimits(**args[role+"_limits"])
        args["material_"+role] = Material(**config.get("material_"+role, {}))
    return config, args


def calculation_inputs(case):
    config, args = prepare_case(case)
    return {"analysis": normalized_config(config), "strength": json_value(args),
            "instant": case["instant"]}


def calculate_case(case):
    """可编辑几何、材料、载荷及修形；不接收参照输出参与数学求值。"""
    config, args = prepare_case(case)
    ltca = analyze_config(config)
    strength = evaluate_iso6336_explicit_factors(ltca.pair, **args)
    return ltca, strength, config


def _field(value, name):
    for part in name.split("."):
        value = value[part]
    return value


def _safe_reference_path(directory, name):
    path = (directory/name).resolve()
    if directory not in path.parents:
        raise ValueError("参照关联文件必须位于该案例分发目录内")
    return path


def _protect_output(output, paths):
    """保护直接输入及自定义参照关联的全部文件，须在任何结果写入前执行。"""
    for path in paths:
        if output == path or path in output.parents or output in path.parents:
            raise ValueError("输出目录不能覆盖案例输入、MASTA对标数据或其父子目录")


def _write_csv(path, rows):
    import csv
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def comparison_rows(case, reference, ltca, strength):
    """共同量及应力组合用显式容差；LTCA差异不预设精度通过。"""
    label = {"case_id": case["case_id"], "label": case["label"]}
    values, rows = common_values(ltca), []
    for field, original in reference["geometry_values"].items():
        row = comparison_row(label, field, values[field], original, "common")
        row["group"] = "geometry"
        rows.append(row)
    for field, original in reference["strength_values"].items():
        actual = _field(strength, field)
        unit = ("sqrt(Pa)" if field == "elastic_factor_sqrt_pa" else
                "Pa" if field.endswith("_pa") else "N" if field.endswith("_n") else "1")
        atol, rtol = {"Pa": 1e-3, "sqrt(Pa)": 1e-7, "N": 1e-8, "1": 1e-10}[unit], 1e-8
        delta = abs(actual-original)
        rows.append({"case_id": case["case_id"], "label": case["label"], "group": "strength",
            "field": field, "unit": unit, "project": actual, "reference": original,
            "signed_difference": actual-original, "absolute_error": delta,
            "relative_error": delta/abs(original) if original else None,
            "absolute_tolerance": atol, "relative_tolerance": rtol, "matched": delta <= atol+rtol*abs(original)})
    for field, original in reference["ltca_values"].items():
        rows.append(comparison_row(label, field, ltca.summary[field], original, "ltca"))
    return rows


def run(case_path, output, *, reference_path=None, plots=False):
    """完整求解及有源数据的对比，输入变化不会继承原工况通过状态。"""
    case_path, output = Path(case_path).resolve(), Path(output).resolve()
    reference_path = Path(reference_path or PROJECT/"examples/helical_18_77_reference.json").resolve()
    directory = reference_path.parent
    protected = [case_path, reference_path, directory/"masta_reference"]
    _protect_output(output, protected)
    case, reference = read(case_path), read(reference_path)
    baseline_path = _safe_reference_path(directory, reference["case_file"])
    _protect_output(output, [baseline_path,
        *[_safe_reference_path(directory, name) for name in reference["files"]],
        _safe_reference_path(directory, reference["contact_metadata"]),
        _safe_reference_path(directory, reference["contact_arrays"])])
    baseline = read(baseline_path)
    manifest_path = directory/"helical_18_77_manifest.json"
    manifest = read(manifest_path) if manifest_path.is_file() else {}
    bound_files = manifest.get("files", {})
    binding = (bound_files.get(reference_path.name) == sha(reference_path)
        and bound_files.get(baseline_path.name) == sha(baseline_path)
        and all(sha(_safe_reference_path(directory, name)) == digest for name, digest in reference["files"].items()))
    contact_path = _safe_reference_path(directory, reference["contact_metadata"])
    arrays_path = _safe_reference_path(directory, reference["contact_arrays"])
    metadata = read(contact_path)
    if sha(arrays_path) != metadata["arrays_sha256"]:
        raise ValueError("原接触数组与已绑定元数据哈希不符，拒绝生成接触斑点比较")
    with np.load(arrays_path, allow_pickle=False) as source:
        arrays = {key: source[key].copy() for key in source.files}
    ltca, strength, config = calculate_case(case)
    same = calculation_inputs(case) == calculation_inputs(baseline)
    count_matches = len(ltca.phase) == len(arrays["rotation_phase"])
    phase_error = (float(np.max(np.abs((ltca.phase-arrays["rotation_phase"]+.5) % 1-.5)))
                   if count_matches else None)
    geometry_same = (normalized_config(config)["pair"] ==
                     calculation_inputs(baseline)["analysis"]["pair"])
    registered = count_matches and phase_error < 1e-10 and geometry_same
    if (not ltca.summary["all_converged"] or ltca.summary["max_force_residual_n"] > 1e-7
            or ltca.summary["max_gap_residual_m"] > 3e-12):
        raise RuntimeError("当前案例的接触平衡或协调残差未达到声明阈值")
    rows = comparison_rows(case, reference, ltca, strength)
    output.mkdir(parents=True, exist_ok=True)
    export_result(ltca, output/"ltca")
    write_json(output/"case_inputs.json", case)
    write_json(output/"strength.json", strength)
    write_json(output/"comparison.json", rows)
    _write_csv(output/"comparison.csv", rows)
    grouped = {group: [row for row in rows if row["group"] == group] for group in ("geometry", "strength")}
    group_pass = {group: bool(items) and binding and same and all(row["matched"] is True for row in items)
                  for group, items in grouped.items()}
    summary = {"case_id": case["case_id"], "label": case["label"],
        "checked_at_utc": datetime.now(timezone.utc).isoformat(), "own_values_recomputed": True,
        "reference_own_values_used": False, "masta_runtime_required": False,
        "reference_binding_verified": binding, "same_inputs_as_reference": same,
        "geometry_comparison_passed": group_pass["geometry"], "strength_comparison_passed": group_pass["strength"],
        "geometry_fields": sum(r["group"] == "geometry" for r in rows),
        "strength_fields": sum(r["group"] == "strength" for r in rows),
        "phase_geometry_registered": registered, "maximum_phase_residual_cycles": phase_error,
        "geometry_inputs_match_reference": geometry_same,
        "phase_registration_basis": metadata["phase_registration"],
        "ltca_response_accuracy_verified": False, "physical_validation": False,
        "independent_strength_factor_prediction": False,
        "loading": {"nominal_pinion_torque_n_m": strength["nominal_pinion_torque_n_m"],
            "nominal_tangential_force_n": strength["nominal_tangential_force_n"],
            "ltca_torque_n_m": config["torque_n_m"], "ltca_target_normal_force_n": ltca.summary["normal_force_n"],
            "ltca_torque_multiplier": case["loading"]["ltca_torque_multiplier"]},
        "ltca_summary": ltca.summary, "strength": strength, "instant": case["instant"],
        "limitations": reference["limitations"],
        "kernel_sha256": {p.name: sha(p) for p in sorted((PROJECT/"gear_contact").glob("*.py"))},
        "validation_tool_sha256": sha(Path(__file__)),
        "input_sha256": sha(case_path), "reference_sha256": sha(reference_path),
        "contact_arrays_sha256": sha(arrays_path)}
    if plots:
        if not registered:
            raise ValueError("当前起采相位或相位数不同，不能生成同相位接触对比；请按参照几何配准")
        from case_contact_plotting import plot_case_contact
        from case_rating_plotting import plot_case_ratings
        summary["rating_figures"] = plot_case_ratings(rows, output)
        selected_metadata = deepcopy(metadata)
        selected_metadata["suggested_instant"] = deepcopy(case["instant"])
        summary["contact_figures"] = plot_case_contact(ltca, selected_metadata, arrays, output)
        if config.get("extended_contact_settings") is not None:
            from case_sensitivity import run_case_sensitivity
            summary["sensitivity"] = run_case_sensitivity(case, output, plots=True)
            if reference.get("native_contact_metadata"):
                from case_native_contact_plotting import plot_native_contact_comparison
                native = read(_safe_reference_path(directory, reference["native_contact_metadata"]))
                summary["native_contact_figures"] = plot_native_contact_comparison(ltca, native,
                    _safe_reference_path(directory, reference["native_contact_image"]),
                    _safe_reference_path(directory, reference["native_contact_arrays"]), output)
    write_json(output/"workflow_summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", nargs="?", type=Path, default=PROJECT/"examples/helical_18_77_case.json")
    parser.add_argument("--reference", type=Path, help="MASTA对标数据，不作为本项目计算输入")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plots", action="store_true", help="生成几何、强度、周期和接触斑点PNG/SVG")
    args = parser.parse_args()
    summary = run(args.case, args.output, reference_path=args.reference, plots=args.plots)
    print(json.dumps({key: summary[key] for key in ("case_id", "geometry_comparison_passed",
        "strength_comparison_passed", "phase_geometry_registered", "ltca_response_accuracy_verified")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
