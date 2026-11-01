"""可移植输入、原始离散结果和SI单位契约的JSON/NPZ/CSV导出。"""
from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from .analysis import transmission_error_spectrum


def json_value(value):
    """保留数组、单位和结构；没有定义的非有限相位输出JSON null。"""
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"对象 {type(value).__name__} 没有可移植JSON表示")


def write_json(path, value):
    Path(path).write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2,
                                    allow_nan=False)+"\n", encoding="utf-8")


def export_result(result, output_dir):
    """每个相位的变长接触数组以offsets无损展开，禁止使用pickle对象数组。"""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    offsets = np.r_[0, np.cumsum([len(r.solution.forces_n) for r in result.rotations])]
    series = {"phase": result.phase, "normal_approach_m": result.normal_approach_m,
              "angular_te_rad": result.angular_te_rad,
              "mesh_stiffness_n_per_m": result.mesh_stiffness_n_per_m,
              "cell_offsets": offsets}
    cell_names = ("face_coordinate_m", "roll_pinion_m", "roll_wheel_m", "widths_m", "tooth_pair_index",
                  "curvature_per_m", "profile_factor_pinion", "profile_factor_wheel",
                  "lead_factor_pinion", "lead_factor_wheel")
    for key in cell_names:
        series[key] = np.concatenate([np.asarray(getattr(r.grid, key)) for r in result.rotations])
    for key in ("gaps_m", "line_load_n_per_m", "pressure_pa", "half_width_m"):
        series[key] = np.concatenate([np.asarray(getattr(r, key)) for r in result.rotations])
    for key in ("forces_n", "residual_gaps_m", "active_mask"):
        series[key] = np.concatenate([np.asarray(getattr(r.solution, key)) for r in result.rotations])
    extended = [r.contact_geometry_data for r in result.rotations]
    if any(value is not None for value in extended):
        if not all(value is not None for value in extended):
            raise ValueError("同一周期不能混合缺少延伸几何的相位")
        for name, attribute in (("eligible", "eligible"), ("extended_pinion", "centre_extended_pinion"),
                                ("extended_wheel", "centre_extended_wheel"), ("additional_gap_m", "gap_correction_m")):
            series[name] = np.concatenate([getattr(value, attribute) for value in extended])
        rotation_lines, line_cells, line_endpoints, endpoint_phases = [0], [0], [0], [0]
        for k, value in enumerate(extended):
            rotation_lines.append(rotation_lines[-1]+len(value.line_offsets)-1)
            line_cells.extend(offsets[k]+value.line_offsets[1:])
            line_endpoints.extend(endpoint_phases[-1]+value.endpoints.line_offsets[1:])
            endpoint_phases.append(endpoint_phases[-1]+len(value.endpoints.face_coordinate_m))
        series.update(rotation_line_offsets=np.asarray(rotation_lines), line_cell_offsets=np.asarray(line_cells),
                      line_endpoint_offsets=np.asarray(line_endpoints), rotation_endpoint_offsets=np.asarray(endpoint_phases))
        for name in ("face_coordinate_m", "roll_pinion_m", "roll_wheel_m", "tooth_pair_index", "curvature_per_m",
                     "radius_pinion_m", "radius_wheel_m", "extended_pinion", "extended_wheel"):
            series["endpoint_"+name] = np.concatenate([getattr(value.endpoints, name) for value in extended])
        for name in ("pressure_pa", "half_width_m", "line_load_n_per_m"):
            series["endpoint_"+name] = np.concatenate([r.endpoint_pressure[name] for r in result.rotations])
    np.savez_compressed(output/"arrays.npz", **series)
    write_json(output/"input.json", {"pair": result.pair, **result.inputs})
    write_json(output/"summary.json", result.summary)
    write_json(output/"diagnostics.json", {
        "contact": [r.solution.diagnostics for r in result.rotations],
        "elasticity": [r.elastic_diagnostics for r in result.rotations]})
    spectrum = transmission_error_spectrum(result.normal_approach_m,
        phase_origin=float(result.inputs.get("phase_offset", 0.0)))
    spectrum["sampling_phase_origin"] = float(result.inputs.get("phase_offset", 0.0))
    write_json(output/"te_spectrum.json", spectrum)
    with (output/"cycle.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        names = ("phase", "normal_approach_m", "angular_te_rad", "mesh_stiffness_n_per_m")
        writer.writerow(names)
        writer.writerows(zip(*(series[name] for name in names)))
    with (output/"contact_cells.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        names = ("phase", "tooth_pair_index", "face_coordinate_m", "roll_pinion_m", "roll_wheel_m", "widths_m",
                 "gaps_m", "forces_n", "line_load_n_per_m", "pressure_pa", "half_width_m",
                 "residual_gaps_m", "active_mask")
        additional_names = (("eligible", "extended_pinion", "extended_wheel", "additional_gap_m")
                            if extended and extended[0] is not None else ())
        writer.writerow((*names, *additional_names))
        phase = np.repeat(result.phase, np.diff(offsets))
        writer.writerows(zip(phase, *(series[name] for name in names[1:]),
                             *(series[name] for name in additional_names)))
    if extended and extended[0] is not None:
        with (output/"endpoint_cells.csv").open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            names = ("tooth_pair_index", "face_coordinate_m", "roll_pinion_m", "roll_wheel_m", "curvature_per_m",
                     "pressure_pa", "half_width_m", "line_load_n_per_m", "extended_pinion", "extended_wheel")
            writer.writerow(("phase", *names))
            phase = np.repeat(result.phase, np.diff(series["rotation_endpoint_offsets"]))
            writer.writerows(zip(phase, *(series["endpoint_"+name] for name in names)))
    files = ("arrays.npz", "input.json", "summary.json", "diagnostics.json", "te_spectrum.json",
             "cycle.csv", "contact_cells.csv")
    if extended and extended[0] is not None:
        files += ("endpoint_cells.csv",)
    manifest = {
        "module": "gear_contact", "version": "0.2.0", "implementation": "新写物理模型计算内核",
        "units": {"length": "m", "angle": "rad", "force": "N", "torque": "N m",
                  "pressure": "Pa", "compliance": "m/N", "normal_stiffness": "N/m"},
        "array_layout": "cell_offsets[k]:cell_offsets[k+1] 对应phase[k]的原始接触单元",
        "line_load_concentration": "每相位max(F_i/w_i)除以总力/理论接触线总长度；非标准评级系数",
        "te_reference": "法向接近量保留输入间隙基准；角度TE参考小轮虚功一致法向力臂",
        "files": {name: hashlib.sha256((output/name).read_bytes()).hexdigest() for name in files}}
    write_json(output/"manifest.json", manifest)
    return manifest
