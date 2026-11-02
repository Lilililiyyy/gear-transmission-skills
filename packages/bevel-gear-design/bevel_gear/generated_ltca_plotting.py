"""新写的实际展成FE承载结果绘图与源数据导出；不运行几何或接触求解。

压力仅为实际未变形小轮FE三角片上的f/A，面积取三维节点，
有符号扭矩与从动弹性角不取绝对值。没有外部LTCA参照或Hertz峰值认证。
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


def _finite(value, name):
    if isinstance(value, bool):
        raise ValueError(name + "须为有限实数")
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(name + "须为有限实数") from error
    if not math.isfinite(value):
        raise ValueError(name + "须为有限实数")
    return value


def _close(first, second, name, *, relative=1e-9, absolute=1e-12):
    if not math.isclose(first, second, rel_tol=relative, abs_tol=absolute):
        raise ValueError(name + "与原始数据或实际节点不一致")


def _metadata(metadata):
    if not isinstance(metadata, dict):
        raise ValueError("须提供显式制造/齿坯/材料与网格身份元数据")
    for member in ("pinion", "wheel"):
        values = metadata["case"]["members"][member]
        material = values["material"]
        if _finite(material["young_modulus_pa"], "弹性模量Pa") <= 0 or \
           not -1 < _finite(material["poisson_ratio"], "泊松比") < .5:
            raise ValueError("材料须满足三维各向同性线弹性条件")
        if _finite(values["blank"]["bore_diameter_m"], "孔径m") <= 0:
            raise ValueError("孔壁固定模型须显式给出正孔径")
    return metadata


def _result(result):
    if result.get("converged") is not True or result.get("equilibrium", {}).get("converged") is not True:
        raise ValueError("未收敛结果不能进入成功承载图库")
    if result.get("preparation", {}).get("full_gears") is not True:
        raise ValueError("本图入口要求实际整轮FE模型，不能将扇区结果标为整轮")
    torque = _finite(result["pinion_torque_nm"], "主动有符号扭矩Nm")
    target = _finite(result["wheel_target_torque_nm"], "从动目标有符号扭矩Nm")
    theta = _finite(result["elastic_follower_lag_rad"], "从动有符号弹性角rad")
    stiffness = _finite(result["incremental_mesh_stiffness_nm_per_rad"], "从动轴增量刚度Nm/rad")
    _finite(result["driving_angle_rad"], "主动姿态rad")
    angle_parts = [_finite(result[k], k) for k in ("unloaded_tca_te_rad", "mesh_unloaded_correction_rad", "total_follower_te_rad")]
    _close(angle_parts[2], angle_parts[0]+angle_parts[1]+theta, "总TE角分量", absolute=1e-13)
    if torque == 0 or target*torque <= 0 or target*theta <= 0 or stiffness <= 0:
        raise ValueError("扭矩/弹性角符号或增量刚度不满足本承载分支的被动弹性关系")
    direction = 1 if torque > 0 else -1
    if result.get("torque_direction") != direction:
        raise ValueError("扭矩方向标记与原有符号载荷不一致")
    contacts = result["contacts"]
    if len(contacts) != result["candidate_count"] or not contacts:
        raise ValueError("候选计数与接触数组不一致")
    active, normal_force, first_torque, second_torque = [], 0., 0., 0.
    for index, contact in enumerate(contacts):
        force = _finite(contact["force_n"], "非负法向力N")
        area = _finite(contact["area_m2"], "实际面片面积m²")
        pressure = _finite(contact["mean_pressure_pa"], "面片平均压力Pa")
        if force < 0 or area <= 0 or pressure < 0:
            raise ValueError("法向力/均压须非负且实际面片面积须为正")
        _close(pressure, force/area, "f/A平均压力", absolute=1e-5)
        first, second = (_finite(contact[k], "有符号轴力臂m") for k in ("pinion_lever_m", "wheel_lever_m"))
        if direction*first <= 0 or direction*second <= 0:
            raise ValueError("工作侧有符号力臂不支持该扭矩方向，不能对力臂取绝对值")
        _finite(contact["signed_gap_m"], "原始有符号间隙m")
        remaining = _finite(contact["remaining_gap_m"], "剩余间隙m")
        gap_tolerance = result["equilibrium"].get("diagnostics", {}).get("gap_acceptance_m", 1e-8)
        if remaining < -gap_tolerance:
            raise ValueError("剩余负间隙超过解的接受容差")
        if contact["active"]:
            if force <= 0:
                raise ValueError("活跃面片必须具有正法向力")
            active.append(index)
            if abs(remaining) > gap_tolerance:
                raise ValueError("活跃面片剩余间隙超过接受容差")
        elif force > result["equilibrium"].get("diagnostics", {}).get("force_acceptance_n", 1e-8):
            raise ValueError("非活跃标记遗漏了具有正法向力的承载面片")
        normal_force += force
        first_torque += force*first
        second_torque += force*second
    if len(active) != result["active_contact_count"] or not active:
        raise ValueError("活跃接触计数不一致或没有承载面片")
    if set(active) != set(result["equilibrium"]["active_indices"]):
        raise ValueError("接触活跃标记与平衡活跃集不一致")
    _close(normal_force, result["total_normal_force_n"], "法向合力N")
    _close(first_torque, result["pinion_torque_reconstructed_nm"], "主动轴重构扭矩Nm")
    _close(second_torque, result["wheel_torque_reconstructed_nm"], "从动轴重构扭矩Nm")
    _close(second_torque, target, "从动KKT扭矩平衡", absolute=1e-8)
    _close(max(c["mean_pressure_pa"] for c in contacts), result["maximum_cell_mean_pressure_pa"], "均压最大值Pa")
    # 主动轴没有第二个独立KKT约束，残差保留，不假设它为零。
    residual = (first_torque-torque)/abs(torque)
    _close(residual, result["pinion_torque_balance_relative_residual"], "主动重构相对残差", absolute=1e-10)
    return active


def _geometry_signature(result):
    fields = ("first_node_indices", "second_node_indices", "first_weights", "second_weights",
              "area_m2", "signed_gap_m", "pinion_lever_m", "wheel_lever_m", "tooth_offset")
    vectors = ("first_force_normal_local", "second_force_normal_local", "normal_world",
               "first_point_world_m", "second_point_world_m")
    payload = {"contacts": [{**{k: c[k] for k in fields}, **{k: c.get(k) for k in vectors}}
                             for c in result["contacts"]],
        "preparation": result["preparation"], "driving_angle_rad": result["driving_angle_rad"],
        "condensation": result.get("condensation"),
        "unloaded_tca_te_rad": result["unloaded_tca_te_rad"],
        "mesh_unloaded_correction_rad": result["mesh_unloaded_correction_rad"]}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _csv(path, rows):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return str(Path(path).resolve())


def _write(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return str(Path(path).resolve())


def export_generated_ltca_source_data(result, nodes_m, output_dir, *, metadata):
    """导出全部候选与实际三顶点；压力用三维面积，绘图投影固定为局部z-y。"""
    _metadata(metadata)
    active = _result(result)
    if len(nodes_m) != result["preparation"]["members"]["pinion"]["nodes"]:
        raise ValueError("小轮节点数与该结果准备网格不一致")
    rows, vertices = [], []
    for index, contact in enumerate(result["contacts"]):
        indices = contact["first_node_indices"]
        if len(indices) != 3 or len(set(indices)) != 3 or any(
                isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(nodes_m) for i in indices):
            raise ValueError("三角片须引用实际小轮网格中三个不同的零基节点")
        points = [tuple(_finite(v, "FE节点坐标m") for v in nodes_m[i]) for i in indices]
        if any(len(p) != 3 for p in points):
            raise ValueError("FE节点须为三维m坐标")
        a, b = ([points[j][k]-points[0][k] for k in range(3)] for j in (1, 2))
        cross = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        area = .5*math.sqrt(sum(v*v for v in cross))
        if area <= 0:
            raise ValueError("实际FE三角片退化，不能绘制为承载面积")
        _close(area, contact["area_m2"], "三维节点三角片面积", relative=1e-8, absolute=1e-18)
        row = {"contact_index": index, "tooth_offset": contact["tooth_offset"], "active": contact["active"],
            "force_n": contact["force_n"], "area_m2": contact["area_m2"], "mean_pressure_pa": contact["mean_pressure_pa"],
            "mean_pressure_mpa": contact["mean_pressure_pa"]*1e-6, "signed_gap_m": contact["signed_gap_m"],
            "remaining_gap_m": contact["remaining_gap_m"], "pinion_lever_m": contact["pinion_lever_m"],
            "wheel_lever_m": contact["wheel_lever_m"]}
        for prefix in ("first", "second"):
            for j in range(3):
                row[f"{prefix}_node_{j}"] = contact[f"{prefix}_node_indices"][j]
                row[f"{prefix}_weight_{j}"] = contact[f"{prefix}_weights"][j]
        rows.append(row)
        for corner, (node, point) in enumerate(zip(indices, points)):
            vertices.append({"contact_index": index, "corner": corner, "node_index": node,
                "local_x_m": point[0], "local_y_m": point[1], "local_z_m": point[2],
                "plot_z_mm": point[2]*1000, "plot_y_mm": point[1]*1000})
    summary_keys = ("driving_angle_rad", "pinion_torque_nm", "wheel_target_torque_nm",
        "wheel_torque_reconstructed_nm", "pinion_torque_reconstructed_nm", "pinion_torque_balance_relative_residual",
        "unloaded_tca_te_rad", "mesh_unloaded_correction_rad", "elastic_follower_lag_rad", "total_follower_te_rad",
        "incremental_mesh_stiffness_nm_per_rad", "candidate_count", "active_contact_count", "total_normal_force_n",
        "maximum_cell_mean_pressure_pa")
    summary = {k: result[k] for k in summary_keys}
    summary.update(converged=True, pressure_peak_validated=result.get("pressure_peak_validated", False),
        discretization_converged=result.get("discretization_converged"),
        wheel_relative_balance_residual=result["equilibrium"]["relative_balance_residual"],
        complementarity_residual_nm=result["equilibrium"]["complementarity_residual_nm"])
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    files = [_csv(output/"actual_generated_ltca_contacts.csv", rows),
             _csv(output/"actual_generated_ltca_facet_vertices.csv", vertices),
             _csv(output/"actual_generated_ltca_summary.csv", [summary])]
    data = {"metadata": metadata, "result_summary": summary, "geometry_signature": _geometry_signature(result),
        "active_contact_indices": active, "coordinate_definition": "实际未变形小轮FE制造局部z-y投影；原xyz为m，显示投影为mm；节点为零基",
        "pressure_definition": "法向力f除以实际三维三角片面积A；二维投影面积不参与计算，非Hertz峰值",
        "external_reference_comparison": False, "fields_source": "同一收敛解的contacts、equilibrium和原FE节点",
        "scope": "实际展成整轮线弹性FE接触；固定孔壁、根底连接及有限齿对域，未验证网格收敛"}
    files.append(_write(output/"actual_generated_ltca_source.json", data))
    return {"files": files, "summary": summary, "active_contact_indices": active}


def plot_generated_ltca(result, nodes_m, output_dir, *, metadata, publication=True):
    """只画真实活跃三角片均压及齿对法向合力，保留主动轴重构残差。"""
    if not isinstance(publication, bool): raise TypeError("publication须为布尔值")
    exported = export_generated_ltca_source_data(result, nodes_m, output_dir, metadata=metadata)
    from .plotting import _pyplot, _panel, _size
    from .benchmark_plotting import _save_figure
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    import numpy as np
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=_size(183, 88), gridspec_kw={"width_ratios": [2.1, 1]})
    fig.subplots_adjust(left=.11, right=.98, bottom=.24, top=.88, wspace=.65)
    active = [result["contacts"][i] for i in exported["active_contact_indices"]]
    polygons = [[(float(nodes_m[i][2])*1000, float(nodes_m[i][1])*1000) for i in c["first_node_indices"]] for c in active]
    pressures = [c["mean_pressure_pa"]*1e-6 for c in active]
    from .halfspace_plotting import pressure_colour_scale
    normalization,pressure_cmap,_=pressure_colour_scale([value*1e6 for value in pressures])
    collection = PolyCollection(polygons, array=np.array(pressures), cmap=pressure_cmap,
        norm=normalization, edgecolors="#394953", linewidths=.45)
    axes[0].add_collection(collection); axes[0].autoscale_view(); axes[0].margins(.12)
    axes[0].set(xlabel="小轮局部 z / mm", ylabel="小轮局部 y / mm")
    axes[0].set_aspect("equal", adjustable="datalim")
    colorbar = fig.colorbar(collection, ax=axes[0], fraction=.06, pad=.04)
    colorbar.set_label("面片平均压力 / MPa", fontsize=7)
    _panel(axes[0], "a", "实际活跃FE三角片：z-y投影")
    shares = {}
    for contact in result["contacts"]:
        key = contact["tooth_offset"]; shares[key] = shares.get(key, 0.)+contact["force_n"]
    offsets = sorted(shares)
    axes[1].bar([str(k) for k in offsets], [shares[k] for k in offsets], color="#2166AC", width=.6)
    axes[1].set_ylim(0, max(shares.values())*1.8)
    axes[1].set(xlabel="齿对偏移", ylabel="法向合力 / N")
    axes[1].grid(axis="y", linewidth=.4, alpha=.5); axes[1].set_axisbelow(True)
    _panel(axes[1], "b", "齿对合力")
    axes[1].text(.02, .98, f"从动相对平衡残差\n{result['equilibrium']['relative_balance_residual']:.2e}\n"
        f"互补残差\n{result['equilibrium']['complementarity_residual_nm']:.2e} N·m",
        transform=axes[1].transAxes, va="top", fontsize=7)
    fig.text(.5, .025, f"T1={result['pinion_torque_nm']:g} N·m；{len(active)}/{result['candidate_count']}个活跃；"
        f"主动重构相对残差{result['pinion_torque_balance_relative_residual']*100:+.4f}%\n"
        "未变形真实FE面片上的f/A，无平滑；未验证压力网格收敛或Hertz峰值，非外部对标",
        ha="center", va="bottom", fontsize=7, color="#52616B")
    files = _save_figure(fig, output_dir, "actual_generated_ltca", publication, "实际展成整轮FE面片平均压力与齿对合力")
    plt.close(fig)
    return {**exported, "figures": files}


def export_generated_ltca_load_sweep(results, output_dir, *, metadata):
    """导出同准备网格/相位的原signed载荷扫描，原顺序和三个角分量分别保留。"""
    _metadata(metadata)
    results = tuple(results)
    if len(results) < 2: raise ValueError("载荷扫描至少需要两个真实载荷解")
    for result in results: _result(result)
    signatures = {_geometry_signature(r) for r in results}
    torques = [r["pinion_torque_nm"] for r in results]
    if len(signatures) != 1 or len(set(torques)) != len(torques) or len({1 if t>0 else -1 for t in torques}) != 1:
        raise ValueError("载荷扫描要求同准备网格/相位/材料/接触几何、同一方向且不重复的signed扭矩")
    rows = []
    keys = ("driving_angle_rad", "pinion_torque_nm", "wheel_target_torque_nm", "wheel_torque_reconstructed_nm",
        "pinion_torque_reconstructed_nm", "pinion_torque_balance_relative_residual", "unloaded_tca_te_rad",
        "mesh_unloaded_correction_rad", "elastic_follower_lag_rad", "total_follower_te_rad",
        "incremental_mesh_stiffness_nm_per_rad", "total_normal_force_n", "candidate_count", "active_contact_count",
        "maximum_cell_mean_pressure_pa")
    for index, result in enumerate(results):
        rows.append({"source_index": index, **{k: result[k] for k in keys},
            "elastic_follower_lag_arcsec": math.degrees(result["elastic_follower_lag_rad"])*3600,
            "wheel_axis_stiffness_knm_per_rad": result["incremental_mesh_stiffness_nm_per_rad"]*.001,
            "wheel_relative_balance_residual": result["equilibrium"]["relative_balance_residual"],
            "complementarity_residual_nm": result["equilibrium"]["complementarity_residual_nm"],
            "pressure_peak_validated": result.get("pressure_peak_validated", False)})
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    files = [_csv(output/"actual_generated_ltca_load_sweep.csv", rows)]
    files.append(_write(output/"actual_generated_ltca_load_sweep_source.json", {"metadata": metadata,
        "geometry_signature": next(iter(signatures)), "source_count": len(rows),
        "plot_order_source_indices": sorted(range(len(rows)), key=lambda i: rows[i]["pinion_torque_nm"]),
        "definitions": {"angle": "从动signed弹性角，rad换角秒；不叠加TCA或网格初相位",
            "stiffness": "从动轴实际活跃集ell2.T Cactive^-1 ell2，Nm/rad换kNm/rad，非T1/theta割线"},
        "maximum_absolute_pinion_torque_relative_residual": max(abs(r["pinion_torque_balance_relative_residual"]) for r in results),
        "external_reference_comparison": False, "interpolation": "仅原载荷点折线，无补点或网格收敛外推"}))
    return {"files": files, "rows": rows}


def plot_generated_ltca_load_sweep(results, output_dir, *, metadata, publication=True):
    """画signed主动扭矩、signed从动弹性角与正的从动轴切线刚度。"""
    if not isinstance(publication, bool): raise TypeError("publication须为布尔值")
    results = tuple(results)
    exported = export_generated_ltca_load_sweep(results, output_dir, metadata=metadata)
    from .plotting import _pyplot, _panel, _size
    from .benchmark_plotting import _save_figure
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=_size(183, 80))
    fig.subplots_adjust(left=.12, right=.98, bottom=.25, top=.88, wspace=.45)
    rows = sorted(exported["rows"], key=lambda r: r["pinion_torque_nm"])
    x = [r["pinion_torque_nm"] for r in rows]
    for i, key in enumerate(("elastic_follower_lag_arcsec", "wheel_axis_stiffness_knm_per_rad")):
        axes[i].plot(x, [r[key] for r in rows], "o-", color=("#2166AC", "#B2182B")[i], markersize=3.5)
        axes[i].set_xticks(x); axes[i].set_xlabel("主动有符号扭矩 T1 / (N·m)")
        axes[i].grid(linewidth=.4, alpha=.5)
    axes[0].set_ylabel("从动有符号弹性滞后 / arcsec")
    axes[1].set_ylabel("从动轴切线刚度 / (kN·m/rad)")
    _panel(axes[0], "a", "有符号弹性响应")
    _panel(axes[1], "b", "实际活跃集增量刚度")
    residual = max(abs(r["pinion_torque_balance_relative_residual"]) for r in results)*100
    fig.text(.5, .025, f"{len(rows)}个真实载荷点；同整轮FE/工作面/相位；主动重构最大相对残差{residual:.3f}%\n"
        "TCA、网格初相位与弹性角分列源数据；未验证网格收敛，无外部LTCA对标",
        ha="center", va="bottom", fontsize=7, color="#52616B")
    files = _save_figure(fig, output_dir, "actual_generated_ltca_load_sweep", publication, "实际展成整轮FE有符号载荷扫描")
    plt.close(fig)
    return {**exported, "figures": files}
