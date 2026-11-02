"""显式双侧数据的参数化对标绘图与源数据导出。

本模块为新写数据与绘图实现，不计算几何、强度或接触结果。
调用者负责确认共同定义与工况；未知结果不得以常数或自身结果代替参照。
原始SI量、显示单位变换和误差均保留，零参照的相对误差不定义。
"""
from __future__ import annotations

import copy
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Mapping, Sequence
from .publication_style import BLUE, RED


_DISPLAY = {"m": ("mm", 1000.0), "rad": ("deg", 180.0 / math.pi),
            "Pa": ("MPa", 1e-6), "N": ("N", 1.0), "N*m": ("N·m", 1.0),
            "1": ("1", 1.0), "N*m/rad": ("N·m/rad", 1.0),
            "sqrt(Pa)": ("√Pa", 1.0), "kg/m": ("kg/m", 1.0),
            "rad/s": ("rad/s", 1.0), "N/m2": ("N/m²", 1.0)}


@dataclass(frozen=True)
class BenchmarkMetric:
    """一个同定义标量；reference/project与容差均采用unit指定的原单位。"""

    field: str
    label: str
    unit: str
    reference: float
    project: float
    group: str
    absolute_tolerance: float | None = None
    relative_tolerance: float | None = None
    note: str = ""


@dataclass(frozen=True)
class BenchmarkCase:
    """一个已配准工况；method和coverage必须明确比较方法与覆盖边界。"""

    case_id: str
    label: str
    method: str
    coverage: str
    metrics: tuple[BenchmarkMetric, ...]
    reference_name: str = "MASTA"
    project_name: str = "本项目"
    registration: dict | None = None


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name}须为有限实数")
    return float(value)


def _validated(cases):
    values = tuple(cases)
    if not values or any(not isinstance(c, BenchmarkCase) for c in values):
        raise ValueError("至少需要一个BenchmarkCase")
    ids = set()
    names = {(c.reference_name, c.project_name) for c in values}
    if len(names) != 1:
        raise ValueError("同一比较导出须保持双方名称一致")
    for case in values:
        if not all(isinstance(v, str) and v.strip() for v in
                   (case.case_id, case.label, case.method, case.coverage,
                    case.reference_name, case.project_name)):
            raise ValueError("工况名称、方法、覆盖边界和双方名称不得为空")
        if case.case_id in ids or not case.metrics:
            raise ValueError("工况标识须唯一且至少含一个共同字段")
        ids.add(case.case_id)
        fields = set()
        for metric in case.metrics:
            if not isinstance(metric, BenchmarkMetric):
                raise TypeError("metrics须含BenchmarkMetric")
            if metric.unit not in _DISPLAY:
                raise ValueError(f"未支持的原单位：{metric.unit}")
            if not metric.field or not metric.label or not metric.group or metric.field in fields:
                raise ValueError("字段标识/标签/分组须非空且同工况字段唯一")
            fields.add(metric.field)
            _finite(metric.reference, "参照值")
            _finite(metric.project, "本项目值")
            for tolerance in (metric.absolute_tolerance, metric.relative_tolerance):
                if tolerance is not None and _finite(tolerance, "容差") < 0:
                    raise ValueError("容差不得为负")
    return values


def benchmark_rows(cases: Sequence[BenchmarkCase]) -> list[dict]:
    """生成独立源数据行，不写文件；相对差为(project-reference)/abs(reference)。"""
    rows = []
    for case in _validated(cases):
        for metric in case.metrics:
            display_unit, scale = _DISPLAY[metric.unit]
            difference = metric.project - metric.reference
            relative = difference / abs(metric.reference) if metric.reference != 0 else None
            tolerances = [v for v in (metric.absolute_tolerance,
                          None if metric.relative_tolerance is None else
                          metric.relative_tolerance * abs(metric.reference)) if v is not None]
            tolerance = max(tolerances) if tolerances else None
            rows.append({
                "case_id": case.case_id, "case_label": case.label, "method": case.method,
                "coverage": case.coverage, "reference_name": case.reference_name,
                "project_name": case.project_name, "field": metric.field,
                "label": metric.label, "group": metric.group, "unit": metric.unit,
                "reference": metric.reference, "project": metric.project,
                "signed_difference": difference, "absolute_difference": abs(difference),
                "relative_difference": relative, "absolute_tolerance": metric.absolute_tolerance,
                "relative_tolerance": metric.relative_tolerance,
                "effective_absolute_tolerance": tolerance,
                "matched": None if tolerance is None else abs(difference) <= tolerance,
                "display_unit": display_unit, "display_scale": scale,
                "reference_display": metric.reference * scale,
                "project_display": metric.project * scale,
                "signed_difference_display": difference * scale, "note": metric.note,
            })
    return rows


def export_benchmark_source_data(cases, output_dir, *, stem="benchmark") -> dict:
    """导出完整比较CSV与可重建输入JSON；不隐去失败或未设置容差的字段。"""
    cases = _validated(cases)
    rows = benchmark_rows(cases)
    if not isinstance(stem, str) or not stem or Path(stem).name != stem:
        raise ValueError("stem须为无目录的文件名前缀")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    csv_path, input_path = output / f"{stem}.csv", output / f"{stem}_inputs.json"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    data = {"schema": "bevel_benchmark_v1", "cases": [asdict(c) for c in cases],
            "scope": "显式双侧标量比较；覆盖范围由各工况给定，非完整标准认证",
            "relative_difference_definition": "(project-reference)/abs(reference)，零参照留空",
            "row_count": len(rows), "matched_count": sum(r["matched"] is True for r in rows),
            "unassessed_count": sum(r["matched"] is None for r in rows)}
    input_path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8")
    return {"csv": str(csv_path.resolve()), "inputs": str(input_path.resolve()),
            "row_count": len(rows)}


def _save_figure(fig, output, stem, publication, description):
    """保存固定物理画布与实际数据说明，不作tight裁剪或附加合成工况描述。"""
    files = []
    for extension in ("png", "svg", *(("pdf",) if publication else ())):
        path = Path(output) / f"{stem}.{extension}"
        metadata = ({"Title": description, "Description": description, "Date": None}
                    if extension == "svg" else
                    {"Title": description, "Subject": description, "CreationDate": None}
                    if extension == "pdf" else {"Description": description})
        fig.savefig(path, dpi=600, facecolor="white", metadata=metadata)
        files.append(str(path.resolve()))
    return files


def plot_benchmark_comparison(cases, output_dir, *, panels=None, stem="benchmark",
                              error_panel=False, publication=True, footer="") -> dict:
    """按显式(case_id, group, 标题)分面绘双侧值，可追加同字段相对误差面板。

    每个面板只能含一种原单位；长度、角度、应力分别显示mm、deg、MPa。
    error_panel仅接受一个数值面板，避免不同定义或量纲默默混合。
    """
    cases = _validated(cases)
    rows = benchmark_rows(cases)
    if not isinstance(error_panel, bool) or not isinstance(publication, bool):
        raise TypeError("error_panel与publication须为布尔值")
    if panels is None:
        panels = [(c.case_id, group, c.label) for c in cases
                  for group in dict.fromkeys(m.group for m in c.metrics)]
    panels = tuple(panels)
    if not panels or (error_panel and len(panels) != 1):
        raise ValueError("误差面板要求恰好一个同单位数值面板")
    selected = []
    for case_id, group, heading in panels:
        subset = [r for r in rows if r["case_id"] == case_id and r["group"] == group]
        if not subset or len({r["unit"] for r in subset}) != 1:
            raise ValueError("每个面板须有共同字段且原单位一致")
        selected.append((heading, subset))
    # 全部数据与分面检查通过后才写文件或加载绘图后端。
    exported = export_benchmark_source_data(cases, output_dir, stem=stem)
    from .plotting import _panel, _pyplot, _size
    plt = _pyplot()
    from matplotlib.lines import Line2D
    count = len(selected) + int(error_panel)
    columns = min(2, count)
    nrows = math.ceil(count / columns)
    max_fields = max(len(v) for _, v in selected)
    height = max(79.0, nrows * (24.0 + max_fields * 4.0) + 17.0)
    fig, axes = plt.subplots(nrows, columns, figsize=_size(183, height), squeeze=False)
    fig.subplots_adjust(left=.16, right=.98, bottom=.23 if nrows == 1 else .16,
                        top=.86 if nrows == 1 else .90, hspace=.48, wspace=.53)
    blue, orange = BLUE, RED
    for i, (heading, subset) in enumerate(selected):
        ax = axes.flat[i]
        y = list(range(len(subset)))
        ax.plot([r["reference_display"] for r in subset], [v-.11 for v in y],
                "o", color=blue, markersize=3.4, label=subset[0]["reference_name"])
        ax.plot([r["project_display"] for r in subset], [v+.11 for v in y],
                "s", markerfacecolor="none", markeredgecolor=orange,
                markersize=3.4, label=subset[0]["project_name"])
        ax.set_yticks(y, [r["label"] for r in subset])
        ax.invert_yaxis()
        ax.set_xlabel(f"值 / {subset[0]['display_unit']}")
        ax.grid(axis="x", linewidth=.4, alpha=.5)
        ax.set_axisbelow(True)
        ax.margins(x=.12, y=.09)
        _panel(ax, chr(97+i), heading)
    if error_panel:
        ax = axes.flat[1]
        subset = selected[0][1]
        valid = [(i, r["relative_difference"]) for i, r in enumerate(subset)
                 if r["relative_difference"] is not None]
        # 误差直接取双侧输入，不放大原值差；坐标轴科学计数明确标示误差尺度。
        ax.axvline(0, color="#7A858C", linewidth=.6)
        for i, value in valid:
            ax.plot([0, value], [i, i], color=orange, linewidth=.8)
        ax.plot([v for _, v in valid], [i for i, _ in valid], "o", color=orange, markersize=3.4)
        ax.set_yticks(range(len(subset)), [r["label"] for r in subset])
        ax.invert_yaxis()
        ax.set_xlabel(f"相对差 ({cases[0].project_name}−{cases[0].reference_name}) / |{cases[0].reference_name}|")
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useMathText=False)
        ax.grid(axis="x", linewidth=.4, alpha=.5)
        ax.margins(x=.2, y=.09)
        _panel(ax, "b", "同字段相对差")
    for ax in list(axes.flat)[count:]:
        ax.set_visible(False)
    handles = [Line2D([], [], color=blue, marker="o", linestyle="none", markersize=3.4),
               Line2D([], [], color=orange, marker="s", markerfacecolor="none",
                      linestyle="none", markersize=3.4)]
    fig.legend(handles, [cases[0].reference_name, cases[0].project_name],
               loc="upper center", ncol=2, bbox_to_anchor=(.57, .995))
    if footer:
        fig.text(.5, .025, footer, ha="center", va="bottom", fontsize=7, color="#52616B")
    files = _save_figure(fig, output_dir, stem, publication, "共同定义双侧数值对标；无统计误差线")
    plt.close(fig)
    return {**exported, "figures": files, "width_mm": 183.0, "height_mm": height,
            "panels": [list(p) for p in panels], "error_panel": error_panel}


def make_ltca_benchmark_case(*, case_id, label, quantity, reference_values,
                             project_values, sample_ids, unit, registration,
                             relative_tolerance=None, absolute_tolerance=None) -> BenchmarkCase:
    """验证双方承载数据与配准声明后建立比较工况，不生成或推断缺失数据。

    quantity为pressure、elastic_follower_lag或incremental_mesh_stiffness；
    registration须明确双方数据有效且同几何、载荷、坐标系、相位原点和评价定义。
    还需保存双方逐样本主动角rad、共同相位原点rad、坐标系和载荷说明。
    压力须保存共同接触采样三维坐标m；函数不代替调用者执行配准或插值。
    """
    units = {"pressure": "Pa", "elastic_follower_lag": "rad",
             "incremental_mesh_stiffness": "N*m/rad"}
    if quantity not in units or unit != units[quantity]:
        raise ValueError("LTCA量名或单位不匹配；不得自动换算未知定义")
    if not isinstance(registration, Mapping):
        raise ValueError("须显式提供LTCA双方配准信息")
    required = ["reference_data_validated", "project_data_validated", "same_geometry",
                "same_load", "same_coordinate_frame", "same_phase_origin", "same_quantity_definition"]
    if quantity == "pressure":
        required.append("same_contact_sampling")
    if any(registration.get(k) is not True for k in required):
        raise ValueError("双方结果或几何/载荷/坐标/相位/评价/采样配准尚未确认")
    if not isinstance(registration.get("description"), str) or not registration["description"].strip():
        raise ValueError("配准说明不得为空")
    if reference_values is None or project_values is None or sample_ids is None:
        raise ValueError("缺少显式双侧LTCA数组，拒绝建立定量对标")
    reference, project, ids = tuple(reference_values), tuple(project_values), tuple(sample_ids)
    if not reference or len(reference) != len(project) or len(reference) != len(ids):
        raise ValueError("双方数组与共同采样标识须非空且等长")
    if any(not isinstance(s, str) or not s for s in ids) or len(set(ids)) != len(ids):
        raise ValueError("共同采样标识须为非空且唯一的字符串")
    for key in ("coordinate_frame", "load_description"):
        if not isinstance(registration.get(key), str) or not registration[key].strip():
            raise ValueError("须明确共同坐标系及载荷说明")
    _finite(registration.get("phase_origin_rad"), "共同相位原点rad")
    angles = []
    for key in ("reference_driving_angles_rad", "project_driving_angles_rad"):
        values = registration.get(key)
        if values is None or len(values) != len(ids):
            raise ValueError("须逐样本提供双方主动角rad，长度与比较数组一致")
        angles.append(tuple(_finite(value, "主动角rad") for value in values))
    if angles[0] != angles[1]:
        raise ValueError("双方主动角尚未显式配准，函数不会自动取近邻或插值")
    if quantity == "pressure":
        points = registration.get("common_sample_points_m")
        if points is None or len(points) != len(ids) or any(len(p) != 3 for p in points):
            raise ValueError("压力比较须提供逐样本共同三维坐标m")
        for point in points:
            for value in point:
                _finite(value, "共同采样坐标m")
    result = BenchmarkCase(case_id, label, "已配准的LTCA双侧数据", registration["description"],
                           tuple(BenchmarkMetric(s, s, unit, a, b, quantity,
                           absolute_tolerance, relative_tolerance) for s, a, b in zip(ids, reference, project)),
                           registration=copy.deepcopy(dict(registration)))
    _validated((result,))
    return result


def export_mesh_study_source_data(assessment, output_dir, *, controls, stem="mesh_study") -> dict:
    """导出独立局部网格变化研究，保留未通过门槛状态，不提升为收敛认证。

    assessment为已完成的多级变化评估，controls是各级共同控制参数。
    本函数只核验绘图字段；生成评估的调用者仍需确认所有级的控制参数一致。
    """
    if not isinstance(assessment, Mapping) or not isinstance(controls, Mapping):
        raise ValueError("需要网格评估字典与共同控制参数")
    rows = assessment.get("rows", ())
    if len(rows) < 2 or assessment.get("checked_levels") != len(rows):
        raise ValueError("至少需要两级已评估网格数据")
    if not isinstance(assessment.get("passed"), bool):
        raise ValueError("评估须显式保留通过或失败状态")
    fields = ("contact_size_m", "driving_angle_rad", "pinion_torque_nm",
              "elastic_follower_lag_rad", "incremental_mesh_stiffness_nm_per_rad",
              "maximum_cell_mean_pressure_pa")
    exported_rows = []
    for index, row in enumerate(rows):
        values = {k: _finite(row[k], k) for k in fields}
        if any(values[k] <= 0 for k in ("contact_size_m", "pinion_torque_nm",
                   "incremental_mesh_stiffness_nm_per_rad", "maximum_cell_mean_pressure_pa")):
            raise ValueError("网格尺寸、扭矩、刚度及面片平均压力最大值须为正")
        if values["driving_angle_rad"] != controls.get("driving_angle_rad") or \
           values["pinion_torque_nm"] != controls.get("pinion_torque_nm"):
            raise ValueError("网格研究各级主动角和扭矩必须与共同控制一致")
        exported_rows.append({"level": index, **values,
            "contact_size_mm": values["contact_size_m"] * 1000,
            "elastic_follower_lag_arcsec": values["elastic_follower_lag_rad"] * 180/math.pi * 3600,
            "incremental_mesh_stiffness_knm_per_rad": values["incremental_mesh_stiffness_nm_per_rad"] / 1000,
            "maximum_cell_mean_pressure_gpa": values["maximum_cell_mean_pressure_pa"] / 1e9,
            "change_thresholds_passed": assessment["passed"],
            "pressure_peak_validated": assessment.get("pressure_peak_validated", False)})
    sizes = [r["contact_size_m"] for r in exported_rows]
    if len(set(sizes)) != len(sizes) or any(a <= b for a, b in zip(sizes, sizes[1:])):
        raise ValueError("网格级须按尺寸严格递减排列，不静默排序或删点")
    if not isinstance(stem, str) or not stem or Path(stem).name != stem:
        raise ValueError("stem须为无目录的文件名前缀")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    csv_path, json_path = output / f"{stem}.csv", output / f"{stem}.json"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(exported_rows[0]))
        writer.writeheader()
        writer.writerows(exported_rows)
    data = {"assessment": assessment, "common_controls": controls,
        "scope": "独立局部网格变化研究，不是外部软件对标或物理验证",
        "pressure_definition": "面片平均压力最大值，不是Hertz峰值",
        "plotting": "原级顺序折线，不平滑、不外推、不估计收敛阶"}
    json_path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return {"csv": str(csv_path.resolve()), "inputs": str(json_path.resolve()), "row_count": len(exported_rows)}


def plot_mesh_study(assessment, output_dir, *, controls, stem="mesh_study", publication=True) -> dict:
    """画弹性角、增量刚度和平均压力最大值的三面板原级变化；失败状态正常展示。"""
    if not isinstance(publication, bool):
        raise TypeError("publication须为布尔值")
    exported = export_mesh_study_source_data(assessment, output_dir, controls=controls, stem=stem)
    from .plotting import _panel, _pyplot, _size
    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=_size(183, 82))
    fig.subplots_adjust(left=.09, right=.98, bottom=.25, top=.85, wspace=.63)
    x = [r["contact_size_m"]*1000 for r in assessment["rows"]]
    definitions = [
        ("elastic_follower_lag_rad", 180/math.pi*3600, "弹性从动滞后", "弹性滞后 / arcsec"),
        ("incremental_mesh_stiffness_nm_per_rad", .001, "增量扭转刚度", "刚度 / (kN·m/rad)"),
        ("maximum_cell_mean_pressure_pa", 1e-9, "最大面片平均压力", "平均压力最大值 / GPa"),
    ]
    for index, (key, scale, title, ylabel) in enumerate(definitions):
        ax = axes[index]
        ax.plot(x, [r[key]*scale for r in assessment["rows"]], "o-", color="#2166AC", markersize=3.5)
        ax.set_xticks(x, [f"{v:g}" for v in x])
        ax.set_xlim(max(x)+.05*(max(x)-min(x)), min(x)-.05*(max(x)-min(x)))
        ax.set(xlabel="局部网格尺寸 / mm", ylabel=ylabel)
        ax.grid(linewidth=.4, alpha=.5)
        _panel(ax, chr(97+index), title)
    status = "变化门槛通过（仍非收敛认证）" if assessment["passed"] else "变化门槛未通过"
    fig.text(.5, .97, f"独立局部网格研究：{status}", ha="center", va="top", fontsize=8)
    fig.text(.5, .025, "相同CAD、材料、载荷与连续折线加密控制；仅改变局部尺寸\n"
             "面片平均压力最大值不是Hertz峰值；不作收敛阶外推或外部软件对标",
             ha="center", va="bottom", fontsize=7, color="#52616B")
    files = _save_figure(fig, output_dir, stem, publication, "独立局部网格变化研究；保留原评估状态")
    plt.close(fig)
    return {**exported, "figures": files, "width_mm": 183.0, "height_mm": 82.0,
            "assessment_passed": assessment["passed"], "external_reference_comparison": False}


def export_generated_tca_source_data(results, output_dir, *, metadata, stem="generated_tca") -> dict:
    """导出实际曲面局部相切的已保存结果与制造/初装参数，不运行TCA。

    局部点由world=R@local+t反解为R.T@(world−t)，所有长度原值为m。
    无载TE须等于从动增量角减名义增量角，不减去零主动角的非零值。
    """
    results = tuple(results)
    if len(results) < 2 or any(r.get("converged") is not True or r.get("status") != "tangent" for r in results):
        raise ValueError("至少需要两个已收敛局部相切结果；失败点不能进入成功路径图")
    if not isinstance(metadata, Mapping) or any(key not in metadata for key in
            ("requests", "initial_angles_rad", "shaft_angle_rad", "build_solids")):
        raise ValueError("须保留制造输入、初装角、轴交角与实体检查构建信息")
    rows = []
    for index, result in enumerate(results):
        frames = result["frames"]
        matrix = frames["pinion_rotation_rows"]
        if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
            raise ValueError("小轮旋转矩阵须为3×3按行布局")
        matrix = [[_finite(v, "旋转矩阵元素") for v in row] for row in matrix]
        if any(abs(sum(matrix[k][i]*matrix[k][j] for k in range(3))-(1 if i==j else 0)) > 1e-9
               for i in range(3) for j in range(3)):
            raise ValueError("小轮旋转矩阵须正交，不能将一般矩阵转置当逆矩阵")
        points = [result[k] for k in ("pinion_contact_point_world_m", "wheel_contact_point_world_m")]
        translation = frames["translation_m"]
        if any(len(point) != 3 for point in (*points, translation)):
            raise ValueError("接触点和安装平移须为三维m坐标")
        for point in (*points, translation):
            for v in point: _finite(v, "接触点/平移m")
        local = [sum(matrix[i][j]*(points[0][i]-translation[i]) for i in range(3)) for j in range(3)]
        angle, follower, nominal, te = (_finite(result[k], k) for k in
            ("driving_angle_rad", "follower_increment_angle_rad", "nominal_follower_increment_angle_rad", "unloaded_te_rad"))
        if not math.isclose(te, follower-nominal, rel_tol=1e-10, abs_tol=1e-14):
            raise ValueError("无载TE与从动增量角减名义增量角不一致")
        row = {"sample_index": index, "driving_angle_rad": angle, "driving_angle_deg": math.degrees(angle),
            "follower_increment_angle_rad": follower, "nominal_follower_increment_angle_rad": nominal,
            "unloaded_te_rad": te, "unloaded_te_arcsec": math.degrees(te)*3600,
            "position_residual_m": _finite(result["position_residual_m"], "位置残差m"),
            "normal_cross_residual": _finite(result["normal_cross_residual"], "法向叉积残差"),
            "opposing_normal_dot": _finite(result["opposing_normal_dot"], "法向点积"),
            "solid_interference_checked": result.get("solid_interference_checked", False)}
        for member, point in zip(("pinion", "wheel"), points):
            for axis, value in zip("xyz", point): row[f"{member}_world_{axis}_m"] = value
            uv = result[f"{member}_uv"]
            if len(uv) != 2 or any(not 0 < _finite(v, "工作面归一参数") < 1 for v in uv):
                raise ValueError("成功路径要求所选工作面的内部归一参数")
            row[f"{member}_u"], row[f"{member}_v"] = uv
            material = result.get("normal_material_checks", {}).get(member, {})
            row[f"{member}_material_normal_classified"] = material.get("classified")
            row[f"{member}_parameter_normal_to_outward"] = material.get("parameter_normal_to_outward")
        for axis, value in zip("xyz", local): row[f"pinion_local_{axis}_m"] = value
        row["pinion_local_axial_mm"] = local[2]*1000
        row["pinion_local_radial_mm"] = math.hypot(local[0], local[1])*1000
        rows.append(row)
    if any(a["driving_angle_rad"] >= b["driving_angle_rad"] for a,b in zip(rows,rows[1:])):
        raise ValueError("主动角须严格递增，函数不静默排序或删除样本")
    if not isinstance(stem, str) or not stem or Path(stem).name != stem:
        raise ValueError("stem须为无目录的文件名前缀")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    csv_path, json_path = output/f"{stem}.csv", output/f"{stem}.json"
    with csv_path.open("w",encoding="utf-8-sig",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=tuple(rows[0]));writer.writeheader();writer.writerows(rows)
    data = {"results": results, "metadata": metadata,
        "coordinate_definition": "world=R@local+t；小轮制造局部轴向z/径向hypot(x,y)的子午投影，m换mm",
        "scope": "所选实际工作面无载局部相切；非外部对标、非承载或全轮无侵入认证",
        "interpolation": "只连接原样本，不补点、不强制零TE或周期闭合"}
    json_path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    return {"csv": str(csv_path.resolve()), "inputs": str(json_path.resolve()), "row_count": len(rows)}


def plot_generated_tca(results, output_dir, *, metadata, stem="generated_tca", publication=True) -> dict:
    """画已保存实际工作面的局部接触路径与无载TE；不调用制造或接触求解器。"""
    if not isinstance(publication,bool): raise TypeError("publication须为布尔值")
    results = tuple(results)
    exported = export_generated_tca_source_data(results,output_dir,metadata=metadata,stem=stem)
    with Path(exported["csv"]).open(encoding="utf-8-sig",newline="") as stream:
        rows = list(csv.DictReader(stream))
    from .plotting import _panel,_pyplot,_size
    plt = _pyplot()
    fig,axes = plt.subplots(1,2,figsize=_size(183,86))
    fig.subplots_adjust(left=.12,right=.98,bottom=.24,top=.88,wspace=.48)
    x,y = ([float(r[k]) for r in rows] for k in ("pinion_local_axial_mm","pinion_local_radial_mm"))
    axes[0].plot(x,y,"o-",color="#2166AC",markersize=3.5)
    for i,row in enumerate(rows):
        axes[0].annotate(f"{float(row['driving_angle_deg']):g}°",(x[i],y[i]),
            xytext=(5,3),textcoords="offset points",fontsize=7)
    axes[0].set(xlabel="小轮局部轴向 z / mm",ylabel="小轮局部径向 ρ / mm")
    axes[0].set_aspect("equal",adjustable="datalim")
    axes[0].margins(.18)
    _panel(axes[0],"a","接触路径：小轮子午投影")
    axes[1].plot([float(r["driving_angle_deg"]) for r in rows],
        [float(r["unloaded_te_arcsec"]) for r in rows],"o-",color="#B2182B",markersize=3.5)
    axes[1].axhline(0,color="#82909A",linewidth=.5)
    axes[1].set(xlabel="主动增量角 / deg",ylabel="无载从动TE / arcsec")
    _panel(axes[1],"b","相切相位相对名义齿数比")
    for ax in axes: ax.grid(linewidth=.4,alpha=.5)
    material_checked = all(r.get("normal_material_checks",{}).get(m,{}).get("classified") is True
                           for r in results for m in ("pinion","wheel"))
    interference_checked = all(r.get("solid_interference_checked") is True for r in results)
    normal_note = "已核对局部材料外法向" if material_checked else "未核对局部材料外法向"
    interference_note = "原解已做实体干涉检查" if interference_checked else "尚未全轮干涉检查"
    fig.text(.5,.025,f"{len(rows)}个实算姿态；{normal_note}；{interference_note}\n"
        "实际NURBS工作面无载局部相切，非正间隙相位、承载或外部对标",ha="center",va="bottom",fontsize=7,color="#52616B")
    files = _save_figure(fig,output_dir,stem,publication,"实际工作面局部无载相切：接触路径与TE")
    plt.close(fig)
    return {**exported,"figures":files,"width_mm":183.0,"height_mm":86.0,
        "external_reference_comparison":False,"solid_interference_checked":interference_checked,
        "local_material_normals_checked":material_checked}
