"""同一物理内核上的参数扫描、多工况比较及 Pareto 方案筛选。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping, Sequence

import numpy as np

from .analysis import AnalysisResult, analyze_pair
from .surfaces import SurfaceModification


_PRESSURE_DEFINITIONS = {
    "center": "全部求解相位中条带中心的局部Hertz峰压最大值",
    "endpoint": "全部求解相位中接触线端点的局部Hertz峰压最大值",
}


def _validate_pressure_metric(pressure_metric):
    """只允许明确的中心或端点定义，避免未知选项落回默认压力。"""
    if not isinstance(pressure_metric, str) or pressure_metric not in _PRESSURE_DEFINITIONS:
        raise ValueError("pressure_metric必须为center或endpoint")
    return pressure_metric


def _selected_pressure_pa(result, pressure_metric):
    """读取所选压力；端点必须保留全部相位的有限原数组，不能用摘要补齐。"""
    if pressure_metric == "center":
        pressure = float(result.summary["max_pressure_pa"])
        if not np.isfinite(pressure) or pressure < 0:
            raise ValueError("中心最大压力必须为非负有限值，单位Pa")
        return pressure
    phases = np.asarray(result.phase)
    if phases.ndim != 1 or not phases.size or len(result.rotations) != phases.size:
        raise ValueError("endpoint压力指标需要全部求解相位的端点压力数组")
    maxima = []
    for index, rotation in enumerate(result.rotations):
        endpoint = rotation.endpoint_pressure
        if not isinstance(endpoint, Mapping) or "pressure_pa" not in endpoint:
            raise ValueError(f"endpoint压力指标缺少相位{index}的端点压力数组")
        try:
            pressure = np.asarray(endpoint["pressure_pa"], dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(f"相位{index}的端点压力必须为非负有限一维数组，单位Pa") from error
        if (pressure.ndim != 1 or not pressure.size or np.any(~np.isfinite(pressure))
                or np.any(pressure < 0)):
            raise ValueError(f"相位{index}的端点压力必须为非负有限非空一维数组，单位Pa")
        maxima.append(float(np.max(pressure)))
    return max(maxima)


@dataclass(frozen=True)
class OperatingPoint:
    """单个工况；失准为重叠齿宽右端减左端的法向间隙m。"""
    name: str
    torque_n_m: float
    misalignment_m: float = 0.0

    def __post_init__(self):
        if not self.name or not isinstance(self.name, str):
            raise ValueError("工况需要非空名称")
        if not np.isfinite(self.torque_n_m) or self.torque_n_m < 0 or not np.isfinite(self.misalignment_m):
            raise ValueError("工况扭矩幅值必须非负有限，失准必须有限")


@dataclass
class DesignStudy:
    """方案→工况结果；行表、最坏压力和Pareto统一采用pressure_metric。"""
    results: dict[str, dict[str, AnalysisResult]]
    rows: list[dict]
    design_summary: list[dict]
    pareto_names: list[str]
    pressure_metric: str = "center"

    def to_dict(self):
        """返回可写JSON的比较摘要，完整离散数组由export_result另存。"""
        metric = _validate_pressure_metric(self.pressure_metric)
        return {"rows": self.rows, "design_summary": self.design_summary,
                "pareto_names": self.pareto_names,
                "pressure_metric": metric,
                "pressure_definition": _PRESSURE_DEFINITIONS[metric],
                "units": {"max_pressure_mpa": "MPa", "worst_pressure_mpa": "MPa",
                          "pressure_limit_pa": "Pa", "te_peak_to_peak_um": "μm",
                          "worst_te_peak_to_peak_um": "μm"},
                "objectives": ["最坏工况法向TE峰峰值",
                               "最坏工况的" + _PRESSURE_DEFINITIONS[metric]],
                "search_scope": "在明确给定的有限候选与工况集合中比较，不宣称连续全局最优"}


def _pareto_indices(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return []
    return [i for i, value in enumerate(values)
            if not np.any(np.all(values <= value, axis=1) & np.any(values < value, axis=1))]


def compare_designs(pair, designs: Mapping[str, object], *,
                    operating_points: Sequence[OperatingPoint] | None = None,
                    pressure_limit_pa: float | None = None,
                    pressure_metric: str = "center", **analysis_kwargs) -> DesignStudy:
    """所有方案在同一工况集合比较；designs定义小轮受载面的修形。

    压力上限由使用者根据实际材料、寿命与可靠度要求给定，此处不猜测
    标准许用应力。默认只按TE和压力两个目标做非支配筛选，不暗含权重。
    pressure_metric为center时保留历史中心峰压；endpoint采用每个求解
    相位的端点原数组最大值，缺失或非有限即报错，不回退到中心值。
    """
    pressure_metric = _validate_pressure_metric(pressure_metric)
    if not designs or any(not isinstance(name, str) or not name for name in designs):
        raise ValueError("需要至少一个具名修形方案")
    if "modification_pinion" in analysis_kwargs:
        raise ValueError("小轮修形由designs指定，不能在分析参数中重复给定")
    if pressure_limit_pa is not None and (not np.isfinite(pressure_limit_pa) or pressure_limit_pa <= 0):
        raise ValueError("压力约束应为正有限值，单位Pa")
    kwargs = dict(analysis_kwargs)
    if operating_points is None:
        if "torque_n_m" not in kwargs:
            raise ValueError("未给工况集合时需要torque_n_m")
        operating_points = [OperatingPoint("给定工况", kwargs.pop("torque_n_m"), kwargs.pop("misalignment_m", 0.0))]
    elif "torque_n_m" in kwargs or "misalignment_m" in kwargs:
        raise ValueError("工况集合已经包含载荷和失准，不能重复给定")
    operating_points = tuple(operating_points)
    if not operating_points or len({p.name for p in operating_points}) != len(operating_points):
        raise ValueError("工况集合应非空且名称不能重复")
    results, rows, summaries = {}, [], []
    for name, modification in designs.items():
        if not hasattr(modification, "evaluate"):
            raise ValueError("修形方案必须提供evaluate(profile_factor,lead_factor)")
        results[name] = {}
        for point in operating_points:
            result = analyze_pair(pair, torque_n_m=point.torque_n_m, misalignment_m=point.misalignment_m,
                                  modification_pinion=modification, **kwargs)
            results[name][point.name] = result
            s = result.summary
            pressure = _selected_pressure_pa(result, pressure_metric)
            rows.append({"design": name, "operating_point": point.name,
                         "torque_n_m": point.torque_n_m, "misalignment_m": point.misalignment_m,
                         "te_peak_to_peak_um": s["te_peak_to_peak_m"]*1e6,
                         "max_pressure_mpa": pressure*1e-6,
                         "pressure_metric": pressure_metric,
                         "mean_stiffness_n_per_m": s["mean_mesh_stiffness_n_per_m"],
                         "load_centroid_offset_um": s["max_abs_load_centroid_offset_m"]*1e6,
                         "contact_fraction": s["mean_contact_fraction"],
                         "force_residual_n": s["max_force_residual_n"],
                         "gap_residual_m": s["max_gap_residual_m"]})
        selected = [r for r in rows if r["design"] == name]
        worst_pressure = max(r["max_pressure_mpa"] for r in selected)
        feasible = pressure_limit_pa is None or worst_pressure*1e6 <= pressure_limit_pa
        summaries.append({"design": name,
                          "worst_te_peak_to_peak_um": max(r["te_peak_to_peak_um"] for r in selected),
                          "worst_pressure_mpa": worst_pressure,
                          "pressure_metric": pressure_metric,
                          "worst_load_centroid_offset_um": max(r["load_centroid_offset_um"] for r in selected),
                          "feasible": bool(feasible), "pareto": False})
    feasible_indices = [i for i, s in enumerate(summaries) if s["feasible"]]
    objectives = [[summaries[i]["worst_te_peak_to_peak_um"], summaries[i]["worst_pressure_mpa"]]
                  for i in feasible_indices]
    for local_index in _pareto_indices(objectives):
        summaries[feasible_indices[local_index]]["pareto"] = True
    pareto = [s["design"] for s in summaries if s["pareto"]]
    return DesignStudy(results, rows, summaries, pareto, pressure_metric=pressure_metric)


def robust_design_search(pair, candidates, operating_points, *, pressure_limit_pa=None,
                         pressure_metric="center", **analysis_kwargs):
    """以给定多工况的最坏TE和压力评价有限候选集，保留完整可复核结果。"""
    return compare_designs(pair, candidates, operating_points=operating_points,
                           pressure_limit_pa=pressure_limit_pa, pressure_metric=pressure_metric,
                           **analysis_kwargs)


def sweep_parameter(pair, parameter: str, values, *, base_modification=SurfaceModification(),
                    operating_points=None, pressure_metric="center", **analysis_kwargs):
    """扫描SurfaceModification中的一个具名物理参数，单位与字段一致。"""
    pressure_metric = _validate_pressure_metric(pressure_metric)
    if parameter not in base_modification.__dataclass_fields__:
        raise ValueError(f"未知修形参数：{parameter}")
    values = tuple(values)
    if not values or not all(np.isfinite(v) for v in values) or len(set(values)) != len(values):
        raise ValueError("扫描数值应为不重复的非空有限序列")
    designs = {f"{parameter}={value:.12g}": replace(base_modification, **{parameter: value}) for value in values}
    return compare_designs(pair, designs, operating_points=operating_points,
                           pressure_metric=pressure_metric, **analysis_kwargs)
