"""新写的参数修形及测量曲面；法向正值表示材料去除，单位m。"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


def _factors(profile_factor, lead_factor):
    p, l = np.broadcast_arrays(np.asarray(profile_factor, dtype=float), np.asarray(lead_factor, dtype=float))
    if np.any(~np.isfinite(p)) or np.any(~np.isfinite(l)) or np.any((p < 0) | (p > 1)) or np.any((l < 0) | (l > 1)):
        raise ValueError("齿廓、齿向因子必须在[0,1]内且为有限数")
    return p, l


@dataclass(frozen=True)
class SurfaceModification:
    """在单轮完整齿廓/齿宽上定义法向去除量。

    profile/lead_crown_m 是相应边缘相对中心的去除量；slope_m 是
    两端去除差，符号决定方向，并以较低端去除0为基准。bias_m 是
    两组对角去除差，非零时以较低对角去除0为基准。顶/根/端部
    修形为二次过渡。参数去除量始终非负，不隐含两轮分摊或建议值。
    """
    profile_crown_m: float = 0.0
    lead_crown_m: float = 0.0
    profile_slope_m: float = 0.0
    lead_slope_m: float = 0.0
    bias_m: float = 0.0
    tip_relief_m: float = 0.0
    tip_start: float = 0.75
    root_relief_m: float = 0.0
    root_end: float = 0.25
    end_relief_m: float = 0.0
    end_fraction: float = 0.1

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            value = getattr(self, name)
            if not np.isfinite(value):
                raise ValueError(f"{name} 必须为有限数")
        for name in ("profile_crown_m", "lead_crown_m", "tip_relief_m", "root_relief_m", "end_relief_m"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} 必须为非负材料去除量，单位m")
        if not 0 <= self.tip_start < 1:
            raise ValueError("齿顶修形起点tip_start必须在[0,1)内")
        if not 0 < self.root_end <= 1:
            raise ValueError("齿根修形终点root_end必须在(0,1]内")
        if not 0 < self.end_fraction <= 0.5:
            raise ValueError("齿端修形长度比例end_fraction必须在(0,0.5]内")

    def evaluate(self, profile_factor, lead_factor):
        p, l = _factors(profile_factor, lead_factor)
        removal = self.profile_crown_m * (2 * p - 1) ** 2 + self.lead_crown_m * (2 * l - 1) ** 2
        removal += self.profile_slope_m * (p - 0.5) + abs(self.profile_slope_m) / 2
        removal += self.lead_slope_m * (l - 0.5) + abs(self.lead_slope_m) / 2
        removal += (self.bias_m * (2 * p - 1) * (2 * l - 1) + abs(self.bias_m)) / 2
        removal += self.tip_relief_m * (np.maximum(p - self.tip_start, 0) / (1 - self.tip_start)) ** 2
        removal += self.root_relief_m * (np.maximum(self.root_end - p, 0) / self.root_end) ** 2
        removal += self.end_relief_m * (
            (np.maximum(self.end_fraction - l, 0) / self.end_fraction) ** 2
            + (np.maximum(l - (1 - self.end_fraction), 0) / self.end_fraction) ** 2
        )
        return removal


def _coordinates(values, name: str) -> np.ndarray:
    data = np.array(values, dtype=float, copy=True)
    if data.ndim != 1 or data.size < 2 or np.any(~np.isfinite(data)) or np.any(np.diff(data) <= 0):
        raise ValueError(f"{name} 必须为至少两个严格递增有限坐标")
    if data[0] < 0 or data[-1] > 1:
        raise ValueError(f"{name} 必须使用单轮完整范围内的[0,1]归一化坐标")
    data.setflags(write=False)
    return data


def _check_bounds(bounds: str) -> None:
    if bounds not in ("error", "clip", "extrapolate"):
        raise ValueError("越界策略只能为error、clip或extrapolate")


def _prepare_query(query: np.ndarray, coordinates: np.ndarray, bounds: str) -> np.ndarray:
    if bounds == "error" and np.any((query < coordinates[0]) | (query > coordinates[-1])):
        raise ValueError("求值坐标超出测量范围；请显式选择clip或extrapolate策略")
    return np.clip(query, coordinates[0], coordinates[-1]) if bounds == "clip" else query


@dataclass(frozen=True)
class MeasuredCurve:
    """明确齿廓或齿向坐标的一维测量偏差，采用分段线性插值。

    removal_m 以调用方给定的理想面为零基准；正值去除、负值凸起。
    数据覆盖不足时默认报错，不自动外推。extrapolate只使用端部斜率。
    """
    coordinates: object
    removal_m: object
    axis: Literal["profile", "lead"] = "profile"
    bounds: Literal["error", "clip", "extrapolate"] = "error"

    def __post_init__(self) -> None:
        c = _coordinates(self.coordinates, "测量曲线坐标")
        values = np.array(self.removal_m, dtype=float, copy=True)
        if values.shape != c.shape or np.any(~np.isfinite(values)):
            raise ValueError("测量曲线偏差必须与坐标同长且有限，单位m")
        if self.axis not in ("profile", "lead"):
            raise ValueError("测量曲线轴只能为profile或lead")
        _check_bounds(self.bounds)
        values.setflags(write=False)
        object.__setattr__(self, "coordinates", c)
        object.__setattr__(self, "removal_m", values)

    def evaluate(self, profile_factor, lead_factor):
        p, l = _factors(profile_factor, lead_factor)
        q = _prepare_query(p if self.axis == "profile" else l, self.coordinates, self.bounds)
        index = np.clip(np.searchsorted(self.coordinates, q, side="right") - 1, 0, len(self.coordinates) - 2)
        ratio = (q - self.coordinates[index]) / (self.coordinates[index + 1] - self.coordinates[index])
        return self.removal_m[index] * (1 - ratio) + self.removal_m[index + 1] * ratio


@dataclass(frozen=True)
class MeasuredSurface:
    """矩形测量网格的双线性法向偏差；矩阵shape=(齿廓点数,齿向点数)。

    两坐标轴均以各自单轮完整范围归一化，不能使用当前接触区域重新
    缩放；正值材料去除、负值凸起。越界策略与MeasuredCurve一致。
    """
    profile_coordinates: object
    lead_coordinates: object
    removal_m: object
    bounds: Literal["error", "clip", "extrapolate"] = "error"

    def __post_init__(self) -> None:
        p = _coordinates(self.profile_coordinates, "测量曲面齿廓坐标")
        l = _coordinates(self.lead_coordinates, "测量曲面齿向坐标")
        values = np.array(self.removal_m, dtype=float, copy=True)
        if values.shape != (len(p), len(l)) or np.any(~np.isfinite(values)):
            raise ValueError("测量曲面偏差矩阵必须为(齿廓点数,齿向点数)且有限，单位m")
        _check_bounds(self.bounds)
        values.setflags(write=False)
        object.__setattr__(self, "profile_coordinates", p)
        object.__setattr__(self, "lead_coordinates", l)
        object.__setattr__(self, "removal_m", values)

    def evaluate(self, profile_factor, lead_factor):
        p, l = _factors(profile_factor, lead_factor)
        p = _prepare_query(p, self.profile_coordinates, self.bounds)
        l = _prepare_query(l, self.lead_coordinates, self.bounds)
        i = np.clip(np.searchsorted(self.profile_coordinates, p, side="right") - 1, 0, len(self.profile_coordinates) - 2)
        j = np.clip(np.searchsorted(self.lead_coordinates, l, side="right") - 1, 0, len(self.lead_coordinates) - 2)
        u = (p - self.profile_coordinates[i]) / (self.profile_coordinates[i + 1] - self.profile_coordinates[i])
        v = (l - self.lead_coordinates[j]) / (self.lead_coordinates[j + 1] - self.lead_coordinates[j])
        z = self.removal_m
        return (1 - u) * ((1 - v) * z[i, j] + v * z[i, j + 1]) + u * ((1 - v) * z[i + 1, j] + v * z[i + 1, j + 1])


@dataclass(frozen=True)
class CompositeModification:
    """相同理想面基准下逐项相加，可组合参数曲面与测量残差。"""
    components: tuple

    def __post_init__(self) -> None:
        components = tuple(self.components)
        if not all(callable(getattr(item, "evaluate", None)) for item in components):
            raise ValueError("组合修形的每一项都必须提供evaluate方法")
        object.__setattr__(self, "components", components)

    def evaluate(self, profile_factor, lead_factor):
        p, l = _factors(profile_factor, lead_factor)
        answer = np.zeros_like(p)
        for item in self.components:
            answer = answer + item.evaluate(p, l)
        return answer
