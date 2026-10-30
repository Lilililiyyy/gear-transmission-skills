"""新写条带接触端点压力后处理；中心与端点是两套独立采样布局。

从N个中心线载荷映射至N+1个端点，再使用端点本身的法截面曲率计算
Hertz峰压与半宽。该附加后处理不修改中心压力、接触力或求解结果。
"""
from __future__ import annotations

import numpy as np

from .elasticity import Material
from .hertz import hertz_line_contact


def endpoint_line_loads(widths_m, line_loads_n_per_m):
    """中心线载荷N/m转为N+1端点值，首末线性外推并逐值截至非负。

    中心坐标为累计线长减半个本单元线长。中间端点评价位置取相邻
    中心的中点；非均匀条带时，它不等于累计单元边界。单条带返回
    两个相同端点，不猜测载荷梯度。长度单位m，不修改输入数组。
    """
    widths = np.asarray(widths_m, dtype=float)
    loads = np.asarray(line_loads_n_per_m, dtype=float)
    if (widths.ndim != 1 or not len(widths) or loads.shape != widths.shape
            or np.any(~np.isfinite(widths)) or np.any(widths <= 0)
            or np.any(~np.isfinite(loads)) or np.any(loads < 0)):
        raise ValueError("端点映射需要同形的非空正有限条带长度和非负有限中心线载荷")
    if len(widths) == 1:
        return np.repeat(loads, 2)
    centres = np.cumsum(widths)-widths/2

    def value(i, j, position):
        return loads[i]+(position-centres[i])*(loads[j]-loads[i])/(centres[j]-centres[i])

    result = np.empty(len(widths)+1)
    result[0] = value(0, 1, 0.)
    result[1:-1] = (loads[:-1]+loads[1:])/2
    result[-1] = value(len(widths)-1, len(widths)-2, float(np.sum(widths)))
    return np.maximum(result, 0.)


def _offsets(value, expected_length, *, name):
    offsets = np.asarray(value)
    if (offsets.ndim != 1 or len(offsets) < 2 or not np.issubdtype(offsets.dtype, np.integer)
            or offsets[0] != 0 or offsets[-1] != expected_length or np.any(np.diff(offsets) <= 0)):
        raise ValueError(f"{name}须为从0开始、严格递增并覆盖完整数组的整数偏移")
    return offsets


def extended_endpoint_pressure(extended, forces_n, material_pinion=Material(), material_wheel=Material()):
    """按延伸接触网格的原线布局生成附加端点q、p0及Hertz半宽数组。

    extended提供完整grid、中心line_offsets和独立endpoints布局。q≤0
    的端点压力及半宽为0，保留其曲率原值；不把无载奇异曲率替换成
    人造正数。q>0时曲率必须有限且为正，否则明确拒绝该压力结果。
    """
    width = np.asarray(extended.grid.widths_m, dtype=float)
    forces = np.asarray(forces_n, dtype=float)
    curvature = np.asarray(extended.endpoints.curvature_per_m, dtype=float)
    if (width.ndim != 1 or not len(width) or forces.shape != width.shape
            or np.any(~np.isfinite(width)) or np.any(width <= 0)
            or np.any(~np.isfinite(forces)) or np.any(forces < 0)
            or curvature.ndim != 1):
        raise ValueError("端点压力需要完整正条带长度、同形非负有限接触力及一维端点曲率")
    centres = _offsets(extended.line_offsets, len(width), name="中心线偏移")
    endpoints = _offsets(extended.endpoints.line_offsets, len(curvature), name="端点线偏移")
    if len(centres) != len(endpoints) or np.any(np.diff(endpoints) != np.diff(centres)+1):
        raise ValueError("每条接触线须保留N个中心与N+1个端点，不能拼接或省略端点")
    loads = np.empty(curvature.shape)
    for line in range(len(centres)-1):
        cs = slice(centres[line], centres[line+1])
        es = slice(endpoints[line], endpoints[line+1])
        loads[es] = endpoint_line_loads(width[cs], forces[cs]/width[cs])
    pressure = np.zeros_like(loads)
    half_width = np.zeros_like(loads)
    loaded = loads > 0
    if np.any(loaded):
        if np.any(~np.isfinite(curvature[loaded])) or np.any(curvature[loaded] <= 0):
            raise ValueError("正载端点的Hertz组合曲率必须有限且为正；不能以绝对值或微小半径修正")
        result = hertz_line_contact(loads[loaded], curvature[loaded], material_pinion, material_wheel)
        pressure[loaded] = result.maximum_pressure_pa
        half_width[loaded] = result.half_width_m
    return {"line_load_n_per_m": loads, "pressure_pa": pressure, "half_width_m": half_width}
