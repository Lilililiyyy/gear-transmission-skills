"""Hertz 线接触及具有有限参照深度的平面应变局部压陷。

标准有效模量为 1/E*=(1−ν₁²)/E₁+(1−ν₂²)/E₂。
压力半宽和最大压力使用经典 Hertz 线接触式。
有限深度的局部位移式可见 Greenwood 等的线接触研究：
https://www.sciencedirect.com/science/article/abs/pii/S0043164810004163
本实现的有限参照长度、解析导数和互补能是显式计算契约；不能将
二维线接触位移定义为相对无穷远的绝对位移。它与结构齿梁能分开。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .elasticity import Material


@dataclass(frozen=True)
class HertzLineResult:
    """半宽 m、中心最大压力 Pa 和标准有效模量 Pa。"""

    half_width_m: np.ndarray
    maximum_pressure_pa: np.ndarray
    effective_modulus_pa: float


def effective_modulus_pa(material_pinion=Material(), material_wheel=Material()):
    """两个各向同性半空间的标准有效模量 E*，单位 Pa。"""
    return 1 / (material_pinion.plane_strain_compliance_per_pa
                + material_wheel.plane_strain_compliance_per_pa)


def hertz_line_contact(line_load_n_per_m, curvature_per_m,
                       material_pinion=Material(), material_wheel=Material()):
    """平行局部圆柱的弹性、无摩擦、非共形线接触，允许数组广播。

    q 为法向线载荷 N/m；κ 为法截面组合曲率 1/m。外啮合 κ=1/R₁+1/R₂，
    内啮合 κ=1/R₁−1/R₂；调用者必须传正的有符号组合结果，不能取绝对值。
    b=√[4q/(πE*κ)]，p₀=2q/(πb)。q=0 时两项均为零。
    此函数只有组合曲率，不能检验两轮各自半径和有限齿体的半空间条件；
    局部位移函数会进一步检验 b/R 和 b/D。
    """
    q, curvature = np.broadcast_arrays(np.asarray(line_load_n_per_m, dtype=float),
                                       np.asarray(curvature_per_m, dtype=float))
    if not np.all(np.isfinite(q)) or np.any(q < 0):
        raise ValueError("法向线载荷必须为有限非负数，单位 N/m")
    if not np.all(np.isfinite(curvature)) or np.any(curvature <= 0):
        raise ValueError("Hertz 线接触要求有限正组合曲率，单位 1/m")
    modulus = effective_modulus_pa(material_pinion, material_wheel)
    half_width = np.sqrt(4 * q / (np.pi * modulus * curvature))
    pressure = np.divide(2 * q, np.pi * half_width,
                         out=np.zeros_like(q), where=half_width > 0)
    return HertzLineResult(half_width, pressure, modulus)


def _local_terms(forces_n, widths_m, radius_pinion_m, radius_wheel_m,
                  material_pinion, material_wheel, curvature_per_m,
                  reference_lengths_m, maximum_half_width_ratio):
    force, width, radius_a, radius_b = np.broadcast_arrays(
        np.asarray(forces_n, dtype=float), np.asarray(widths_m, dtype=float),
        np.asarray(radius_pinion_m, dtype=float), np.asarray(radius_wheel_m, dtype=float))
    if not np.all(np.isfinite(force)) or np.any(force < 0):
        raise ValueError("局部接触力必须为有限非负数，单位 N")
    for name, value in (("积分单元长度", width), ("小轮曲率半径", radius_a),
                        ("大轮曲率半径", radius_b)):
        if not np.all(np.isfinite(value)) or np.any(value <= 0):
            raise ValueError(f"{name}必须为有限正数，单位 m")
    if not np.isfinite(maximum_half_width_ratio) or not 0 < maximum_half_width_ratio <= 0.2:
        raise ValueError("接触半宽适用域比例必须在 (0, 0.2] 内")
    curvature = 1 / radius_a + 1 / radius_b if curvature_per_m is None else curvature_per_m
    contact = hertz_line_contact(force / width, curvature, material_pinion, material_wheel)
    half_width = contact.half_width_m
    if half_width.shape != force.shape:
        raise ValueError("组合曲率广播后必须保持接触单元数组形状")
    if reference_lengths_m is None:
        depth_a, depth_b = radius_a, radius_b
    else:
        if len(reference_lengths_m) != 2:
            raise ValueError("有限参照深度需要分别提供两轮的长度")
        depth_a = np.broadcast_to(np.asarray(reference_lengths_m[0], dtype=float), force.shape)
        depth_b = np.broadcast_to(np.asarray(reference_lengths_m[1], dtype=float), force.shape)
        if any(not np.all(np.isfinite(d)) or np.any(d <= 0) for d in (depth_a, depth_b)):
            raise ValueError("有限参照深度必须为有限正数，单位 m")
    length_limit = np.minimum(np.minimum(radius_a, radius_b), np.minimum(depth_a, depth_b))
    if np.any(half_width > maximum_half_width_ratio * length_limit):
        ratio = float(np.max(half_width / length_limit))
        raise ValueError(f"接触半宽超出局部半空间适用域：最大 b/min(R,D)={ratio:.6g}，"
                         f"限值为 {maximum_half_width_ratio:.6g}")
    c_a = material_pinion.plane_strain_compliance_per_pa
    c_b = material_wheel.plane_strain_compliance_per_pa
    # 零力分支不参与对数。由 F ln(F)→0 显式给出零位移/零能量极限。
    safe_b = np.where(force > 0, half_width, 1.0)
    nu_a, nu_b = material_pinion.poisson_ratio, material_wheel.poisson_ratio
    logarithm = (c_a * (np.log(2 * depth_a / safe_b) - nu_a / (2 * (1 - nu_a)))
                 + c_b * (np.log(2 * depth_b / safe_b) - nu_b / (2 * (1 - nu_b))))
    return force, width, logarithm, c_a + c_b


def local_line_approach(forces_n, widths_m, radius_pinion_m, radius_wheel_m,
                         material_pinion=Material(), material_wheel=Material(), *,
                         curvature_per_m=None, reference_lengths_m=None,
                         maximum_half_width_ratio=0.1):
    """返回局部法向压陷 m 及其对单元力的解析导数 m/N。

    q=F/w，cᵢ=(1−νᵢ²)/Eᵢ，L=Σcᵢ[ln(2Dᵢ/b)−νᵢ/(2(1−νᵢ))]。
    h=2qL/π，dh/dF=2[L−Σcᵢ/2]/(πw)。Dᵢ 是接触点下有限深度处
    的位移参照；默认两轮曲率半径，调用者可提供适合有限齿体的深度。
    这是局部半空间近似，不是整个有限圆柱或内齿圈的整体压缩解。

    内啮合必须显式传正的 curvature_per_m；缺省只代表外啮合。
    b/min(Rᵢ,Dᵢ) 默认限制为 0.1，超域明确报错。线接触的沿线端部
    效应没有建模，widths 是积分权重，不能以网格单元宽度替代实体尺寸。
    F=0 的极限为 h=0、dh/dF=+Inf；接触求解应在正力点计算活动集切线。
    结构柔度不得再次包含本局部接触压陷。
    """
    force, width, logarithm, total_c = _local_terms(
        forces_n, widths_m, radius_pinion_m, radius_wheel_m,
        material_pinion, material_wheel, curvature_per_m, reference_lengths_m,
        maximum_half_width_ratio)
    approach = np.where(force > 0, 2 * force * logarithm / (np.pi * width), 0.0)
    derivative = np.where(force > 0, 2 * (logarithm - total_c / 2) / (np.pi * width), np.inf)
    return approach, derivative


def local_line_complementary_energy(forces_n, widths_m, radius_pinion_m,
                                     radius_wheel_m, material_pinion=Material(),
                                     material_wheel=Material(), *, curvature_per_m=None,
                                     reference_lengths_m=None,
                                     maximum_half_width_ratio=0.1):
    """返回逐接触单元的互补能 ∫₀ᶠh(f)df，单位 J。

    因 b∝√F，解析积分为 F²[L(F)+Σcᵢ/4]/(πw)。
    对 F 求导严格得到 local_line_approach；零力极限为零。
    各单元能量相加后可用于法向总力约束下的凸接触优化。
    """
    force, width, logarithm, total_c = _local_terms(
        forces_n, widths_m, radius_pinion_m, radius_wheel_m,
        material_pinion, material_wheel, curvature_per_m, reference_lengths_m,
        maximum_half_width_ratio)
    return np.where(force > 0, force**2 * (logarithm + total_c / 4) / (np.pi * width), 0.0)
