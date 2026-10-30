"""新写ISO 6336单齿刚度与Love条带耦合模型，内部采用m、N、Pa、rad。

ISO等效单齿刚度已代表综合接触变形，本模式不能再叠加独立Hertz
压陷。Hertz公式仍可用于求解后的半宽和压力。Love修正按接触线去
普通均值，等效算子一般不对称；本模块不将其伪装成互易弹性矩阵。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .elasticity import Material
from .geometry import GearPair
from .solver import ContactSolution, ContactConvergenceError, solve_contact


@dataclass(frozen=True)
class ISOContactSettings:
    """ISO刚度输入与耦合求解设置。

    nominal_pinion_torque_n_m是参考圆名义扭矩；application_factor只
    进入低载修正的Ft/b，不能用LTCA评级后的目标扭矩替代名义扭矩。
    material_reference_youngs_modulus_pa显式选择207或206 GPa等材料
    基准。未给定根高系数时从当前实体根圆、参考圆和变位反算，不假定
    内外齿均为1.25。耦合松弛和迭代保护为本模块新增求解策略。
    """
    nominal_pinion_torque_n_m: float
    application_factor: float = 1.0
    material_reference_youngs_modulus_pa: float = 207e9
    dedendum_coeff_pinion: float | None = None
    dedendum_coeff_wheel: float | None = None
    coupling_relaxation: float = 1.0
    maximum_coupling_iterations: int = 500
    coupling_tolerance_m: float = 2e-12

    def __post_init__(self):
        if (isinstance(self.nominal_pinion_torque_n_m, bool)
                or not math.isfinite(self.nominal_pinion_torque_n_m)
                or self.nominal_pinion_torque_n_m < 0):
            raise ValueError("ISO名义扭矩必须为非负有限数，单位N·m")
        for name in ("application_factor", "material_reference_youngs_modulus_pa", "coupling_tolerance_m"):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name}必须为正有限数")
        for name in ("dedendum_coeff_pinion", "dedendum_coeff_wheel"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name}须为正有限根高系数，或省略以由实际几何计算")
        if (isinstance(self.coupling_relaxation, bool) or not math.isfinite(self.coupling_relaxation)
                or not 0 < self.coupling_relaxation <= 1):
            raise ValueError("耦合松弛系数必须在(0,1]内")
        if (isinstance(self.maximum_coupling_iterations, bool)
                or not isinstance(self.maximum_coupling_iterations, int)
                or self.maximum_coupling_iterations < 1):
            raise ValueError("耦合迭代上限必须为正整数")


def iso_single_tooth_stiffness(pair: GearPair, material_pinion: Material,
                               material_wheel: Material, settings: ISOContactSettings) -> dict:
    """计算单位齿宽单齿刚度c′，单位N/m²，包含名义载荷的1/4次修正。

    有符号当量齿数zv=z/(cos²βb*cosβ)，内齿为负。经验式的倒刚度
    系数乘1e-9转换为m²/N；100 N/mm对应100000 N/m。零名义载荷
    采用低载修正1，保留这一明确分支，不把零载经验刚度解释为实验值。
    """
    p, w = pair.pinion, pair.wheel
    beta = abs(p.helix_angle_rad)
    denominator = math.cos(p.base_helix_angle_rad)**2*math.cos(beta)
    za = p.teeth/denominator
    zb = (-1 if w.internal else 1)*w.teeth/denominator
    xa, xb = p.profile_shift, w.profile_shift
    if not (za > 0 and za <= abs(zb)):
        za, zb, xa, xb = zb, za, xb, xa
    if za <= 0 or zb == 0:
        raise ValueError("ISO刚度需要至少一个正的外齿当量齿数")
    reciprocal = (0.04723+0.15551/za-0.00635*xa-0.11654*xa/za
                  -0.00193*xb+0.00529*xa**2+0.00182*xb**2)
    if zb >= 0:
        reciprocal += 0.25791/zb-0.24188*xb/zb
    reciprocal *= 1e-9
    if not math.isfinite(reciprocal) or reciprocal <= 0:
        raise ValueError("当前参数使ISO经验倒刚度非正，不能用于接触计算")

    def dedendum(gear, override):
        if override is not None:
            return override
        sign = -1 if gear.internal else 1
        value = sign*(gear.reference_radius_m-gear.root_radius_m)/gear.normal_module_m+gear.profile_shift
        if not math.isfinite(value) or value <= 0:
            raise ValueError("由实际根圆反算的ISO根高系数非正，请明确给定适用值")
        return value

    ha = dedendum(p, settings.dedendum_coeff_pinion)
    hb = dedendum(w, settings.dedendum_coeff_wheel)
    rack_a = (1+.5*(1.2-ha))*(1-.02*(20-math.degrees(p.pressure_angle_rad)))
    rack_b = (1+.5*(1.2-hb))*(1-.02*(20-math.degrees(w.pressure_angle_rad)))
    rack_factor = .5*(rack_a+rack_b)
    ea, eb = material_pinion.youngs_modulus_pa, material_wheel.youngs_modulus_pa
    material_factor = (2*ea*eb/(ea+eb))/settings.material_reference_youngs_modulus_pa
    ft = settings.nominal_pinion_torque_n_m/p.reference_radius_m
    load_per_width = ft*settings.application_factor/pair.face_overlap_m
    low_load_factor = (load_per_width/100000.)**.25 if 0 < load_per_width < 100000. else 1.
    stiffness = 1/reciprocal*.8*rack_factor*math.cos(beta)*material_factor*low_load_factor
    if not math.isfinite(stiffness) or stiffness <= 0:
        raise ValueError("ISO单齿刚度非正或非有限，请检查参数适用域")
    return {"model": "ISO6336单齿等效刚度", "single_stiffness_n_per_m2": stiffness,
        "signed_virtual_teeth_pinion": p.teeth/denominator,
        "signed_virtual_teeth_wheel": (-1 if w.internal else 1)*w.teeth/denominator,
        "reciprocal_stiffness_m2_per_n": reciprocal, "rack_factor": rack_factor,
        "dedendum_coeff_pinion": ha, "dedendum_coeff_wheel": hb,
        "material_factor": material_factor,
        "material_reference_youngs_modulus_pa": settings.material_reference_youngs_modulus_pa,
        "nominal_tangential_force_n": ft, "load_per_width_n_per_m": load_per_width,
        "low_load_factor": low_load_factor, "local_contact_included": True}


def _layout(grid):
    width = np.asarray(grid.widths_m, dtype=float)
    ids = np.asarray(grid.tooth_pair_index)
    if (width.ndim != 1 or not len(width) or ids.shape != width.shape
            or np.any(~np.isfinite(width)) or np.any(width <= 0)):
        raise ValueError("ISO接触网格必须保留同形的正有限条带长度和齿对索引")
    return width, ids


def build_iso_compliance(pair, grid, material_pinion, material_wheel, settings):
    """返回ISO法向对角柔度m/N及诊断；本矩阵不含Love非对称修正。

    width为实际接触线长度，不能用轴向宽度替代；法向条带刚度为
    c′*width/cosβb，故对角柔度为cosβb/(c′*width)。
    """
    width, _ = _layout(grid)
    details = iso_single_tooth_stiffness(pair, material_pinion, material_wheel, settings)
    diagonal = math.cos(pair.base_helix_angle_rad)/(details["single_stiffness_n_per_m2"]*width)
    return np.diag(diagonal), details


def love_coupling_matrix(grid, material_pinion):
    """构造逐接触线去普通均值的Love位移算子，单位m/N。

    非对角影响为2*(1-νp²)/(π*Ep)*log((2d+1)/(2d-1))/width_i；
    只使用小轮材料。均值扣除针对每条线的全部条带，不只针对活动点。
    因此算子一般不对称；不能直接传入要求弹性互易性的solve_contact。
    """
    width, ids = _layout(grid)
    operator = np.zeros((len(width), len(width)))
    scale = 2*(1-material_pinion.poisson_ratio**2)/(math.pi*material_pinion.youngs_modulus_pa)
    for line in np.unique(ids):
        indices = np.flatnonzero(ids == line)
        distance = np.abs(np.arange(len(indices))[:, None]-np.arange(len(indices))[None, :])
        influence = np.zeros(distance.shape, dtype=float)
        mask = distance > 0
        influence[mask] = scale*np.log((2*distance[mask]+1)/(2*distance[mask]-1))
        block = influence/width[indices, None]
        block -= block.mean(axis=0, keepdims=True)
        operator[np.ix_(indices, indices)] = block
    return operator


def _solve_coupled(gaps, diagonal, operator, force, settings):
    if force == 0:
        approach = float(np.min(gaps))
        minimum = gaps == approach
        # 右切线仍需施加单边约束；只对最小间隙分支求单位力线性化。
        right = _solve_coupled(np.zeros(np.count_nonzero(minimum)), diagonal[minimum],
            operator[np.ix_(minimum, minimum)], 1., settings)
        return ContactSolution(np.zeros_like(gaps), approach, gaps-approach,
            np.zeros_like(gaps, dtype=bool), right.tangent_stiffness_n_per_m,
            {"converged": True, "iterations": 0, "force_balance_n": 0.,
             "maximum_penetration_m": 0., "active_gap_error_m": 0., "complementarity_n_m": 0.,
             "zero_load": True, "gap_tolerance_m": settings.coupling_tolerance_m,
             "coupling_iterations": 0, "model_energy_available": False})
    displacement = np.zeros_like(gaps)
    previous = None
    matrix = np.diag(diagonal)
    deformation_scale = max(float(np.median(diagonal))*force/len(gaps), 1e-9)
    inner_relative_tolerance = max(np.finfo(float).eps,
        min(2e-9, settings.coupling_tolerance_m/(10*deformation_scale)))
    for iteration in range(1, settings.maximum_coupling_iterations+1):
        frozen = solve_contact(gaps+displacement, matrix, force, initial_forces_n=previous,
            relative_tolerance=inner_relative_tolerance,
            absolute_gap_tolerance_m=settings.coupling_tolerance_m/10)
        forces = frozen.forces_n
        target = operator@forces
        positive = forces > force*1e-13
        deformation = gaps+diagonal*forces+target
        approach = float(np.mean(deformation[positive]))
        residual = deformation-approach
        penetration = float(max(0., -np.min(residual)))
        active_error = float(np.max(np.abs(residual[positive])))
        if max(penetration, active_error) <= settings.coupling_tolerance_m:
            break
        displacement += settings.coupling_relaxation*(target-displacement)
        previous = forces
    else:
        raise ContactConvergenceError("ISO/Love耦合达到迭代上限，未满足接触协调，不返回有效结果")
    active_matrix = (matrix+operator)[np.ix_(positive, positive)]
    try:
        tangent = float(np.sum(np.linalg.solve(active_matrix, np.ones(np.count_nonzero(positive)))))
    except np.linalg.LinAlgError as error:
        raise ContactConvergenceError("ISO/Love活动分支切线矩阵奇异") from error
    if not math.isfinite(tangent) or tangent <= 0:
        raise ContactConvergenceError("ISO/Love活动分支的总法向切线非正")
    diagnostics = {"converged": True, "iterations": iteration, "coupling_iterations": iteration,
        "force_balance_n": float(np.sum(forces)-force), "maximum_penetration_m": penetration,
        "active_gap_error_m": active_error, "complementarity_n_m": float(np.max(np.abs(forces*residual))),
        "active_cells": int(np.count_nonzero(positive)), "gap_tolerance_m": settings.coupling_tolerance_m,
        "model_energy_available": False,
        "coupling_displacement_peak_to_peak_m": float(np.ptp(target))}
    return ContactSolution(forces, approach, residual, positive, tangent, diagnostics)


def solve_iso_contact(pair, grid, gaps_m, normal_force_n, material_pinion, material_wheel,
                      settings, *, coupled=True, eligible=None):
    """求ISO单边条带接触；Love模式用冻结位移外迭代，不对算子强制对称化。

    返回(ContactSolution,弹性诊断)。切线是当前活动分支对总法向载荷的
    导数；未进行额外有限差分扰动，不把数值收敛声明为物理精度。
    eligible可排除几何上不允许接触的候选。Love仍先按完整条带布局
    构造，再限制到允许集合；输出保留完整数组及排除点的真实间隙，
    接触残差验收只针对允许集合，不将排除点间隙伪置零。
    """
    gaps = np.asarray(gaps_m, dtype=float)
    width, _ = _layout(grid)
    if (gaps.shape != width.shape or np.any(~np.isfinite(gaps))
            or not math.isfinite(normal_force_n) or normal_force_n < 0):
        raise ValueError("ISO接触间隙须同形且有限，总法向力须非负有限")
    allowed = np.ones(width.shape, dtype=bool) if eligible is None else np.asarray(eligible, dtype=bool)
    if allowed.shape != width.shape or not np.any(allowed):
        raise ValueError("ISO接触允许标志须与完整条带同形，且至少有一个允许候选")
    compliance, details = build_iso_compliance(pair, grid, material_pinion, material_wheel, settings)
    diagonal = np.diag(compliance)
    if coupled:
        operator = love_coupling_matrix(grid, material_pinion)
        solution = _solve_coupled(gaps[allowed], diagonal[allowed], operator[np.ix_(allowed, allowed)],
                                  normal_force_n, settings)
        scale = float(np.max(np.abs(compliance+operator)))
        details = dict(details, model="ISO6336单齿刚度与去均值Love耦合",
            influence_operator_reciprocal=bool(np.allclose(operator, operator.T, rtol=1e-10, atol=scale*1e-13)),
            influence_operator_asymmetry_m_per_n=float(np.max(np.abs(operator-operator.T))),
            coupling_reference_material="pinion", coupling_iterations=solution.diagnostics["coupling_iterations"])
    else:
        operator = np.zeros_like(compliance)
        solution = solve_contact(gaps[allowed], compliance[np.ix_(allowed, allowed)], normal_force_n)
        details = dict(details, model="ISO6336单齿刚度，无齿宽Love耦合",
            influence_operator_reciprocal=True, coupling_iterations=0)
    if not np.all(allowed):
        force_full = np.zeros_like(gaps)
        force_full[allowed] = solution.forces_n
        active_full = np.zeros_like(allowed)
        active_full[allowed] = solution.active_mask
        residual_full = gaps+diagonal*force_full+operator@force_full-solution.approach_m
        diagnostics = dict(solution.diagnostics,
            excluded_candidates=int(np.count_nonzero(~allowed)),
            maximum_excluded_penetration_m=float(max(0., -np.min(residual_full[~allowed]))),
            residual_check_scope="eligible候选；排除点的真实间隙仍保留在完整数组")
        solution = ContactSolution(force_full, solution.approach_m, residual_full, active_full,
                                   solution.tangent_stiffness_n_per_m, diagnostics)
    details["excluded_candidates"] = int(np.count_nonzero(~allowed))
    return solution, details
