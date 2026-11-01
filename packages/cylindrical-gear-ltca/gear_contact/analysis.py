"""齿面几何、应变能柔度、单边接触与压力后处理的统一计算入口。"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable

import numpy as np

from .geometry import GearPair, contact_grid
from .surfaces import SurfaceModification
from .elasticity import Material, build_structural_compliance
from .hertz import hertz_line_contact, local_line_approach, local_line_complementary_energy
from .solver import ContactSolution, solve_contact
from .iso_contact import ISOContactSettings, solve_iso_contact
from .extended_contact import ExtendedContactSettings, extended_contact_grid


@dataclass
class RotationResult:
    """一个啮合位置的原始离散接触结果，不包含展示插值。"""
    grid: object
    gaps_m: np.ndarray
    solution: ContactSolution
    line_load_n_per_m: np.ndarray
    pressure_pa: np.ndarray
    half_width_m: np.ndarray
    elastic_diagnostics: dict
    contact_geometry_data: object | None = None
    endpoint_pressure: dict | None = None


@dataclass
class AnalysisResult:
    """一个齿轮副、一个工况的单周期结果。"""
    pair: GearPair
    phase: np.ndarray
    normal_approach_m: np.ndarray
    angular_te_rad: np.ndarray
    mesh_stiffness_n_per_m: np.ndarray
    rotations: list[RotationResult]
    summary: dict
    inputs: dict


def transmission_error_spectrum(values, count: int = 10, *, phase_origin: float = 0.0) -> dict:
    """均匀单周期采样的实数傅里叶谱。

    第一个样本位于phase_origin，返回相位仍参考物理相位0，约定
    x(t)=mean+ΣA cos(2πnt+φ)。仅适用于均匀完整周期采样。
    正常正频率幅值乘2，Nyquist项不乘2。近零幅值相位为NaN，避免
    将浮点噪声报告为确定相位；均值单独返回，频率单位为啮合阶次。
    """
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or a.size < 2 or not np.all(np.isfinite(a)):
        raise ValueError("频谱需要至少两个有限、均匀单周期样本")
    if not isinstance(count, int) or count < 1:
        raise ValueError("谐波阶数必须为正整数")
    if not np.isfinite(phase_origin):
        raise ValueError("频谱的起采物理相位必须为有限数")
    coefficients = np.fft.rfft(a)/a.size
    last = min(count, a.size//2)
    orders = np.arange(1, last+1)
    amplitudes = 2*np.abs(coefficients[orders])
    if a.size % 2 == 0 and last == a.size//2:
        amplitudes[-1] /= 2
    selected = coefficients[orders]
    if phase_origin != 0:
        selected = selected * np.exp(-2j * np.pi * orders * phase_origin)
    phase = np.angle(selected)
    threshold = max(float(np.max(np.abs(a)))*1e-12, np.finfo(float).tiny)
    phase[amplitudes < threshold] = np.nan
    return {"mean": float(np.mean(a)), "order": orders,
            "amplitude": amplitudes, "phase_rad": phase}


class _LocalLaw:
    def __init__(self, grid, pinion_material, wheel_material, reference_depths=None):
        self.args = (grid.widths_m, grid.radius_pinion_m, grid.radius_wheel_m,
                     pinion_material, wheel_material)
        self.kwargs = {"curvature_per_m": grid.curvature_per_m,
                       "reference_lengths_m": reference_depths}

    def __call__(self, forces):
        return local_line_approach(forces, *self.args, **self.kwargs)

    def energy(self, forces):
        return local_line_complementary_energy(forces, *self.args, **self.kwargs)


def analyze_pair(pair: GearPair, *, torque_n_m: float,
                 material_pinion: Material = Material(), material_wheel: Material = Material(),
                 modification_pinion=SurfaceModification(), modification_wheel=SurfaceModification(),
                 misalignment_m: float = 0.0, phases: int = 32, strips: int = 21,
                 phase_offset: float = 0.0,
                 elastic_mode: str = "coupled", local_contact: bool = True,
                 iso_settings: ISOContactSettings | dict | None = None,
                 extended_contact_settings: ExtendedContactSettings | dict | None = None,
                 local_reference_depths_m=None,
                 compliance_provider: Callable | None = None,
                 boundary_gap: Callable | None = None) -> AnalysisResult:
    """求一个均匀采样啮合周期的静态接触响应。

    torque_n_m 为已选齿面的非负载荷幅值；正反受载面的修形由调用方明确
    选择。misalignment_m 是有效重叠齿宽右端减左端的法向间隙。
    可用 compliance_provider(pair,grid) 接入自行生成的结构柔度矩阵，
    该矩阵不能再包含本函数 local_contact 所计算的局部压陷。
    boundary_gap(phase,grid) 可提供附加法向边界间隙，长度单位m。
    local_reference_depths_m=(D小轮,D大轮) 可显式给出局部位移参照深度；
    省略时每个单元以相应曲率半径为有限参照，不能解释为无限远位移。
    phase_offset 为均匀完整周期的起采相位，单位为基节周期；仅通过
    接触线几何配准确定，不允许通过响应曲线拟合来最小化对标差异。
    elastic_mode为iso/iso_love时采用ISO综合齿对刚度，iso_settings
    必须提供名义载荷，且local_contact必须False，避免重复计入压陷。
    """
    if not np.isfinite(torque_n_m) or torque_n_m < 0 or not np.isfinite(misalignment_m):
        raise ValueError("扭矩幅值必须非负有限，失准间隙必须有限")
    if isinstance(phases, bool) or not isinstance(phases, int) or phases < 2:
        raise ValueError("每周期相位数必须为至少2的整数")
    if isinstance(strips, bool) or not isinstance(strips, int) or strips < 2:
        raise ValueError("每条接触线的条带数必须为至少2的整数")
    if isinstance(phase_offset, bool) or not np.isfinite(phase_offset):
        raise ValueError("起采相位必须为有限实数，以基节周期计")
    phase_offset = float(phase_offset) % 1.0
    iso_mode = elastic_mode in ("iso", "iso_love")
    if iso_mode:
        if local_contact:
            raise ValueError("ISO刚度已包含综合接触变形，必须显式设local_contact=False")
        if compliance_provider is not None:
            raise ValueError("ISO模式不能同时提供外部结构柔度")
        if iso_settings is None:
            raise ValueError("ISO模式需要iso_settings明确给出名义扭矩，不能使用评级后目标代替")
        if isinstance(iso_settings, dict):
            iso_settings = ISOContactSettings(**iso_settings)
        if not isinstance(iso_settings, ISOContactSettings):
            raise ValueError("iso_settings必须为ISOContactSettings或对应字段对象")
    elif iso_settings is not None:
        raise ValueError("iso_settings只适用于iso或iso_love模式")
    if extended_contact_settings is not None:
        if not iso_mode:
            raise ValueError("齿顶延伸目前只与ISO模式组合，不能用于未定义延伸齿梁截面的模型")
        if isinstance(extended_contact_settings, dict):
            extended_contact_settings = ExtendedContactSettings(**extended_contact_settings)
        if not isinstance(extended_contact_settings, ExtendedContactSettings):
            raise ValueError("extended_contact_settings必须为ExtendedContactSettings或对应字段对象")
    force = pair.normal_force_from_torque(torque_n_m)
    phase_array = (np.arange(phases, dtype=float)/phases + phase_offset) % 1.0
    rotations = []
    zlo = max(-pair.pinion.face_width_m/2, pair.axial_offset_m-pair.wheel.face_width_m/2)
    zhi = min(pair.pinion.face_width_m/2, pair.axial_offset_m+pair.wheel.face_width_m/2)
    zmid = (zlo+zhi)/2
    centroids, concentrations, contact_fractions = [], [], []
    for phase in phase_array:
        extended = (extended_contact_grid(pair, float(phase), strips=strips, settings=extended_contact_settings)
                    if extended_contact_settings is not None else None)
        grid = extended.grid if extended is not None else contact_grid(pair, float(phase), strips=strips)
        if not iso_mode and (np.any(~np.isfinite(grid.curvature_per_m)) or np.any(grid.curvature_per_m <= 0)
                or np.any(grid.radius_pinion_m <= 0) or np.any(grid.radius_wheel_m <= 0)):
            raise ValueError("接触网格含基圆尖点或非正组合曲率，局部线接触近似不适用；"
                             "请明确给定实际有效根形径并检查啮合干涉，不能用微小半径替代")
        removal_p = modification_pinion.evaluate(grid.profile_factor_pinion, grid.lead_factor_pinion)
        removal_w = modification_wheel.evaluate(grid.profile_factor_wheel, grid.lead_factor_wheel)
        gaps = np.asarray(removal_p+removal_w + misalignment_m*(grid.face_coordinate_m-zmid)/(zhi-zlo), dtype=float)
        if extended is not None:
            gaps = gaps+extended.gap_correction_m
        if boundary_gap is not None:
            additional = np.asarray(boundary_gap(float(phase), grid), dtype=float)
            if additional.shape != gaps.shape or not np.all(np.isfinite(additional)):
                raise ValueError("边界间隙回调必须返回与接触网格同形的有限米制数组")
            gaps = gaps+additional
        if iso_mode:
            solution, elastic_diagnostics = solve_iso_contact(pair, grid, gaps, force,
                material_pinion, material_wheel, iso_settings, coupled=elastic_mode == "iso_love",
                eligible=extended.eligible if extended is not None else None)
        elif compliance_provider is None:
            compliance, elastic_diagnostics = build_structural_compliance(
                pair, grid, material_pinion, material_wheel, mode=elastic_mode)
        else:
            compliance = compliance_provider(pair, grid)
            elastic_diagnostics = {"model": "调用方给定的结构柔度矩阵", "local_contact_included": False}
        if not iso_mode:
            local_law = _LocalLaw(grid, material_pinion, material_wheel, local_reference_depths_m) if local_contact else None
            solution = solve_contact(gaps, compliance, force, local_law=local_law)
        line_load = solution.forces_n/grid.widths_m
        if not iso_mode:
            pressure = hertz_line_contact(line_load, grid.curvature_per_m, material_pinion, material_wheel)
            pressure_values = np.asarray(pressure.maximum_pressure_pa)
            half_width_values = np.asarray(pressure.half_width_m)
        else:
            loaded = solution.forces_n > 0
            if np.any(loaded & (~np.isfinite(grid.curvature_per_m) | (grid.curvature_per_m <= 0))):
                raise ValueError("承载接触点曲率无效，不能计算Hertz压力")
            pressure_values, half_width_values = np.zeros_like(line_load), np.zeros_like(line_load)
            if np.any(loaded):
                pressure = hertz_line_contact(line_load[loaded], grid.curvature_per_m[loaded], material_pinion, material_wheel)
                pressure_values[loaded], half_width_values[loaded] = pressure.maximum_pressure_pa, pressure.half_width_m
        endpoint_pressure = None
        if extended is not None:
            from .contact_pressure import extended_endpoint_pressure
            endpoint_pressure = extended_endpoint_pressure(extended, solution.forces_n, material_pinion, material_wheel)
        rotations.append(RotationResult(grid, gaps, solution, line_load,
                                        pressure_values, half_width_values, elastic_diagnostics, extended, endpoint_pressure))
        centroids.append(float(np.dot(solution.forces_n, grid.face_coordinate_m-zmid)/force) if force > 0 else 0.0)
        mean_load = force/float(np.sum(grid.widths_m))
        concentrations.append(float(np.max(line_load)/mean_load) if force > 0 else 0.0)
        contact_fractions.append(float(np.sum(grid.widths_m[solution.active_mask])/np.sum(grid.widths_m)))
    approach = np.array([r.solution.approach_m for r in rotations])
    angular = np.array([pair.normal_approach_to_angular_te(value) for value in approach])
    stiffness = np.array([r.solution.tangent_stiffness_n_per_m for r in rotations])
    summary = {
        "torque_n_m": float(torque_n_m), "normal_force_n": float(force),
        "mean_normal_approach_m": float(np.mean(approach)),
        "te_peak_to_peak_m": float(np.ptp(approach)),
        "te_rms_ac_m": float(np.sqrt(np.mean((approach-np.mean(approach))**2))),
        "angular_te_peak_to_peak_rad": float(np.ptp(angular)),
        "max_pressure_pa": max(float(np.max(r.pressure_pa)) for r in rotations),
        "max_line_load_n_per_m": max(float(np.max(r.line_load_n_per_m)) for r in rotations),
        "mean_mesh_stiffness_n_per_m": float(np.mean(stiffness)),
        "mesh_stiffness_fluctuation_n_per_m": float(np.ptp(stiffness)),
        "max_line_load_concentration_factor": max(concentrations),
        "max_abs_load_centroid_offset_m": max(abs(value) for value in centroids),
        "mean_contact_fraction": float(np.mean(contact_fractions)),
        "max_force_residual_n": max(abs(r.solution.diagnostics["force_balance_n"]) for r in rotations),
        "max_gap_residual_m": max(max(r.solution.diagnostics["maximum_penetration_m"],
                                        r.solution.diagnostics["active_gap_error_m"]) for r in rotations),
        "all_converged": all(r.solution.diagnostics["converged"] for r in rotations),
        "elastic_mode": elastic_mode if compliance_provider is None else "provided",
        "local_contact": bool(local_contact),
        "extended_contact": extended_contact_settings is not None,
        "transverse_contact_ratio": float(pair.transverse_contact_ratio),
        "overlap_ratio": float(pair.overlap_ratio),
        "model_scope": ("ISO综合齿对刚度、可选Love条带影响；无摩擦、准静态，压力采用局部线接触近似"
                        if iso_mode else "固定齿根弹性降阶、无摩擦、小变形、准静态；压力采用局部线接触近似"),
    }
    if extended_contact_settings is not None:
        summary["max_endpoint_pressure_pa"] = max(float(np.max(r.endpoint_pressure["pressure_pa"])) for r in rotations)
    inputs = {"torque_n_m": float(torque_n_m), "misalignment_m": float(misalignment_m),
              "phase_offset": phase_offset,
              "phases": phases, "strips": strips, "elastic_mode": elastic_mode,
              "local_contact": bool(local_contact), "material_pinion": material_pinion,
              "local_reference_depths_m": local_reference_depths_m,
              "material_wheel": material_wheel, "modification_pinion": modification_pinion,
              "modification_wheel": modification_wheel,
              "external_compliance_provider": compliance_provider is not None,
              "external_boundary_gap": boundary_gap is not None}
    if iso_settings is not None:
        inputs["iso_settings"] = asdict(iso_settings)
    if extended_contact_settings is not None:
        inputs["extended_contact_settings"] = asdict(extended_contact_settings)
    return AnalysisResult(pair, phase_array, approach, angular, stiffness, rotations, summary, inputs)
