"""新写的ISO 6336:2006显式因子名义强度链，单位m、N、Pa、rad。

几何、弹性系数和名义力由当前输入求得；齿形、应力修正、载荷分配等
评级因子须由调用方显式提供。本模块不独立预测这些因子，不生成齿根
过渡，不等同完整ISO评级软件。材料/寿命应力限为可选外部输入；标准
接触应力是评级标量，不能作为LTCA空间压力或接触斑点。
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import math

from .elasticity import Material
from .geometry import GearPair


def _positive_values(record):
    """因子与应力限均须为正有限数，不以缺省猜测代替缺失的输入。"""
    for field in fields(record):
        value = getattr(record, field.name)
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{field.name} 必须为正有限数")


@dataclass(frozen=True)
class MeshRatingFactors:
    """评级公共因子，全部无量纲；不由LTCA峰值压力反标。

    additional_load_factor保留来源额外乘子，不能未经核验将其解释为
    行星均载系数。无载荷放大的基础名义计算可显式采用相应因子1。
    """
    application_factor: float
    dynamic_factor: float
    additional_load_factor: float
    face_contact_factor: float
    transverse_contact_factor: float
    face_bending_factor: float
    transverse_bending_factor: float
    contact_ratio_factor: float
    helix_contact_factor: float
    helix_bending_factor: float
    root_contact_ratio_factor: float
    deep_tooth_factor: float

    def __post_init__(self):
        _positive_values(self)


@dataclass(frozen=True)
class GearRatingFactors:
    """单轮评级输入；bF为齿根应力有效宽度m，其余无量纲。

    YF与YS应来自对应制造齿根几何或经验证的标准因子计算，不以齿数
    粗略替代；ZB/ZD保留各轮单对齿接触因子，YB为轮缘厚度因子。
    """
    form_factor: float
    stress_correction_factor: float
    rim_factor: float
    single_pair_contact_factor: float
    root_stress_face_width_m: float

    def __post_init__(self):
        _positive_values(self)


@dataclass(frozen=True)
class StressLimits:
    """给定材料/寿命/可靠度下的sigmaHG、sigmaFG，单位Pa。

    这里接收的是尚未除以最低设计安全系数的应力限；如果调用方只有
    sigmaHP或sigmaFP，必须明确恢复其安全系数基准后再使用，不能混用。
    本模块不独立建立材料S-N曲线或重新推导寿命修正系数。
    """
    contact_limit_pa: float
    bending_limit_pa: float

    def __post_init__(self):
        _positive_values(self)


def evaluate_iso6336_explicit_factors(
    pair: GearPair, *, nominal_pinion_torque_n_m: float, face_width_m: float,
    mesh_factors: MeshRatingFactors, pinion_factors: GearRatingFactors,
    wheel_factors: GearRatingFactors, material_pinion: Material = Material(),
    material_wheel: Material = Material(), pinion_limits: StressLimits | None = None,
    wheel_limits: StressLimits | None = None,
) -> dict:
    """由名义扭矩和显式因子计算两轮2006评级应力及可选安全系数。

    Ft=2*Tnom/d1，在参考圆定义，不能以LTCA评级放大后的目标扭矩代入。
    sigmaH0=ZH*ZE*Zepsilon*Zbeta*sqrt(Ft/(b*d1)*(u+1)/u)，内啮合u取负。
    sigmaH=sigmaH0*ZB/D*sqrt(KA*KV*Kextra*KHbeta*KHalpha)。
    sigmaF0=Ft/(bF*mn)*YF*YS*Ybeta*Yepsilon*YB*YDT；
    sigmaF=sigmaF0*KA*KV*Kextra*KFbeta*KFalpha。

    该接口保留显式来源因子，只复算应力组合链；改变几何时调用方应
    重新计算依赖几何的外部因子，不能把固定因子比较当作完整重评级。
    零载时应力为0，安全系数返回None，避免输出无穷JSON或虚假有限值。
    """
    if (isinstance(nominal_pinion_torque_n_m, bool) or not math.isfinite(nominal_pinion_torque_n_m)
            or nominal_pinion_torque_n_m < 0):
        raise ValueError("名义小轮扭矩幅值必须为非负有限数，单位N·m")
    if isinstance(face_width_m, bool) or not math.isfinite(face_width_m) or face_width_m <= 0:
        raise ValueError("评级接触齿宽必须为正有限数，单位m")
    if face_width_m > pair.face_overlap_m * (1+1e-10):
        raise ValueError("评级接触齿宽不能大于当前几何的重叠齿宽")
    p, w, f = pair.pinion, pair.wheel, mesh_factors
    d1 = 2*p.reference_radius_m
    u = (-1 if w.internal else 1)*w.teeth/p.teeth
    ratio_term = (u+1)/u
    if ratio_term <= 0:
        raise ValueError("有符号齿数比不能形成正的接触评级曲率项")
    ft = 2*nominal_pinion_torque_n_m/d1
    alpha_t = p.transverse_pressure_angle_rad
    alpha_wt = pair.working_pressure_angle_rad
    zone_factor = math.sqrt(2*math.cos(p.base_helix_angle_rad) /
                            (math.cos(alpha_t)**2*math.tan(alpha_wt)))
    elastic_factor = math.sqrt(1/(math.pi*((1-material_pinion.poisson_ratio**2)/material_pinion.youngs_modulus_pa
                                +(1-material_wheel.poisson_ratio**2)/material_wheel.youngs_modulus_pa)))
    contact_nominal = (zone_factor*elastic_factor*f.contact_ratio_factor*f.helix_contact_factor
                       * math.sqrt(ft/(face_width_m*d1)*ratio_term))
    common_load = f.application_factor*f.dynamic_factor*f.additional_load_factor
    contact_load = math.sqrt(common_load*f.face_contact_factor*f.transverse_contact_factor)
    bending_load = common_load*f.face_bending_factor*f.transverse_bending_factor
    gears = {}
    for role, factors, limits in (("pinion", pinion_factors, pinion_limits),
                                  ("wheel", wheel_factors, wheel_limits)):
        root_nominal = (ft/(factors.root_stress_face_width_m*p.normal_module_m)
                        *factors.form_factor*factors.stress_correction_factor*f.helix_bending_factor
                        *f.root_contact_ratio_factor*factors.rim_factor*f.deep_tooth_factor)
        contact = contact_nominal*factors.single_pair_contact_factor*contact_load
        bending = root_nominal*bending_load
        gears[role] = {"nominal_root_stress_pa": root_nominal, "contact_stress_pa": contact,
            "root_stress_pa": bending,
            "contact_safety_factor": limits.contact_limit_pa/contact if limits is not None and contact > 0 else None,
            "bending_safety_factor": limits.bending_limit_pa/bending if limits is not None and bending > 0 else None,
            "stress_limits": asdict(limits) if limits else None, "factors": asdict(factors)}
    return {"method": "ISO6336_2006_explicit_factors", "implementation": "新写显式因子名义评级链",
        "scope": "由当前几何、材料及名义载荷复算应力组合；外部评级因子和材料/寿命应力限不在本模块预测范围",
        "nominal_pinion_torque_n_m": float(nominal_pinion_torque_n_m),
        "nominal_tangential_force_n": ft, "face_width_m": float(face_width_m),
        "zone_factor": zone_factor, "elastic_factor_sqrt_pa": elastic_factor,
        "nominal_contact_stress_pa": contact_nominal, "signed_ratio": u,
        "mesh_factors": asdict(f), "pinion": gears["pinion"], "wheel": gears["wheel"],
        "physical_validation": False}
