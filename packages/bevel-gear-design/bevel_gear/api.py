"""锥齿轮几何、啮合力和常幅疲劳评级的具名工程接口。

长度m、应力Pa、角度rad、扭矩N·m；转速rpm及持续时间h在字段名中
显式标注。毫米/MPa换算仅在底层评级边界完成。
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
import math

from . import factors
from .forces import calculate_mesh_forces
from .geometry import GeometryInput, GeometryResult, calculate_geometry
from .rating import SNCurve, evaluate_rating_mm


@dataclass(frozen=True)
class FatigueCurve:
    """调用方提供的三点S-N数据；应力Pa，循环数无量纲，不内置材料库。"""
    cycles: tuple[float, float, float]
    stresses_pa: tuple[float, float, float]

    def __post_init__(self):
        validated = SNCurve(tuple(self.cycles), tuple(float(s) / 1e6 for s in self.stresses_pa))
        object.__setattr__(self, "cycles", validated.cycles)
        object.__setattr__(self, "stresses_pa", tuple(s * 1e6 for s in validated.stresses_mpa))

    def to_mpa(self) -> SNCurve:
        return SNCurve(self.cycles, tuple(value / 1e6 for value in self.stresses_pa))


@dataclass(frozen=True)
class SteelMaterial:
    """单轮钢材与表面输入。曲线和极限须由使用者提供，不代表推荐牌号。

    category可选case_hardened、through_hardened、normalized、nitrided、
    induction_hardened。Ra/Rz为m；类别用于缺口、尺寸、磨合和表面系数。
    """
    category: str
    hardness_hb: float
    yield_strength_pa: float
    contact_limit_pa: float
    bending_limit_pa: float
    ra_m: float
    rz_m: float
    contact_curve: FatigueCurve
    bending_curve: FatigueCurve

    def __post_init__(self):
        object.__setattr__(self, "category", factors.material_category(self.category))
        for key in ("hardness_hb", "yield_strength_pa", "contact_limit_pa", "bending_limit_pa", "ra_m", "rz_m"):
            factors.positive(getattr(self, key), key)
        for key in ("contact_curve", "bending_curve"):
            curve = getattr(self, key)
            if isinstance(curve, dict):
                curve = FatigueCurve(**curve)
                object.__setattr__(self, key, curve)
            if not isinstance(curve, FatigueCurve):
                raise TypeError(f"{key}必须为FatigueCurve或具名曲线对象")


@dataclass(frozen=True)
class LoadCase:
    """单一常幅工况；零扭矩可独立计算几何/力，但不进入本疲劳评级。"""
    pinion_torque_n_m: float
    pinion_speed_rpm: float
    duration_h: float
    flank: str = "drive"
    spiral_hand: int = 1

    def __post_init__(self):
        factors.positive(self.pinion_torque_n_m, "小轮常幅扭矩N·m")
        factors.positive(self.pinion_speed_rpm, "小轮转速rpm")
        factors.nonnegative(self.duration_h, "持续时间h")
        if self.flank not in ("drive", "coast"):
            raise ValueError("工作齿面须为drive或coast")
        if isinstance(self.spiral_hand, bool) or self.spiral_hand not in (-1, 1):
            raise ValueError("旋向标志须为+1或-1，按局部力方向约定解释")


@dataclass(frozen=True)
class RatingSettings:
    """明确的B1/B2评级边界，刀尖半径是刀具输入而非齿根实际圆角。

    required_*是本次验收目标，安全系数=容量/实际应力；目标许用应力
    =容量/目标安全系数。默认弯曲目标为直齿1.5、弧齿1.3，仅为本模型
    默认策略，项目验收应显式给定目标。钢材默认弹性系数为189800√Pa。
    """
    pinion_cutter_tip_radius_m: float
    wheel_cutter_tip_radius_m: float
    viscosity_40_m2_s: float
    quality_grade: float = 8.0
    dynamic_method: str = "C"
    transverse_method: str = "B"
    wheel_cut_method: str = "generated"
    application_factor: float | None = None
    face_load_factor: float | None = None
    required_contact_safety: float = 1.0
    required_bending_safety: float | None = None
    fatigue_curve_safety: float = 1.0
    elasticity_factor_sqrt_pa: float = 189800.0
    pinion_protuberance_m: float = 0.0
    wheel_protuberance_m: float = 0.0
    driver_duty: str = "uniform"
    driven_duty: str = "uniform"
    assembly: str = "full_load_matched"
    support: str = "both_supported"
    root_method: str = "B1"
    static_load_distribution: bool = False
    b2_iteration_limit: int = 10000

    def __post_init__(self):
        for key in ("pinion_cutter_tip_radius_m", "wheel_cutter_tip_radius_m", "viscosity_40_m2_s",
                    "required_contact_safety", "fatigue_curve_safety", "elasticity_factor_sqrt_pa"):
            factors.positive(getattr(self, key), key)
        for key in ("application_factor", "face_load_factor", "required_bending_safety"):
            if getattr(self, key) is not None:
                factors.positive(getattr(self, key), key)
        for key in ("pinion_protuberance_m", "wheel_protuberance_m"):
            factors.nonnegative(getattr(self, key), key)
        if not math.isfinite(self.quality_grade):
            raise ValueError("精度等级必须为有限数")
        if self.dynamic_method not in ("B", "C") or self.transverse_method not in ("B", "C"):
            raise NotImplementedError("动载及横向分配只开放B/C方法")
        if self.wheel_cut_method not in ("generated", "form_cut"):
            raise NotImplementedError("大轮加工只开放展成或成形分支")
        if self.root_method not in ("B1","B2"):
            raise NotImplementedError("齿根评级方法须为B1或B2")
        if not isinstance(self.static_load_distribution,bool):
            raise ValueError("static_load_distribution须为布尔值")
        if self.static_load_distribution and self.root_method != "B2":
            raise ValueError("静载接触分配标志只适用于B2，不是独立静强度校核")
        if isinstance(self.b2_iteration_limit,bool) or not isinstance(self.b2_iteration_limit,int) or not 1 <= self.b2_iteration_limit <= 1000000:
            raise ValueError("B2截面最大迭代次数须为1至1000000的整数")


def _rating_inputs(geometry: GeometryResult, load: LoadCase, pinion: SteelMaterial,
                   wheel: SteelMaterial, settings: RatingSettings) -> dict:
    """单一换算边界：SI尺寸/应力→mm/MPa，粗糙度→µm，黏度→mm²/s。"""
    data = geometry.inputs
    values = {
        "z_1": data.pinion_teeth, "z_2": data.wheel_teeth,
        "sigma": data.shaft_angle_rad, "d_e2": data.wheel_outer_pitch_diameter_m * 1e3,
        "b_2": data.face_width_m * 1e3, "beta_m1": data.spiral_angle_rad, "beta_m2": data.spiral_angle_rad,
        "alpha_nD": data.pressure_angle_rad, "alpha_nC": data.pressure_angle_rad,
        "k_hap": data.addendum_coefficient, "k_hfp": data.dedendum_coefficient,
        "x_hm1": data.pinion_profile_shift_coefficient, "x_smn": data.tooth_thickness_modification_coefficient,
        "r_c0": (data.cutter_radius_m or 0.0) * 1e3,
        "geartype": 0 if geometry.diagnostics["tooth_trace"] == "straight" else (1 if data.spiral_angle_rad == 0 else 2),
        "methodB": 1 if settings.root_method == "B1" else 0, "isstaticload": int(settings.static_load_distribution),
        "T_1": load.pinion_torque_n_m, "n_1": load.pinion_speed_rpm, "duration": load.duration_h,
        "Kvmethod": {"B": 1, "C": 0}[settings.dynamic_method],
        "KHalphamethod": {"B": 1, "C": 0}[settings.transverse_method],
        "wheelcutmethod": 2 if settings.wheel_cut_method == "generated" else 0,
        "isograde": settings.quality_grade, "v40": settings.viscosity_40_m2_s * 1e6,
        "rho_a01": settings.pinion_cutter_tip_radius_m * 1e3,
        "rho_a02": settings.wheel_cutter_tip_radius_m * 1e3,
        "s_pr1": settings.pinion_protuberance_m * 1e3,
        "s_pr2": settings.wheel_protuberance_m * 1e3,
        "S_dutycucle": settings.fatigue_curve_safety,
        "elasticity_factor_sqrt_mpa": settings.elasticity_factor_sqrt_pa / 1e3,
        "driver_duty": settings.driver_duty, "driven_duty": settings.driven_duty,
        "assembly": settings.assembly, "support": settings.support,
    }
    if settings.application_factor is not None:
        values["factor_KA"] = settings.application_factor
    if settings.face_load_factor is not None:
        values["factor_KHbeta"] = settings.face_load_factor
    surface = {"case_hardened": 1, "induction_hardened": 2, "nitrided": 3,
               "through_hardened": 2, "normalized": 0}
    for index, material, prefix in ((1, pinion, "pinion"), (2, wheel, "wheel")):
        values.update({f"material_{prefix}": material.category, f"HB_{index}": material.hardness_hb,
                       f"sigma_S{index}": material.yield_strength_pa / 1e6,
                       f"sigma_Hlim{index}": material.contact_limit_pa / 1e6,
                       f"sigma_Flim{index}": material.bending_limit_pa / 1e6,
                       f"R_a{index}": material.ra_m * 1e6, f"R_z{index}": material.rz_m * 1e6,
                       f"surfacetreatment{index}": surface[material.category],
                       f"contact_sn_{prefix}": material.contact_curve.to_mpa(),
                       f"bending_sn_{prefix}": material.bending_curve.to_mpa()})
    return values


def _life_record(raw: dict, index: int, kind: str) -> dict:
    life = raw[f"life{kind}{index}_cycles"]
    state = "out_of_curve" if life is None else ("endurance_assumption" if math.isinf(life) else "finite")
    return {"cycles": raw[f"cycles{index}"], "life_cycles": life,
            "damage": raw[f"damage{kind}{index}"], "state": state}


def analyze_pair(geometry: GeometryInput, load: LoadCase, material_pinion: SteelMaterial,
                 material_wheel: SteelMaterial, rating: RatingSettings) -> dict:
    """参数化求一对零偏置外啮合锥齿轮；输入对象保持不变。"""
    for obj, expected in ((geometry, GeometryInput), (load, LoadCase),
                          (material_pinion, SteelMaterial), (material_wheel, SteelMaterial),
                          (rating, RatingSettings)):
        if not isinstance(obj, expected):
            raise TypeError(f"输入必须为{expected.__name__}")
    if geometry.backlash_m != 0:
        raise NotImplementedError("几何入口可计算非零侧隙；当前集成评级尚未验证侧隙减薄对两轮根截面的影响，要求backlash_m=0")
    result_geometry = calculate_geometry(geometry)
    result_forces = calculate_mesh_forces(result_geometry, load.pinion_torque_n_m,
                                         flank=load.flank, spiral_hand=load.spiral_hand)
    inputs_mm = _rating_inputs(result_geometry, load, material_pinion, material_wheel, rating)
    if rating.root_method == "B2":
        from .rating_b2 import evaluate_b2_rating_mm
        raw = evaluate_b2_rating_mm(inputs_mm,iteration_limit=rating.b2_iteration_limit)
    else:
        raw = evaluate_rating_mm(inputs_mm)
    contact_target = rating.required_contact_safety
    bending_target = rating.required_bending_safety
    if bending_target is None:
        bending_target = 1.5 if result_geometry.diagnostics["tooth_trace"] == "straight" else 1.3
    output = {"contact_stress_pa": raw["sigma_H"] * 1e6,
              "required_contact_safety": contact_target, "required_bending_safety": bending_target,
              "coefficients": {key: value for key, value in raw.items()
                               if key != "Z_E" and (rating.root_method == "B2" or key not in ("Y_A1", "Y_A2", "Y_J1", "Y_J2"))
                               and (key.startswith(("K_", "Y_", "Z_")) or key == "ZLZvZR")},
              "elasticity_factor_sqrt_pa": raw["Z_E"] * 1e3,
              "life": {},
              "rating_geometry": {"transverse_contact_ratio": raw["geometry"]["epsilon_valpha"],
                                  "overlap_contact_ratio": raw["geometry"]["epsilon_vbeta"],
                                  "resultant_contact_ratio": raw["geometry"]["epsilon_vgamma"]},
              "diagnostics": raw["diagnostics"]}
    output["root_method"] = rating.root_method
    if rating.root_method == "B2":
        output["static_load_distribution"] = rating.static_load_distribution
        output["b2_load_sharing_factor"] = raw["b2_load_sharing_factor"]
        output["b2_contact_ratio_factor"] = raw["b2_contact_ratio_factor"]
        output["root_sections"] = {}
        for index,member in ((1,"pinion"),(2,"wheel")):
            section = raw[f"root_section{index}"]
            output["root_sections"][member] = {key.replace("_mm","_m"):value/1000 if key.endswith("_mm") else value for key,value in section.items()}
    for index, prefix in ((1, "pinion"), (2, "wheel")):
        output.update({f"{prefix}_bending_stress_pa": raw[f"sigma_F{index}"] * 1e6,
                       f"{prefix}_contact_capacity_pa": raw[f"sigma_HP{index}"] * 1e6,
                       f"{prefix}_bending_capacity_pa": raw[f"sigma_FP{index}"] * 1e6,
                       f"{prefix}_allowable_contact_stress_pa": raw[f"sigma_HP{index}"] * 1e6 / contact_target,
                       f"{prefix}_allowable_bending_stress_pa": raw[f"sigma_FP{index}"] * 1e6 / bending_target,
                       f"{prefix}_contact_safety": raw[f"S_H{index}"],
                       f"{prefix}_bending_safety": raw[f"S_F{index}"]})
        output["life"][prefix] = {"contact": _life_record(raw, index, "H"),
                                  "bending": _life_record(raw, index, "F")}
    output["meets_targets"] = all(output[f"{role}_contact_safety"] >= contact_target
                                   and output[f"{role}_bending_safety"] >= bending_target
                                   for role in ("pinion", "wheel"))
    diagnostics = list(raw["diagnostics"])
    diagnostics.append("几何重合度和评级当量重合度采用不同评价模型，分别保留；均不表示齿面LTCA载荷分布")
    return {"inputs": {"geometry": asdict(geometry), "load": asdict(load),
                       "materials": {"pinion": asdict(material_pinion), "wheel": asdict(material_wheel)},
                       "rating": asdict(rating)},
            "geometry": asdict(result_geometry), "forces": asdict(result_forces),
            "rating": output, "diagnostics": diagnostics,
            "model": "zero_offset_bevel_constant_amplitude_"+rating.root_method.lower()}


def analyze_design(config: dict) -> dict:
    """读取严格具名的JSON兼容配置；不填补缺失材料或疲劳曲线。"""
    if not isinstance(config, dict):
        raise TypeError("设计配置必须为对象")
    allowed = {"geometry", "load", "materials", "rating", "name", "description", "standard"}
    if set(config) - allowed:
        raise ValueError(f"未知顶层配置字段：{sorted(set(config) - allowed)}")
    standard = config.get("standard","ISO10300")
    if standard == "AGMA2003C10":
        return analyze_agma_design(config)
    if standard != "ISO10300":
        raise NotImplementedError("评级路线须为ISO10300或AGMA2003C10")
    materials = config["materials"]
    result = analyze_pair(GeometryInput(**config["geometry"]), LoadCase(**config["load"]),
                          SteelMaterial(**materials["pinion"]), SteelMaterial(**materials["wheel"]),
                          RatingSettings(**config["rating"]))
    result["name"] = config.get("name", "锥齿轮设计")
    result["description"] = config.get("description", "")
    return result


def analyze_agma_pair(geometry: GeometryInput, load: LoadCase, material_pinion: SteelMaterial,
                      material_wheel: SteelMaterial, settings, *, required_contact_safety=1.,
                      required_bending_safety=None) -> dict:
    """AGMA完整I/J路线的统一SI门面，输出结构可接现有图表和谱分析。

    算法、插值和系数按所选AGMA路线保持独立；额定功率逐轮输出，不合并
    为未验证的整机容量。最高点以上的一次失效假设在life.assumption标出。
    """
    from .agma_rating import AGMARatingSettings,evaluate_agma_rating_si
    if not isinstance(geometry,GeometryInput) or not isinstance(load,LoadCase):
        raise TypeError("AGMA门面需要GeometryInput和LoadCase")
    if not all(isinstance(m,SteelMaterial) for m in (material_pinion,material_wheel)) or not isinstance(settings,AGMARatingSettings):
        raise TypeError("AGMA门面需要SteelMaterial和AGMARatingSettings")
    g=calculate_geometry(geometry)
    f=calculate_mesh_forces(g,load.pinion_torque_n_m,flank=load.flank,spiral_hand=load.spiral_hand)
    si=evaluate_agma_rating_si(g,load,material_pinion,material_wheel,settings)
    raw=si["raw_mm_mpa"]
    contact_target=factors.positive(required_contact_safety,"接触目标安全系数")
    bend_target=required_bending_safety
    if bend_target is None:
        bend_target=1.5 if g.diagnostics["tooth_trace"]=="straight" else 1.3
    bend_target=factors.positive(bend_target,"弯曲目标安全系数")
    output={"standard":"AGMA2003C10","root_method":"AGMA2003C10","contact_stress_pa":si["contact_stress_pa"],
            "required_contact_safety":contact_target,"required_bending_safety":bend_target,
            "coefficients":si["coefficients"],"elasticity_factor_sqrt_pa":raw["C_p"]*1000,
            "rated_contact_power_w":si["rated_contact_power_w"],"rated_bending_power_w":si["rated_bending_power_w"],
            "damage_interpolation":settings.damage_interpolation,"hardness_ratio_policy":settings.hardness_ratio_policy,
            "life":{},"diagnostics":list(si["diagnostics"])}
    for index,member,material,letter in ((0,"pinion",material_pinion,"P"),(1,"wheel",material_wheel,"G")):
        output.update({f"{member}_bending_stress_pa":si["bending_stress_pa"][index],
                       f"{member}_contact_capacity_pa":si["contact_capacity_pa"][index],
                       f"{member}_bending_capacity_pa":si["bending_capacity_pa"][index],
                       f"{member}_allowable_contact_stress_pa":si["contact_capacity_pa"][index]/contact_target,
                       f"{member}_allowable_bending_stress_pa":si["bending_capacity_pa"][index]/bend_target,
                       f"{member}_contact_safety":si["contact_safety"][index],f"{member}_bending_safety":si["bending_safety"][index]})
        life={}
        for mode,symbol,curve in (("contact","H",material.contact_curve),("bending","F",material.bending_curve)):
            nfail=si["life_cycles"][symbol][index]
            state="out_of_curve" if nfail is None else ("endurance_assumption" if math.isinf(nfail) else "finite")
            equivalent=(raw["s_c"]/raw["C_H"]*raw["K_T"]*raw["C_R"] if symbol=="H" else raw[f"s_t{letter}"]*raw["K_T"]*raw["K_R"])
            above=equivalent*1e6*settings.curve_safety_factor>curve.stresses_pa[0]
            life[mode]={"cycles":raw[f"cycles{index+1}"],"life_cycles":nfail,"damage":si["damage"][symbol][index],
                        "state":state,"assumption":"one_cycle_above_curve" if above and settings.damage_interpolation=="source_semilog" else None}
        output["life"][member]=life
    output["meets_targets"]=all(output[f"{member}_contact_safety"]>=contact_target and output[f"{member}_bending_safety"]>=bend_target for member in ("pinion","wheel"))
    rating_input={**asdict(settings),"required_contact_safety":contact_target,"required_bending_safety":bend_target}
    return {"inputs":{"standard":"AGMA2003C10","geometry":asdict(geometry),"load":asdict(load),
                      "materials":{"pinion":asdict(material_pinion),"wheel":asdict(material_wheel)},"rating":rating_input},
            "geometry":asdict(g),"forces":asdict(f),"rating":output,"diagnostics":list(si["diagnostics"]),
            "model":"zero_offset_bevel_agma2003_constant_amplitude"}


def analyze_agma_design(config: dict) -> dict:
    """读取AGMA具名配置，沿用几何/材料/载荷契约并单独解析评级设置。"""
    from .agma_rating import AGMARatingSettings
    allowed={"standard","geometry","load","materials","rating","name","description"}
    if not isinstance(config,dict) or set(config)-allowed:
        raise ValueError("AGMA配置包含未知字段或不是字典")
    if config.get("standard","AGMA2003C10")!="AGMA2003C10":
        raise ValueError("AGMA入口的standard必须为AGMA2003C10")
    values=dict(config["rating"])
    contact=values.pop("required_contact_safety",1.)
    bending=values.pop("required_bending_safety",None)
    materials=config["materials"]
    result=analyze_agma_pair(GeometryInput(**config["geometry"]),LoadCase(**config["load"]),
             SteelMaterial(**materials["pinion"]),SteelMaterial(**materials["wheel"]),AGMARatingSettings(**values),
             required_contact_safety=contact,required_bending_safety=bending)
    result["name"]=config.get("name","AGMA锥齿轮设计")
    result["description"]=config.get("description","")
    return result


def sweep_parameter(config: dict, parameter: str, values) -> list[dict]:
    """逐点完整重算具名参数；失败点保留诊断，不用插值或旧应力缩放替代。"""
    parts = parameter.split(".")
    if len(parts) != 2 or parts[0] not in ("geometry", "load", "rating"):
        raise ValueError("扫描参数应为geometry/load/rating中的一个具名字段")
    if parts[1] not in config[parts[0]]:
        raise ValueError("扫描字段必须显式存在于配置")
    rows = []
    for value in values:
        point = deepcopy(config)
        point[parts[0]][parts[1]] = value
        try:
            result = analyze_design(point)
            rows.append({"parameter": parameter, "value": value, "status": "ok", "result": result})
        except (ValueError, NotImplementedError) as error:
            rows.append({"parameter": parameter, "value": value, "status": "invalid", "diagnostic": str(error)})
    return rows
