"""零偏置锥齿轮、钢材、常幅载荷的 B1 齿根与点蚀评级。

公开顶层接口应使用 SI；本文件的 evaluate_rating_mm 是明确采用 mm、
MPa、N·m、rpm、h、rad、µm 的底层数值入口。ISO 10300 仅表示选用的
计算路线，不表示覆盖所有版本、材料、制造方法或完整标准认证。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from . import factors


@dataclass(frozen=True)
class SNCurve:
    """三点常幅 S-N 曲线：N 严增、应力 MPa 严减，末点为耐久阈值。

    点间双对数插值；低于第一循环数时强度取第一点。低于末点应力时
    本模型取无限寿命。超出最高应力的失效寿命没有数据，不推算单次
    断裂或虚构曲线延长线。输入曲线不在计算中修改。
    """
    cycles: tuple[float, float, float]
    stresses_mpa: tuple[float, float, float]

    def __post_init__(self):
        if len(self.cycles) != 3 or len(self.stresses_mpa) != 3:
            raise ValueError("S-N 曲线必须恰好有三个循环数/应力点")
        cycles = tuple(factors.positive(v, "S-N 循环数") for v in self.cycles)
        stress = tuple(factors.positive(v, "S-N 应力 MPa") for v in self.stresses_mpa)
        if not cycles[0] < cycles[1] < cycles[2] or not stress[0] > stress[1] > stress[2]:
            raise ValueError("S-N 循环数须严格增加，应力须严格减少")
        object.__setattr__(self, "cycles", cycles)
        object.__setattr__(self, "stresses_mpa", stress)

    def stress_at_cycles(self, cycles):
        cycles = factors.nonnegative(cycles, "累计循环数")
        if cycles <= self.cycles[0]:
            return self.stresses_mpa[0]
        for index in range(2):
            if cycles < self.cycles[index + 1]:
                fraction = math.log(cycles / self.cycles[index]) / math.log(self.cycles[index + 1] / self.cycles[index])
                return self.stresses_mpa[index] * (self.stresses_mpa[index + 1] / self.stresses_mpa[index])**fraction
        return self.stresses_mpa[-1]

    def cycles_at_stress(self, stress_mpa):
        stress = factors.nonnegative(stress_mpa, "等效应力 MPa")
        if stress <= self.stresses_mpa[-1]:
            return math.inf
        if stress > self.stresses_mpa[0]:
            return None
        for index in range(2):
            if stress > self.stresses_mpa[index + 1]:
                fraction = math.log(stress / self.stresses_mpa[index]) / math.log(self.stresses_mpa[index + 1] / self.stresses_mpa[index])
                return self.cycles[index] * (self.cycles[index + 1] / self.cycles[index])**fraction
        return self.cycles[-1]

    def damage(self, cycles, equivalent_stress_mpa, *, curve_safety_factor=1.0):
        safety = factors.positive(curve_safety_factor, "寿命曲线安全系数")
        life = self.cycles_at_stress(equivalent_stress_mpa * safety)
        if life is None:
            return None, None
        return cycles / life, life


def _curve(values, prefix, kind):
    descriptive = {("pinion", "H"): "contact_sn_pinion", ("wheel", "H"): "contact_sn_wheel",
                   ("pinion", "F"): "bending_sn_pinion", ("wheel", "F"): "bending_sn_wheel"}
    key = descriptive[prefix, kind]
    if key in values:
        curve = values[key]
        if isinstance(curve, SNCurve):
            return curve
        return SNCurve(tuple(row[0] for row in curve), tuple(row[1] for row in curve))
    required = [f"{prefix}N{kind}{i}" for i in range(1, 4)] + [f"{prefix}sigma{kind}{i}" for i in range(1, 4)]
    if any(name not in values for name in required):
        raise ValueError(f"缺少完整三点 S-N 曲线：{key}")
    return SNCurve(tuple(float(values[f"{prefix}N{kind}{i}"]) for i in range(1, 4)),
                   tuple(float(values[f"{prefix}sigma{kind}{i}"]) for i in range(1, 4)))


def _virtual_geometry(values):
    """从外端尺寸建立评级的当量圆柱几何，所有直径/路径长度为 mm。"""
    z1, z2 = values["z_1"], values["z_2"]
    if any(not isinstance(z, int) or isinstance(z, bool) or z < 3 for z in (z1, z2)):
        raise ValueError("两轮齿数必须为不小于 3 的整数")
    ratio = z2 / z1
    shaft = float(values["sigma"])
    if not 0 < shaft < math.pi:
        raise ValueError("轴交角必须在 (0,π) 内")
    delta1 = math.atan2(math.sin(shaft), math.cos(shaft) + ratio)
    delta2 = shaft - delta1
    if any(not 1e-9 < d < math.pi / 2 - 1e-9 for d in (delta1, delta2)):
        raise NotImplementedError("评级仅支持两节锥角均为锐角的零偏置外啮合锥齿轮")
    outer_wheel = factors.positive(values["d_e2"], "大轮外端节圆直径 mm")
    width = factors.positive(values["b_2"], "齿宽 mm")
    outer_cone = outer_wheel / (2 * math.sin(delta2))
    if width >= outer_cone:
        raise ValueError("齿宽必须小于外端锥距")
    mean_cone = outer_cone - width / 2
    outer_pinion = outer_wheel * math.sin(delta1) / math.sin(delta2)
    mean_pinion = outer_pinion - width * math.sin(delta1)
    mean_wheel = outer_wheel - width * math.sin(delta2)
    beta = float(values["beta_m2"])
    if not 0 <= beta < math.pi / 2:
        raise ValueError("评级中的平均螺旋角使用 [0,π/2) 内的大小值")
    normal = float(values["alpha_nD"])
    if not 0 < normal < math.pi / 3:
        raise ValueError("法向压力角必须在 (0,π/3) 内")
    outer_module = outer_wheel / z2
    transverse_module = mean_cone / outer_cone * outer_module
    normal_module = transverse_module * math.cos(beta)
    addendum_coeff = factors.positive(values["k_hap"], "齿顶高系数")
    dedendum_coeff = factors.positive(values["k_hfp"], "齿根高系数")
    shift = float(values["x_hm1"])
    thickness_shift = float(values["x_smn"])
    addendum_p = normal_module * (addendum_coeff + shift)
    addendum_w = normal_module * (addendum_coeff - shift)
    if min(addendum_p, addendum_w) <= 0:
        raise ValueError("变位后两轮中点齿顶高必须均为正")
    virtual_teeth_p = z1 / math.cos(delta1)
    virtual_teeth_w = z2 / math.cos(delta2)
    virtual_ratio = virtual_teeth_w / virtual_teeth_p
    pitch_p = mean_pinion / math.cos(delta1)
    pitch_w = mean_wheel / math.cos(delta2)
    virtual_center = (pitch_p + pitch_w) / 2
    tip_p, tip_w = pitch_p + 2 * addendum_p, pitch_w + 2 * addendum_w
    transverse_pressure = math.atan(math.tan(normal) / math.cos(beta))
    base_p = pitch_p * math.cos(transverse_pressure)
    base_w = pitch_w * math.cos(transverse_pressure)
    base_helix = math.asin(math.sin(beta) * math.cos(normal))
    base_pitch = transverse_module * math.pi * math.cos(transverse_pressure)
    path = (math.sqrt(tip_p**2 - base_p**2) + math.sqrt(tip_w**2 - base_w**2)) / 2 - virtual_center * math.sin(transverse_pressure)
    transverse_ratio = factors.positive(path / base_pitch, "端面重合度")
    overlap_ratio = width * math.sin(beta) / (normal_module * math.pi)
    if transverse_ratio <= 1:
        raise ValueError("本评级路线要求端面重合度大于 1")
    normal_teeth_p = virtual_teeth_p / (math.cos(base_helix)**2 * math.cos(beta))
    normal_teeth_w = virtual_ratio * normal_teeth_p
    normal_pitch_p = pitch_p / math.cos(base_helix)**2
    normal_pitch_w = virtual_ratio * normal_pitch_p
    return {
        "ratio": ratio, "delta_1": delta1, "delta_2": delta2, "d_e1": outer_pinion,
        "d_m1": mean_pinion, "d_m2": mean_wheel, "R_e": outer_cone, "R_m": mean_cone,
        "m_et": outer_module, "m_mt": transverse_module, "m_mn": normal_module,
        "d_v1": pitch_p, "d_v2": pitch_w, "u_v": virtual_ratio, "a_v": virtual_center,
        "d_va1": tip_p, "d_va2": tip_w, "d_vb1": base_p, "d_vb2": base_w,
        "alpha_vt": transverse_pressure, "beta_vb": base_helix, "p_et": base_pitch,
        "g_valpha": path, "epsilon_valpha": transverse_ratio, "epsilon_vbeta": overlap_ratio,
        "epsilon_vgamma": math.hypot(transverse_ratio, overlap_ratio),
        "z_vn1": normal_teeth_p, "z_vn2": normal_teeth_w,
        "d_van1": normal_pitch_p + tip_p - pitch_p,
        "d_van2": normal_pitch_w + tip_w - pitch_w,
        "d_vbn1": normal_pitch_p * math.cos(normal),
        "d_vbn2": normal_pitch_w * math.cos(normal),
        "thickness_shift": thickness_shift, "profile_shift": shift,
        "normal_pressure_angle_rad": normal, "spiral_angle_rad": beta,
        "face_width_mm": width, "cutter_height_coefficient": dedendum_coeff,
    }


def evaluate_rating_mm(inputs):
    """计算明确单位字典的常幅钢材 B1 全链，返回应力/安全系数/寿命诊断。

    长度输入 d_e2/b_2/rho_a0*/s_pr*/r_c0 为 mm，角度 sigma/alpha_n*/
    beta_m* 为 rad，T_1 为 N·m，n_1 为 rpm，duration 为 h。Ra/Rz 与
    返回 f_pt/y_alpha 为 µm，sigma* 应力输入输出为 MPa。标准工程顶层
    应使用 SI 类型化门面完成换算，不能将此字典直接作为 SI 接口。

    必须提供四条三点 S-N 曲线或其 N/sigma 字段；本模块没有材料数据库。
    KA/KHβ 是外部工况/安装边界，未指定时按明确的均匀载荷及满载调整
    默认策略计算；其余已覆盖系数由几何、制造、载荷与材料自动计算。
    不支持整体预设系数模式，不接收未定义材料/B2/偏置分支。
    """
    values = dict(inputs)
    required = ("z_1", "z_2", "sigma", "d_e2", "b_2", "beta_m1", "beta_m2", "alpha_nD",
                "k_hap", "k_hfp", "x_hm1", "x_smn", "T_1", "n_1", "duration", "r_c0",
                "rho_a01", "rho_a02", "s_pr1", "s_pr2", "sigma_Hlim1", "sigma_Hlim2",
                "sigma_Flim1", "sigma_Flim2", "sigma_S1", "sigma_S2", "HB_1", "HB_2",
                "R_a1", "R_a2", "R_z1", "R_z2", "v40", "isograde", "surfacetreatment1", "surfacetreatment2")
    missing = [name for name in required if name not in values]
    if missing:
        raise ValueError("缺少 B1 评级输入字段：" + ", ".join(missing))
    if values.get("userdefined", False):
        raise NotImplementedError("不支持绕过自动系数链的整体预设系数模式")
    if values.get("methodB", 1) != 1:
        raise NotImplementedError("当前只实现 B1 齿根计算，B2 尚未开放")
    if values.get("isstaticload", 0):
        raise NotImplementedError("此入口为常幅疲劳载荷，静强度单独评级尚未开放")
    torque = factors.positive(values["T_1"], "常幅扭矩 N·m；零载不属于本疲劳评级入口")
    speed = factors.positive(values["n_1"], "转速 rpm")
    hours = factors.nonnegative(values["duration"], "持续时间 h")
    if not 0 <= float(values["beta_m1"]) < math.pi / 2:
        raise ValueError("小轮平均螺旋角须为 [0,π/2) 内的大小值")
    g = _virtual_geometry(values)
    gear_type = values.get("geartype", 0 if g["spiral_angle_rad"] == 0 else 2)
    if gear_type not in (0, 1, 2):
        raise NotImplementedError("B1支持直齿=0、零度弧齿=1和螺旋锥齿=2")
    straight = gear_type == 0
    if (gear_type < 2) != (g["spiral_angle_rad"] == 0) or (gear_type < 2 and float(values["beta_m1"]) != 0):
        raise ValueError("直齿/零度弧齿要求两轮平均螺旋角为零，螺旋类别要求正平均螺旋角")
    material_p = factors.material_category(values.get("material_pinion", values.get("abbreviation1", "case_hardened")))
    material_w = factors.material_category(values.get("material_wheel", values.get("abbreviation2", "case_hardened")))
    limit_h = [factors.positive(values[f"sigma_Hlim{i}"], "接触疲劳极限 MPa") for i in (1, 2)]
    limit_f = [factors.positive(values[f"sigma_Flim{i}"], "弯曲疲劳极限 MPa") for i in (1, 2)]
    curves_h = [_curve(values, prefix, "H") for prefix in ("pinion", "wheel")]
    curves_f = [_curve(values, prefix, "F") for prefix in ("pinion", "wheel")]
    width, module, beta = g["face_width_mm"], g["m_mn"], g["spiral_angle_rad"]
    tangential = 2000 * torque / g["d_m1"]
    velocity = g["d_m1"] * speed / 19098
    effective_width = 0.85 * width
    out = {"geometry": g, "F_mt1": tangential,
           "F_mt2": tangential * math.cos(float(values["beta_m2"])) / math.cos(float(values["beta_m1"])),
           "delta_1": g["delta_1"], "delta_2": g["delta_2"], "d_m1": g["d_m1"],
           "diagnostics": [], "units": {"stress": "MPa", "geometry": "mm", "roughness": "µm"}}
    for suffix, sign in (("D", 1), ("C", -1)):
        angle = float(values.get(f"alpha_n{suffix}", g["normal_pressure_angle_rad"]))
        if not 0 < angle < math.pi / 3:
            raise ValueError("两齿面的法向压力角须在 (0,π/3) 内")
        for index, wheel_sign in ((1, sign), (2, -sign)):
            beta_i, delta = float(values[f"beta_m{index}"]), g[f"delta_{index}"]
            force_i = out[f"F_mt{index}"]
            out[f"F_ax{index}{suffix}"] = force_i * (math.tan(angle) * math.sin(delta) / math.cos(beta_i) + wheel_sign * math.tan(beta_i) * math.cos(delta))
            out[f"F_rad{index}{suffix}"] = force_i * (math.tan(angle) * math.cos(delta) / math.cos(beta_i) - wheel_sign * math.tan(beta_i) * math.sin(delta))

    ratio_a, ratio_b, ratio_total = g["epsilon_valpha"], g["epsilon_vbeta"], g["epsilon_vgamma"]
    line_h = factors.contact_line_geometry(g["g_valpha"], g["beta_vb"], width, g["p_et"], ratio_a, ratio_b)
    line_f = factors.contact_line_geometry(g["g_valpha"], g["beta_vb"], width, g["p_et"], ratio_a, ratio_b, for_bending=True)
    out["Z_H"] = 2 * math.sqrt(math.cos(g["beta_vb"]) / math.sin(2 * g["alpha_vt"]))
    evaluation_p = 2 + (ratio_a - 2) * min(ratio_b, 1)
    evaluation_w = 2 * ratio_a - 2 + (2 - ratio_a) * min(ratio_b, 1)
    if ratio_b > 1:
        evaluation_p = evaluation_w = ratio_a
    curvature_a = math.sqrt(g["d_va1"]**2 / g["d_vb1"]**2 - 1) - evaluation_p * math.pi / (values["z_1"] / math.cos(g["delta_1"]))
    curvature_b = math.sqrt(g["d_va2"]**2 / g["d_vb2"]**2 - 1) - evaluation_w * math.pi / (values["z_2"] / math.cos(g["delta_2"]))
    out["Z_MB"] = math.tan(g["alpha_vt"]) / math.sqrt(factors.positive(curvature_a, "小轮点蚀评价曲率") * factors.positive(curvature_b, "大轮点蚀评价曲率"))
    out.update(Z_E=factors.positive(values.get("elasticity_factor_sqrt_mpa", 189.8), "钢材弹性系数 √MPa"),
               Z_beta=math.sqrt(math.cos(beta)), Z_K=0.8, Z_X1=1.0, Z_X2=1.0)
    out["Z_LS"] = factors.load_sharing_factor(ratio_total, ratio_b, g["p_et"], g["beta_vb"], line_h)
    out["ZLZvZR"] = factors.lubricant_film_factor(*limit_h, values["v40"], velocity,
        values["R_a1"], values["R_a2"], g["a_v"], g["alpha_vt"], g["beta_vb"], g["u_v"])
    out["Z_W"] = factors.hardness_ratio_factor(values["HB_1"], values["HB_2"])
    out["K_A"] = factors.positive(values.get("factor_KA", factors.application_factor(
        values.get("driver_duty", "uniform"), values.get("driven_duty", "uniform"))), "应用系数 KA")
    out["K_Hbeta"] = factors.positive(values.get("factor_KHbeta", factors.face_load_factor(
        0.85, values.get("assembly", "full_load_matched"), values.get("support", "both_supported"))), "齿向系数 KHβ")
    out["f_pt"] = factors.single_pitch_deviation(module, g["d_m2"], values["isograde"])
    out["y_alpha"] = factors.running_in_allowance(material_p, limit_h[0], material_w, limit_h[1], velocity, out["f_pt"])
    dynamic_method = values.get("Kvmethod", 0)
    if dynamic_method == 1:
        out["K_v"], out["dynamic_diagnostics"] = factors.dynamic_factor_b(speed, g["d_m1"], values["z_1"],
            g["normal_pressure_angle_rad"], g["ratio"], tangential, width, effective_width, out["K_A"],
            ratio_total, out["f_pt"], out["y_alpha"])
    elif dynamic_method == 0:
        out["K_v"] = factors.dynamic_factor_c(module, g["d_m1"], g["d_e1"], values["isograde"], velocity)
    else:
        raise NotImplementedError("动载系数只支持方法 B=1 或 C=0")
    out["dynamicFactorErrorCode"] = 0
    # 零度弧齿仍走非直齿的横向分配/目标分支；这里只取β→0的齿向修正极限。
    out["K_Fbeta"] = out["K_Hbeta"] if beta == 0 else factors.bending_face_load_factor(straight, beta, g["R_m"], values["r_c0"], out["K_Hbeta"])
    out["Y_epsilon"] = factors.contact_ratio_factor(ratio_a, ratio_b)
    transverse_method = values.get("KHalphamethod", 1)
    if transverse_method == 1:
        out["K_Halpha"], out["K_Falpha"] = factors.transverse_load_factors_b(width, ratio_a, ratio_b,
            ratio_total, tangential, out["K_A"], out["K_v"], out["K_Hbeta"], out["Z_LS"],
            out["Y_epsilon"], out["y_alpha"], out["f_pt"])
    elif transverse_method == 0:
        factor = max(factors.transverse_load_factor_c(values[f"surfacetreatment{i}"], straight,
            values["isograde"], tangential, out["K_A"], effective_width, ratio_a) for i in (1, 2))
        out["K_Halpha"] = out["K_Falpha"] = factor
    else:
        raise NotImplementedError("横向分配只支持方法 B=1 或 C=0")
    out["Y_K"] = factors.bevel_root_factor(width, line_f.projected_middle_length_mm)
    out["Y_LS"] = factors.load_sharing_factor(ratio_total, ratio_b, g["p_et"], g["beta_vb"], line_f)**2
    wheel_cut = values.get("wheelcutmethod", 2)
    if wheel_cut not in (0, 1, 2):
        raise NotImplementedError("大轮加工方式只支持成形=0/1 或展成=2")
    for index, category in ((1, material_p), (2, material_w)):
        direction = 1 if index == 1 else -1
        shape, correction, section = factors.root_form_factors(g["cutter_height_coefficient"], module,
            direction * g["thickness_shift"], direction * g["profile_shift"], g["normal_pressure_angle_rad"],
            values[f"rho_a0{index}"], values[f"s_pr{index}"], g[f"z_vn{index}"],
            g[f"d_van{index}"], g[f"d_vbn{index}"], generated=index == 1 or wheel_cut == 2)
        out[f"Y_Fa{index}"], out[f"Y_Sa{index}"] = shape, correction
        out[f"Y_deltarelT{index}"] = factors.relative_sensitivity_factor(category, values[f"sigma_S{index}"], section)
        out[f"Y_RrelT{index}"] = factors.relative_surface_factor(values[f"surfacetreatment{index}"], values[f"R_z{index}"])
        out[f"Y_X{index}"] = factors.size_factor(module, category)
        out[f"sigma_F0{index}"] = tangential / (width * module) * shape * correction * out["Y_epsilon"] * out["Y_K"] * out["Y_LS"]
        out[f"sigma_F{index}"] = out[f"sigma_F0{index}"] * out["K_A"] * out["K_v"] * out["K_Fbeta"] * out["K_Falpha"]
        out[f"root_section{index}"] = {"thickness_mm": section.thickness_mm, "fillet_radius_mm": section.fillet_radius_mm,
                                     "iterations": section.iterations}
        out[f"Y_A{index}"] = out[f"Y_J{index}"] = 0.0
    nominal_h = math.sqrt(tangential / (g["d_v1"] * factors.positive(line_h.middle_length_mm, "点蚀中心接触线 mm")) * (g["u_v"] + 1) / g["u_v"])
    out["sigma_H"] = nominal_h * math.prod(out[key] for key in ("Z_MB", "Z_H", "Z_E", "Z_LS", "Z_beta", "Z_K")) * math.sqrt(out["K_A"] * out["K_v"] * out["K_Halpha"] * out["K_Hbeta"])
    out["S_Hmin"], out["S_Fmin"] = 1.0, 1.5 if straight else 1.3
    cycles = (speed * 60 * hours, speed * 60 * hours / g["ratio"])
    curve_safety = values.get("S_dutycucle", 1.0)
    for index in (1, 2):
        row = index - 1
        out[f"Z_NT{index}"] = curves_h[row].stress_at_cycles(cycles[row]) / limit_h[row]
        out[f"Y_NT{index}"] = curves_f[row].stress_at_cycles(cycles[row]) / limit_f[row]
        out[f"sigma_HP{index}"] = limit_h[row] * out[f"Z_NT{index}"] * out[f"Z_X{index}"] * out["ZLZvZR"] * out["Z_W"]
        root_multiplier = 2 * out[f"Y_deltarelT{index}"] * out[f"Y_RrelT{index}"] * out[f"Y_X{index}"]
        out[f"sigma_FP{index}"] = limit_f[row] * out[f"Y_NT{index}"] * root_multiplier
        out[f"S_H{index}"] = out[f"sigma_HP{index}"] / out["sigma_H"]
        out[f"S_F{index}"] = out[f"sigma_FP{index}"] / out[f"sigma_F{index}"]
        contact_equivalent = out["sigma_H"] / (out[f"Z_X{index}"] * out["ZLZvZR"] * out["Z_W"])
        bending_equivalent = out[f"sigma_F{index}"] / root_multiplier
        out[f"sigmaF{index}"] = bending_equivalent
        out[f"damageH{index}"], out[f"lifeH{index}_cycles"] = curves_h[row].damage(cycles[row], contact_equivalent, curve_safety_factor=curve_safety)
        out[f"damageF{index}"], out[f"lifeF{index}_cycles"] = curves_f[row].damage(cycles[row], bending_equivalent, curve_safety_factor=curve_safety)
        out[f"cycles{index}"] = cycles[row]
        if out[f"damageH{index}"] is None or out[f"damageF{index}"] is None:
            out["diagnostics"].append(f"第 {index} 轮等效应力高于给定 S-N 最高应力，寿命/损伤不作外推")
    out["contact_line_pitting_mm"] = (line_h.tip_length_mm, line_h.middle_length_mm, line_h.root_length_mm)
    out["contact_line_bending_mm"] = (line_f.tip_length_mm, line_f.middle_length_mm, line_f.root_length_mm)
    out["diagnostics"].append("钢材弹性系数默认 189.8 √MPa；KA/KHβ 为工况与安装边界；曲线末点按耐久阈值处理")
    if gear_type == 1:
        out["diagnostics"].append("零度弧齿齿向修正取零角极限；横向分配C及默认弯曲目标按弧齿分支")
    return out
