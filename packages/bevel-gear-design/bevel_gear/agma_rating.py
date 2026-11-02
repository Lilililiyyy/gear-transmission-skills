"""AGMA 2003C10 路线的零偏置钢材锥齿轮疲劳计算。

独立计算解析接触几何I、两轮危险截面J、载荷/温度/硬度系数、工作
应力、三点曲线容量与常幅损伤。标准名称表示计算路线，不代表完整
标准认证。材料数据、安装/载荷边界必须显式给定，没有材料数据库。
底层单位mm、MPa、N·m、rpm、h、rad；额定功率kW，Ra为µm。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .factors import positive, nonnegative, material_category
from .rating import _curve
from .rating_b2 import solve_critical_section


def _sqrt(value, name):
    if not math.isfinite(value) or value < -1e-10:
        raise ValueError(f"{name}不是有效非负几何量")
    return math.sqrt(max(value, 0.0))


def _sharing(chord, pitch, offset, limit):
    """所有实际存在接触弦的三次权重，不截成预设接触线数量。"""
    total = chord**3
    for sign in (1, -1):
        for number in range(1, limit + 1):
            radicand = chord**2 - 4 * number * pitch * (number * pitch + sign * 2 * offset)
            if radicand <= 0:
                break
            total += radicand**1.5
        else:
            raise ArithmeticError("AGMA接触弦求和达到次数上限")
    return chord**3 / positive(total, "AGMA接触弦三次权重 mm³")


def _geometry(values, limit):
    """使用显式外端/中点几何输入形成I/J所需的虚拟法截面。"""
    for key in ("z_1", "z_2"):
        if not isinstance(values[key], int) or isinstance(values[key], bool) or values[key] < 3:
            raise ValueError("AGMA两轮齿数必须为不小于3的整数")
    beta, alpha = float(values["beta_m2"]), float(values["alpha_nD"])
    if not 0 <= beta < math.pi / 2 or not 0 < alpha < math.pi / 3:
        raise ValueError("AGMA螺旋角须在[0,π/2)、压力角在(0,π/3)内")
    outer = positive(values["A_o"], "外端锥距 mm")
    width = positive(values["b_2"], "齿宽 mm")
    if width >= outer:
        raise ValueError("AGMA齿宽必须小于外端锥距")
    mean = outer - width / 2
    diameters = [positive(values[f"d_e{i}"], "外端节圆直径 mm") for i in (1, 2)]
    angles = [float(values[f"delta{i}"]) for i in (1, 2)]
    if any(not 0 < angle < math.pi / 2 for angle in angles):
        raise NotImplementedError("AGMA入口仅支持两节锥均为锐角的零偏置外啮合")
    pitch = positive(values["P_d"], "外端径节 1/mm")
    mean_pitch = outer / mean * pitch
    base_pitch = mean / outer * math.pi / pitch * math.cos(beta) * math.cos(alpha)
    normal_pitch = base_pitch / (math.cos(alpha)**2 * (math.cos(beta)**2 + math.tan(alpha)**2))
    virtual_radii = [diameter / (2 * math.cos(angle)) * mean / outer for diameter, angle in zip(diameters, angles)]
    normal_radii = [radius / math.cos(beta)**2 for radius in virtual_radii]
    bases = [radius * math.cos(alpha) for radius in normal_radii]
    mean_addenda = [positive(values[f"hae{i}"], "外端齿顶高 mm") - width / 2 * math.tan(values[f"deltaa{i}"] - angles[i - 1]) for i in (1, 2)]
    tips = [radius + addendum for radius, addendum in zip(normal_radii, mean_addenda)]
    rolls = [_sqrt(tip**2 - base**2, "法向齿顶展长平方") for tip, base in zip(tips, bases)]
    paths = [roll - radius * math.sin(alpha) for roll, radius in zip(rolls, normal_radii)]
    path = positive(sum(paths), "AGMA法向接触路径 mm")
    transverse = path / normal_pitch
    factor = width / outer * (2 - width / outer) / (2 * (1 - width / outer))
    overlap = (factor * math.tan(beta) - (factor * math.tan(beta))**3 / 3) * outer * pitch / math.pi
    if beta <= math.pi / 18:
        overlap = 0.0
    total = math.hypot(transverse, overlap)
    cosine = math.cos(alpha) * math.sqrt(math.cos(beta)**2 + math.tan(alpha)**2)
    if not -1 - 1e-14 <= cosine <= 1 + 1e-14:
        raise ValueError("AGMA基螺旋角余弦超域")
    base_helix = math.acos(min(1.0, max(-1.0, cosine)))
    squared_span = path**2 * math.cos(base_helix)**4 + (width * math.sin(base_helix))**2
    asymmetry = (values["z_2"] - values["z_1"]) / (3.2 * values["z_2"] + 4 * values["z_1"])
    return dict(beta=beta, alpha=alpha, outer=outer, mean=mean, width=width,
                diameters=diameters, angles=angles, pitch=pitch, mean_pitch=mean_pitch,
                base_pitch=base_pitch, virtual_radii=virtual_radii, normal_radii=normal_radii,
                bases=bases, tips=tips, rolls=rolls, paths=paths, path=path,
                transverse_contact_ratio=transverse, overlap_contact_ratio=overlap,
                total_contact_ratio=total, base_helix=base_helix, squared_span=squared_span,
                asymmetry=asymmetry)


def _pitting_geometry(g, limit):
    """接触曲率与实际弦权重，并保留步长逐级缩小的I评价搜索。"""
    beta, alpha, width, path = (g[key] for key in ("beta", "alpha", "width", "path"))
    cosine = math.cos(g["base_helix"])
    arms = [radius * math.sin(alpha) / cosine**2 for radius in g["virtual_radii"]]
    ratio = 2 / g["total_contact_ratio"] if g["total_contact_ratio"] < 2 else 1.0
    evaluations = 0

    def at(offset):
        nonlocal evaluations
        evaluations += 1
        if evaluations > limit:
            raise ArithmeticError("AGMA接触几何I搜索未在次数上限内收敛")
        chord = _sqrt(g["squared_span"] - 4 * offset**2, "I接触弦长平方")
        position = (path / 2 + path**2 * offset * cosine**2 / g["squared_span"]
                    + width * path * chord * g["asymmetry"] * math.sin(beta) / g["squared_span"] - g["paths"][1])
        first, second = arms[0] + position, arms[1] - position
        curvature_radius = positive(first, "I小轮曲率半径 mm") * positive(second, "I大轮曲率半径 mm") / (first + second)
        length = width * path * chord * cosine / g["squared_span"]
        sharing = _sharing(chord, g["base_pitch"], offset, limit)
        return positive(length * curvature_radius * math.cos(beta) * math.cos(alpha) * g["pitch"]
                        / (width * g["diameters"][0] * ratio * sharing * g["mean_pitch"]), "AGMA几何I")

    if beta <= math.pi / 18:
        offset = path / 2 - g["base_pitch"]
        result = at(offset)
        return result, {"evaluations": evaluations, "last_evaluated_offset_mm": offset, "search_step_mm": 0.0}
    offset, result, step = 0.0, 0.0, -math.sqrt(g["squared_span"])
    while True:
        step /= 10
        while True:
            previous_result, previous_offset = result, offset
            offset += step
            result = at(offset)
            if result - previous_result >= 0:
                break
        evaluated_offset = offset
        offset = previous_offset
        # 搜索输出是最后实际评价点；不把回退后的未评价位置冒充返回点。
        if abs(result - previous_result) <= 0.001:
            return result, {"evaluations": evaluations, "last_evaluated_offset_mm": evaluated_offset,
                            "restored_offset_mm": offset, "search_step_mm": step,
                            "last_factor_increment": result - previous_result}


def _bending_geometry(values, g, static, limit):
    beta, alpha, width, path = (g[key] for key in ("beta", "alpha", "width", "path"))
    cb, sb = math.cos(g["base_helix"]), math.sin(g["base_helix"])
    span = math.sqrt(g["squared_span"])
    offset = g["base_pitch"] - span / 2 if g["total_contact_ratio"] <= 2 else 0.0
    chord = _sqrt(g["squared_span"] - 4 * offset**2, "J接触弦长平方")
    if beta < math.pi / 18 and static:
        offset, chord = span / 2, 0.0
    positions, shifts = [], []
    for sign in (1, -1):
        if beta < math.pi / 18:
            positions.append(path / 2 + path**2 * offset / g["squared_span"])
            shifts.append(width * path * chord * g["asymmetry"] / g["squared_span"])
        else:
            positions.append(path / 2 + (path**2 * offset * cb**2 + sign * width * path * chord * g["asymmetry"] * sb) / g["squared_span"])
            shifts.append((width * path * chord * g["asymmetry"] * cb**2 - sign * width**2 * offset * sb) / g["squared_span"])
    sharing = 1.0 if beta < math.pi / 18 and static else _sharing(chord, g["base_pitch"], offset, limit)
    ratio = 1.0 if static or g["total_contact_ratio"] >= 2 else 2 / g["total_contact_ratio"]
    center_width = width * path * chord * cb**2 / g["squared_span"]
    sections, factors = [], []
    for index in (0, 1):
        radius, base = g["normal_radii"][index], g["bases"][index]
        pressure = math.atan((positions[index] + sum(g["normal_radii"]) * math.sin(alpha) - g["rolls"][1 - index]) / base)
        thickness = positive(values[f"smn{index + 1}"], "中点法向弧齿厚 mm")
        tooth_angle = thickness / (2 * radius) - math.tan(pressure) + pressure + math.tan(alpha) - alpha
        load_angle = pressure - tooth_angle
        load_height = base / math.cos(load_angle) - radius
        actual_radius = g["virtual_radii"][index] * (g["mean"] + shifts[index]) / g["mean"] + load_height
        taper = values["theta_a2"] if index == 0 else values["theta_f2"]
        root_height = values[f"hfe{index + 1}"] - width / 2 * math.tan(taper)
        section = solve_critical_section(2 * radius, root_height, values[f"rho_a0{index + 1}"], thickness,
            alpha, load_height, 1 / g["pitch"], iteration_limit=limit,
            residual_tolerance=0.001, initial_step_factor=0.127)
        half, lever = section.half_thickness_mm, section.lever_mm
        base_form = 2 / 3 * g["pitch"] / positive(lever / half**2 - math.tan(load_angle) / (3 * half), "AGMA截面形状表达式 1/mm")
        degrees = math.degrees(alpha)
        correction = (0.3254545 - 0.0072727 * degrees
                      + (2 * half / section.fillet_radius_mm)**(0.3318182 - 0.0090909 * degrees)
                      * (2 * half / lever)**(0.2681818 + 0.0090909 * degrees))
        total_width = (width - center_width) / math.cos(beta)
        left, right = total_width / 2 + shifts[index] / math.cos(beta), total_width / 2 - shifts[index] / math.cos(beta)
        effective_left = left if left > 0 and right > 0 else (total_width if left > 0 and right < 0 else 0.0)
        effective_right = right if left > 0 and right > 0 else (total_width if right > 0 and left < 0 else 0.0)
        effective = lever * math.cos(beta) * (math.atan(effective_left / lever) + math.atan(effective_right / lever)) + center_width
        factor = (base_form / positive(correction, "AGMA应力修正") / (sharing * ratio)
                  * actual_radius / g["virtual_radii"][index] * effective / width * g["pitch"] / g["mean_pitch"])
        factors.append(positive(factor, "AGMA几何J"))
        sections.append({"thickness_mm": 2 * half, "lever_mm": lever, "fillet_radius_mm": section.fillet_radius_mm,
                         "effective_width_mm": effective, "iterations": section.iterations, "residual": section.residual})
    return factors, sections, {"load_sharing_factor": sharing, "contact_ratio_factor": ratio}


def _damage(curve, cycles, stress, safety, policy):
    """兼容半对数损伤或显式双对数改进；不修改调用方曲线。"""
    stress *= safety
    if policy == "loglog":
        life = curve.cycles_at_stress(stress)
    elif stress > curve.stresses_mpa[0]:
        life = 1.0
    elif stress <= curve.stresses_mpa[-1]:
        life = math.inf
    else:
        for index in (0, 1):
            if stress > curve.stresses_mpa[index + 1]:
                high, low = curve.stresses_mpa[index:index + 2]
                first, second = curve.cycles[index:index + 2]
                life = first * math.exp((stress - high) * math.log(first / second) / (high - low))
                break
    return (None, None) if life is None else (cycles / life, life)


def evaluate_agma_rating_mm(inputs, *, static_load_distribution=None,
                            damage_interpolation="source_semilog", hardness_ratio_policy="floating",
                            iteration_limit=10000):
    """完整解析I/J与常幅疲劳链；长度mm，应力MPa，额定功率kW。

    source_semilog显式保留容量双对数而损伤半对数的策略；最高应力点
    以上该策略取失效寿命1循环。loglog使用一致双对数逆解，最高点以上
    返回未知寿命。曲线安全系数只作用损伤，不改容量或输入曲线。
    hardness_ratio_policy='integer_compatibility'仅恢复硬度比修正中的
    齿数比取整；默认floating使用真实比值。静载标志只改变J载荷分配。
    """
    values = dict(inputs)
    if values.get("offset_m", 0) != 0:
        raise NotImplementedError("AGMA2003评级入口尚未支持偏置")
    if values.get("userdefined", False):
        raise NotImplementedError("AGMA不开放跳过完整几何/系数链的预置系数模式")
    if damage_interpolation not in ("source_semilog", "loglog"):
        raise ValueError("损伤插值仅支持source_semilog或loglog")
    if hardness_ratio_policy not in ("floating", "integer_compatibility"):
        raise ValueError("硬度齿数比策略仅支持floating或integer_compatibility")
    if not isinstance(iteration_limit, int) or isinstance(iteration_limit, bool) or not 1 <= iteration_limit <= 1000000:
        raise ValueError("AGMA迭代上限须为1…1000000的整数")
    required = ("z_1", "z_2", "beta_m2", "alpha_nD", "A_o", "b_2", "d_e1", "d_e2", "delta1", "delta2",
                "deltaa1", "deltaa2", "hae1", "hae2", "hfe1", "hfe2", "P_d", "smn1", "smn2",
                "rho_a01", "rho_a02", "theta_a2", "theta_f2", "r_c0", "T_P", "n_P", "duration",
                "HB_1", "HB_2", "R_a1", "THETAoil", "isograde", "factor_KA",
                "sigma_Hlim1", "sigma_Hlim2", "sigma_Flim1", "sigma_Flim2")
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError("缺少AGMA评级字段：" + ", ".join(missing))
    static = values.get("isstaticload", 0) if static_load_distribution is None else static_load_distribution
    if static not in (0, 1, False, True):
        raise ValueError("AGMA静载分配标志须为布尔值")
    torque, speed = positive(values["T_P"], "扭矩 N·m"), positive(values["n_P"], "转速 rpm")
    hours = nonnegative(values["duration"], "时长 h")
    for index in (1, 2):
        positive(values[f"HB_{index}"], "布氏硬度")
        positive(values[f"sigma_Hlim{index}"], "接触极限 MPa")
        positive(values[f"sigma_Flim{index}"], "弯曲极限 MPa")
        name = "pinion" if index == 1 else "wheel"
        material_category(values.get(f"material_{name}", "case_hardened"))
    roughness = nonnegative(values["R_a1"], "小轮Ra µm")
    cutter = nonnegative(values["r_c0"], "刀盘半径 mm")
    temperature = float(values["THETAoil"])
    if not math.isfinite(temperature) or temperature <= -273:
        raise ValueError("油温须为大于−273的有限摄氏温度")
    grade = values["isograde"]
    if not isinstance(grade, int) or isinstance(grade, bool) or not 5 <= grade <= 12:
        raise ValueError("本AGMA动载路线支持整数精度等级5…12")
    g = _geometry(values, iteration_limit)
    if cutter == 0 and g["beta"] != 0:
        raise ValueError("非零螺旋角须给定正刀盘半径mm")
    out = {"geometry": g, "units": {"geometry": "mm", "stress": "MPa", "rated_power": "kW"},
           "static_load_distribution": bool(static), "damage_interpolation": damage_interpolation,
           "hardness_ratio_policy": hardness_ratio_policy, "S_Hmin": 1.0, "S_Fmin": 1.0,
           "n": values["z_1"], "N": values["z_2"], "diagnostics": []}
    out["Ifactor"], out["I_search"] = _pitting_geometry(g, iteration_limit)
    factors, sections, out["J_distribution"] = _bending_geometry(values, g, bool(static), iteration_limit)
    out["Jfactor_P"], out["Jfactor_G"] = factors
    out["root_section1"], out["root_section2"] = sections
    out["transversereContactRatioSmall"] = g["beta"] <= math.pi / 18 and g["transverse_contact_ratio"] < 1
    quality = 17 - grade
    exponent = 0.25 * (12 - quality)**0.667
    velocity_constant = 50 + 56 * (1 - exponent)
    velocity = 5.236e-5 * values["d_e1"] * speed
    limit_velocity = (velocity_constant + quality - 3)**2 / 200
    if velocity > limit_velocity:
        raise ValueError(f"AGMA节线速度{velocity:.6g}m/s超过动载式上限{limit_velocity:.6g}m/s")
    out["pitchLineVelocityIsOutoftheExtent"] = False
    out["K_v"] = (velocity_constant / (velocity_constant + math.sqrt(200 * velocity)))**(-exponent)
    out["K_O"] = positive(values["factor_KA"], "工况系数KO")
    width, pitch = g["width"], g["pitch"]
    out["C_s"] = 0.5 if width < 12.7 else (0.00492 * width + 0.4375 if width <= 79.8 else 0.83)
    out["K_s"] = 0.5 if 1 / pitch < 1.6 else 0.4867 + 0.008399 / pitch
    crowning = values.get("iscrowning", 1)
    mounting = values.get("mounting", 0)
    if crowning not in (0, 1) or mounting not in (0, 1, 2):
        raise ValueError("修鼓标志须为0/1，安装类别须为0/1/2")
    if values.get("K_m_defined", False):
        if "K_m" not in values:
            raise ValueError("显式载荷分布模式缺少K_m")
        out["K_m"] = positive(values.get("K_m"), "显式载荷分布系数Km")
    else:
        out["K_mb"] = (1.0, 1.1, 1.25)[mounting]
        out["K_m"] = (out["K_mb"] + 5.6e-6 * width**2) * (2 if crowning == 0 else 1)
    if values.get("crowningFactorDefined", False) and "C_xc" not in values:
        raise ValueError("显式修鼓模式缺少C_xc")
    out["C_xc"] = positive(values.get("C_xc"), "显式修鼓系数Cxc") if values.get("crowningFactorDefined", False) else (2.0 if crowning == 0 else 1.5)
    out["K_x"] = (1.0 if g["beta"] == 0 else min(1.15, max(1.0, 0.211 * (cutter / g["mean"])**(0.279 / math.log10(math.sin(g["beta"]))) + 0.789)))
    hardness_ratio = values["HB_1"] / values["HB_2"]
    ratio = values["z_2"] / values["z_1"] if hardness_ratio_policy == "floating" else values["z_2"] // values["z_1"]
    out["C_H"] = (1.0 if hardness_ratio < 1.2 else (1 + (0.00898 * hardness_ratio - 0.00829) * (ratio - 1) if hardness_ratio <= 1.7 else 1 + 0.00075 * math.exp(-0.52 * roughness) * (450 - values["HB_2"])))
    positive(out["C_H"], "AGMA硬度系数CH")
    out.update(K_T=(273 + temperature) / 393 if temperature > 120 else 1.0, C_R=1.0, K_R=1.0, C_p=190.0)
    out["s_c"] = out["C_p"] * math.sqrt(2000 * torque / (width * values["d_e1"]**2 * out["Ifactor"]) * math.prod(out[k] for k in ("K_O", "K_v", "K_m", "C_s", "C_xc")))
    cycles = speed * 60 * hours
    out["cycles1"], out["cycles2"] = cycles, cycles * values["z_1"] / values["z_2"]
    safety = positive(values.get("S_dutycucle", 1), "损伤曲线安全系数")
    for index, side in ((1, "P"), (2, "G")):
        h, f = _curve(values, "pinion" if index == 1 else "wheel", "H"), _curve(values, "pinion" if index == 1 else "wheel", "F")
        count = out[f"cycles{index}"]
        contact_capacity, root_capacity = h.stress_at_cycles(count), f.stress_at_cycles(count)
        out[f"s_ac{side}"], out[f"s_at{side}"] = values[f"sigma_Hlim{index}"], values[f"sigma_Flim{index}"]
        out[f"C_L{side}"], out[f"K_L{side}"] = contact_capacity / values[f"sigma_Hlim{index}"], root_capacity / values[f"sigma_Flim{index}"]
        out[f"s_wc{side}"] = contact_capacity * out["C_H"] / (out["K_T"] * out["C_R"])
        out[f"S_H{side}"] = out[f"s_wc{side}"] / out["s_c"]
        out[f"s_t{side}"] = (2000 * torque / (width * values["d_e1"]) * pitch * math.prod(out[k] for k in ("K_O", "K_v", "K_s", "K_m")) / (out["K_x"] * out[f"Jfactor_{side}"]))
        out[f"s_wt{side}"] = root_capacity / (out["K_T"] * out["K_R"])
        out[f"S_F{side}"] = out[f"s_wt{side}"] / out[f"s_t{side}"]
        member_speed = speed if index == 1 else speed * values["d_e1"] / values["d_e2"]
        diameter = values[f"d_e{index}"]
        out[f"P_ac{side}"] = (member_speed * width / 19100000 * out["Ifactor"] / math.prod(out[k] for k in ("K_v", "K_m", "K_O", "C_s", "C_xc")) * (contact_capacity * diameter * out["C_H"] / (out["C_p"] * out["K_T"] * out["C_R"]))**2)
        out[f"P_at{side}"] = (member_speed * width / 19100000 * out[f"Jfactor_{side}"] * out["K_x"] / math.prod(out[k] for k in ("K_s", "K_m", "K_O", "K_v")) * root_capacity * diameter / pitch / (out["K_T"] * out["K_R"]))
        out[f"damageH{side}"], out[f"lifeH{side}_cycles"] = _damage(h, count, out["s_c"] / out["C_H"] * out["K_T"] * out["C_R"], safety, damage_interpolation)
        out[f"damageF{side}"], out[f"lifeF{side}_cycles"] = _damage(f, count, out[f"s_t{side}"] * out["K_T"] * out["K_R"], safety, damage_interpolation)
    out["diagnostics"].append("静载标志只改变齿根接触分配；本输出仍为常幅疲劳，不含屈服/塑性校核")
    if damage_interpolation == "source_semilog":
        out["diagnostics"].append("容量采用双对数插值，损伤显式采用半对数；最高S-N应力以上取1循环失效")
    out["diagnostics"].append("两轮额定接触功率分别沿用该路线公式；大轮值含额外直径比，不能直接当作整机额定功率")
    return out


@dataclass(frozen=True)
class AGMARatingSettings:
    """独立AGMA工艺/工况输入；长度m，油温℃，系数无量纲。"""
    pinion_cutter_tip_radius_m: float
    wheel_cutter_tip_radius_m: float
    application_factor: float
    quality_grade: int = 8
    mounting: int = 0
    crowned: bool = True
    oil_temperature_c: float = 80.0
    curve_safety_factor: float = 1.0
    face_load_factor: float | None = None
    crowning_factor: float | None = None
    static_load_distribution: bool = False
    damage_interpolation: str = "source_semilog"
    hardness_ratio_policy: str = "floating"

    def __post_init__(self):
        for key in ("pinion_cutter_tip_radius_m", "wheel_cutter_tip_radius_m", "application_factor", "curve_safety_factor"):
            positive(getattr(self, key), key)
        for key in ("face_load_factor", "crowning_factor"):
            if getattr(self, key) is not None:
                positive(getattr(self, key), key)
        for key in ("crowned", "static_load_distribution"):
            if not isinstance(getattr(self, key), bool):
                raise ValueError(f"{key}须为布尔值")
        if self.damage_interpolation not in ("source_semilog", "loglog"):
            raise ValueError("损伤策略须为source_semilog或loglog")
        if self.hardness_ratio_policy not in ("floating", "integer_compatibility"):
            raise ValueError("硬度比策略须为floating或integer_compatibility")


def agma_inputs_from_si(geometry, load, pinion_material, wheel_material, settings):
    """接收现有具名几何/载荷/材料；只作单位与语义映射，不套用ISO系数。"""
    if geometry.inputs.offset_m != 0:
        raise NotImplementedError("AGMA2003入口只接收零偏置几何")
    if geometry.inputs.backlash_m != 0:
        raise NotImplementedError("AGMA集成暂不接受非零侧隙，底层显式齿厚字典可单独验证")
    if geometry.inputs.cutter_radius_m is None and geometry.inputs.spiral_angle_rad != 0:
        raise ValueError("AGMA齿线修正必须显式给定刀盘半径")
    data = geometry.inputs
    values = dict(z_1=data.pinion_teeth, z_2=data.wheel_teeth,
                  beta_m1=data.spiral_angle_rad, beta_m2=data.spiral_angle_rad,
                  alpha_nD=data.pressure_angle_rad, alpha_nC=data.pressure_angle_rad,
                  A_o=geometry.outer_cone_distance_m * 1e3, b_2=data.face_width_m * 1e3,
                  P_d=1 / (geometry.outer_transverse_module_m * 1e3),
                  r_c0=(data.cutter_radius_m or 0.0) * 1e3, THETAoil=settings.oil_temperature_c,
                  theta_a2=geometry.wheel.addendum_taper_angle_rad,
                  theta_f2=geometry.wheel.dedendum_taper_angle_rad,
                  T_P=load.pinion_torque_n_m, n_P=load.pinion_speed_rpm, duration=load.duration_h,
                  isograde=settings.quality_grade, mounting=settings.mounting, iscrowning=int(settings.crowned),
                  factor_KA=settings.application_factor, S_dutycucle=settings.curve_safety_factor,
                  isstaticload=int(settings.static_load_distribution), K_m_defined=settings.face_load_factor is not None,
                  crowningFactorDefined=settings.crowning_factor is not None)
    if settings.face_load_factor is not None:
        values["K_m"] = settings.face_load_factor
    if settings.crowning_factor is not None:
        values["C_xc"] = settings.crowning_factor
    for index, member, material in ((1, geometry.pinion, pinion_material), (2, geometry.wheel, wheel_material)):
        if material.category not in ("case_hardened", "through_hardened", "normalized", "nitrided", "induction_hardened"):
            raise NotImplementedError("AGMA固定弹性常数入口仅接收明确钢材数据")
        values.update({f"d_e{index}": member.outer_pitch_diameter_m * 1e3,
                       f"delta{index}": member.pitch_cone_angle_rad, f"deltaa{index}": member.tip_cone_angle_rad,
                       f"hae{index}": member.outer_addendum_m * 1e3, f"hfe{index}": member.outer_dedendum_m * 1e3,
                       f"smn{index}": member.mean_normal_tooth_thickness_m * 1e3,
                       f"rho_a0{index}": getattr(settings, "pinion_cutter_tip_radius_m" if index == 1 else "wheel_cutter_tip_radius_m") * 1e3,
                       f"HB_{index}": material.hardness_hb, f"sigma_Hlim{index}": material.contact_limit_pa / 1e6,
                       f"sigma_Flim{index}": material.bending_limit_pa / 1e6})
        prefix = "pinion" if index == 1 else "wheel"
        values[f"contact_sn_{prefix}"] = material.contact_curve.to_mpa()
        values[f"bending_sn_{prefix}"] = material.bending_curve.to_mpa()
    values["R_a1"] = pinion_material.ra_m * 1e6
    return values


def evaluate_agma_rating_si(geometry, load, pinion_material, wheel_material, settings, *, iteration_limit=10000):
    """SI适配：应力/容量Pa，功率W；原系数独立返回，不与ISO因子混用。"""
    values = agma_inputs_from_si(geometry, load, pinion_material, wheel_material, settings)
    raw = evaluate_agma_rating_mm(values, damage_interpolation=settings.damage_interpolation,
                                 hardness_ratio_policy=settings.hardness_ratio_policy, iteration_limit=iteration_limit)
    return {"contact_stress_pa": raw["s_c"] * 1e6,
            "elasticity_factor_sqrt_pa": raw["C_p"] * 1e3,
            "bending_stress_pa": tuple(raw[f"s_t{side}"] * 1e6 for side in ("P", "G")),
            "contact_capacity_pa": tuple(raw[f"s_wc{side}"] * 1e6 for side in ("P", "G")),
            "bending_capacity_pa": tuple(raw[f"s_wt{side}"] * 1e6 for side in ("P", "G")),
            "contact_safety": tuple(raw[f"S_H{side}"] for side in ("P", "G")),
            "bending_safety": tuple(raw[f"S_F{side}"] for side in ("P", "G")),
            "rated_contact_power_w": tuple(raw[f"P_ac{side}"] * 1e3 for side in ("P", "G")),
            "rated_bending_power_w": tuple(raw[f"P_at{side}"] * 1e3 for side in ("P", "G")),
            "life_cycles": {mode: tuple(raw[f"life{mode}{side}_cycles"] for side in ("P", "G")) for mode in ("H", "F")},
            "damage": {mode: tuple(raw[f"damage{mode}{side}"] for side in ("P", "G")) for mode in ("H", "F")},
            "coefficients": {key: raw[key] for key in ("Ifactor", "Jfactor_P", "Jfactor_G", "K_O", "K_v", "C_s", "K_s", "K_m", "C_xc", "K_x", "C_H", "K_T", "C_R", "K_R")},
            "diagnostics": tuple(raw["diagnostics"]), "raw_mm_mpa": raw}
