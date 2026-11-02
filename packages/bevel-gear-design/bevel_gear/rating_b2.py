"""零偏置钢材锥齿轮的 B2 齿根评级与静载接触分配。

长度 mm、应力 MPa、角度 rad、扭矩 N·m、转速 rpm、时间 h。
B2 从危险截面、有效载荷宽度与 Y_A/Y_J 计算两轮齿根应力；不以
B1 的齿形系数代替。静载标志仅改变接触分配，不构成屈服/塑性校核。
两轮危险截面均用有符号残差割线；失败明确报错，不伪装成收敛。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from . import factors as f
from .rating import _curve, _virtual_geometry


@dataclass(frozen=True)
class CriticalSection:
    """危险截面：半齿厚、力臂、截面倾角及根圆角，长度 mm。"""
    half_thickness_mm: float
    lever_mm: float
    section_angle_rad: float
    parameter_mm: float
    fillet_radius_mm: float
    iterations: int
    residual: float


def solve_critical_section(pitch_diameter_mm, root_height_mm, cutter_radius_mm,
                            mean_tooth_thickness_mm, pressure_angle_rad,
                            load_height_mm, outer_module_mm, *, iteration_limit=10000,
                            residual_tolerance=1e-5, initial_step_factor=0.005):
    """解 s/(tanθ·h)−2=0；默认采用0.005模数探步及1e−5残差。

    两轮均保留残差符号。上限保护、零割线分母、非有限值和无效危险
    截面报中文异常；返回的是最后实际评价点，不返回未评价的下一步。
    首步/模数和残差容差可显式指定，用于其他解析评级路线。
    """
    diameter = f.positive(pitch_diameter_mm, "虚拟法向节圆直径 mm")
    radius = f.positive(cutter_radius_mm, "刀尖圆角半径 mm")
    height = f.positive(root_height_mm, "中点齿根高 mm")
    tooth = f.positive(mean_tooth_thickness_mm, "中点法向齿厚 mm")
    module = f.positive(outer_module_mm, "外端模数 mm")
    tolerance = f.positive(residual_tolerance, "危险截面残差容差")
    step = f.positive(initial_step_factor, "危险截面首探步/模数") * module
    pressure_angle_rad = float(pressure_angle_rad)
    load_height_mm = float(load_height_mm)
    if not 0 < pressure_angle_rad < math.pi / 3:
        raise ValueError("危险截面压力角须在(0,π/3)内")
    if not math.isfinite(load_height_mm):
        raise ValueError("加载点相对节圆高度必须为有限数，单位mm")
    if not isinstance(iteration_limit, int) or isinstance(iteration_limit, bool) or not 1 <= iteration_limit <= 1000000:
        raise ValueError("危险截面迭代上限须为1…1000000的整数")
    radial = height - radius
    offset = tooth / 2 + height * math.tan(pressure_angle_rad) + radius * (1 - math.sin(pressure_angle_rad)) / math.cos(pressure_angle_rad)
    parameter = offset + radial
    previous_parameter = previous_residual = None
    for iteration in range(1, iteration_limit + 1):
        roll_angle = 2 * parameter / diameter
        distance = parameter - offset
        radial_part = radial * math.cos(roll_angle) - distance * math.sin(roll_angle)
        tangent_part = radial * math.sin(roll_angle) + distance * math.cos(roll_angle)
        if tangent_part == 0:
            raise ArithmeticError("危险截面切向坐标为零，截面倾角未定义")
        section_angle = math.atan(radial_part / tangent_part)
        half_thickness = diameter / 2 * math.sin(roll_angle) - radius * math.cos(section_angle) - tangent_part
        lever = (load_height_mm + diameter / 2 * (1 - math.cos(roll_angle))
                 + radius * math.sin(section_angle) + radial_part)
        denominator = math.tan(section_angle) * lever
        if denominator == 0:
            raise ArithmeticError("危险截面残差分母为零")
        residual = half_thickness / denominator - 2
        if not all(math.isfinite(v) for v in (parameter, residual, half_thickness, lever)):
            raise ArithmeticError("危险截面迭代产生非有限数")
        if abs(residual) <= tolerance:
            fillet_denominator = diameter / 2 + height - radius
            fillet = radial**2 / f.positive(fillet_denominator, "根圆角计算分母 mm") + radius
            return CriticalSection(f.positive(half_thickness, "危险截面半齿厚 mm"),
                f.positive(lever, "危险截面力臂 mm"), section_angle, parameter,
                f.positive(fillet, "根圆角半径 mm"), iteration, residual)
        if previous_parameter is None:
            trial = parameter + step
        else:
            difference = residual - previous_residual
            if difference == 0:
                raise ArithmeticError("危险截面割线残差不再变化，不能报告收敛")
            trial = previous_parameter - previous_residual * (parameter - previous_parameter) / difference
        previous_parameter, previous_residual = parameter, residual
        parameter = trial
    raise ArithmeticError(f"危险截面在{iteration_limit}次内未达到{tolerance:g}残差，末残差={residual:.6g}")


def relative_sensitivity_b2(section):
    """B2相对缺口敏感性：截面齿厚/(2根圆角)≥1.5时为1，否则0.95。"""
    ratio = section.thickness_mm / (2 * section.fillet_radius_mm)
    return 1.0 if ratio >= 1.5 else 0.95


def _common_state(values):
    """复用几何与系数函数形成方法无关的点蚀/载荷状态，不计算B1齿根。"""
    required = ("z_1", "z_2", "sigma", "d_e2", "b_2", "beta_m1", "beta_m2", "alpha_nD",
                "k_hap", "k_hfp", "x_hm1", "x_smn", "T_1", "n_1", "duration", "r_c0",
                "rho_a01", "rho_a02", "s_pr1", "s_pr2", "sigma_Hlim1", "sigma_Hlim2",
                "sigma_Flim1", "sigma_Flim2", "HB_1", "HB_2", "R_a1", "R_a2", "v40", "isograde")
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError("缺少B2评级字段：" + ", ".join(missing))
    if values.get("userdefined", False):
        raise NotImplementedError("B2不支持绕过完整系数链的预设系数模式")
    geometry = _virtual_geometry(values)
    beta = geometry["spiral_angle_rad"]
    gear_type = values.get("geartype", 0 if beta == 0 else 2)
    if gear_type not in (0, 1, 2) or (gear_type < 2) != (beta == 0):
        raise ValueError("B2支持零平均螺旋角的直齿/零度齿及正角度螺旋锥齿")
    materials = tuple(f.material_category(values.get(f"material_{name}", values.get(f"abbreviation{i}", "case_hardened")))
                      for i, name in ((1, "pinion"), (2, "wheel")))
    torque = f.positive(values["T_1"], "B2常幅扭矩 N·m")
    speed = f.positive(values["n_1"], "转速 rpm")
    hours = f.nonnegative(values["duration"], "持续时间 h")
    width, module = geometry["face_width_mm"], geometry["m_mn"]
    tangential = 2000 * torque / geometry["d_m1"]
    velocity = geometry["d_m1"] * speed / 19098
    out = {"geometry": geometry, "delta_1": geometry["delta_1"], "delta_2": geometry["delta_2"],
           "d_m1": geometry["d_m1"], "F_mt1": tangential, "diagnostics": [],
           "root_method": "B2", "units": {"stress": "MPa", "geometry": "mm", "roughness": "µm"}}
    for index in (1, 2):
        angle = float(values[f"beta_m{index}"])
        if not 0 <= angle < math.pi / 2:
            raise ValueError("B2平均螺旋角须取[0,π/2)内大小值")
    out["F_mt2"] = tangential * math.cos(float(values["beta_m2"])) / math.cos(float(values["beta_m1"]))
    for flank, sign in (("D", 1), ("C", -1)):
        pressure = float(values.get(f"alpha_n{flank}", geometry["normal_pressure_angle_rad"]))
        if not 0 < pressure < math.pi / 3:
            raise ValueError("法向压力角须在(0,π/3)内")
        for index, orientation in ((1, sign), (2, -sign)):
            spiral, delta = float(values[f"beta_m{index}"]), geometry[f"delta_{index}"]
            out[f"F_ax{index}{flank}"] = out[f"F_mt{index}"] * (math.tan(pressure) * math.sin(delta) / math.cos(spiral) + orientation * math.tan(spiral) * math.cos(delta))
            out[f"F_rad{index}{flank}"] = out[f"F_mt{index}"] * (math.tan(pressure) * math.cos(delta) / math.cos(spiral) - orientation * math.tan(spiral) * math.sin(delta))
    epsilon, overlap, total = (geometry[key] for key in ("epsilon_valpha", "epsilon_vbeta", "epsilon_vgamma"))
    line = f.contact_line_geometry(geometry["g_valpha"], geometry["beta_vb"], width, geometry["p_et"], epsilon, overlap)
    cosine = math.cos(geometry["beta_vb"])
    out["Z_H"] = 2 * math.sqrt(cosine / math.sin(2 * geometry["alpha_vt"]))
    point_a = epsilon if overlap > 1 else 2 + (epsilon - 2) * overlap
    point_b = epsilon if overlap > 1 else 2 * epsilon - 2 + (2 - epsilon) * overlap
    terms = [math.sqrt(geometry[f"d_va{i}"]**2 / geometry[f"d_vb{i}"]**2 - 1)
             - point * math.pi / (values[f"z_{i}"] / math.cos(geometry[f"delta_{i}"]))
             for i, point in ((1, point_a), (2, point_b))]
    out["Z_MB"] = math.tan(geometry["alpha_vt"]) / math.sqrt(f.positive(terms[0], "小轮点蚀曲率") * f.positive(terms[1], "大轮点蚀曲率"))
    out.update(Z_E=f.positive(values.get("elasticity_factor_sqrt_mpa", 189.8), "钢材弹性系数 √MPa"),
               Z_beta=math.sqrt(math.cos(beta)), Z_K=0.8, Z_X1=1.0, Z_X2=1.0)
    limits = [f.positive(values[f"sigma_Hlim{i}"], "接触疲劳极限 MPa") for i in (1, 2)]
    out["Z_LS"] = f.load_sharing_factor(total, overlap, geometry["p_et"], geometry["beta_vb"], line)
    out["ZLZvZR"] = f.lubricant_film_factor(*limits, values["v40"], velocity, values["R_a1"], values["R_a2"], geometry["a_v"], geometry["alpha_vt"], geometry["beta_vb"], geometry["u_v"])
    out["Z_W"] = f.hardness_ratio_factor(values["HB_1"], values["HB_2"])
    out["K_A"] = f.positive(values["factor_KA"] if "factor_KA" in values else f.application_factor(values.get("driver_duty", "uniform"), values.get("driven_duty", "uniform")), "KA")
    out["K_Hbeta"] = f.positive(values["factor_KHbeta"] if "factor_KHbeta" in values else f.face_load_factor(0.85, values.get("assembly", "full_load_matched"), values.get("support", "both_supported")), "KHβ")
    out["f_pt"] = f.single_pitch_deviation(module, geometry["d_m2"], values["isograde"])
    out["y_alpha"] = f.running_in_allowance(materials[0], limits[0], materials[1], limits[1], velocity, out["f_pt"])
    method = values.get("Kvmethod", 0)
    if method == 1:
        out["K_v"], out["dynamic_diagnostics"] = f.dynamic_factor_b(speed, geometry["d_m1"], values["z_1"], geometry["normal_pressure_angle_rad"], geometry["ratio"], tangential, width, 0.85 * width, out["K_A"], total, out["f_pt"], out["y_alpha"])
    elif method == 0:
        out["K_v"] = f.dynamic_factor_c(module, geometry["d_m1"], geometry["d_e1"], values["isograde"], velocity)
    else:
        raise NotImplementedError("B2动载仅支持B=1或C=0")
    out["dynamicFactorErrorCode"] = 0
    # 零度齿的sinβ→0极限使修正项趋于1；此处避免直接log(0)。
    out["K_Fbeta"] = out["K_Hbeta"] if beta == 0 else f.bending_face_load_factor(False, beta, geometry["R_m"], values["r_c0"], out["K_Hbeta"])
    out["Y_epsilon"] = f.contact_ratio_factor(epsilon, overlap)
    method = values.get("KHalphamethod", 1)
    if method == 1:
        out["K_Halpha"], out["K_Falpha"] = f.transverse_load_factors_b(width, epsilon, overlap, total, tangential, out["K_A"], out["K_v"], out["K_Hbeta"], out["Z_LS"], out["Y_epsilon"], out["y_alpha"], out["f_pt"])
    elif method == 0:
        missing = [f"surfacetreatment{i}" for i in (1, 2) if f"surfacetreatment{i}" not in values]
        if missing:
            raise ValueError("横向分配C缺少热处理类别：" + ", ".join(missing))
        out["K_Halpha"] = out["K_Falpha"] = max(f.transverse_load_factor_c(values[f"surfacetreatment{i}"], gear_type == 0, values["isograde"], tangential, out["K_A"], 0.85 * width, epsilon) for i in (1, 2))
    else:
        raise NotImplementedError("B2横向分配仅支持B=1或C=0")
    nominal = math.sqrt(tangential / (geometry["d_v1"] * f.positive(line.middle_length_mm, "点蚀接触线 mm")) * (geometry["u_v"] + 1) / geometry["u_v"])
    out["sigma_H"] = nominal * math.prod(out[key] for key in ("Z_MB", "Z_H", "Z_E", "Z_LS", "Z_beta", "Z_K")) * math.sqrt(out["K_A"] * out["K_v"] * out["K_Halpha"] * out["K_Hbeta"])
    out["contact_line_pitting_mm"] = (line.tip_length_mm, line.middle_length_mm, line.root_length_mm)
    return out, materials, gear_type


def _positive_sqrt(value, name):
    """只容忍量级内的舍入负数，不把真实非法几何裁成零。"""
    if value < -1e-10:
        raise ValueError(f"{name}的平方根参数为负")
    return math.sqrt(max(value, 0))


def _b2_sections(values, geometry, gear_type, static, iteration_limit):
    beta, alpha = geometry["spiral_angle_rad"], geometry["normal_pressure_angle_rad"]
    base_helix, width = geometry["beta_vb"], geometry["face_width_mm"]
    cosine, sine = math.cos(base_helix), math.sin(base_helix)
    pitches = [geometry[f"d_v{i}"] / cosine**2 for i in (1, 2)]
    normal_center = sum(pitches) / 2
    path = sum(_positive_sqrt(geometry[f"d_van{i}"]**2 - geometry[f"d_vbn{i}"]**2, "法向接触路径") for i in (1, 2)) / 2 - normal_center * math.sin(alpha)
    path = f.positive(path, "法向接触路径 mm")
    span = math.hypot(path * cosine**2, width * sine)
    offset = geometry["p_et"] * cosine - span / 2 if geometry["epsilon_vgamma"] <= 2 else 0.0
    if static and gear_type < 2:
        offset = span / 2
    chord = _positive_sqrt(span**2 - 4 * offset**2, "有效接触弦长")
    tooth_difference = (values["z_2"] - values["z_1"]) / (3.2 * values["z_2"] + 4 * values["z_1"])
    positions, shifts = [], []
    for direction in (1, -1):
        if gear_type < 2:
            positions.append(path / 2 + path**2 * offset / span**2)
            shifts.append(width * path * chord * tooth_difference / span**2)
        else:
            positions.append(path / 2 + (path**2 * offset * cosine**2 + direction * width * path * chord * tooth_difference * sine) / span**2)
            shifts.append((width * path * chord * tooth_difference * cosine**2 - direction * width**2 * offset * sine) / span**2)
    sections, load_angles, radii = [], [], []
    for index, direction in ((1, 1), (2, -1)):
        opposite = 3 - index
        diameter, pitch = geometry[f"d_vbn{index}"], pitches[index - 1]
        pressure = math.atan((positions[index - 1] + normal_center * math.sin(alpha)
                             - _positive_sqrt(geometry[f"d_van{opposite}"]**2 - geometry[f"d_vbn{opposite}"]**2, "对轮齿顶展长") / 2) / (diameter / 2))
        thickness = (math.pi / 2 + 2 * direction * (geometry["thickness_shift"] + geometry["profile_shift"] * math.tan(alpha))) * geometry["m_mn"]
        tooth_angle = thickness / pitch - math.tan(pressure) + pressure + math.tan(alpha) - alpha
        load_angle = pressure - tooth_angle
        load_height = (diameter / math.cos(load_angle) - pitch) / 2
        actual_radius = geometry[f"d_v{index}"] / 2 * (geometry["R_m"] + shifts[index - 1]) / geometry["R_m"] + load_height
        root_height = geometry["m_mn"] * (values["k_hfp"] - direction * geometry["profile_shift"])
        section = solve_critical_section(pitch, root_height, values[f"rho_a0{index}"], thickness,
            alpha, load_height, geometry["m_et"], iteration_limit=iteration_limit)
        sections.append(section)
        load_angles.append(load_angle)
        radii.append(f.positive(actual_radius, "B2实际加载半径 mm"))
    # 多条接触弦的三次权重；静载直/零度齿存在弦长零极限，直接取其定义值1。
    pitch_projection = geometry["p_et"] * cosine
    cube_sum = chord**3
    for multiplier in (4, 8):
        for direction in (1, -1):
            radicand = chord**2 - multiplier * pitch_projection * (multiplier / 4 * pitch_projection + direction * 2 * offset)
            if radicand > 0:
                cube_sum += radicand**1.5
    sharing = 1.0 if static and gear_type < 2 else chord**3 / f.positive(cube_sum, "B2接触弦权重 mm³")
    ratio_factor = 1.0 if static or geometry["epsilon_vgamma"] >= 2 else 2 / geometry["epsilon_vgamma"]
    central_width = width * path * chord * cosine**2 / span**2
    actual_widths = []
    for shift, section in zip(shifts, sections):
        total = (width - central_width) / math.cos(beta)
        left = total / 2 + shift / math.cos(beta)
        right = total / 2 - shift / math.cos(beta)
        if left > 0 and right > 0:
            pass
        elif left > 0 and right < 0:
            left, right = total, 0.0
        elif left < 0 and right > 0:
            left, right = 0.0, total
        else:
            left = right = 0.0
        actual_widths.append(f.positive(section.lever_mm * math.cos(beta)
            * (math.atan(left / section.lever_mm) + math.atan(right / section.lever_mm))
            + central_width, "B2有效载荷宽度 mm"))
    return sections, load_angles, radii, actual_widths, sharing, ratio_factor


def evaluate_b2_rating_mm(inputs, *, static_load_distribution=None, iteration_limit=10000):
    """完整B2钢材常幅数值入口，输出结构与B1的应力/容量/寿命字段一致。

    输入与B1底层字典单位相同，不原地修改。static_load_distribution
    只改变B2接触载荷分配；后续仍用转速、持续时间和S-N疲劳容量。
    不输出虚假的静屈服安全系数。B1专用Y_Fa/Y_Sa/Y_K/Y_LS字段不适用，
    因而不以零填充；B2直接返回实际计算的Y_A/Y_J及截面诊断。
    """
    values = dict(inputs)
    flag = values.get("isstaticload", 0) if static_load_distribution is None else static_load_distribution
    if flag not in (False, True, 0, 1):
        raise ValueError("静载分配标志必须为布尔值")
    static = bool(flag)
    out, materials, gear_type = _common_state(values)
    geometry = out["geometry"]
    sections, angles, radii, widths, sharing, ratio_factor = _b2_sections(
        values, geometry, gear_type, static, iteration_limit)
    pressure_deg = math.degrees(geometry["normal_pressure_angle_rad"])
    constant = 0.3254545 - 0.0072727 * pressure_deg
    notch_exponent = 0.3318183 - 0.0090909 * pressure_deg
    slender_exponent = 0.2681818 + 0.0090909 * pressure_deg
    wheel_cut = values.get("wheelcutmethod", 2)
    if wheel_cut not in (0, 1, 2):
        raise NotImplementedError("B2缺口敏感性仅支持成形或展成刀具")
    out["static_load_distribution"] = static
    out["b2_load_sharing_factor"] = sharing
    out["b2_contact_ratio_factor"] = ratio_factor
    out["S_Hmin"], out["S_Fmin"] = 1.0, 1.5 if gear_type == 0 else 1.3
    cycles = values["n_1"] * 60 * values["duration"]
    out["cycles1"], out["cycles2"] = cycles, cycles / geometry["ratio"]
    for index, direction in ((1, 1), (2, -1)):
        section = sections[index - 1]
        half, lever = section.half_thickness_mm, section.lever_mm
        base_form = 2 / (3 * geometry["m_et"] * f.positive(lever / half**2 - math.tan(angles[index - 1]) / (3 * half), "B2齿形表达式 1/mm"))
        correction = constant + (2 * half / section.fillet_radius_mm)**notch_exponent * (2 * half / lever)**slender_exponent
        form = base_form / f.positive(correction, "B2应力修正")
        out[f"Y_A{index}"] = correction / f.positive(2.3 * (1 - half / (3 * lever) * math.tan(geometry["normal_pressure_angle_rad"])), "B2应力系数分母")
        out[f"Y_J{index}"] = form / (sharing * ratio_factor) * 2 * radii[index - 1] / geometry[f"d_v{index}"] * widths[index - 1] / geometry["face_width_mm"] * geometry["m_mt"] / geometry["m_et"]
        out[f"sigma_F0{index}"] = out["F_mt1"] / geometry["face_width_mm"] * geometry["m_mt"] / geometry["m_et"]**2 * out[f"Y_A{index}"] / f.positive(out[f"Y_J{index}"], "B2几何系数YJ")
        out[f"sigma_F{index}"] = out[f"sigma_F0{index}"] * out["K_A"] * out["K_v"] * out["K_Fbeta"] * out["K_Falpha"]
        # B2缺口评价明确按给定刀具加工方式，两轮分别传相反变位/齿厚修正。
        sensitivity_section = f.root_section(geometry["cutter_height_coefficient"], geometry["m_mn"],
            direction * geometry["thickness_shift"], direction * geometry["profile_shift"],
            geometry["normal_pressure_angle_rad"], values[f"rho_a0{index}"], values[f"s_pr{index}"],
            geometry[f"z_vn{index}"], generated=wheel_cut == 2)
        out[f"Y_deltarelT{index}"] = relative_sensitivity_b2(sensitivity_section)
        out[f"Y_RrelT{index}"] = 1.0
        out[f"Y_X{index}"] = f.size_factor(geometry["m_mn"], materials[index - 1])
        out[f"root_section{index}"] = {"thickness_mm": 2 * half, "fillet_radius_mm": section.fillet_radius_mm,
            "lever_mm": lever, "effective_width_mm": widths[index - 1], "section_angle_rad": section.section_angle_rad,
            "iterations": section.iterations, "residual": section.residual}
        prefix = "pinion" if index == 1 else "wheel"
        curve_h, curve_f = _curve(values, prefix, "H"), _curve(values, prefix, "F")
        cycle_count = out[f"cycles{index}"]
        out[f"Z_NT{index}"] = curve_h.stress_at_cycles(cycle_count) / f.positive(values[f"sigma_Hlim{index}"], "接触疲劳极限 MPa")
        out[f"Y_NT{index}"] = curve_f.stress_at_cycles(cycle_count) / f.positive(values[f"sigma_Flim{index}"], "弯曲疲劳极限 MPa")
        out[f"sigma_HP{index}"] = curve_h.stress_at_cycles(cycle_count) * out[f"Z_X{index}"] * out["ZLZvZR"] * out["Z_W"]
        multiplier = 2 * out[f"Y_deltarelT{index}"] * out[f"Y_RrelT{index}"] * out[f"Y_X{index}"]
        out[f"sigma_FP{index}"] = curve_f.stress_at_cycles(cycle_count) * multiplier
        out[f"S_H{index}"] = out[f"sigma_HP{index}"] / out["sigma_H"]
        out[f"S_F{index}"] = out[f"sigma_FP{index}"] / out[f"sigma_F{index}"]
        out[f"sigmaF{index}"] = out[f"sigma_F{index}"] / multiplier
        safety = values.get("S_dutycucle", 1)
        out[f"damageH{index}"], out[f"lifeH{index}_cycles"] = curve_h.damage(cycle_count, out["sigma_H"] / (out[f"Z_X{index}"] * out["ZLZvZR"] * out["Z_W"]), curve_safety_factor=safety)
        out[f"damageF{index}"], out[f"lifeF{index}_cycles"] = curve_f.damage(cycle_count, out[f"sigmaF{index}"], curve_safety_factor=safety)
        if out[f"damageH{index}"] is None or out[f"damageF{index}"] is None:
            out["diagnostics"].append(f"第{index}轮等效应力超过三点S-N上界，寿命不作外推")
    line_f = f.contact_line_geometry(geometry["g_valpha"], geometry["beta_vb"], geometry["face_width_mm"], geometry["p_et"], geometry["epsilon_valpha"], geometry["epsilon_vbeta"], for_bending=True)
    out["contact_line_bending_mm"] = (line_f.tip_length_mm, line_f.middle_length_mm, line_f.root_length_mm)
    out["not_applicable_factors"] = ("Y_Fa1", "Y_Fa2", "Y_Sa1", "Y_Sa2", "Y_K", "Y_LS")
    out["diagnostics"].append("B2两轮使用有符号割线及1e−5危险截面残差；静载标志仅为接触分配，S-N容量仍为疲劳模型")
    return out
