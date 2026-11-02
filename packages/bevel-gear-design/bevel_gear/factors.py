"""ISO 10300 常幅钢材评级所用的明确计算分支。

本模块为重新组织的新写数学实现；仅覆盖说明中的 B1 齿根及 B/C
载荷系数，不代表完整标准认证。底层长度 mm、角度 rad、应力 MPa，
粗糙度与齿距偏差为 µm，运动黏度 mm²/s，线速度 m/s。
"""
from __future__ import annotations

from dataclasses import dataclass
import math


def positive(value, name):
    """检查具有明确物理单位的有限正输入。"""
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name}必须为有限正数")
    return value


def nonnegative(value, name):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name}必须为有限非负数")
    return value


def material_category(value):
    """接受显式钢材类别及通用冶金缩写，未知类别不猜测。"""
    names = {
        "case_hardened": "case_hardened", "Eh": "case_hardened",
        "through_hardened": "through_hardened", "V": "through_hardened",
        "normalized": "normalized", "St": "normalized",
        "nitrided": "nitrided", "NT": "nitrided", "NV": "nitrided",
        "induction_hardened": "induction_hardened", "IF": "induction_hardened",
    }
    if value not in names:
        raise NotImplementedError(f"尚未实现钢材类别 {value!r} 的完整系数分支")
    return names[value]


@dataclass(frozen=True)
class ContactLine:
    """三个评价位置的接触线长度及位置，全部为 mm。"""
    tip_length_mm: float
    middle_length_mm: float
    root_length_mm: float
    tip_position_mm: float
    middle_position_mm: float
    root_position_mm: float
    projected_middle_length_mm: float


def contact_line_geometry(path_length_mm, base_helix_angle_rad, face_width_mm,
                           base_pitch_mm, transverse_ratio, overlap_ratio, *,
                           for_bending=False):
    """由路径/齿宽椭圆式求三个接触线截长；齿根与点蚀评价中心符号相反。"""
    path = positive(path_length_mm, "啮合路径长度 mm")
    width = positive(face_width_mm, "齿宽 mm")
    pitch = positive(base_pitch_mm, "基节 mm")
    alpha_ratio = positive(transverse_ratio, "端面重合度")
    overlap = nonnegative(overlap_ratio, "轴向重合度")
    cosine, sine = math.cos(base_helix_angle_rad), math.sin(base_helix_angle_rad)
    if cosine <= 0:
        raise ValueError("基圆螺旋角余弦必须为正")
    center = (pitch - pitch * alpha_ratio / 2) * cosine * max(1 - overlap, 0)
    if not for_bending:
        center = -center
    positions = (center + pitch * cosine, center, center - pitch * cosine)
    ellipse = path**2 * cosine**2 + width**2 * sine**2
    lengths = tuple(width * path * math.sqrt(max(ellipse - 4 * x*x, 0)) / ellipse
                    for x in positions)
    return ContactLine(*lengths, *positions, lengths[1] * cosine)


def load_sharing_factor(resultant_ratio, overlap_ratio, base_pitch_mm,
                         base_helix_angle_rad, line: ContactLine):
    """求点蚀载荷分担因子 Z_LS；齿根 Y_LS 为它的平方。"""
    ratio = positive(resultant_ratio, "合成重合度")
    if ratio <= 2:
        return 1.0
    if overlap_ratio >= 1:
        return (1 + 2 * (1 - (2 / ratio)**1.5) * math.sqrt(1 - 4 / ratio**2))**-0.5
    scale = ratio * base_pitch_mm * math.cos(base_helix_angle_rad) / 2
    weights = [(1 - abs(position / scale)**1.5) * length * math.pi / 4
               for position, length in zip(
                   (line.tip_position_mm, line.middle_position_mm, line.root_position_mm),
                   (line.tip_length_mm, line.middle_length_mm, line.root_length_mm))]
    return math.sqrt(positive(weights[1] / positive(sum(weights), "接触线权重和"), "中心接触权重"))


def contact_ratio_factor(transverse_ratio, overlap_ratio):
    """B1 齿根重合度系数；高轴向重合度采用 0.625 下界。"""
    transverse_ratio = positive(transverse_ratio, "端面重合度")
    overlap_ratio = nonnegative(overlap_ratio, "轴向重合度")
    value = 0.25 + 0.75 / transverse_ratio
    if overlap_ratio >= 1:
        return 0.625
    return max(0.625, value - overlap_ratio * (0.75 / transverse_ratio - 0.375))


def bevel_root_factor(face_width_mm, projected_middle_length_mm):
    ratio = positive(projected_middle_length_mm, "投影接触线长度 mm") / positive(face_width_mm, "齿宽 mm")
    return (0.5 + 0.5 * ratio)**2 / ratio


def lubricant_film_factor(limit_pinion_mpa, limit_wheel_mpa, viscosity_40_cst,
                           pitch_speed_m_s, ra_pinion_um, ra_wheel_um,
                           virtual_center_mm, transverse_pressure_angle_rad,
                           base_helix_angle_rad, virtual_ratio):
    """方法 B 的黏度、速度与粗糙度联合影响系数。两轮 Ra 均须为正。"""
    limit = min(max(min(limit_pinion_mpa, limit_wheel_mpa), 850), 1200)
    viscosity = positive(viscosity_40_cst, "40°C 运动黏度 mm²/s")
    velocity = positive(pitch_speed_m_s, "节线速度 m/s")
    roughness = (positive(ra_pinion_um, "小轮 Ra µm") + positive(ra_wheel_um, "大轮 Ra µm")) * 3
    baseline_l = 0.83 + 0.08 * (limit - 850) / 350
    viscosity_factor = baseline_l + 4 * (1 - baseline_l) / (1.2 + 134 / viscosity)**2
    baseline_v = 0.85 + 0.08 * (limit - 850) / 350
    velocity_factor = baseline_v + 2 * (1 - baseline_v) / math.sqrt(0.8 + 32 / velocity)
    radius = (virtual_center_mm * math.sin(transverse_pressure_angle_rad)
              / math.cos(base_helix_angle_rad) * virtual_ratio / (1 + virtual_ratio)**2)
    relative_roughness = roughness * (10 / positive(radius, "联合曲率半径 mm"))**(1 / 3)
    roughness_factor = (3 / relative_roughness)**(0.12 + (1000 - limit) / 5000)
    return viscosity_factor * velocity_factor * roughness_factor


def hardness_ratio_factor(hardness_pinion_hb, hardness_wheel_hb):
    """硬度匹配影响；同硬度为 1，异硬度按较软轮硬度插值。"""
    a = positive(hardness_pinion_hb, "小轮 HB")
    b = positive(hardness_wheel_hb, "大轮 HB")
    return 1.0 if a == b else min(1.2, max(1.0, 1.2 - (min(a, b) - 130) / 1700))


@dataclass(frozen=True)
class RootSection:
    """危险截面的齿厚、圆角半径及展成参数；长度 mm。"""
    thickness_mm: float
    fillet_radius_mm: float
    auxiliary_angle_rad: float | None
    cutter_height_mm: float
    cutter_space_mm: float
    relative_height: float
    iterations: int


def root_section(cutter_height_coefficient, module_mm, thickness_shift,
                   profile_shift, pressure_angle_rad, cutter_tip_radius_mm,
                   protuberance_mm, virtual_teeth, *, generated=True):
    """B1 危险截面；展成分支显式限制固定点迭代并检查收敛。"""
    module = positive(module_mm, "中点法向模数 mm")
    cutter_radius = positive(cutter_tip_radius_mm, "刀具齿顶圆角半径 mm")
    protuberance = nonnegative(protuberance_mm, "刀具凸角量 mm")
    cutter_height = positive(cutter_height_coefficient, "刀具齿顶高系数") * module
    cutter_space = ((math.pi / 4 - thickness_shift) * module
                    - cutter_height * math.tan(pressure_angle_rad)
                    - (cutter_radius * (1 - math.sin(pressure_angle_rad)) - protuberance)
                    / math.cos(pressure_angle_rad))
    relative_height = cutter_radius / module - cutter_height / module + profile_shift
    if not generated:
        thickness = module * math.pi - 2 * cutter_space - 2 * cutter_radius * math.cos(math.pi / 6)
        return RootSection(positive(thickness, "成形危险截面齿厚 mm"), cutter_radius,
                           None, cutter_height, cutter_space, relative_height, 0)
    teeth = positive(virtual_teeth, "虚拟法向齿数")
    offset = 2 / teeth * (math.pi / 2 - cutter_space / module) - math.pi / 3
    angle = math.pi / 6
    for iteration in range(1, 1001):
        previous = angle
        angle = 2 * relative_height / teeth * math.tan(previous) - offset
        if not math.isfinite(angle) or abs(angle) >= math.pi / 2:
            raise ArithmeticError("展成危险截面迭代离开有效角度域")
        if abs(angle - previous) <= 1e-5:
            # 危险截面按满足停止判据的前一评价点计算，误差界由步长限定。
            angle = previous
            break
    else:
        raise ArithmeticError("展成危险截面固定点迭代未收敛")
    cosine = math.cos(angle)
    thickness = module * (teeth * math.sin(math.pi / 3 - angle)
                          + math.sqrt(3) * (relative_height / cosine - cutter_radius / module))
    radius = cutter_radius + 2 * module * relative_height**2 / (
        cosine * (teeth * cosine**2 - 2 * relative_height))
    return RootSection(positive(thickness, "展成危险截面齿厚 mm"),
                       positive(radius, "危险截面曲率半径 mm"), angle,
                       cutter_height, cutter_space, relative_height, iteration)


def root_form_factors(cutter_height_coefficient, module_mm, thickness_shift,
                        profile_shift, pressure_angle_rad, cutter_tip_radius_mm,
                        protuberance_mm, virtual_teeth, virtual_tip_diameter_mm,
                        virtual_base_diameter_mm, *, generated=True):
    """返回 B1 齿形系数 Y_Fa、应力修正 Y_Sa 及危险截面诊断。"""
    section = root_section(cutter_height_coefficient, module_mm, thickness_shift,
                           profile_shift, pressure_angle_rad, cutter_tip_radius_mm,
                           protuberance_mm, virtual_teeth, generated=generated)
    module, normal = module_mm, pressure_angle_rad
    if generated:
        diameter_ratio = virtual_base_diameter_mm / virtual_tip_diameter_mm
        if not 0 < diameter_ratio < 1:
            raise ValueError("虚拟齿顶圆必须位于基圆以外")
        tip_pressure = math.acos(diameter_ratio)
        half_tooth_angle = ((math.pi / 2 + 2 * (profile_shift * math.tan(normal) + thickness_shift))
                            / virtual_teeth + math.tan(normal) - normal
                            - math.tan(tip_pressure) + tip_pressure)
        load_angle = tip_pressure - half_tooth_angle
        root_angle = section.auxiliary_angle_rad
        lever = 0.5 * ((math.cos(half_tooth_angle) - math.sin(half_tooth_angle)
                       * math.tan(load_angle)) * virtual_tip_diameter_mm
                      - virtual_teeth * math.cos(math.pi / 3 - root_angle) * module
                      - section.relative_height / math.cos(root_angle) * module
                      + cutter_tip_radius_mm)
        shape = (6 * lever * math.cos(load_angle) * module
                 / (section.thickness_mm**2 * math.cos(normal)))
    else:
        lever = (section.cutter_height_mm - cutter_tip_radius_mm / 2 + module
                 - (math.pi / 4 + thickness_shift - math.tan(normal)) * module * math.tan(normal))
        shape = 6 * lever * module / section.thickness_mm**2
    positive(lever, "弯曲力臂 mm")
    slenderness = section.thickness_mm / lever
    notch = min(max(section.thickness_mm / (2 * section.fillet_radius_mm), 1), 7.9999)
    correction = (1.2 + 0.13 * slenderness) * notch**(1 / (1.21 + 2.3 / slenderness))
    return positive(shape, "齿形系数"), positive(correction, "应力修正系数"), section


def relative_sensitivity_factor(category, proof_strength_mpa, section: RootSection):
    """B1 相对缺口敏感性；仅接受已明确的钢材与证明强度表点。"""
    category = material_category(category)
    if category == "case_hardened":
        coefficient = 0.003
    elif category == "nitrided":
        coefficient = 0.1005
    elif category == "normalized":
        table = {300: 0.0833, 400: 0.0445}
        if proof_strength_mpa not in table:
            raise NotImplementedError("正火钢的证明强度仅实现 300/400 MPa 表点")
        coefficient = table[proof_strength_mpa]
    else:
        table = {500: 0.0281, 600: 0.0194, 800: 0.0064, 1000: 0.0014}
        if proof_strength_mpa not in table:
            raise NotImplementedError("所选钢材的证明强度仅实现 500/600/800/1000 MPa 表点")
        coefficient = table[proof_strength_mpa]
    notch = max(section.thickness_mm / (2 * section.fillet_radius_mm), 1)
    return (1 + math.sqrt(coefficient * 0.2 * (1 + 2 * notch))) / (1 + math.sqrt(coefficient * 1.2))


def relative_surface_factor(heat_treatment, roughness_rz_um):
    """B1 相对齿根表面状态；0=正火、1=渗碳、2=调质、3=氮化。"""
    rz = nonnegative(roughness_rz_um, "齿根 Rz µm")
    if rz > 40:
        raise NotImplementedError("齿根表面状态式仅覆盖 Rz≤40 µm")
    if heat_treatment in (1, 2):
        return 1.12 if rz < 1 else 1.674 - 0.529 * (rz + 1)**0.1
    if heat_treatment == 3:
        return 1.025 if rz < 1 else 4.299 - 3.259 * (rz + 1)**0.005
    if heat_treatment != 0:
        raise NotImplementedError("未定义的钢材表面热处理类别")
    return 1.07 if rz < 1 else 5.306 - 4.203 * (rz + 1)**0.01


def size_factor(module_mm, category):
    module = positive(module_mm, "中点法向模数 mm")
    category = material_category(category)
    if category in ("case_hardened", "nitrided", "induction_hardened"):
        return min(1.0, max(0.8, 1.05 - 0.01 * module))
    return min(1.0, max(0.85, 1.03 - 0.006 * module))


def application_factor(driver="uniform", driven="uniform"):
    """按主动/从动机械的冲击等级选择应用系数。"""
    levels = ("uniform", "light_shocks", "medium_shocks", "heavy_shocks")
    if driver not in levels or driven not in levels:
        raise ValueError("冲击等级必须为 uniform/light_shocks/medium_shocks/heavy_shocks")
    table = ((1, 1.25, 1.5, 1.75), (1.1, 1.35, 1.6, 1.85),
             (1.25, 1.5, 1.75, 2), (1.5, 1.75, 2, 2.25))
    return table[levels.index(driver)][levels.index(driven)]


def face_load_factor(effective_width_ratio, assembly="full_load_matched", support="both_supported"):
    """安装调整与悬臂支承产生的齿向载荷系数 K_Hβ（方法 C）。"""
    ratio = positive(effective_width_ratio, "有效齿宽比")
    if ratio > 1:
        raise ValueError("有效齿宽比不得超过 1")
    if assembly not in ("full_load_matched", "light_load_matched", "unmatched"):
        raise ValueError("未知的安装接触调整方式")
    supports = {"both_supported": 0, "one_cantilever": 1, "both_cantilever": 2}
    if support not in supports:
        raise ValueError("未知的齿轮支承方式")
    multiplier = 1.0 if assembly == "full_load_matched" else (
        (1.05, 1.1, 1.25) if assembly == "light_load_matched" else (1.2, 1.32, 1.5))[supports[support]]
    return 1.5 * multiplier * max(1, 0.85 / ratio)


def bending_face_load_factor(straight, spiral_angle_rad, mean_cone_distance_mm,
                              cutter_radius_mm, contact_face_factor):
    if straight:
        return positive(contact_face_factor, "K_Hβ")
    sine = math.sin(spiral_angle_rad)
    if not 0 < sine < 1:
        raise ValueError("螺旋锥齿齿根齿向系数要求 0<βm<π/2")
    exponent = 0.279 / math.log10(sine)
    ratio = positive(cutter_radius_mm, "刀盘半径 mm") / positive(mean_cone_distance_mm, "中点锥距 mm")
    correction = min(1.15, max(1.0, 0.211 * ratio**exponent + 0.789))
    return contact_face_factor / correction


def single_pitch_deviation(module_mm, reference_diameter_mm, grade):
    """按规定尺寸段的几何均值及等级系数求单节距偏差，返回 µm。"""
    module = positive(module_mm, "法向模数 mm")
    diameter = positive(reference_diameter_mm, "参考直径 mm")
    if not 0.5 <= module <= 70 or not 5 <= diameter <= 10000:
        raise NotImplementedError("齿距偏差分段仅覆盖 mn=0.5…70 mm、d=5…10000 mm")
    if not math.isfinite(grade) or not 1 <= grade <= 12:
        raise ValueError("精度等级必须在 1…12 内")
    diameter_limits = (5, 20, 50, 125, 280, 560, 1000, 1600, 2500, 4000, 6000, 8000, 10000)
    diameter_mid = next(math.sqrt(a * b) for a, b in zip(diameter_limits[:-1], diameter_limits[1:]) if diameter <= b)
    module_limits = (0.5, 2, 3.5, 6, 10, 16, 25, 40, 70)
    module_mid = math.sqrt(40 * 70)
    for low, high in zip(module_limits[:-1], module_limits[1:]):
        if module < high:
            module_mid = math.sqrt(low * high)
            break
    deviation = (0.3 * (module_mid + 0.4 * math.sqrt(diameter_mid)) + 4) * 2**((grade - 5) / 2)
    if deviation < 5:
        return round(deviation, 1)
    if deviation < 10:
        if int(deviation) == int(deviation + 0.75):
            return float(int(deviation))
        if int(deviation) != int(deviation + 0.25):
            return float(round(deviation))
        return int(deviation) + 0.5
    return float(round(deviation))


def running_in_allowance(category_pinion, limit_pinion_mpa, category_wheel,
                           limit_wheel_mpa, pitch_speed_m_s, pitch_deviation_um):
    """两轮跑合量取平均；调质钢与其余已支持钢材分开。"""
    def one(category, limit):
        if material_category(category) == "through_hardened":
            value = 160 / positive(limit, "接触疲劳极限 MPa") * pitch_deviation_um
            if pitch_speed_m_s > 10:
                return min(value, 6400 / limit)
            if pitch_speed_m_s > 5:
                return min(value, 12800 / limit)
            return value
        return min(0.075 * pitch_deviation_um, 3)
    return (one(category_pinion, limit_pinion_mpa) + one(category_wheel, limit_wheel_mpa)) / 2


def dynamic_factor_c(module_mm, mean_diameter_mm, outer_diameter_mm, grade, pitch_speed_m_s):
    """动载方法 C；对适用域外输入报错，避免返回零系数或失效错误码。"""
    if not 1.25 <= module_mm <= 50 or not 1 <= grade <= 9:
        raise NotImplementedError("动载方法 C 仅覆盖 mn=1.25…50 mm、精度等级 1…9")
    if grade < 5:
        return 1 + 0.02 * grade
    exponent = 0.25 * (grade - 5)**0.0667
    constant = 50 + 56 * (1 - exponent)
    equivalent_speed = nonnegative(pitch_speed_m_s, "节线速度 m/s") * positive(outer_diameter_mm, "外端直径 mm") / positive(mean_diameter_mm, "中点直径 mm")
    return (constant / (constant + math.sqrt(200 * equivalent_speed)))**(-exponent)


def dynamic_factor_b(speed_rpm, mean_diameter_mm, teeth, normal_pressure_angle_rad,
                       ratio, tangential_force_n, face_width_mm, effective_width_mm,
                       application, resultant_ratio, pitch_deviation_um, running_in_um):
    """动载方法 B 的亚临界、共振和超临界分段；返回 Kv 与转速比诊断。"""
    force = positive(tangential_force_n, "切向力 N")
    loading = min(1, force * application / effective_width_mm / 100)
    width_factor = min(1, effective_width_mm / (0.85 * face_width_mm))
    mesh = 20 * loading * width_factor
    inertia = 9.825e-7 * mean_diameter_mm**2 / math.cos(normal_pressure_angle_rad)**2 * ratio**2 / (1 + ratio**2)
    resonance_rpm = 30000 / (math.pi * teeth) * math.sqrt(mesh / inertia)
    speed_ratio = speed_rpm / resonance_rpm
    if resultant_ratio <= 2:
        a, b, c, d, e, f = 0.32, 0.34, 0.23, 0.9, 0.47, 0.47
    else:
        a, b, c, d, e, f = (0.32, 0.57 / (resultant_ratio - 0.3),
                            0.096 / (resultant_ratio - 1.56),
                            (0.57 - 0.05 * resultant_ratio) / (resultant_ratio - 1.44),
                            0.47, 0.12 / (resultant_ratio - 1.74))
    high_baseline = (0.75 if resultant_ratio <= 1.5 else 1 if resultant_ratio > 2.5
                     else 0.125 * math.sin(math.pi * (resultant_ratio - 2)) + 0.875)
    profile_load = face_width_mm * (pitch_deviation_um - running_in_um) / (force * application)
    resonance = profile_load * (a + b) + d + 1
    supercritical = profile_load * (e + f) + high_baseline
    if speed_ratio <= 0.75:
        kv = speed_ratio * (profile_load * 14 * loading * width_factor * (a + b) + c) + 1
    elif speed_ratio <= 1.25:
        kv = resonance
    elif speed_ratio < 1.5:
        kv = supercritical + (resonance - supercritical) * (1.5 - speed_ratio) / 0.25
    else:
        kv = supercritical
    return positive(kv, "动载系数 Kv"), {"resonance_speed_rpm": resonance_rpm, "speed_ratio": speed_ratio}


def transverse_load_factors_b(face_width_mm, transverse_ratio, overlap_ratio,
                               resultant_ratio, tangential_force_n, application,
                               dynamic, contact_face, sharing, root_contact_ratio,
                               running_in_um, pitch_deviation_um):
    load = positive(tangential_force_n * application * dynamic * contact_face / face_width_mm, "单位齿宽计算载荷 N/mm")
    error_term = 20 * (pitch_deviation_um - running_in_um) / load
    if resultant_ratio <= 2:
        contact = 0.5 * resultant_ratio * (0.9 + 0.4 * error_term)
    else:
        contact = 0.9 + 0.4 * math.sqrt(2 * (resultant_ratio - 1) / resultant_ratio) * error_term
    value = max(contact, 1)
    return (min(value, resultant_ratio / (transverse_ratio * sharing**2)),
            min(value, resultant_ratio / (transverse_ratio * root_contact_ratio)))


def transverse_load_factor_c(heat_treatment, straight, grade, tangential_force_n,
                              application, effective_width_mm, transverse_ratio):
    if heat_treatment not in (0, 1, 2, 3):
        raise NotImplementedError("未知的表面热处理类别")
    high_load = tangential_force_n * application / effective_width_mm >= 100
    if straight:
        if not high_load:
            return 1.2
        low = 7 if heat_treatment else 8
        return 1.0 if grade <= low else 1.1 if grade < low + 2 else 1.2
    if not high_load:
        return 1.4
    low = 6 if heat_treatment else 7
    return 1.0 if grade <= low else 1.1 if grade < low + 2 else 1.2 if grade < low + 3 else 1.4
