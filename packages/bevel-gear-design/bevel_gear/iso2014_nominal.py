"""ISO10300:2014 B1 的独立名义应力及方法B载荷修正链，全部采用SI。

覆盖零偏置、等齿宽、对称压力角、面铣、展成齿根；不含容量、寿命或
齿面接触求解。有效齿宽比例是测量/接触模式边界，默认.85仅为设计估计。
本模块保留独立的版本语义，不改变其他评级方法，也不代表标准认证。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .factors import positive, nonnegative, material_category
from .geometry import GeometryInput, GeometryResult, calculate_geometry


@dataclass(frozen=True)
class NominalMaterial:
    """名义链所需材料；E为Pa、密度kg/m³，接触应力数仅用于磨合量。"""
    young_modulus_pa: float
    poisson_ratio: float
    density_kg_m3: float
    heat_treatment: str = "case_hardened"
    contact_stress_number_pa: float | None = None

    def __post_init__(self):
        positive(self.young_modulus_pa, "弹性模量 Pa")
        positive(self.density_kg_m3, "密度 kg/m³")
        if not math.isfinite(self.poisson_ratio) or not -1 < self.poisson_ratio < .5:
            raise ValueError("泊松比必须在(-1,.5)内")
        category = material_category(self.heat_treatment)
        if category in ("normalized", "through_hardened"):
            if self.contact_stress_number_pa is None:
                raise ValueError("调质/正火材料须给定磨合量所需的允许接触应力数 Pa")
            positive(self.contact_stress_number_pa, "磨合量所需允许接触应力数 Pa")
        elif self.contact_stress_number_pa is not None:
            positive(self.contact_stress_number_pa, "允许接触应力数 Pa")


def _running_in(material, pitch_deviation_m, speed_m_s):
    """磨合量按各轮材料分别计算，最后由调用方取算术平均。"""
    category = material_category(material.heat_treatment)
    if category in ("case_hardened", "nitrided", "induction_hardened"):
        return min(.075 * pitch_deviation_m, 3e-6)
    limit_mpa = material.contact_stress_number_pa / 1e6
    result = 160 * pitch_deviation_m / limit_mpa
    if 5 < speed_m_s <= 10:
        result = min(result, 12800 / limit_mpa * 1e-6)
    elif speed_m_s > 10:
        result = min(result, 6400 / limit_mpa * 1e-6)
    return result


def _root_factors(module, height_factor, thickness_shift, profile_shift, alpha,
                  cutter_radius, protuberance, teeth, tip_diameter, base_diameter,
                  opposite_cutter_radius=None):
    """展成根部标量方程；G<0的正锐角根唯一，采用有界二分而非粗步固定点。"""
    height = height_factor * module
    space = ((math.pi / 4 - thickness_shift) * module - height * math.tan(alpha)
             - (cutter_radius * (1 - math.sin(alpha)) - protuberance) / math.cos(alpha))
    relative_height = cutter_radius / module - height_factor + profile_shift
    offset = 2 / teeth * (math.pi / 2 - space / module) - math.pi / 3
    if relative_height >= 0:
        raise NotImplementedError("当前展成齿根求根分支要求G<0；其他多根分支未验证")
    def residual(angle):
        return 2 * relative_height / teeth * math.tan(angle) - offset - angle
    lower, upper = 0., math.pi / 2 - 1e-12
    if residual(lower) <= 0 or residual(upper) >= 0:
        raise ValueError("展成齿根方程在已覆盖正锐角区间内无根")
    for iterations in range(1, 201):
        angle = (lower + upper) / 2
        if residual(angle) > 0:
            lower = angle
        else:
            upper = angle
        if upper - lower <= 1e-13:
            break
    else:
        raise ArithmeticError("展成齿根二分求解未收敛")
    angle = (lower + upper) / 2
    cosine = math.cos(angle)
    thickness = module * (teeth * math.sin(math.pi / 3 - angle)
                         + math.sqrt(3) * (relative_height / cosine - cutter_radius / module))
    radius = cutter_radius + 2 * relative_height**2 * module / (
        cosine * (teeth * cosine**2 - 2 * relative_height))
    if not 0 < base_diameter / tip_diameter < 1:
        raise ValueError("虚拟法向齿顶圆必须位于基圆以外")
    tip_pressure = math.acos(base_diameter / tip_diameter)
    half_angle = ((math.pi / 2 + 2 * (profile_shift * math.tan(alpha) + thickness_shift)) / teeth
                  + math.tan(alpha) - alpha - math.tan(tip_pressure) + tip_pressure)
    load_angle = tip_pressure - half_angle
    lever = .5 * ((math.cos(half_angle) - math.sin(half_angle) * math.tan(load_angle)) * tip_diameter
                  - module * teeth * math.cos(math.pi / 3 - angle)
                  - module * relative_height / cosine + cutter_radius)
    loaded_thickness = thickness
    opposite_thickness = thickness
    if opposite_cutter_radius is not None and opposite_cutter_radius != cutter_radius:
        opposite = _root_factors(module, height_factor, thickness_shift, profile_shift, alpha,
            opposite_cutter_radius, protuberance, teeth, tip_diameter, base_diameter)
        opposite_thickness = opposite['thickness_m']
        thickness = .5 * (loaded_thickness + opposite_thickness)
    # 两侧根形分别求解，危险截面弦齿厚取均值；载荷臂/受载角和根圆半径
    # 仍取当前工作侧。不能以一个平均刀尖半径先替代两侧几何再求根。
    for value, name in ((thickness, "危险截面齿厚 m"), (radius, "危险截面圆角 m"), (lever, "弯曲力臂 m")):
        positive(value, name)
    shape = 6 * lever * module * math.cos(load_angle) / (thickness**2 * math.cos(alpha))
    slenderness, notch = thickness / lever, thickness / (2 * radius)
    correction = (1.2 + .13 * slenderness) * notch**(1 / (1.21 + 2.3 / slenderness))
    return {"Y_Fa": positive(shape, "齿形系数"), "Y_Sa": positive(correction, "应力修正系数"),
            "thickness_m": thickness, "fillet_radius_m": radius, "lever_arm_m": lever,
            "theta_rad": angle, "theta_residual": residual(angle), "iterations": iterations,
            "loaded_side_thickness_m": loaded_thickness, "opposite_side_thickness_m": opposite_thickness,
            "loaded_cutter_tip_radius_m": cutter_radius,
            "opposite_cutter_tip_radius_m": cutter_radius if opposite_cutter_radius is None else opposite_cutter_radius,
            "load_pressure_angle_rad": load_angle,
            "E_m": space, "G": relative_height, "H": offset, "L_a": slenderness, "q_s": notch}


def _contact_lines(path, base_pitch, base_helix, width, effective_width, alpha_ratio,
                   overlap_ratio, crown_exponent, *, bending):
    """顺序固定为齿顶/中间/齿根；返回三线的实际长度、位置和权重面积。"""
    cb, sb, tb = math.cos(base_helix), math.sin(base_helix), math.tan(base_helix)
    center = -(base_pitch - .5 * base_pitch * alpha_ratio) * cb * max(1 - overlap_ratio, 0)
    if bending:
        center = -center
    positions = (center + base_pitch * cb, center, center - base_pitch * cb)
    # 相交轴、零偏置模型的投影辅助角为零。
    maximum = .5 * (path + effective_width * tb) * cb
    lengths, corrections, areas = [], [], []
    for position in positions:
        if base_helix == 0:
            theoretical = effective_width
        elif abs(position) > maximum:
            theoretical = 0.
        else:
            intercept = position * cb + tb * (position * sb + effective_width / 2)
            start = min(effective_width, max(0., (intercept + path / 2) / tb))
            end = min(effective_width, max(0., (intercept - path / 2) / tb))
            theoretical = math.hypot(start - end, (-start * tb + intercept) - (-end * tb + intercept))
        if abs(position) > maximum:
            correction = 1.
        else:
            correction = math.sqrt(max(0., 1 - (position / maximum)**2)) * (1 - math.sqrt(effective_width / width))
        length = theoretical * (1 - correction)
        area = .25 * math.pi * (1 - (abs(position) / maximum)**crown_exponent) * length
        lengths.append(length)
        corrections.append(correction)
        areas.append(area)
    total = positive(sum(areas), "三条接触线权重面积之和 m")
    return {"order": ("tip", "middle", "root"), "position_m": tuple(positions),
            "length_m": tuple(lengths), "length_correction": tuple(corrections),
            "weighted_area_m": tuple(areas), "maximum_distance_m": maximum,
            "middle_load_fraction": positive(areas[1] / total, "中线载荷分担比例")}


def evaluate_iso2014_nominal(geometry: GeometryInput | GeometryResult, *, torque_nm: float,
        pinion_speed_rad_s: float, pinion: NominalMaterial, wheel: NominalMaterial,
        pinion_pitch_deviation_m: float, wheel_pitch_deviation_m: float,
        pinion_cutter_tip_radius_m: float, wheel_cutter_tip_radius_m: float,
        pinion_opposite_cutter_tip_radius_m: float | None = None,
        wheel_opposite_cutter_tip_radius_m: float | None = None,
        pinion_protuberance_m: float = 0., wheel_protuberance_m: float = 0.,
        application_factor: float = 1., face_load_factor: float = 1.,
        effective_face_width_ratio: float = .85, profile_crowning: str = "low",
        dynamic_method: str = "B", manufacturing_method: str = "face_milling",
        wheel_generated: bool = True, working_flank: str = "drive") -> dict:
    """计算H0/F0与载荷修正后的H/F，返回Pa及全部具名中间量。

    齿距偏差为每轮实际最大单齿距偏差m，不把不同质量标准的等级整数混用。
    KA/KHβ是工况和安装边界；有效齿宽比例来自测量/接触模式估计，不是LTCA。
    方法C、准双曲、非展成/面滚切不在本入口覆盖范围。正扭矩且非零转速；
    正负转向只通过角速度绝对值进入动态公式，不推导左右面非对称制造差别。
    """
    if dynamic_method != "B":
        raise NotImplementedError("2014名义链当前仅实现动态方法B；方法C未开放")
    if manufacturing_method != "face_milling" or not wheel_generated:
        raise NotImplementedError("当前仅覆盖面铣及两轮展成齿根")
    if working_flank not in ('drive', 'coast'):
        raise ValueError('标准受载面须为drive或coast；四刀尖须按该工作面显式映射')
    g = calculate_geometry(geometry) if isinstance(geometry, GeometryInput) else geometry
    if not isinstance(g, GeometryResult):
        raise TypeError("geometry须为GeometryInput或GeometryResult")
    i = g.inputs
    if i.offset_m != 0 or g.diagnostics['tooth_trace'] not in ('straight', 'circular_arc'):
        raise NotImplementedError("当前仅覆盖零偏置的直齿或圆弧齿线")
    if profile_crowning not in ("low", "high"):
        raise ValueError("齿廓鼓形类别须为low或high")
    if not math.isfinite(effective_face_width_ratio) or not 0 < effective_face_width_ratio <= 1:
        raise ValueError("有效齿宽比例须在(0,1]内")
    torque = positive(torque_nm, "小轮扭矩 N·m")
    speed = abs(float(pinion_speed_rad_s))
    positive(speed, "小轮角速度 rad/s")
    ka, khb = positive(application_factor, "KA"), positive(face_load_factor, "KHβ")
    fp = max(nonnegative(pinion_pitch_deviation_m, "小轮单齿距偏差 m"),
             nonnegative(wheel_pitch_deviation_m, "大轮单齿距偏差 m"))
    radii = [positive(pinion_cutter_tip_radius_m, "小轮刀尖圆角 m"), positive(wheel_cutter_tip_radius_m, "大轮刀尖圆角 m")]
    opposite_radii = [loaded if other is None else positive(other, '非受载侧刀尖圆角m') for loaded, other in
                      zip(radii, (pinion_opposite_cutter_tip_radius_m, wheel_opposite_cutter_tip_radius_m))]
    protuberances = [nonnegative(pinion_protuberance_m, "小轮凸角量 m"), nonnegative(wheel_protuberance_m, "大轮凸角量 m")]
    beta, alpha = i.spiral_angle_rad, i.pressure_angle_rad
    module, width = g.mean_normal_module_m, i.face_width_m
    effective_width = effective_face_width_ratio * width
    transverse_alpha = math.atan(math.tan(alpha) / math.cos(beta))
    base_helix = math.asin(math.sin(beta) * math.cos(alpha))
    pitch, tip, base, teeth, normal_tip, normal_base = [], [], [], [], [], []
    dimensions = (g.pinion, g.wheel)
    for gear in dimensions:
        dv = gear.mean_pitch_diameter_m / math.cos(gear.pitch_cone_angle_rad)
        da = dv + 2 * gear.mean_addendum_m
        db = dv * math.cos(transverse_alpha)
        pitch.append(dv); tip.append(da); base.append(db)
        teeth.append(dv / (module * math.cos(base_helix)**2))
        normal_tip.append(dv / math.cos(base_helix)**2 + da - dv)
        normal_base.append(dv / math.cos(base_helix)**2 * math.cos(alpha))
    base_pitch = math.pi * module * math.cos(transverse_alpha) / math.cos(beta)
    path = sum(.5 * (math.sqrt(positive(da**2 - db**2, "齿顶至基圆平方差 m²"))
                     - dv * math.sin(transverse_alpha)) for dv, da, db in zip(pitch, tip, base))
    ea = path / base_pitch
    if ea <= 1:
        raise ValueError("2014名义链要求虚拟端面重合度大于1")
    eb = math.sin(beta) * effective_width / (math.pi * module)
    ev = ea + eb
    exponent = 3. if profile_crowning == "low" else 1.5
    contact = _contact_lines(path, base_pitch, base_helix, width, effective_width, ea, eb, exponent, bending=False)
    bending = _contact_lines(path, base_pitch, base_helix, width, effective_width, ea, eb, exponent, bending=True)
    fraction = contact["middle_load_fraction"]
    zls, yls = math.sqrt(fraction), fraction
    if eb >= 1:
        factors_k = (ea, ea)
    elif eb > 0:
        factors_k = (2 + (ea - 2) * eb, 2 * ea - 2 + (2 - ea) * eb)
    else:
        factors_k = (2., 2 * (ea - 1))
    mid_terms = [math.sqrt((da / db)**2 - 1) - k * math.pi / (dv * math.cos(beta) / module)
                 for dv, da, db, k in zip(pitch, tip, base, factors_k)]
    if min(mid_terms) <= 0:
        raise ValueError("中区系数超出本名义模型的实数有效域")
    zmb = math.tan(transverse_alpha) / math.sqrt(mid_terms[0] * mid_terms[1])
    psi = math.atan(math.tan(beta) * math.sin(alpha))
    inverse_back_radius = sum(1 / (g.mean_cone_distance_m * math.tan(x.pitch_cone_angle_rad)) for x in dimensions)
    radius_normal = math.sin(alpha) / (math.cos(beta)**2 * inverse_back_radius)
    radius_relative = abs(radius_normal) * math.cos(psi)**2
    tangential = 2 * torque / g.pinion.mean_pitch_diameter_m
    normal_force = tangential / (math.cos(alpha) * math.cos(beta))
    ze = math.sqrt(1 / (math.pi * ((1 - pinion.poisson_ratio**2) / pinion.young_modulus_pa
                                  + (1 - wheel.poisson_ratio**2) / wheel.young_modulus_pa)))
    zk = .85
    h0 = zmb * zls * ze * zk * math.sqrt(normal_force / (contact["length_m"][1] * radius_relative))
    ye = .625 if eb > 1 else max(.625, .25 + .75 / ea - eb * (.75 / ea - .375))
    ba = width / math.cos(beta)
    height = sum(x.mean_addendum_m + x.mean_dedendum_m for x in dimensions) / 2
    ratio = ba / height
    a = -.0182 * ratio**2 + .4736 * ratio - .32
    b = -.0032 * ratio**2 + .0526 * ratio + .712
    c = -.005 * ratio**2 + .085 * ratio + .54
    positive(c, "YBS分母辅助系数")
    lbb = bending["length_m"][1] * math.cos(base_helix) / math.cos(beta)
    ybs = 1 + a / c * (lbb / ba - 1.05 * b)**2
    # 中点法向总侧隙j_n由两轮各分一半，转换为弧齿厚减薄量j_n/(2cosα)。
    # s_mn=πm_n/2+2m_n(x_sm+x_hm tanα)，故两轮x_sm均减j_n/(4m_n cosα)。
    # 大轮原有的理论齿厚变位为-xs；侧隙减薄不随主从轮反号。
    backlash_shift = i.backlash_m / (4 * module * math.cos(alpha))
    thickness_shifts = tuple(sign * i.tooth_thickness_modification_coefficient - backlash_shift for sign in (1, -1))
    roots = [_root_factors(module, i.dedendum_coefficient, shift,
                           sign * i.pinion_profile_shift_coefficient, alpha, radius, protuberance,
                           z, da, db, opposite_cutter_radius=other)
             for sign, shift, radius, other, protuberance, z, da, db in zip((1, -1), thickness_shifts, radii, opposite_radii, protuberances, teeth, normal_tip, normal_base)]
    f0 = [tangential / (width * module) * r["Y_Fa"] * r["Y_Sa"] * ye * ybs * yls for r in roots]
    pitch_speed = g.pinion.mean_pitch_diameter_m * speed / 2
    running = [_running_in(material, fp, pitch_speed) for material in (pinion, wheel)]
    allowance = sum(running) / 2
    effective_pitch = fp - allowance
    cf = min(1., tangential * ka / (effective_width * 100000.))
    mesh_stiffness, single_stiffness = 20e9 * cf, 14e9 * cf
    masses = [material.density_kg_m3 * math.pi / 8 * (gear.mean_pitch_diameter_m / math.cos(alpha))**2
              for material, gear in zip((pinion, wheel), dimensions)]
    reduced_mass = masses[0] * masses[1] / sum(masses)
    resonance = math.sqrt(mesh_stiffness / reduced_mass) / i.pinion_teeth
    relative_speed = speed / resonance
    cv1, cv5 = .32, .47
    cv2, cv3, cv4, cv6 = ((.34, .23, .9, .47) if ev <= 2 else
                         (.57 / (ev - .3), .096 / (ev - 1.56), (.57 - .05 * ev) / (ev - 1.44), .12 / (ev - 1.74)))
    cv7 = .75 if ev <= 1.5 else (.125 * math.sin(math.pi * (ev - 2)) + .875 if ev <= 2.5 else 1.)
    bp = width * effective_pitch * single_stiffness / (ka * tangential)
    subcritical = 1 + (bp * (cv1 + cv2) + cv3) * relative_speed
    main_resonance = 1 + bp * (cv1 + cv2) + cv4
    supercritical = 1 + bp * (cv5 + cv6) + cv7
    intermediate = supercritical + (main_resonance - supercritical) / .25 * (1.5 - relative_speed)
    if relative_speed <= .75:
        sector, kv = "subcritical", subcritical
    elif relative_speed <= 1.25:
        sector, kv = "main_resonance", main_resonance
    elif relative_speed <= 1.5:
        sector, kv = "intermediate", intermediate
    else:
        sector, kv = "supercritical", supercritical
    kv = max(1., kv)
    ratio_load = mesh_stiffness * effective_pitch * width / (tangential * ka * kv * khb)
    raw_transverse = (ev / 2 * (.9 + .4 * ratio_load) if ev <= 2 else
                      .9 + .4 * math.sqrt(2 * (ev - 1) / ev) * ratio_load)
    upper_h, upper_f = ev / (ea * zls**2), ev / (ea * yls)
    kha = 1. if raw_transverse < 1 else min(raw_transverse, upper_h)
    kfa = 1. if raw_transverse < 1 else min(raw_transverse, upper_f)
    if beta == 0:
        kf0 = 1.
    else:
        cutter = positive(i.cutter_radius_m, "面铣刀盘半径 m")
        kf0 = 1. if cutter > g.mean_cone_distance_m else min(1.15, max(1.,
            .211 * (cutter / g.mean_cone_distance_m)**(.279 / math.log10(math.sin(beta))) + .789))
    kfb = khb / kf0
    stress = {"H0": h0, "F01": f0[0], "F02": f0[1],
              "H": h0 * math.sqrt(ka * kv * khb * kha),
              "F1": f0[0] * ka * kv * kfb * kfa, "F2": f0[1] * ka * kv * kfb * kfa}
    if not all(math.isfinite(value) and value > 0 for value in stress.values()):
        raise ArithmeticError("名义链出现非有限或非正应力，工况未通过")
    return {"stress_pa": stress,
            "coefficients": {"Z_K": zk, "Z_E_sqrt_pa": ze, "Z_MB": zmb, "Z_LS": zls, "Y_LS": yls,
                "Y_epsilon": ye, "Y_BS": ybs, "a_BS": a, "b_BS": b, "c_BS": c,
                "K_A": ka, "K_v": kv, "K_Hbeta": khb, "K_Fbeta": kfb, "K_F0": kf0,
                "K_Halpha": kha, "K_Falpha": kfa, "K_alpha_raw": raw_transverse,
                "K_Halpha_upper": upper_h, "K_Falpha_upper": upper_f},
            "geometry": {"module_m": module, "virtual_diameter_m": tuple(pitch), "virtual_tip_m": tuple(tip),
                "virtual_base_m": tuple(base), "virtual_normal_teeth": tuple(teeth),
                "base_helix_angle_rad": base_helix, "transverse_pressure_angle_rad": transverse_alpha,
                "path_length_m": path, "base_pitch_m": base_pitch, "face_width_m": width,
                "effective_face_width_m": effective_width, "epsilon_alpha": ea, "epsilon_beta": eb,
                "epsilon_virtual": ev, "relative_curvature_m": radius_relative,
                "normal_curvature_m": radius_normal if working_flank == 'drive' else -radius_normal,
                "developed_width_m": ba, "mean_tooth_depth_m": height, "covered_width_m": lbb},
            "contact_lines": {"contact": contact, "root": bending}, "root_sections": tuple(roots),
            "dynamic": {"max_pitch_deviation_m": fp, "running_in_per_gear_m": tuple(running),
                "running_in_m": allowance, "effective_pitch_deviation_m": effective_pitch,
                "mass_per_width_kg_m": tuple(masses), "reduced_mass_kg_m": reduced_mass,
                "stiffness_correction": cf, "mesh_stiffness_n_m2": mesh_stiffness,
                "single_stiffness_n_m2": single_stiffness, "resonance_rad_s": resonance,
                "relative_speed": relative_speed, "sector": sector, "pitch_influence": bp,
                "cv": (cv1, cv2, cv3, cv4, cv5, cv6, cv7),
                "subcritical": subcritical, "main_resonance": main_resonance,
                "intermediate": intermediate, "supercritical": supercritical},
            "forces_n": {"tangential": tangential, "normal": normal_force},
            "diagnostics": {"standard_scope": "ISO10300:2014 B1名义应力与方法B载荷链",
                "working_flank": working_flank,
                "normal_curvature_definition": "normal_curvature_m为带工作面法向符号的法截面相对曲率半径，单位m；relative_curvature_m为Hertz式使用的正相对曲率半径",
                "mean_normal_backlash_m": i.backlash_m,
                "root_thickness_shifts_with_backlash": thickness_shifts,
                "backlash_root_model": "两轮分别按中点法向总侧隙减薄，不将大轮减薄反号；危险截面仍为标准等效展成根",
                "root_solver": "G<0正锐角唯一根的有界二分，角度区间公差1e-13rad",
                "effective_width_boundary": "显式比例边界；默认.85为设计估计，不是计算得到的接触图",
                "capacity_life_included": False, "surface_contact_solver_included": False}}
