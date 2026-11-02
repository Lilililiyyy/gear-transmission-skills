"""新写的零偏置锥齿轮宏观几何；长度m，角度rad。

采用共同锥顶的外啮合节锥、对称法向压力角及中点齿高/齿厚参数。
圆弧齿线只描述理想刀盘轨迹；几何重合度及根切标记是设计近似，
不等于三维齿面接触、制造验收或标准符合性认证。
"""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class GeometryInput:
    """独立几何输入。

    backlash_m 是中点沿齿面法向测量的总侧隙，按两轮各半分摊。
    profile_shift只改变中点齿高和齿厚分配，不改变共同节锥。
    tooth_trace='auto'按中点螺旋角选直线或圆弧；零度圆弧须显式指定。
    齿高锥度可选standard、uniform、constant_slot、limited_slot。
    根部圆角半径是显式工艺输入，本模块不推定缺失值。
    """
    pinion_teeth: int
    wheel_teeth: int
    wheel_outer_pitch_diameter_m: float
    face_width_m: float
    shaft_angle_rad: float = math.pi / 2
    spiral_angle_rad: float = 0.0
    pressure_angle_rad: float = math.radians(20)
    addendum_coefficient: float = 1.0
    dedendum_coefficient: float = 1.25
    pinion_profile_shift_coefficient: float = 0.0
    tooth_thickness_modification_coefficient: float = 0.0
    backlash_m: float = 0.0
    cutter_radius_m: float | None = None
    depth_taper: str = "standard"
    tooth_trace: str = "auto"
    offset_m: float = 0.0
    root_fillet_radius_pinion_m: float | None = None
    root_fillet_radius_wheel_m: float | None = None
    cutter_group_count: int = 0
    wheel_addendum_taper_angle_rad: float | None = None
    wheel_dedendum_taper_angle_rad: float | None = None

    def __post_init__(self) -> None:
        for name in ("pinion_teeth", "wheel_teeth"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 3:
                raise ValueError(f"{name}必须为不小于3的整数")
        optional = ("cutter_radius_m", "root_fillet_radius_pinion_m", "root_fillet_radius_wheel_m")
        optional_angles = ("wheel_addendum_taper_angle_rad", "wheel_dedendum_taper_angle_rad")
        for name, value in self.__dict__.items():
            if name in ("depth_taper", "tooth_trace", "pinion_teeth", "wheel_teeth") or (name in optional+optional_angles and value is None):
                continue
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{name}必须为有限数")
        if self.offset_m != 0:
            raise NotImplementedError("本项目仅支持零偏置直齿和圆弧齿")
        if self.wheel_outer_pitch_diameter_m <= 0 or self.face_width_m <= 0:
            raise ValueError("外端节圆直径和齿宽必须为正，单位m")
        if not 0 < self.shaft_angle_rad < math.pi:
            raise ValueError("轴交角必须在0和π之间，单位rad")
        if not 0 <= self.spiral_angle_rad < math.pi / 2:
            raise ValueError("中点螺旋角幅值必须在[0,π/2)内，旋向在啮合力接口单独指定")
        if not 0 < self.pressure_angle_rad < math.pi / 2:
            raise ValueError("对称法向压力角必须在0和π/2之间")
        if self.addendum_coefficient <= 0 or self.dedendum_coefficient <= self.addendum_coefficient:
            raise ValueError("齿顶高系数须为正，齿根高系数须大于齿顶高系数以保持正顶隙")
        if self.backlash_m < 0:
            raise ValueError("中点法向总侧隙不能为负，单位m")
        for name in optional:
            value = getattr(self, name)
            if value is not None and value <= 0:
                raise ValueError(f"{name}必须为正长度，单位m")
        if self.depth_taper not in ("standard", "uniform", "constant_slot", "limited_slot", "custom"):
            raise ValueError("未知齿高锥度；可选standard、uniform、constant_slot、limited_slot或custom")
        angles = tuple(getattr(self,key) for key in optional_angles)
        if self.depth_taper == "custom":
            if any(value is None or not 0 <= value < math.pi/2 for value in angles):
                raise ValueError("custom锥度须显式给齐大轮齿顶/齿根锥度角，范围[0,π/2)")
        elif any(value is not None for value in angles):
            raise ValueError("给定锥度角时须选择custom，不能静默忽略输入")
        if self.tooth_trace not in ("auto", "straight", "circular_arc"):
            raise ValueError("齿线类型须为auto、straight或circular_arc")
        if self.tooth_trace == "straight" and self.spiral_angle_rad != 0:
            raise ValueError("直齿线不能同时指定非零中点螺旋角")
        if isinstance(self.cutter_group_count,bool) or not isinstance(self.cutter_group_count,int) or self.cutter_group_count < 0:
            raise ValueError("刀盘分组数须为非负整数")
        if self.cutter_group_count:
            raise ValueError("本项目直齿和圆弧齿接口不使用刀盘分组数，协议字段须为0")


@dataclass(frozen=True)
class GearDimensions:
    """单轮具名尺寸，轴向坐标从共同锥顶向齿轮外端取正。"""
    teeth: int
    pitch_cone_angle_rad: float
    tip_cone_angle_rad: float
    root_cone_angle_rad: float
    addendum_taper_angle_rad: float
    dedendum_taper_angle_rad: float
    outer_pitch_diameter_m: float
    mean_pitch_diameter_m: float
    inner_pitch_diameter_m: float
    outer_tip_diameter_m: float
    outer_root_diameter_m: float
    inner_tip_diameter_m: float
    inner_root_diameter_m: float
    outer_addendum_m: float
    mean_addendum_m: float
    inner_addendum_m: float
    outer_dedendum_m: float
    mean_dedendum_m: float
    inner_dedendum_m: float
    mean_normal_tooth_thickness_m: float
    mean_transverse_tooth_thickness_m: float
    mean_normal_chordal_tooth_thickness_m: float
    mean_normal_chordal_addendum_m: float
    virtual_normal_pitch_diameter_m: float
    outer_spiral_angle_rad: float
    mean_spiral_angle_rad: float
    inner_spiral_angle_rad: float
    outer_axial_pitch_coordinate_m: float
    mean_axial_pitch_coordinate_m: float
    inner_axial_pitch_coordinate_m: float
    outer_tip_axial_coordinate_m: float
    inner_tip_axial_coordinate_m: float
    outer_root_axial_coordinate_m: float
    inner_root_axial_coordinate_m: float
    root_fillet_radius_m: float | None


@dataclass(frozen=True)
class GeometryResult:
    """成对宏观几何；resultant为两重合度的欧氏合成，不能当LTCA。"""
    inputs: GeometryInput
    pinion: GearDimensions
    wheel: GearDimensions
    ratio: float
    outer_cone_distance_m: float
    mean_cone_distance_m: float
    inner_cone_distance_m: float
    outer_transverse_module_m: float
    mean_transverse_module_m: float
    mean_normal_module_m: float
    normal_clearance_m: float
    mean_whole_tooth_depth_m: float
    mean_working_tooth_depth_m: float
    mean_normal_backlash_m: float
    outer_normal_backlash_m: float
    mean_transverse_backlash_m: float
    outer_transverse_backlash_m: float
    transverse_contact_ratio: float
    overlap_contact_ratio: float
    resultant_contact_ratio: float
    diagnostics: dict

    @property
    def profile_contact_ratio(self) -> float:
        """当量齿轮法截面近似的齿廓重合度。"""
        return self.transverse_contact_ratio


def _asin_domain(value: float, description: str) -> float:
    if abs(value) > 1 + 2e-13:
        raise ValueError(f"{description}超出反正弦实数域，请检查齿宽、刀盘半径或螺旋角")
    return math.asin(min(1.0, max(-1.0, value)))


def _taper_angles(inputs: GeometryInput, mean_distance: float, outer_distance: float,
                  normal_module: float, outer_module: float) -> tuple[float, float]:
    """返回大轮齿顶/齿根锥度；零偏置时小轮取配对的根/顶锥度。"""
    ha, hf = inputs.addendum_coefficient, inputs.dedendum_coefficient
    x = inputs.pinion_profile_shift_coefficient
    hf1, hf2 = normal_module * (hf - x), normal_module * (hf + x)
    standard_top, standard_root = math.atan(hf1 / mean_distance), math.atan(hf2 / mean_distance)
    if inputs.depth_taper == "custom":
        return inputs.wheel_addendum_taper_angle_rad, inputs.wheel_dedendum_taper_angle_rad
    if inputs.depth_taper == "uniform":
        return 0.0, 0.0
    if inputs.depth_taper == "standard":
        return standard_top, standard_root
    if inputs.cutter_radius_m is None:
        raise ValueError("齿槽锥度法必须显式给定刀盘半径cutter_radius_m")
    beta, alpha = inputs.spiral_angle_rad, inputs.pressure_angle_rad
    total = outer_module / (outer_distance * math.tan(alpha) * math.cos(beta)) * (
        1 - mean_distance * math.sin(beta) / inputs.cutter_radius_m
    ) * math.pi / 2
    if total <= 0:
        raise ValueError("给定刀盘半径导致非正齿高锥度，当前工艺近似不适用")
    if inputs.depth_taper == "limited_slot":
        total = min(total, 1.3 * (standard_top + standard_root))
    top = total * (ha - x) / (2 * ha)
    return top, total - top


def _gear_dimensions(inputs: GeometryInput, teeth: int, pitch_angle: float,
                     top_taper: float, root_taper: float, mean_addendum: float,
                     mean_dedendum: float, normal_thickness: float,
                     outer_distance: float, mean_distance: float, inner_distance: float,
                     beta_outer: float, beta_inner: float, fillet: float | None) -> GearDimensions:
    half_width = inputs.face_width_m / 2
    ha_outer = mean_addendum + half_width * math.tan(top_taper)
    ha_inner = mean_addendum - half_width * math.tan(top_taper)
    hf_outer = mean_dedendum + half_width * math.tan(root_taper)
    hf_inner = mean_dedendum - half_width * math.tan(root_taper)
    if min(ha_outer, ha_inner, hf_outer, hf_inner) <= 0:
        raise ValueError("齿高锥度使内/外端齿顶高或齿根高非正，请减小齿宽或变位")
    tip_angle, root_angle = pitch_angle + top_taper, pitch_angle - root_taper
    if not 0 < root_angle <= pitch_angle <= tip_angle < math.pi / 2:
        raise ValueError("齿顶/齿根锥角超出有效外啮合锥体范围，当前尺寸不可用")
    sine, cosine = math.sin(pitch_angle), math.cos(pitch_angle)
    dp_outer, dp_mean, dp_inner = (2 * radius * sine for radius in (outer_distance, mean_distance, inner_distance))
    tip_outer, tip_inner = dp_outer + 2 * ha_outer * cosine, dp_inner + 2 * ha_inner * cosine
    root_outer, root_inner = dp_outer - 2 * hf_outer * cosine, dp_inner - 2 * hf_inner * cosine
    if min(root_outer, root_inner) <= 0:
        raise ValueError("实际根圆直径非正，请检查齿数、齿高及齿宽")
    beta, alpha = inputs.spiral_angle_rad, inputs.pressure_angle_rad
    normal_diameter = dp_mean / ((1 - math.sin(beta)**2 * math.cos(alpha)**2) * math.cos(beta) * cosine)
    chord = normal_diameter * math.sin(normal_thickness / normal_diameter)
    chord_addendum = mean_addendum + normal_diameter * cosine / 2 * (1 - math.cos(normal_thickness / normal_diameter))
    return GearDimensions(
        teeth, pitch_angle, tip_angle, root_angle, top_taper, root_taper,
        dp_outer, dp_mean, dp_inner, tip_outer, root_outer, tip_inner, root_inner,
        ha_outer, mean_addendum, ha_inner, hf_outer, mean_dedendum, hf_inner,
        normal_thickness, normal_thickness / math.cos(beta), chord, chord_addendum,
        normal_diameter, beta_outer, beta, beta_inner,
        outer_distance * cosine, mean_distance * cosine, inner_distance * cosine,
        outer_distance * cosine - ha_outer * sine, inner_distance * cosine - ha_inner * sine,
        outer_distance * cosine + hf_outer * sine, inner_distance * cosine + hf_inner * sine,
        fillet,
    )


def calculate_geometry(inputs: GeometryInput) -> GeometryResult:
    """计算零偏置宏观尺寸，不依赖外部组件或商业运行时。"""
    if not isinstance(inputs, GeometryInput):
        raise TypeError("几何入口需要GeometryInput参数")
    u = inputs.wheel_teeth / inputs.pinion_teeth
    sigma, beta, alpha = inputs.shaft_angle_rad, inputs.spiral_angle_rad, inputs.pressure_angle_rad
    delta_p = math.atan2(math.sin(sigma), u + math.cos(sigma))
    delta_w = sigma - delta_p
    if min(delta_p, delta_w) <= 1e-10 or max(delta_p, delta_w) >= math.pi / 2 - 1e-10:
        raise ValueError("两节锥角必须严格处于0和π/2之间，冠轮/反锥及退化节锥不在当前模型内")
    re = inputs.wheel_outer_pitch_diameter_m / (2 * math.sin(delta_w))
    ri, rm = re - inputs.face_width_m, re - inputs.face_width_m / 2
    if ri <= 0:
        raise ValueError("齿宽必须小于外锥距，内锥距不能为零或负")
    me = inputs.wheel_outer_pitch_diameter_m / inputs.wheel_teeth
    mt = 2 * rm * math.sin(delta_w) / inputs.wheel_teeth
    mn = mt * math.cos(beta)
    ha, hf, x = inputs.addendum_coefficient, inputs.dedendum_coefficient, inputs.pinion_profile_shift_coefficient
    if min(ha + x, ha - x, hf - x, hf + x) <= 0:
        raise ValueError("变位使某轮中点齿顶高或齿根高非正")
    hm_ap, hm_aw, hm_fp, hm_fw = mn * (ha + x), mn * (ha - x), mn * (hf - x), mn * (hf + x)
    theta_aw, theta_fw = _taper_angles(inputs, rm, re, mn, me)
    trace = inputs.tooth_trace
    if trace == "auto":
        trace = "straight" if beta == 0 else "circular_arc"
    if trace == "circular_arc":
        if inputs.cutter_radius_m is None:
            raise ValueError("圆弧齿线须显式提供刀盘半径cutter_radius_m，零螺旋角圆弧也不例外")
        rc = inputs.cutter_radius_m
        beta_outer = _asin_domain((2 * rm * rc * math.sin(beta) - rm**2 + re**2) / (2 * re * rc), "外端螺旋角")
        beta_inner = _asin_domain((2 * rm * rc * math.sin(beta) - rm**2 + ri**2) / (2 * ri * rc), "内端螺旋角")
    else:
        beta_outer = beta_inner = 0.0
    xs, backlash = inputs.tooth_thickness_modification_coefficient, inputs.backlash_m
    thickness_p = math.pi * mn / 2 + 2 * mn * (xs + x * math.tan(alpha)) - backlash / (2 * math.cos(alpha))
    thickness_w = math.pi * mn / 2 - 2 * mn * (xs + x * math.tan(alpha)) - backlash / (2 * math.cos(alpha))
    if min(thickness_p, thickness_w) <= 0:
        raise ValueError("侧隙、变位或齿厚修正使中点法向齿厚非正")
    p = _gear_dimensions(inputs, inputs.pinion_teeth, delta_p, theta_fw, theta_aw,
                         hm_ap, hm_fp, thickness_p, re, rm, ri, beta_outer, beta_inner,
                         inputs.root_fillet_radius_pinion_m)
    w = _gear_dimensions(inputs, inputs.wheel_teeth, delta_w, theta_aw, theta_fw,
                         hm_aw, hm_fw, thickness_w, re, rm, ri, beta_outer, beta_inner,
                         inputs.root_fillet_radius_wheel_m)
    # 当量齿轮的法截面近似：基节投影和齿顶至节点的两段展长分别计算。
    virtual_p = p.mean_pitch_diameter_m / (2 * math.cos(delta_p) * math.cos(beta)**2)
    virtual_w = w.mean_pitch_diameter_m / (2 * math.cos(delta_w) * math.cos(beta)**2)
    path = sum(math.sqrt((radius + addendum)**2 - (radius * math.cos(alpha))**2) - radius * math.sin(alpha)
               for radius, addendum in ((virtual_p, hm_ap), (virtual_w, hm_aw)))
    effective_base_pitch = math.pi * mt * math.cos(beta) / (math.cos(alpha) * (math.cos(beta)**2 + math.tan(alpha)**2))
    epsilon_a = path / effective_base_pitch
    width_ratio = inputs.face_width_m / re
    q = width_ratio * (2 - width_ratio) / (2 * (1 - width_ratio)) * math.tan(beta)
    if abs(q) >= 1:
        raise ValueError("齿宽/螺旋角组合超出重合度级数近似的|q|<1适用域")
    epsilon_b = re / (math.pi * me) * (q - q**3 / 3)
    # 根切只给保守近似标记，不能代替刀具包络、过渡曲线和实体干涉检查。
    inner_module = 2 * ri * math.sin(delta_w) * math.cos(beta_inner) / inputs.wheel_teeth
    def virtual_teeth(teeth, delta):
        return teeth / ((1 - math.sin(beta_inner)**2 * math.cos(alpha)**2) * math.cos(beta_inner) * math.cos(delta))
    threshold_p = 1.1 * (ha + (ri - rm) * math.tan(theta_aw) / mn) - virtual_teeth(inputs.pinion_teeth, delta_p) * inner_module * math.sin(alpha)**2 / (2 * mn)
    threshold_w = -(1.1 * (ha + (ri - rm) * math.tan(theta_fw) / mn) - virtual_teeth(inputs.wheel_teeth, delta_w) * inner_module * math.sin(alpha)**2 / (2 * mn))
    warnings = []
    if width_ratio > 1 / 3:
        warnings.append("齿宽超过外锥距的1/3，须另行复核齿高、工艺与接触条件")
    if epsilon_a < 1:
        warnings.append("当量齿廓重合度小于1，不能保证连续啮合")
    if x < threshold_p:
        warnings.append("小轮触发内端近似根切指标，须做刀具包络校核")
    if x > threshold_w:
        warnings.append("大轮触发内端近似根切指标，须做刀具包络校核")
    if abs(q) > 0.5:
        warnings.append("重叠度的三次级数参数超过0.5，请检查误差上界")
    outer_backlash = backlash * re / rm * math.cos(beta_outer) / math.cos(beta)
    return GeometryResult(
        inputs, p, w, u, re, rm, ri, me, mt, mn, mn * (hf - ha), mn * (ha + hf),
        2 * mn * ha, backlash, outer_backlash, backlash / (math.cos(alpha) * math.cos(beta)),
        outer_backlash / (math.cos(alpha) * math.cos(beta_outer)), epsilon_a, epsilon_b,
        math.hypot(epsilon_a, epsilon_b),
        {"tooth_trace": trace, "warnings": tuple(warnings),
         "pinion_undercut_indicator": x < threshold_p, "wheel_undercut_indicator": x > threshold_w,
         "pinion_undercut_threshold": threshold_p, "wheel_undercut_threshold": threshold_w,
         "width_to_outer_cone_ratio": width_ratio, "overlap_series_parameter": q,
         "overlap_series_vs_atan_error_bound": re / (math.pi * me) * abs(q)**5 / 5,
         "overlap_error_bound_scope": "只界定三次级数与arctan(q)的差，不是实际齿面重叠度的误差界",
         "contact_ratio_model": "当量齿轮近似；重叠度采用三次级数；合成为欧氏范数",
         "manufacturing_validation": "未进行刀具包络、真实齿面或制造工艺验证"},
    )


@dataclass(frozen=True)
class ManufacturingSuggestion:
    value_m: float
    out_of_range: bool
    description: str


def suggest_cutter_edge_radius(point_width_m: float, role: str) -> ManufacturingSuggestion:
    """独立的经验刀尖圆角表，不自动进入几何或强度输入，不作工艺认证。"""
    if not math.isfinite(point_width_m) or point_width_m < 0:
        raise ValueError("刀尖点宽必须为非负有限长度，单位m")
    width_inch = point_width_m / 0.0254
    if role == "wheel":
        upper = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.38]
        radii = [0, 0.01, 0.015, 0.02, 0.025, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165]
        outside = 0.17
    elif role == "pinion":
        upper = [0.02, 0.03, 0.04, 0.05, 0.10, 0.13, 0.16, 0.17, 0.19, 0.21, 0.30, 0.38]
        radii = [0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.095, 0.11]
        outside = 0.11
    else:
        raise ValueError("角色只能为pinion或wheel")
    for bound, value in zip(upper, radii):
        if width_inch <= bound:
            return ManufacturingSuggestion(value * 0.0254, width_inch <= upper[0], "经验刀尖圆角参考；必须另行验证刀具与齿根")
    return ManufacturingSuggestion(outside * 0.0254, True, "超出经验刀尖圆角表范围，不作为认证推荐")


def suggest_blade_point_width(point_width_m: float) -> ManufacturingSuggestion:
    """经验刀尖宽度阶梯表，输入和返回均为m。"""
    if not math.isfinite(point_width_m) or point_width_m < 0:
        raise ValueError("点宽必须为非负有限长度，单位m")
    upper = [0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.13, 0.16, 0.19, 0.21, 0.29, 0.38]
    values = [0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.065, 0.08, 0.11, 0.125, 0.15, 0.20]
    width_inch = point_width_m / 0.0254
    for bound, value in zip(upper, values):
        if width_inch <= bound:
            return ManufacturingSuggestion(value * 0.0254, width_inch <= upper[0], "经验刀尖宽度参考，未进行制造验证")
    return ManufacturingSuggestion(0.205 * 0.0254, True, "超出经验刀尖宽度表范围")


def suggested_normal_backlash(outer_transverse_module_m: float, accuracy_grade: int) -> ManufacturingSuggestion:
    """外端法向侧隙的经验表；不自动覆盖调用方给定的中点侧隙。"""
    module_mm = outer_transverse_module_m * 1000
    if not math.isfinite(module_mm) or not 1 <= module_mm <= 25.5:
        raise ValueError("侧隙经验表只覆盖1至25.5mm的外端模数")
    if isinstance(accuracy_grade, bool) or not isinstance(accuracy_grade, int) or accuracy_grade < 1:
        raise ValueError("精度等级必须为正整数；该经验表只按是否大于7分组")
    upper = [1.5, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 20, math.inf]
    coarse = [0.05, 0.08, 0.10, 0.13, 0.15, 0.20, 0.25, 0.33, 0.41, 0.51, 0.69, 0.81]
    fine = [0.03, 0.05, 0.08, 0.10, 0.13, 0.15, 0.20, 0.25, 0.30, 0.38, 0.51, 0.61]
    values = coarse if accuracy_grade > 7 else fine
    for bound, value in zip(upper, values):
        if module_mm < bound:
            return ManufacturingSuggestion(value / 1000, False, "外端法向侧隙经验参考；不代表标准公差验收")
    raise AssertionError("内部侧隙表未覆盖有效参数")
