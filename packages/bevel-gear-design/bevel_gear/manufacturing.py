"""新写的零偏置锥齿轮经验刀具设计链，公开长度m、角度rad。

输出点宽、槽宽、刀尖／顶侧／退刀圆角、经验侧隙和近似量具尺寸，
不表示实体刀具包络、机床设置、完整齿根过渡或制造认证。
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import math
from pathlib import Path

from .geometry import GeometryInput, ManufacturingSuggestion, calculate_geometry
from .io import json_safe, write_json


def _length(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name}必须为有限实数")
    if value <= 0 if positive else value < 0:
        raise ValueError(f"{name}必须为{'正' if positive else '非负'}数")
    return float(value)


@dataclass(frozen=True)
class CutterOverrides:
    """显式刀具值；每项m。小轮槽宽的手动值不改写自动圆角建议的基础宽度。"""
    wheel_point_width_m: float
    pinion_slot_width_m: float
    pinion_tip_radius_m: float
    pinion_top_side_radius_m: float
    pinion_return_radius_m: float
    wheel_tip_radius_m: float

    def __post_init__(self):
        for key, value in self.__dict__.items():
            object.__setattr__(self, key, _length(value, key))


@dataclass(frozen=True)
class ManufacturingInput:
    """用工作齿高、顶隙比和齿高分配定义宏观几何及刀具经验链。

    working_depth_factor=k1、clearance_ratio=k2、wheel_addendum_fraction=c1、
    thickness_factor=k3均无量纲；外端端面模数与所有尺寸用m。
    backlash_quality_grade只用于本经验侧隙表，不等同于其它精度系数入口。
    """
    pinion_teeth: int
    wheel_teeth: int
    outer_transverse_module_m: float
    face_width_m: float
    shaft_angle_rad: float = math.pi / 2
    spiral_angle_rad: float = 0.
    pressure_angle_rad: float = math.radians(20)
    cutter_radius_m: float | None = None
    working_depth_factor: float = 2.
    clearance_ratio: float = .125
    wheel_addendum_fraction: float = .4
    thickness_factor: float = 0.
    backlash_quality_grade: int = 8
    depth_taper: str = "standard"
    tooth_trace: str = "auto"
    cutting_process: str = "single_index"
    cutter_group_count: int = 0
    wheel_addendum_taper_angle_rad: float | None = None
    wheel_dedendum_taper_angle_rad: float | None = None
    manual_cutter: CutterOverrides | None = None

    def __post_init__(self):
        for key in ("outer_transverse_module_m", "face_width_m", "working_depth_factor", "clearance_ratio"):
            object.__setattr__(self, key, _length(getattr(self, key), key, positive=True))
        for key in ("shaft_angle_rad", "spiral_angle_rad", "pressure_angle_rad"):
            object.__setattr__(self, key, _length(getattr(self, key), key))
        if self.cutter_radius_m is not None:
            object.__setattr__(self, "cutter_radius_m", _length(self.cutter_radius_m, "刀盘半径m", positive=True))
        if not isinstance(self.wheel_addendum_fraction, (int, float)) or isinstance(self.wheel_addendum_fraction, bool) or not 0 < self.wheel_addendum_fraction < 1:
            raise ValueError("大轮齿顶高占比必须严格处于0和1之间")
        if not isinstance(self.thickness_factor, (int, float)) or isinstance(self.thickness_factor, bool) or not math.isfinite(self.thickness_factor):
            raise ValueError("齿厚分配系数必须为有限实数")
        if isinstance(self.backlash_quality_grade, bool) or not isinstance(self.backlash_quality_grade, int) or self.backlash_quality_grade < 1:
            raise ValueError("经验侧隙等级必须为正整数")
        if self.cutting_process != "single_index":
            raise NotImplementedError("制造辅助只支持single_index单次分度过程")
        if isinstance(self.cutter_group_count, bool) or not isinstance(self.cutter_group_count, int) or self.cutter_group_count < 0:
            raise ValueError("刀组数必须为非负整数，单次分度应为0")
        if self.tooth_trace == "continuous_index":
            raise NotImplementedError("制造辅助不提供连续分度齿线")
        if self.cutter_group_count != 0:
            raise ValueError("单次分度的cutter_group_count必须为0")
        if self.depth_taper not in ("standard", "uniform", "constant_slot", "limited_slot", "custom"):
            raise NotImplementedError("未知的制造齿高锥度")
        if self.depth_taper == "custom":
            for key in ("wheel_addendum_taper_angle_rad", "wheel_dedendum_taper_angle_rad"):
                if getattr(self, key) is None:
                    raise ValueError("自定义齿高锥度必须显式提供大轮齿顶／齿根锥度角rad")
                object.__setattr__(self, key, _length(getattr(self, key), key))
        elif self.wheel_addendum_taper_angle_rad is not None or self.wheel_dedendum_taper_angle_rad is not None:
            raise ValueError("显式齿高锥度角只能与depth_taper='custom'配合")
        if isinstance(self.manual_cutter, dict):
            object.__setattr__(self, "manual_cutter", CutterOverrides(**self.manual_cutter))
        if self.manual_cutter is not None and not isinstance(self.manual_cutter, CutterOverrides):
            raise TypeError("手动刀具必须为CutterOverrides或对应字典")
        manufacturing_input_to_geometry(self)


def manufacturing_input_to_geometry(inputs: ManufacturingInput) -> GeometryInput:
    """四个制造系数显式映射到现有几何接口，侧隙先取0保留名义齿厚。"""
    if not isinstance(inputs, ManufacturingInput):
        raise TypeError("几何映射需要ManufacturingInput")
    k1, k2, fraction = inputs.working_depth_factor, inputs.clearance_ratio, inputs.wheel_addendum_fraction
    extensions = {}
    trace = inputs.tooth_trace
    if inputs.depth_taper == "custom":
        extensions.update(wheel_addendum_taper_angle_rad=inputs.wheel_addendum_taper_angle_rad,
                          wheel_dedendum_taper_angle_rad=inputs.wheel_dedendum_taper_angle_rad)
    return GeometryInput(
        inputs.pinion_teeth, inputs.wheel_teeth, inputs.wheel_teeth * inputs.outer_transverse_module_m,
        inputs.face_width_m, shaft_angle_rad=inputs.shaft_angle_rad,
        spiral_angle_rad=inputs.spiral_angle_rad, pressure_angle_rad=inputs.pressure_angle_rad,
        addendum_coefficient=k1 / 2, dedendum_coefficient=k1 * (k2 + .5),
        pinion_profile_shift_coefficient=k1 * (.5 - fraction),
        tooth_thickness_modification_coefficient=inputs.thickness_factor / 2,
        cutter_radius_m=inputs.cutter_radius_m, depth_taper=inputs.depth_taper, tooth_trace=trace,
        **extensions,
    )


def empirical_backlash(outer_transverse_module_m: float, quality_grade: int) -> ManufacturingSuggestion:
    """本制造链的径节经验侧隙，1≤P_D≤25；超域保留经验分支并标记。

    P_D=0.0254/m；质量等级≤9与>9分两表。不能与其它模数侧隙表混用。
    """
    module = _length(outer_transverse_module_m, "外端端面模数m", positive=True)
    if isinstance(quality_grade, bool) or not isinstance(quality_grade, int) or quality_grade < 1:
        raise ValueError("经验侧隙等级必须为正整数")
    pitch = 25.4 / (module * 1000)
    coarse = quality_grade <= 9
    if pitch < 1:
        inch, outside = (.002 if coarse else .001), True
    elif pitch > 25:
        inch, outside = (.032 if coarse else .024), True
    else:
        boundaries = (16., 10., 8., 6., 5., 4., 3., 2.5, 2., 1.5, 1.25, 0.)
        table = (.002, .003, .004, .005, .006, .008, .010, .013, .016, .020, .027, .032) if coarse else (
                 .001, .002, .003, .004, .005, .006, .008, .010, .012, .015, .020, .024)
        inch = next(value for boundary, value in zip(boundaries, table) if pitch >= boundary)
        outside = False
    return ManufacturingSuggestion(inch * .0254, outside, "径节侧隙经验参考；超域分支不代表有效制造推荐")


def quantize_wheel_point_width(raw_width_m: float) -> float:
    """保留0.254mm步进与截向零后加1；即使恰在整数步也向下一步。"""
    if not math.isfinite(raw_width_m):
        raise ValueError("点宽量化输入必须为有限长度")
    return (math.trunc(raw_width_m * 1000 / .254) + 1) * .254 / 1000


def quantize_pinion_slot_width(raw_width_m: float) -> float:
    """以0.254mm为步长作银行家舍入，不使用四舍五入的远离零规则。"""
    if not math.isfinite(raw_width_m):
        raise ValueError("槽宽量化输入必须为有限长度")
    return round(raw_width_m * 1000 / .254) * .254 / 1000


def _edge_radius(width, role):
    # 表的分支比较保留mm→inch的运算次序；m/0.0254在边界可能差1ULP。
    inch = (width * 1000) / 25.4
    if role == "wheel":
        bounds = (.02, .03, .04, .05, .06, .07, .08, .09, .10, .11, .12, .13, .14,
                  .15, .16, .17, .18, .19, .20, .21, .22, .23, .24, .25, .26, .27,
                  .28, .29, .30, .31, .32, .33, .38)
        radii = (0., .01, .015, .02, .025, .025, .03, .035, .04, .045, .05, .055, .06,
                 .07, .075, .08, .085, .09, .095, .10, .105, .11, .115, .12, .125,
                 .13, .135, .14, .145, .15, .155, .16, .165)
        outside = .17
    else:
        bounds = (.02, .03, .04, .05, .10, .13, .16, .17, .19, .21, .30, .38)
        radii = (0., .01, .015, .02, .025, .03, .04, .05, .06, .075, .095, .11)
        outside = .11
    for bound, radius in zip(bounds, radii):
        if inch <= bound:
            return ManufacturingSuggestion((radius * 25.4) / 1000, inch <= bounds[0], "圆角经验表，未进行刀具包络验证")
    return ManufacturingSuggestion((outside * 25.4) / 1000, True, "圆角经验表超域分支")


def _blade_width(width):
    inch = (width * 1000) / 25.4
    bounds = (.015, .02, .03, .04, .05, .06, .08, .10, .13, .16, .19, .21, .29, .38)
    values = (0., .01, .015, .02, .025, .03, .04, .05, .065, .08, .11, .125, .15, .20)
    for bound, value in zip(bounds, values):
        if inch <= bound:
            return ManufacturingSuggestion((value * 25.4) / 1000, inch <= bounds[0], "刀尖宽度经验表，未进行制造验证")
    return ManufacturingSuggestion((.205 * 25.4) / 1000, True, "刀尖宽度经验表超域分支")


def calculate_manufacturing(inputs: ManufacturingInput | dict, *, strict: bool = True) -> dict:
    """求实算刀具链；strict拒绝不可用尺寸，False保留有限原值和invalid状态。

    strict不夹改圆角或宽度，不能把负值改0后当推荐值。反三角函数及
    正常外啮合锥体域仍由几何层检查；strict=False不是NaN放行模式。
    """
    if isinstance(inputs, dict):
        inputs = ManufacturingInput(**inputs)
    if not isinstance(inputs, ManufacturingInput):
        raise TypeError("制造入口需要ManufacturingInput或具名字典")
    if not isinstance(strict, bool):
        raise TypeError("strict必须为布尔值")
    mapped = manufacturing_input_to_geometry(inputs)
    geometry = calculate_geometry(mapped)
    p, w = geometry.pinion, geometry.wheel
    outer, mean, inner = geometry.outer_cone_distance_m, geometry.mean_cone_distance_m, geometry.inner_cone_distance_m
    beta, alpha, width = inputs.spiral_angle_rad, inputs.pressure_angle_rad, inputs.face_width_m
    beta_outer, beta_inner = p.outer_spiral_angle_rad, p.inner_spiral_angle_rad
    pitch_mean = math.pi * geometry.mean_transverse_module_m
    tan_alpha = math.tan(alpha)
    normal_pitch = pitch_mean * math.cos(beta)
    backlash = empirical_backlash(inputs.outer_transverse_module_m, inputs.backlash_quality_grade)
    jmin, jmax = backlash.value_m, 1.5 * backlash.value_m
    automatic_wheel_raw = normal_pitch - w.mean_normal_tooth_thickness_m - 2 * w.mean_dedendum_m * tan_alpha
    automatic_wheel = quantize_wheel_point_width(automatic_wheel_raw)
    manual = inputs.manual_cutter
    wheel_width = automatic_wheel if manual is None else manual.wheel_point_width_m
    gaps = {"outer": math.pi * inputs.outer_transverse_module_m * math.cos(beta_outer)
                     - 2 * (w.outer_dedendum_m + p.outer_dedendum_m) * tan_alpha - wheel_width + jmin,
            "mean": normal_pitch - 2 * (w.mean_dedendum_m + p.mean_dedendum_m) * tan_alpha - wheel_width + jmin,
            "inner": inner / outer * math.pi * inputs.outer_transverse_module_m * math.cos(beta_inner)
                     - 2 * (w.inner_dedendum_m + p.inner_dedendum_m) * tan_alpha - wheel_width + jmin}
    automatic_slot_raw = min(gaps.values()) - .0007
    automatic_slot = quantize_pinion_slot_width(automatic_slot_raw)
    slot_width = automatic_slot if manual is None else manual.pinion_slot_width_m
    wheel_radius_suggestion = _edge_radius(wheel_width, "wheel")
    return_radius_suggestion = _edge_radius(automatic_slot, "pinion")
    blade_suggestion = _blade_width(automatic_slot + .0007)
    blade_width = blade_suggestion.value_m
    denominator = 1 - math.sin(alpha)
    clearance_bound = (geometry.normal_clearance_m - jmin * inner / (2 * outer * tan_alpha)) / denominator
    blade_bound = (blade_width - .00038) * math.cos(alpha) / denominator
    pressure_term = 1 / math.cos(alpha) - tan_alpha
    excess_inch = (automatic_slot + .0007 - blade_width) / .0254
    radicand = .001 + pressure_term * excess_inch
    empirical_bound = (None if radicand < 0 else .0254 * (
        .002 + pressure_term * excess_inch + .063 * math.sqrt(radicand)) / pressure_term**2)
    automatic_tip = None if empirical_bound is None else min(clearance_bound, blade_bound, empirical_bound)
    tip = automatic_tip if manual is None else manual.pinion_tip_radius_m
    top_side = automatic_tip if manual is None else manual.pinion_top_side_radius_m
    return_radius = return_radius_suggestion.value_m if manual is None else manual.pinion_return_radius_m
    wheel_radius = wheel_radius_suggestion.value_m if manual is None else manual.wheel_tip_radius_m
    backlash_mean = jmin * mean / outer * math.cos(beta_outer) / math.cos(beta)
    chordal_gap = backlash_mean / (2 * math.cos(alpha))
    chord_p = p.mean_normal_tooth_thickness_m - p.mean_normal_tooth_thickness_m**3 / (6 * p.mean_pitch_diameter_m**2) - chordal_gap
    chord_w = w.mean_normal_tooth_thickness_m - w.mean_normal_tooth_thickness_m**3 / (6 * w.mean_pitch_diameter_m**2) - chordal_gap
    transverse_inner_angle = math.atan(tan_alpha / math.cos(beta_inner))
    root_cutting_limit = inner * math.tan(p.pitch_cone_angle_rad) * math.sin(transverse_inner_angle)**2
    root_cutting_indicator = root_cutting_limit <= p.inner_dedendum_m
    offset_wheel_tip = -(mean * math.sin(p.dedendum_taper_angle_rad) - w.mean_addendum_m * math.cos(p.dedendum_taper_angle_rad)) / math.sin(w.tip_cone_angle_rad)
    offset_wheel_root = (mean * math.sin(w.dedendum_taper_angle_rad) - w.mean_dedendum_m * math.cos(w.dedendum_taper_angle_rad)) / math.sin(w.root_cone_angle_rad)
    offset_pinion_tip = (-offset_wheel_root * math.sin(w.root_cone_angle_rad) - geometry.normal_clearance_m) / math.sin(p.tip_cone_angle_rad)
    offset_pinion_root = (-offset_wheel_tip * math.sin(w.tip_cone_angle_rad) - geometry.normal_clearance_m) / math.sin(p.root_cone_angle_rad)
    actual = {"wheel_point_width_m": wheel_width, "pinion_slot_width_m": slot_width,
              "pinion_tip_radius_m": tip, "pinion_top_side_radius_m": top_side,
              "pinion_return_radius_m": return_radius, "wheel_tip_radius_m": wheel_radius}
    unusable = [name for name, value in actual.items() if value is None or not math.isfinite(value)
                or value < 0 or ("width" in name and value == 0)]
    if empirical_bound is None:
        unusable.append("pinion_empirical_radius_domain")
    if min(chord_p, chord_w) <= 0:
        unusable.append("chordal_tooth_thickness_m")
    warnings = list(geometry.diagnostics["warnings"])
    if manual is not None:
        warnings.append("手动小轮槽宽保留为有效值；圆角及刀尖宽度建议仍按自动槽宽计算，二者分别输出")
    if backlash.out_of_range:
        warnings.append("模数超出径节侧隙经验表域；保留经验分支值，不代表有效推荐")
    for title, suggestion in (("大轮圆角", wheel_radius_suggestion), ("小轮退刀圆角", return_radius_suggestion), ("小轮刀尖宽度", blade_suggestion)):
        if suggestion.out_of_range:
            warnings.append(title + "超出经验表范围")
    if empirical_bound is None:
        warnings.append("圆角经验根式超出实数域，自动建议未知；没有将NaN替换成0")
    if unusable:
        warnings.append("不可用经验尺寸：" + ", ".join(unusable))
        if strict:
            raise ValueError(warnings[-1] + "；可用strict=False保留原计算值审查")
    suggested_geometry = asdict(mapped)
    suggested_geometry["backlash_m"] = backlash_mean
    return {"model": "zero_offset_bevel_empirical_cutter_chain", "status": "invalid" if unusable else "ok",
            "inputs": asdict(inputs), "geometry": asdict(geometry), "mapped_geometry_input": asdict(mapped),
            "suggested_geometry_input_with_backlash": suggested_geometry,
            "backlash": {"minimum_outer_m": jmin, "maximum_outer_m": jmax,
                          "minimum_mean_normal_m": backlash_mean, "out_of_range": backlash.out_of_range,
                          "diametral_pitch_per_inch": .0254 / inputs.outer_transverse_module_m},
            "cutter": {**actual, "manual": manual is not None, "automatic_wheel_raw_width_m": automatic_wheel_raw,
                       "automatic_wheel_point_width_m": automatic_wheel,
                       "automatic_pinion_raw_slot_width_m": automatic_slot_raw,
                       "automatic_pinion_slot_width_m": automatic_slot, "pinion_blade_point_width_m": blade_width,
                       "automatic_pinion_tip_radius_m": automatic_tip,
                       "suggested_pinion_return_radius_m": return_radius_suggestion.value_m,
                       "suggested_wheel_tip_radius_m": wheel_radius_suggestion.value_m,
                       "pinion_radius_bounds_m": {"clearance": clearance_bound, "blade_width": blade_bound,
                                                   "empirical": empirical_bound},
                       "section_slot_limits_m": gaps,
                       "out_of_range": {"wheel_radius": wheel_radius_suggestion.out_of_range,
                                        "pinion_return_radius": return_radius_suggestion.out_of_range,
                                        "blade_point_width": blade_suggestion.out_of_range}},
            "inspection": {"pinion_approximate_chordal_thickness_m": chord_p, "wheel_approximate_chordal_thickness_m": chord_w,
                           "pinion_approximate_chordal_addendum_m": p.mean_addendum_m + .25 * p.mean_normal_tooth_thickness_m**2 * math.cos(p.pitch_cone_angle_rad) / p.mean_pitch_diameter_m,
                           "wheel_approximate_chordal_addendum_m": w.mean_addendum_m + .25 * w.mean_normal_tooth_thickness_m**2 * math.cos(w.pitch_cone_angle_rad) / w.mean_pitch_diameter_m},
            "root_cutting": {"pinion_inner_limit_m": root_cutting_limit,
                             "pinion_inner_dedendum_m": p.inner_dedendum_m,
                             "inner_transverse_pressure_angle_rad": transverse_inner_angle,
                             "pinion_undercut_indicator": root_cutting_indicator},
            "cone_reference_offsets_m": {"wheel_tip": offset_wheel_tip, "wheel_root": offset_wheel_root,
                                         "pinion_tip": offset_pinion_tip, "pinion_root": offset_pinion_root},
            "diagnostics": {"warnings": warnings, "unusable_fields": unusable,
                            "geometry_backlash": "geometry保留零侧隙名义齿厚；建议侧隙在独立输入副本中给出，不自动进入强度",
                            "scope": "经验刀具参数与近似量具尺寸，未作实体刀具包络或制造工艺认证"}}


def export_manufacturing_result(result: dict, output_dir) -> list[str]:
    """保存完整JSON及具名制造参数CSV；未知根式保留空值，不补0。"""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    files = [write_json(result, output / "manufacturing_result.json")]
    rows = []
    def flatten(value, path, unit="1"):
        if isinstance(value, dict):
            for key, item in value.items():
                child_unit = "m" if key.endswith("_m") else "rad" if key.endswith("_rad") else unit
                if key == "diametral_pitch_per_inch":
                    child_unit = "1/inch"
                flatten(item, path + "." + key, child_unit)
        else:
            rows.append((path, json_safe(value), "状态" if isinstance(value, bool) else unit))
    for key in ("backlash", "cutter", "inspection", "root_cutting", "cone_reference_offsets_m"):
        flatten(result[key], key, "m" if key.endswith("_m") else "1")
    target = output / "manufacturing_parameters.csv"
    with target.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("quantity", "value", "unit"))
        writer.writerows(rows)
    files.append(target)
    return [str(path.resolve()) for path in files]
