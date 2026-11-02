"""双面刀盘展成五面网格，公开长度m、角度rad，核心仅使用标准库。

两侧有效面和刀尖根过渡都来自逐点包络。槽底采用根直径上连接两侧
根端点的圆弧；这项几何封底约定不表述为槽底的完整刀具包络。
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

from .generated_cutter import CutterMachineSettings, evaluate_generated_point


PATCH_NAMES = ("left_flank", "left_fillet", "root_bottom", "right_fillet", "right_flank")


@dataclass(frozen=True)
class GeneratedGearInput:
    """显式单轮机床和边界，不自动设计配对机床设置或求装配接触。"""
    machine: CutterMachineSettings | dict
    tooth_count: int
    module_m: float
    face_width_m: float
    addendum_m: float
    tip_cone_angle_rad: float
    tip_relief_m: float = 0.
    tip_relief_start: float = .65

    def __post_init__(self):
        if isinstance(self.machine, dict):
            object.__setattr__(self, "machine", CutterMachineSettings(**self.machine))
        if not isinstance(self.machine, CutterMachineSettings):
            raise TypeError("machine须为CutterMachineSettings或具名字典")
        if isinstance(self.tooth_count, bool) or not isinstance(self.tooth_count, int) or self.tooth_count < 3:
            raise ValueError("tooth_count须为不小于3的整数")
        for name in ("module_m", "face_width_m", "addendum_m", "tip_cone_angle_rad"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name}须为正有限实数")
        if self.tip_cone_angle_rad >= math.pi/2:
            raise ValueError("齿顶锥角须为正锐角")
        expected = self.tooth_count*self.module_m/(2*math.sin(self.machine.pitch_cone_angle_rad))
        if not math.isclose(self.machine.reference_cone_distance_m, expected, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("机床参考锥距须与齿数、模数及节锥角一致")
        if self.face_width_m+0.00002 >= self.machine.reference_cone_distance_m:
            raise ValueError("含内端裁切余量的齿宽须小于参考锥距")
        if self.machine.dedendum_m <= 0:
            raise ValueError("生成根过渡须给正齿根高")
        if (isinstance(self.tip_relief_m,bool) or not isinstance(self.tip_relief_m,(int,float)) or
                not math.isfinite(self.tip_relief_m) or self.tip_relief_m < 0):
            raise ValueError('法向齿顶修缘量须为非负有限m')
        if (isinstance(self.tip_relief_start,bool) or not isinstance(self.tip_relief_start,(int,float)) or
                not math.isfinite(self.tip_relief_start) or not 0 <= self.tip_relief_start < 1):
            raise ValueError('法向齿顶修缘起点须为[0,1)内齿廓参数')


def _request(value):
    if isinstance(value, dict):
        value = GeneratedGearInput(**value)
    if not isinstance(value, GeneratedGearInput):
        raise TypeError("展成网格需要GeneratedGearInput或具名字典")
    return value


def _sample_count(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 2:
        raise ValueError(f"{name}须为不小于2的整数")


def _normalize_angle(angle):
    while angle > math.pi:
        angle -= 2*math.pi
    while angle < -math.pi:
        angle += 2*math.pi
    return angle


class _Sections:
    """同一宽向截面复用两侧包络评价，左右高度仍由内核分别求解。"""
    def __init__(self, request):
        self.request = request
        self.cache = {}
        self.reference_mm = request.tooth_count*(request.module_m*1000)*0.5/math.sin(request.machine.pitch_cone_angle_rad)
        self.pitch_diameter_mm = request.tooth_count*(request.module_m*1000)
        self.tip_mm = self.pitch_diameter_mm+2*(request.addendum_m*1000)*math.cos(request.machine.pitch_cone_angle_rad)
        self.root_mm = self.pitch_diameter_mm-2*(request.machine.dedendum_m*1000)*math.cos(request.machine.pitch_cone_angle_rad)
        # 毫米整档的顶径余量和宽端余量属于明确建模参数，不能假定尺度不变。
        self.tip_extension_mm = math.ceil(.01*(request.module_m*1000)*1000)/1000

    def station(self, parameter):
        return (self.reference_mm+.02-(self.request.face_width_m*1000+.04)*parameter)/1000

    def diameter(self, large_diameter_mm, cone_angle, station):
        pitch = self.request.machine.pitch_cone_angle_rad
        distance = (self.reference_mm-station*1000)/math.cos(pitch-cone_angle)
        radius = .5*large_diameter_mm/math.sin(cone_angle)
        return 2*(radius-distance)*math.sin(cone_angle)/1000

    def section(self, station):
        if station in self.cache:
            return self.cache[station]
        request = self.request
        tip = self.diameter(self.tip_mm+2*self.tip_extension_mm, request.tip_cone_angle_rad, station)
        root = self.diameter(self.root_mm, request.machine.machine_root_angle_rad, station)
        if not 0 < root < tip:
            raise ValueError("该宽向截面的根/顶直径须为正并严格递增")
        evaluated = {}
        def evaluate(diameter):
            if diameter not in evaluated:
                evaluated[diameter] = evaluate_generated_point(request.machine, diameter, station)
            return evaluated[diameter]
        transitions = {}
        for flank in ("left", "right"):
            if evaluate(tip)["flanks"][flank]["region"] != "straight_blade":
                raise ValueError("延伸齿顶未落在直刃包络区，请调整机床设置或边界")
            coarse = None
            for index in range(1, 65):
                diameter = tip-(tip-root)*index/64
                if evaluate(diameter)["flanks"][flank]["region"] == "cutter_tip_fillet":
                    coarse = index
                    break
            if coarse is None:
                raise ValueError("实际根径前未找到直刃至刀尖圆角的切换，不能填造根过渡")
            start, end = (coarse-1)*8, coarse*8
            upper = tip-(tip-root)*start/512
            first_fillet = None
            for index in range(start+1, end+1):
                lower = tip-(tip-root)*index/512
                if (evaluate(upper)["flanks"][flank]["region"] == "straight_blade"
                        and evaluate(lower)["flanks"][flank]["region"] == "cutter_tip_fillet"):
                    first_fillet = lower
                    break
                upper = lower
            if first_fillet is None:
                raise ValueError("细分区间未包含有效直刃/圆角切换")
            lower = first_fillet
            for _ in range(48):
                middle = (upper+lower)/2
                if evaluate(middle)["flanks"][flank]["region"] == "straight_blade":
                    upper = middle
                else:
                    lower = middle
            transitions[flank] = {"last_straight_diameter_m": upper, "first_fillet_diameter_m": lower,
                "last_straight": evaluate(upper)["flanks"][flank],
                "first_fillet": evaluate(lower)["flanks"][flank]}
        result = {"station_m": station, "root_diameter_m": root, "extended_tip_diameter_m": tip,
                  "root_cutter_height_upper_m": evaluate(root)["maximum_cutter_height_m"],
                  "transitions": transitions, "evaluate": evaluate}
        self.cache[station] = result
        return result


def _patch_point(section, name, parameter):
    if name == "root_bottom":
        left = section["evaluate"](section["root_diameter_m"])["flanks"]["left"]["point_m"]
        right = section["evaluate"](section["root_diameter_m"])["flanks"]["right"]["point_m"]
        if parameter == 0:
            return left, "root_bottom_connector"
        if parameter == 1:
            return right, "root_bottom_connector"
        angle = math.atan2(left[1], left[0])
        angle += parameter*_normalize_angle(math.atan2(right[1], right[0])-angle)
        radius = section["root_diameter_m"]/2
        return (radius*math.cos(angle), radius*math.sin(angle), left[2]), "root_bottom_connector"
    flank = "left" if name.startswith("left") else "right"
    transition = section["transitions"][flank]
    boundary = transition["last_straight_diameter_m"]
    if name.endswith("_flank"):
        diameter = boundary+(section["extended_tip_diameter_m"]-boundary)*parameter
        value = transition["last_straight"] if parameter == 0 else section["evaluate"](diameter)["flanks"][flank]
    else:
        diameter = section["root_diameter_m"]+(boundary-section["root_diameter_m"])*parameter
        value = transition["last_straight"] if parameter == 1 else section["evaluate"](diameter)["flanks"][flank]
    return value["point_m"], value["region"]


def build_generated_surface_grid(request: GeneratedGearInput | dict, *, width_count=25, profile_count=25) -> dict:
    """构造单齿槽五面开放网格；点序为profile_index*width_count+width_index。

    宽向从外至内均匀取样，侧面从直刃根界至延伸齿顶，圆角从实际根径
    至直刃切换界。共享边复用相同包络评价，槽底明确使用圆弧连接。
    """
    request = _request(request)
    _sample_count(width_count, "width_count")
    _sample_count(profile_count, "profile_count")
    if 5*width_count*profile_count > 2_000_000:
        raise ValueError("五面网格超过两百万点，请降低采样规模")
    sections = _Sections(request)
    patches = {name: {"width_count": width_count, "profile_count": profile_count,
                      "points_m": [], "uv": [], "region_labels": []} for name in PATCH_NAMES}
    boundaries = []
    stations = []
    for width_index in range(width_count):
        u = width_index/(width_count-1)
        section = sections.section(sections.station(u))
        stations.append(section["station_m"])
        boundaries.append({key: value for key, value in section.items() if key != "evaluate"})
    for profile_index in range(profile_count):
        v = profile_index/(profile_count-1)
        for width_index, station in enumerate(stations):
            section = sections.section(station)
            for name, patch in patches.items():
                point, region = _patch_point(section, name, v)
                patch["points_m"].append(point)
                patch["uv"].append((width_index/(width_count-1), v))
                patch["region_labels"].append(region)
    grid = {"inputs": asdict(request), "units": "m", "angle_units": "rad", "patches": patches,
            "stations_m": stations, "sections": boundaries,
            "metadata": {"coordinate_frame": "本轮局部+Z向外；+X径向、+Y切向；不含装配变换",
                "station_definition": "节锥母线截面位置，不是任意齿面点的球面半径",
                "point_order": "profile_index*width_count+width_index", "closed_solid": False,
                "root_transition_model": "两侧均由刀尖圆弧参与的实际展成包络求得",
                "root_bottom_model": "根直径圆上连接两侧根端点的圆弧，不宣称槽底完整刀具包络",
                "width_end_extension_m": .00002, "tip_diameter_extension_per_side_m": sections.tip_extension_mm/1000,
                "negative_root_height_interval_sections": sum(section["root_cutter_height_upper_m"] < 0 for section in boundaries),
                "minimum_root_cutter_height_upper_m": min(section["root_cutter_height_upper_m"] for section in boundaries),
                "root_height_interval_policy": "保留内核默认兼容返回与负区间诊断，不夹零；负区间不作为全物理有效域证明",
                "width_boundary_model": "含外/内端裁切余量的机床截面，尚未按齿坯裁切",
                "limitations": ["单齿槽开放五面，不包含周期加工、轮体实体或齿坯裁切",
                                "显式机床输入不保证配对共轭或制造可行性", "不计算TCA、LTCA或接触压力"]}}
    if request.tip_relief_m:
        from .flank_relief import apply_normal_tip_relief
        grid = apply_normal_tip_relief(grid,request)
    return grid


def build_generated_accuracy_grid(request: GeneratedGearInput | dict) -> dict:
    """独立于均匀插值节点的27个几何精度点，包含复合根过渡检查。"""
    request = _request(request)
    sections = _Sections(request)
    parameters = (.5-.5*math.sqrt(.6), .5, .5+.5*math.sqrt(.6))
    groups = {name: {"width_count": 3, "profile_count": 3, "points_m": []}
              for name in ("left_flank", "root_transition", "right_flank")}
    for v in parameters:
        for u in parameters:
            section = sections.section(sections.station(u))
            for name in ("left_flank", "right_flank"):
                if request.tip_relief_m:
                    from .flank_relief import independent_relieved_point
                    point=independent_relieved_point(request,sections,name,u,v)
                else:
                    point=_patch_point(section,name,v)[0]
                groups[name]["points_m"].append(point)
            if v < 1/3:
                point = _patch_point(section, "left_fillet", 1-3*v)[0]
            elif v > 2/3:
                point = _patch_point(section, "right_fillet", 3*v-2)[0]
            else:
                point = _patch_point(section, "root_bottom", 3*v-1)[0]
            groups["root_transition"]["points_m"].append(point)
    return {"units": "m", "groups": groups, "count": 27,
            "metadata": {"parameters": parameters, "purpose": "均匀网格插值节点之外的点面误差检查"}}
