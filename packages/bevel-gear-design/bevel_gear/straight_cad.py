"""直齿球面渐开线的闭实体工程模型；公开长度m，OCP内部mm。

工作面来自已有解析球面渐开线；工程根延伸为径向直线连接，不是切削
圆角。轴向平面端界与开放齿面模块的球面端界不同。插值/放样为CAD
近似，必须按采样和容差验证；闭实体不代表制造、共轭或承载验证。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .geometry import GeometryResult
from .tooth_surface import SphericalToothSurface, build_spherical_tooth_surface


@dataclass(frozen=True)
class StraightCadSettings:
    """孔径及端界为m；默认端界穿过内、外端节锥站点。

    engineering_root_extension_m强制保留正长度工程根连接；base_margin_rad
    避开基锥曲率奇点。二者属于明确建模假设，不能当作加工圆角半径。
    """
    bore_diameter_m: float
    inner_axial_m: float | None = None
    outer_axial_m: float | None = None
    profile_count: int = 33
    width_sections: int = 2
    engineering_root_extension_m: float = 1e-4
    base_margin_rad: float = 1e-5
    tolerance_m: float = 1e-7

    def __post_init__(self):
        for name in ("bore_diameter_m", "engineering_root_extension_m", "base_margin_rad", "tolerance_m"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name}须为正有限实数")
        for name in ("inner_axial_m", "outer_axial_m"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name}须为正有限轴坐标m或None")
        for name, minimum in (("profile_count", 5), ("width_sections", 2)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= 1001:
                raise ValueError(f"{name}须为{minimum}至1001的整数")
        if self.base_margin_rad >= .05:
            raise ValueError("基锥避让角须小于0.05rad，不能代替大范围齿廓删减")


@dataclass
class StraightCadModel:
    """shape内部mm；face_groups为当前BRep唯一面0起始序号。

    STEP往返后的面编号可能重排；FE约束应使用boundaries的SI几何选择器。
    """
    shape: object
    diagnostics: dict
    face_groups: dict
    boundaries: dict
    surface: SphericalToothSurface
    settings: StraightCadSettings
    tooth_index: int | None


@dataclass
class StraightCadPair:
    """两轮局部闭实体及安装后的复合体；各轮仍为独立实体。"""
    shape: object
    pinion: StraightCadModel
    wheel: StraightCadModel
    placements: dict
    diagnostics: dict


def _surface(value, member):
    if isinstance(value, GeometryResult):
        return build_spherical_tooth_surface(value, member)
    if isinstance(value, SphericalToothSurface):
        return value
    raise TypeError("需要GeometryResult或SphericalToothSurface")


def _ends(surface, settings):
    cosine = math.cos(surface.profile.pitch_cone_angle_rad)
    low = settings.inner_axial_m if settings.inner_axial_m is not None else surface.inner_cone_distance_m * cosine
    high = settings.outer_axial_m if settings.outer_axial_m is not None else surface.outer_cone_distance_m * cosine
    if low >= high or high - low <= 10 * settings.tolerance_m:
        raise ValueError("外端轴坐标须大于内端，且实体厚度须超过十倍建模容差")
    return low, high


def _section(surface, settings, z):
    # 将宏观锥面h=A+B*s与轴向平面z联立；s为节锥母线坐标。
    delta = surface.profile.pitch_cone_angle_rad
    sine, cosine = math.sin(delta), math.cos(delta)
    radii = []
    for height, taper, sign in ((surface.mean_dedendum_m, surface.dedendum_taper_angle_rad, -1),
                                (surface.mean_addendum_m, surface.addendum_taper_angle_rad, 1)):
        b = math.tan(taper)
        a = height - surface.reference_cone_distance_m * b
        denominator = cosine - sign * b * sine
        if denominator <= 0:
            raise ValueError("宏观锥面在该轴向截面不形成有限正半径")
        radius = (z * (sine + sign * b * cosine) + sign * a) / denominator
        radii.append(radius)
    root, tip = radii
    lower = max(math.atan2(root + settings.engineering_root_extension_m, z),
                surface.profile.base_cone_angle_rad + settings.base_margin_rad)
    upper = math.atan2(tip, z)
    if not 0 < root < tip or not lower < upper < math.pi / 2:
        raise ValueError("工程根连接占满有效齿高，或宏观截面不是正锐角锥体")
    if root <= settings.bore_diameter_m / 2 + 10 * settings.tolerance_m:
        raise ValueError("孔径与根圆间必须保留超过十倍容差的正轮体厚度")
    lower_half = surface.profile.half_material_angle(lower)
    upper_half = surface.profile.half_material_angle(upper)
    if not 0 < upper_half < lower_half < math.pi / surface.teeth:
        raise ValueError("有效齿廓尖齿或根侧材料占满齿距，无法生成分离齿槽")
    return {"root_radius_m": root, "tip_radius_m": tip, "lower_polar_rad": lower,
            "upper_polar_rad": upper, "lower_half_angle_rad": lower_half,
            "upper_half_angle_rad": upper_half}


def _ocp():
    try:
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
        from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
        from OCP.GeomAPI import GeomAPI_Interpolate
        from OCP.TColgp import TColgp_HArray1OfPnt
        from OCP.GC import GC_MakeArcOfCircle
        from OCP.gp import gp_Pnt, gp_Ax2, gp_Dir
        from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
        from OCP.TopTools import TopTools_IndexedMapOfShape
        from OCP.TopExp import TopExp
        from OCP.TopAbs import TopAbs_FACE
        from OCP.BRepTools import BRepTools_WireExplorer
    except ImportError as error:
        raise ImportError("闭实体直齿CAD需要可选OCP依赖") from error
    return locals()


def _modified_face_indices(cut,face,mapping):
    """复制布尔历史列表并逐首项读取，避免借用原生列表的迭代器生命周期。"""
    from OCP.TopTools import TopTools_ListOfShape
    history=TopTools_ListOfShape(cut.Modified(face))
    indices=[]
    while not history.IsEmpty():
        index=mapping.FindIndex(history.First())
        if index:indices.append(index-1)
        history.RemoveFirst()
    return indices


def _wire(surface, settings, z, first_tooth, count, ocp, sector):
    section = _section(surface, settings, z)
    pitch = 2 * math.pi / surface.teeth
    root, tip = section["root_radius_m"], section["tip_radius_m"]
    lower, upper = section["lower_polar_rad"], section["upper_polar_rad"]
    half_low, half_tip = section["lower_half_angle_rad"], section["upper_half_angle_rad"]
    low_radius = z * math.tan(lower)
    wire = ocp["BRepBuilderAPI_MakeWire"]()
    groups = []

    def point(radius, angle):
        return ocp["gp_Pnt"](1000 * radius * math.cos(angle), 1000 * radius * math.sin(angle), 1000 * z)

    def add(edge, label):
        wire.Add(edge)
        groups.append((label, edge))

    def line(r1, a1, r2, a2, label):
        add(ocp["BRepBuilderAPI_MakeEdge"](point(r1, a1), point(r2, a2)).Edge(), label)

    def arc(radius, start, stop, label):
        curve = ocp["GC_MakeArcOfCircle"](point(radius, start), point(radius, (start + stop) / 2), point(radius, stop)).Value()
        add(ocp["BRepBuilderAPI_MakeEdge"](curve).Edge(), label)

    def flank(phase, side, reverse=False):
        points = ocp["TColgp_HArray1OfPnt"](1, settings.profile_count)
        for index in range(settings.profile_count):
            fraction = index / (settings.profile_count - 1)
            if reverse:
                fraction = 1 - fraction
            theta = lower + fraction * (upper - lower)
            angle = phase + surface.profile.flank_angle(theta, side)
            points.SetValue(index + 1, point(z * math.tan(theta), angle))
        interpolation = ocp["GeomAPI_Interpolate"](points, False, settings.tolerance_m * 1000)
        interpolation.Perform()
        if not interpolation.IsDone():
            raise ArithmeticError("渐开线截面插值未完成")
        add(ocp["BRepBuilderAPI_MakeEdge"](interpolation.Curve()).Edge(), "involute_" + side)

    for tooth in range(first_tooth, first_tooth + count):
        phase = tooth * pitch
        arc(root, phase - pitch / 2, phase - half_low, "root_land_left")
        line(root, phase - half_low, low_radius, phase - half_low, "engineering_root_left")
        flank(phase, "left")
        arc(tip, phase - half_tip, phase + half_tip, "tip")
        flank(phase, "right", True)
        line(low_radius, phase + half_low, root, phase + half_low, "engineering_root_right")
        arc(root, phase + half_low, phase + pitch / 2, "root_land_right")
    if sector:
        phase = first_tooth * pitch
        bore = settings.bore_diameter_m / 2
        line(root, phase + pitch / 2, bore, phase + pitch / 2, "sector_right")
        arc(bore, phase + pitch / 2, phase - pitch / 2, "bore")
        line(bore, phase - pitch / 2, root, phase - pitch / 2, "sector_left")
    if not wire.IsDone():
        raise ArithmeticError("直齿截面闭线构造失败")
    result = wire.Wire()
    # MakeWire合并公共端点时可能复制边；历史查询须使用成线后的实际拓扑边。
    explorer = ocp["BRepTools_WireExplorer"](result)
    actual = []
    while explorer.More():
        actual.append(explorer.Current())
        explorer.Next()
    if len(actual) != len(groups):
        raise ArithmeticError("截面实际边数与物理区域标签不一致")
    return result, [(label, edge) for (label, _), edge in zip(groups, actual)]


def _build(value, settings, member, tooth_index, full):
    if not isinstance(settings, StraightCadSettings):
        raise TypeError("settings须为StraightCadSettings")
    surface = _surface(value, member)
    if isinstance(tooth_index, bool) or not isinstance(tooth_index, int) or not 0 <= tooth_index < surface.teeth:
        raise ValueError("tooth_index须为齿数范围内的0起始整数")
    if full and surface.teeth * settings.profile_count * settings.width_sections > 1_000_000:
        raise ValueError("整轮采样规模超过一百万点，请降低设置")
    low, high = _ends(surface, settings)
    sections = [_section(surface, settings, low + (high - low) * i / (settings.width_sections - 1))
                for i in range(settings.width_sections)]
    ocp = _ocp()
    loft = ocp["BRepOffsetAPI_ThruSections"](True, True, settings.tolerance_m * 1000)
    loft.CheckCompatibility(False)
    edge_groups = []
    for index in range(settings.width_sections):
        z = low + (high - low) * index / (settings.width_sections - 1)
        wire, groups = _wire(surface, settings, z, 0 if full else tooth_index,
                             surface.teeth if full else 1, ocp, not full)
        loft.AddWire(wire)
        if index < settings.width_sections - 1:
            edge_groups.extend(groups)
    loft.Build()
    if not loft.IsDone():
        raise ArithmeticError("直齿闭实体放样未完成")
    shape = loft.Shape()
    cut = None
    if full:
        margin = max(high - low, .001)
        cylinder = ocp["BRepPrimAPI_MakeCylinder"](
            ocp["gp_Ax2"](ocp["gp_Pnt"](0, 0, (low - margin) * 1000), ocp["gp_Dir"](0, 0, 1)),
            settings.bore_diameter_m * 500, (high - low + 2 * margin) * 1000).Shape()
        cut = ocp["BRepAlgoAPI_Cut"](shape, cylinder)
        cut.Build()
        if not cut.IsDone():
            raise ArithmeticError("整轮通孔布尔切削未完成")
        shape = cut.Shape()
    from .cad_import import _describe, _ocp as import_ocp
    diagnostics = _describe(shape, import_ocp())
    if not diagnostics["closed_solid"] or diagnostics["solids"] != 1 or not diagnostics["volume_m3"] > 0:
        raise ArithmeticError("直齿模型未形成单一有效闭实体，未返回不可靠FE对象")
    mapping = ocp["TopTools_IndexedMapOfShape"]()
    ocp["TopExp"].MapShapes_s(shape, ocp["TopAbs_FACE"], mapping)
    face_groups = {}
    for label, edge in edge_groups:
        face = loft.GeneratedFace(edge)
        if not face.IsNull():
            index = mapping.FindIndex(face)
            if index:
                face_groups.setdefault(label, []).append(index - 1)
            elif cut is not None:
                face_groups.setdefault(label, []).extend(_modified_face_indices(cut,face,mapping))
    for label, face in (("end_inner", loft.FirstShape()), ("end_outer", loft.LastShape())):
        index = mapping.FindIndex(face)
        if index:
            face_groups[label] = [index - 1]
        elif cut is not None:
            face_groups[label] = _modified_face_indices(cut,face,mapping)
    face_groups = {key: sorted(set(indices)) for key, indices in face_groups.items()}
    if full:
        covered = {index for indices in face_groups.values() for index in indices}
        bore_faces = sorted(set(range(mapping.Extent())) - covered)
        if len(bore_faces) != 1:
            raise ArithmeticError("整轮通孔面与放样区域历史不一致，不能可靠标注边界")
        face_groups["bore"] = bore_faces
    pitch = 2 * math.pi / surface.teeth
    boundaries = {"units": "m", "axis_origin_m": [0., 0., 0.], "axis_direction": [0., 0., 1.],
                  "bore": {"type": "cylinder", "radius_m": settings.bore_diameter_m / 2,
                           "axial_min_m": low, "axial_max_m": high,
                           "suggested_constraint": "仅显式选择孔壁固定；端面/扇区侧面默认不固定"},
                  "end_inner": {"type": "plane", "z_m": low}, "end_outer": {"type": "plane", "z_m": high},
                  "sector_angles_rad": None if full else [(tooth_index - .5) * pitch, (tooth_index + .5) * pitch],
                  "section_bounds": [{"z_m": low + (high - low) * i / (settings.width_sections - 1), **section}
                                     for i, section in enumerate(sections)]}
    diagnostics.update(model="straight_spherical_involute_engineering_solid", member=surface.member,
                       object_scope="full_gear" if full else "single_tooth_sector", teeth=surface.teeth,
                       tooth_index=None if full else tooth_index, profile_count=settings.profile_count,
                       width_sections=settings.width_sections, tolerance_m=settings.tolerance_m,
                       engineering_root_extension_m=settings.engineering_root_extension_m,
                       base_margin_rad=settings.base_margin_rad,
                       requested_cutting_root_radius_m=surface.requested_root_fillet_radius_m,
                       end_boundary="axial_planes", root_construction="radial_engineering_extension_without_cutting_fillet",
                       width_interpolation="piecewise_ruled", profile_interpolation="interpolating_cubic_bspline",
                       face_groups_scope="当前BRep唯一面序号；STEP回读须按几何选择器重新识别")
    diagnostics["limitations"] = ["根部为工程径向连接，根应力不能当作真实切削圆角应力",
        "工作侧面由解析点插值并分段直纹放样；采样数不是精度保证",
        "端界为轴向平面，与球面端界开放齿面网格不同",
        "未校核加工包络、安装干涉、配对共轭或实体材料疲劳性能"]
    return StraightCadModel(shape, diagnostics, face_groups, boundaries, surface, settings, None if full else tooth_index)


def build_straight_tooth_sector(value, settings: StraightCadSettings, *, member="pinion", tooth_index=0) -> StraightCadModel:
    """构建一个齿距扇区：完整单齿、工程根区和孔壁之间的闭合轮体。"""
    return _build(value, settings, member, tooth_index, False)


def build_straight_gear(value, settings: StraightCadSettings, *, member="pinion") -> StraightCadModel:
    """整轮周期轮廓一次放样并切通孔，不逐齿布尔融合。"""
    return _build(value, settings, member, 0, True)


def build_straight_pair(geometry: GeometryResult, pinion_settings: StraightCadSettings,
                        wheel_settings: StraightCadSettings, *, full_gears=True,
                        pinion_phase_rad=0., wheel_phase_rad=None) -> StraightCadPair:
    """按共同锥顶和轴交角放置两轮，默认大轮错开半齿距。

    这是理想宏观安装，不优化相位、不消除工程根/端界干涉，也不求TCA。
    full_gears=False只安装每轮0号单齿扇区，便于独立FE算例。
    placements为局部SI坐标到安装SI坐标的4×4刚体矩阵。
    """
    if not isinstance(geometry, GeometryResult):
        raise TypeError("成对安装须提供GeometryResult，以明确共同轴交角")
    if not isinstance(full_gears, bool):
        raise ValueError("full_gears须为布尔值")
    if wheel_phase_rad is None:
        wheel_phase_rad = math.pi / geometry.wheel.teeth
    for value in (pinion_phase_rad, wheel_phase_rad):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("两轮相位须为有限角度rad")
    builder = build_straight_gear if full_gears else build_straight_tooth_sector
    pinion = builder(geometry, pinion_settings, member="pinion")
    wheel = builder(geometry, wheel_settings, member="wheel")
    from OCP.BRep import BRep_Builder
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.TopoDS import TopoDS_Compound
    from OCP.gp import gp_Trsf
    shaft = geometry.inputs.shaft_angle_rad
    c, s = math.cos(shaft), math.sin(shaft)
    bases = {"pinion": ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.)),
             "wheel": ((-c, 0., s), (0., -1., 0.), (s, 0., c))}
    compound, topology = TopoDS_Compound(), BRep_Builder()
    topology.MakeCompound(compound)
    placements = {}
    for name, model, phase in (("pinion", pinion, pinion_phase_rad), ("wheel", wheel, wheel_phase_rad)):
        x, y, z = bases[name]
        cosine, sine = math.cos(phase), math.sin(phase)
        columns = (tuple(cosine * x[i] + sine * y[i] for i in range(3)),
                   tuple(-sine * x[i] + cosine * y[i] for i in range(3)), z)
        matrix = [[columns[j][i] for j in range(3)] + [0.] for i in range(3)] + [[0., 0., 0., 1.]]
        transform = gp_Trsf()
        transform.SetValues(*(value for row in matrix[:3] for value in row))
        topology.Add(compound, BRepBuilderAPI_Transform(model.shape, transform, True).Shape())
        placements[name] = matrix
    from .cad_import import _describe, _ocp as import_ocp
    diagnostics = _describe(compound, import_ocp())
    diagnostics.update(model="straight_spherical_involute_engineering_pair", units="m",
                       pinion_phase_rad=pinion_phase_rad, wheel_phase_rad=wheel_phase_rad,
                       installation_scope="共同锥顶理想安装；不求接触相位或共轭误差",
                       body_overlap_checked=False)
    return StraightCadPair(compound, pinion, wheel, placements, diagnostics)


def sample_straight_contact_points(model: StraightCadModel, *, width_count=9, profile_count=17,
                                   flank="right", tooth_index=0, edge_margin=0.02) -> dict:
    """按本CAD的轴向端界采样解析工作面/外法向，供FE加载点定位。

    返回点不是CAD投影结果；调用者必须独立检查点到CAD工作面的距离。
    edge_margin为两方向归一化避边比例，不在工程根连接上放接触点。
    """
    if not isinstance(model, StraightCadModel):
        raise TypeError("需要StraightCadModel")
    for name, value in (("width_count", width_count), ("profile_count", profile_count)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 2:
            raise ValueError(f"{name}须为不小于2的整数")
    if isinstance(edge_margin, bool) or not isinstance(edge_margin, (int, float)) or not math.isfinite(edge_margin) or not 0 <= edge_margin < .5:
        raise ValueError("edge_margin须在[0,0.5)内")
    surface, settings = model.surface, model.settings
    if flank not in ("left", "right"):
        raise ValueError("工作面须为left或right")
    if isinstance(tooth_index, bool) or not isinstance(tooth_index, int) or not 0 <= tooth_index < surface.teeth:
        raise ValueError("tooth_index须在齿数范围内")
    if model.tooth_index is not None and tooth_index != model.tooth_index:
        raise ValueError("单齿模型只能采样其实际tooth_index")
    low, high = _ends(surface, settings)
    phase = tooth_index * 2 * math.pi / surface.teeth
    cosine, sine = math.cos(phase), math.sin(phase)
    points, normals, uv = [], [], []
    for i in range(width_count):
        u = edge_margin + (1 - 2 * edge_margin) * i / (width_count - 1)
        z = low + u * (high - low)
        section = _section(surface, settings, z)
        for j in range(profile_count):
            v = edge_margin + (1 - 2 * edge_margin) * j / (profile_count - 1)
            theta = section["lower_polar_rad"] + v * (section["upper_polar_rad"] - section["lower_polar_rad"])
            point = surface.profile.evaluate(z / math.cos(theta), theta, flank)
            normal = surface.profile.normal(theta, flank)
            points.append((cosine * point[0] - sine * point[1], sine * point[0] + cosine * point[1], point[2]))
            normals.append((cosine * normal[0] - sine * normal[1], sine * normal[0] + cosine * normal[1], normal[2]))
            uv.append((u, v))
    return {"units": "m", "points_m": points, "normals": normals, "uv": uv,
            "flank": flank, "tooth_index": tooth_index, "shape": [width_count, profile_count],
            "scope": "解析渐开线采样；不是有限元节点或已投影CAD节点"}
