"""实际刀盘展成槽面裁切回转齿坯，生成有孔闭实体；公开SI、OCP内部mm。

侧面与刀尖根过渡沿用逐点展成内核。槽底仍是根圆连接约定；轮体、
轴孔和端界是显式齿坯参数，不宣称全部槽底刀具包络或制造认证。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .generated_surface import (GeneratedGearInput, _request, build_generated_surface_grid,
                                build_generated_accuracy_grid)


@dataclass(frozen=True)
class GeneratedBlankSettings:
    """有孔八点齿坯；支承深度为根高倍数，外端顶径修边沿节锥法向。

    compatible保留既有非正根高度诊断；positive_only遇到非正区间直接
    拒绝，不把高度夹零。布尔容差只用于几何裁切，不作为工艺精度。
    """
    bore_diameter_m: float
    support_depth_factor: float = 2.
    outer_tip_trim_m: float = 1e-5
    boolean_tolerance_m: float = 1e-8
    sewing_tolerance_m: float = 1e-9
    root_height_policy: str = "compatible"
    partition_relative_tolerance: float = 1e-4
    periodic_relative_tolerance: float = 1e-3

    def __post_init__(self):
        for name in ("bore_diameter_m", "boolean_tolerance_m", "sewing_tolerance_m",
                     "partition_relative_tolerance", "periodic_relative_tolerance"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name}须为正有限实数")
        if isinstance(self.support_depth_factor, bool) or not isinstance(self.support_depth_factor, (int, float)) or not math.isfinite(self.support_depth_factor) or self.support_depth_factor <= 1:
            raise ValueError("支承深度倍数须大于1，确保根圆之下有正厚度轮体")
        if isinstance(self.outer_tip_trim_m, bool) or not isinstance(self.outer_tip_trim_m, (int, float)) or not math.isfinite(self.outer_tip_trim_m) or self.outer_tip_trim_m < 0:
            raise ValueError("外端齿顶修边量须为非负有限长度m")
        if self.root_height_policy not in ("compatible", "positive_only"):
            raise ValueError("根高度策略须为compatible或positive_only")
        if self.partition_relative_tolerance > 1e-3 or self.periodic_relative_tolerance > .01:
            raise ValueError("体积验收门槛过宽，不能以放宽容差掩盖裁切失败")


@dataclass
class GeneratedBlank:
    """回转齿坯实体，剖面为SI的(radius,z)八点。"""
    shape: object
    diagnostics: dict
    profile_m: tuple
    request: GeneratedGearInput
    settings: GeneratedBlankSettings


@dataclass
class GeneratedSolidModel:
    """shape是保留材料；slot_shapes为齿坯裁切出的闭合槽腔，不是齿体。

    slot_shapes按体积降序存放，不代表齿号；slot_face_groups才按阵列齿号
    标记。face_groups使用当前BRep唯一面0起始序号，STEP回读后须重识别。
    五面、原网格和齿坯随对象保存，供独立点面检验与分区展示。
    """
    shape: object
    diagnostics: dict
    blank: GeneratedBlank
    shell: object
    grid: dict
    slot_shapes: tuple
    face_groups: dict
    slot_face_groups: dict

    @property
    def request(self):
        """返回原显式机床/齿轮请求；冻结数据类不自动补机床设置。"""
        return self.blank.request

    @property
    def settings(self):
        """返回显式齿坯参数；不与刀尖圆弧半径混用。"""
        return self.blank.settings


def _ocp():
    try:
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Transform
        from OCP.BRepPrimAPI import BRepPrimAPI_MakeRevol
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Splitter
        from OCP.TopTools import TopTools_ListOfShape, TopTools_IndexedMapOfShape
        from OCP.TopExp import TopExp
        from OCP.TopAbs import TopAbs_SOLID, TopAbs_FACE
        from OCP.TopoDS import TopoDS
        from OCP.gp import gp_Pnt, gp_Ax1, gp_Dir, gp_Trsf
        from OCP.BRep import BRep_Tool
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
    except ImportError as error:
        raise ImportError("展成闭实体需要可选OCP依赖") from error
    return locals()


def _volume(shape, ocp):
    properties = ocp["GProp_GProps"]()
    estimate = ocp["BRepGProp"].VolumeProperties_s(shape, properties, 1e-9, True, False)
    return abs(properties.Mass()) * 1e-9, estimate


def _describe(shape, ocp):
    from .cad_import import _describe as describe, _ocp as import_ocp
    result = describe(shape, import_ocp())
    if not result["valid"] or not result["closed_solid"] or result["solids"] != 1:
        raise ArithmeticError("制造几何没有形成单一有效闭实体")
    result["nonadaptive_volume_m3"] = result["volume_m3"]
    result["volume_m3"], result["volume_integration_error_estimate"] = _volume(shape, ocp)
    result["volume_integration"] = "自适应2D Gauss，面相对门槛1e-9；误差估计非全域严格上界"
    if not result["volume_m3"] > 0:
        raise ArithmeticError("制造几何没有正体积")
    return result


def _history_face_indices(splitter, face, mapping):
    """复制布尔历史并及时取面号，避免借用列表迭代器及元素的生命周期问题。"""
    from OCP.TopTools import TopTools_ListOfShape
    index = mapping.FindIndex(face)
    indices = [index-1] if index else []
    for history in (splitter.Modified, splitter.Generated):
        owned = TopTools_ListOfShape(history(face))
        while not owned.IsEmpty():
            # First返回的Shape须在RemoveFirst销毁列表节点前立即使用。
            index = mapping.FindIndex(owned.First())
            if index:
                indices.append(index-1)
            owned.RemoveFirst()
    return sorted(set(indices))


def build_generated_blank(request, settings: GeneratedBlankSettings) -> GeneratedBlank:
    """构建八点剖面回转齿坯，轴线上的两点改为显式孔半径。

    内外端根高之下按支承深度倍数延伸；该轮体构造与展成刀具根圆角
    是不同参数。所有顶/根余量均保留在诊断，不偷偷修改机床设置。
    """
    request = _request(request)
    if isinstance(settings, dict):
        settings = GeneratedBlankSettings(**settings)
    if not isinstance(settings, GeneratedBlankSettings):
        raise TypeError("settings须为GeneratedBlankSettings或具名字典")
    machine = request.machine
    delta, width = machine.pitch_cone_angle_rad, request.face_width_m
    reference = machine.reference_cone_distance_m
    outer_radius, inner_radius = reference * math.sin(delta), (reference - width) * math.sin(delta)
    outer_z, inner_z = reference * math.cos(delta), (reference - width) * math.cos(delta)
    inner_ha = request.addendum_m - math.tan(request.tip_cone_angle_rad - delta) * width
    inner_hf = machine.dedendum_m - math.tan(delta - machine.machine_root_angle_rad) * width
    if min(inner_ha, inner_hf, request.addendum_m - settings.outer_tip_trim_m) <= 0:
        raise ValueError("修边后的内外端齿顶高和齿根高必须为正，未自动外推无效齿坯")
    depth, bore = settings.support_depth_factor, settings.bore_diameter_m / 2
    sine, cosine = math.sin(delta), math.cos(delta)
    profile = ((outer_radius + (request.addendum_m-settings.outer_tip_trim_m)*cosine,
                outer_z-(request.addendum_m-settings.outer_tip_trim_m)*sine),
               (inner_radius+inner_ha*cosine, inner_z-inner_ha*sine),
               (inner_radius-inner_hf*cosine, inner_z+inner_hf*sine),
               (inner_radius-depth*inner_hf*cosine, inner_z+depth*inner_hf*sine),
               (bore, inner_z+depth*inner_hf*sine),
               (bore, outer_z+depth*machine.dedendum_m*sine),
               (outer_radius-depth*machine.dedendum_m*cosine, outer_z+depth*machine.dedendum_m*sine),
               (outer_radius-machine.dedendum_m*cosine, outer_z+machine.dedendum_m*sine))
    if min(profile[i][0] for i in (0, 1, 2, 3, 6, 7)) <= bore + 10 * settings.boolean_tolerance_m:
        raise ValueError("轴孔与根部支承层之间没有足够正壁厚")
    if profile[5][1] - profile[4][1] <= 10 * settings.boolean_tolerance_m:
        raise ValueError("齿坯孔壁轴向长度不为正，当前回转剖面无效")
    ocp = _ocp()
    polygon = ocp["BRepBuilderAPI_MakePolygon"]()
    for radius, z in profile:
        polygon.Add(ocp["gp_Pnt"](radius*1000, 0., z*1000))
    polygon.Close()
    face = ocp["BRepBuilderAPI_MakeFace"](polygon.Wire())
    if not face.IsDone():
        raise ArithmeticError("齿坯剖面建面失败")
    axis = ocp["gp_Ax1"](ocp["gp_Pnt"](0., 0., 0.), ocp["gp_Dir"](0., 0., 1.))
    revolve = ocp["BRepPrimAPI_MakeRevol"](face.Face(), axis, 2*math.pi, True)
    if not revolve.IsDone():
        raise ArithmeticError("齿坯回转未完成")
    shape = revolve.Shape()
    diagnostics = _describe(shape, ocp)
    diagnostics.update(model="generated_eight_point_revolved_blank", bore_radius_m=bore,
        bore_axial_range_m=(profile[4][1], profile[5][1]), profile_m=profile,
        outer_tip_trim_m=settings.outer_tip_trim_m, support_depth_factor=depth)
    diagnostics["limitations"] = ["轮体、轴孔与端界为显式齿坯边界，未建轴、轴承或配合连接",
                                  "齿坯回转不是展成切齿，须实际裁切并核验材料移除"]
    return GeneratedBlank(shape, diagnostics, profile, request, settings)


def _build(request, settings, width_count, profile_count, full):
    blank = build_generated_blank(request, settings)
    request, settings = blank.request, blank.settings
    grid = build_generated_surface_grid(request, width_count=width_count, profile_count=profile_count)
    if settings.root_height_policy == "positive_only" and grid["metadata"]["minimum_root_cutter_height_upper_m"] <= 0:
        raise ValueError("positive_only拒绝非正根高度区间；未夹零或替换既有内核结果")
    from .cad import build_generated_shell
    shell = build_generated_shell(grid, sewing_tolerance_m=settings.sewing_tolerance_m)
    ocp = _ocp()
    original_faces = ocp["TopTools_IndexedMapOfShape"]()
    ocp["TopExp"].MapShapes_s(shell.shape, ocp["TopAbs_FACE"], original_faces)
    labels = {}
    for i in range(1, original_faces.Extent()+1):
        face = ocp["TopoDS"].Face_s(original_faces.FindKey(i))
        surface = ocp["BRep_Tool"].Surface_s(face)
        matches = [name for name, reference in shell.surfaces.items() if surface == reference]
        if len(matches) != 1:
            raise ArithmeticError("五面壳的区域标签与实际曲面句柄不一致")
        labels[matches[0]] = face
    count = request.tooth_count if full else 1
    arguments, tools = ocp["TopTools_ListOfShape"](), ocp["TopTools_ListOfShape"]()
    arguments.Append(blank.shape)
    rotated_faces = {}
    axis = ocp["gp_Ax1"](ocp["gp_Pnt"](0., 0., 0.), ocp["gp_Dir"](0., 0., 1.))
    for k in range(count):
        transform = ocp["gp_Trsf"]()
        transform.SetRotation(axis, k*2*math.pi/request.tooth_count)
        rotated = ocp["BRepBuilderAPI_Transform"](shell.shape, transform, True)
        tools.Append(rotated.Shape())
        rotated_faces[k] = {name: rotated.ModifiedShape(face) for name, face in labels.items()}
    splitter = ocp["BRepAlgoAPI_Splitter"]()
    splitter.SetArguments(arguments)
    splitter.SetTools(tools)
    splitter.SetFuzzyValue(settings.boolean_tolerance_m*1000)
    splitter.Build()
    if not splitter.IsDone():
        raise ArithmeticError("周期展成槽壳拆分齿坯未完成")
    solids = ocp["TopTools_IndexedMapOfShape"]()
    ocp["TopExp"].MapShapes_s(splitter.Shape(), ocp["TopAbs_SOLID"], solids)
    if solids.Extent() != count+1:
        raise ArithmeticError(f"期望材料与{count}个闭合槽腔，实际{solids.Extent()}个实体；拒绝未切透、重复或交叠槽")
    partitions = []
    for i in range(1, solids.Extent()+1):
        shape = solids.FindKey(i)
        diagnostics = _describe(shape, ocp)
        partitions.append((diagnostics["volume_m3"], shape, diagnostics))
    partitions.sort(key=lambda item: item[0], reverse=True)
    kept_volume, shape, diagnostics = partitions[0]
    removed = partitions[1:]
    blank_volume = blank.diagnostics["volume_m3"]
    partition_error = abs(math.fsum(item[0] for item in partitions)/blank_volume - 1)
    removed_volume = math.fsum(item[0] for item in removed)
    spread = (max(item[0] for item in removed)-min(item[0] for item in removed))/(removed_volume/count)
    if partition_error > settings.partition_relative_tolerance:
        raise ArithmeticError("拆分体积不守恒，不能以最大实体掩盖布尔失败")
    if removed_volume <= blank_volume*1e-7 or kept_volume <= blank_volume*.1:
        raise ArithmeticError("齿槽没有移除足够材料或轮体保留过少，拒绝无效裁切")
    if spread > settings.periodic_relative_tolerance:
        raise ArithmeticError("周期齿槽体积不一致，拒绝不完整槽或重复裁切")
    face_map = ocp["TopTools_IndexedMapOfShape"]()
    ocp["TopExp"].MapShapes_s(shape, ocp["TopAbs_FACE"], face_map)
    groups, by_slot = {name: [] for name in labels}, {}
    for k, faces in rotated_faces.items():
        by_slot[k] = {}
        for name, face in faces.items():
            indices = _history_face_indices(splitter, face, face_map)
            if not indices:
                raise ArithmeticError(f"第{k}槽{name}未进入保留材料边界，拒绝不完整展成区域")
            groups[name].extend(indices)
            by_slot[k][name] = indices
    groups = {name: sorted(set(indices)) for name, indices in groups.items()}
    known = {index for indices in groups.values() for index in indices}
    groups["blank"] = sorted(set(range(face_map.Extent()))-known)
    diagnostics.update(model="generated_cutter_revolved_blank_split", object_scope="full_gear" if full else "single_slot_cut",
        tooth_count=request.tooth_count, patterned_slot_count=count, blank_volume_m3=blank_volume,
        removed_volume_m3=removed_volume, removed_slot_volumes_m3=[item[0] for item in removed],
        partition_solid_count=len(partitions), partition_relative_volume_error=partition_error,
        periodic_slot_relative_volume_spread=spread, root_height_policy=settings.root_height_policy,
        minimum_root_cutter_height_upper_m=grid["metadata"]["minimum_root_cutter_height_upper_m"],
        negative_root_height_interval_sections=grid["metadata"]["negative_root_height_interval_sections"],
        root_bottom_model=grid["metadata"]["root_bottom_model"],
        bore=blank.diagnostics["bore_radius_m"], width_count=width_count, profile_count=profile_count,
        boolean_tolerance_m=settings.boolean_tolerance_m, sewing_tolerance_m=settings.sewing_tolerance_m,
        construction="实际展成五面→有孔八点回转齿坯→周期槽壳拆分→逐实体/体积/区域验收",
        face_group_scope="当前BRep唯一面0起始编号；STEP回读不保证编号不变")
    diagnostics["limitations"] = ["工作侧面和刀尖根过渡来自显式机床包络；不自动设计配对机床设置",
        "槽底仍为根径圆连接约定，不宣称槽底全部刀具包络或完整工艺认证",
        "兼容策略保留非正根高度区间，不把布尔有效性当作该区间物理有效证明",
        "轴孔、轮体和端界属于显式齿坯构造，不代表实际毛坯或装夹设计",
        "点插值、布尔与体积验收不是制造精度、配对共轭、TCA或LTCA验证"]
    return GeneratedSolidModel(shape, diagnostics, blank, shell, grid,
                               tuple(item[1] for item in removed), groups, by_slot)


def build_generated_gear(request, settings, *, width_count=25, profile_count=25) -> GeneratedSolidModel:
    """按显式齿数周期裁切整轮；每个槽必须形成独立闭合切除腔。"""
    return _build(request, settings, width_count, profile_count, True)


def build_generated_slot(request, settings, *, width_count=25, profile_count=25) -> GeneratedSolidModel:
    """只裁一个槽，返回带单槽齿坯及一个闭合槽腔，不把槽腔称齿体。"""
    return _build(request, settings, width_count, profile_count, False)


def measure_generated_solid_accuracy(model: GeneratedSolidModel) -> dict:
    """独立27点至实际裁切实体边界距离，含真实侧面和复合根区。

    沿用非插值节点的Gauss位置；有限点检查不构成全域Hausdorff界。
    """
    if not isinstance(model, GeneratedSolidModel):
        raise TypeError("需要GeneratedSolidModel")
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.gp import gp_Pnt
    grid = build_generated_accuracy_grid(model.blank.request)
    groups = {}
    for name, group in grid["groups"].items():
        distances = []
        for point in group["points_m"]:
            vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*(1000*x for x in point))).Vertex()
            query = BRepExtrema_DistShapeShape(vertex, model.shape)
            if not query.IsDone():
                raise ArithmeticError("独立精度点至裁切实体边界查询失败")
            distances.append(query.Value()/1000)
        groups[name] = {"count": len(distances), "maximum_distance_m": max(distances), "distances_m": distances}
    return {"units": "m", "groups": groups, "count": sum(group["count"] for group in groups.values()),
            "maximum_distance_m": max(group["maximum_distance_m"] for group in groups.values()),
            "scope": "指定27个独立点到实际裁切BRep边界；非全域、工艺或接触精度证明"}
