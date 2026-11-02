"""可选OCP曲面后端：真实五面网格建壳、STEP/BREP及三角化。

模块导入仅用标准库；调用CAD功能时才需要OCP。输入/返回坐标m，
CAD内部采用mm。后端不生成缺失圆角、不自动封闭实体、不导入CadQuery。
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path


_PATCHES = ("left_flank", "left_fillet", "root_bottom", "right_fillet", "right_flank")


@dataclass
class GeneratedCadShell:
    """持有可选内核对象及SI诊断；shape是开放曲面壳，不是齿轮实体。"""
    shape: object
    surfaces: dict
    diagnostics: dict


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name}须为正有限实数")


def _point(value):
    if len(value) != 3 or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in value):
        raise ValueError("曲面点须为三个有限SI坐标")
    return tuple(v*1000 for v in value)


def _ocp():
    try:
        from OCP.Approx import Approx_IsoParametric
        from OCP.gp import gp_Pnt
        from OCP.TColgp import TColgp_Array2OfPnt
        from OCP.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.TopExp import TopExp, TopExp_Explorer
        from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID, TopAbs_REVERSED
        from OCP.TopTools import TopTools_IndexedMapOfShape
        from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.Interface import Interface_Static
        from OCP.BRepTools import BRepTools
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.BRep import BRep_Tool
        from OCP.TopLoc import TopLoc_Location
        from OCP.TopoDS import TopoDS
    except ImportError as error:
        raise ImportError("CAD功能需要可选OCP依赖；标准库几何网格可独立使用") from error
    return locals()


def build_generated_shell(grid: dict, *, sewing_tolerance_m=1e-9) -> GeneratedCadShell:
    """以统一均匀参数三次C2插值建五面并缝合，不扩展为封闭实体。

    网格内的根过渡点必须由数学内核先行提供。本函数没有通用圆角
    补面；独立插值策略也不承诺与任意其他建模器的内部算法相同。
    """
    _positive(sewing_tolerance_m, "sewing_tolerance_m")
    if not isinstance(grid, dict) or grid.get("units") != "m" or set(grid.get("patches", {})) != set(_PATCHES):
        raise ValueError("CAD入口须为单位m且含完整五面的展成网格")
    geometry = _ocp()
    surfaces = {}
    sewing = geometry["BRepBuilderAPI_Sewing"](sewing_tolerance_m*1000, True, True, True, False)
    dimensions = None
    for name in _PATCHES:
        patch = grid["patches"][name]
        width, profile = patch["width_count"], patch["profile_count"]
        if any(isinstance(v, bool) or not isinstance(v, int) or v < 4 for v in (width, profile)):
            raise ValueError("三次CAD插值的宽向/齿廓点数均须为不小于4的整数")
        if dimensions is not None and dimensions != (width, profile):
            raise ValueError("五面的宽向/齿廓点数必须相同，以保持共享边参数一致")
        dimensions = width, profile
        if len(patch["points_m"]) != width*profile:
            raise ValueError("曲面点数组长度与网格维度不符")
        points = geometry["TColgp_Array2OfPnt"](1, width, 1, profile)
        for u in range(width):
            for v in range(profile):
                points.SetValue(u+1, v+1, geometry["gp_Pnt"](*_point(patch["points_m"][v*width+u])))
        interpolation = geometry["GeomAPI_PointsToBSplineSurface"]()
        interpolation.Interpolate(points, geometry["Approx_IsoParametric"], False)
        if not interpolation.IsDone():
            raise ValueError(f"{name}的曲面插值失败")
        surface = interpolation.Surface()
        surfaces[name] = surface
        face = geometry["BRepBuilderAPI_MakeFace"](surface, sewing_tolerance_m*1000)
        if not face.IsDone():
            raise ValueError(f"{name}的曲面建面失败")
        sewing.Add(face.Face())
    sewing.Perform()
    shape = sewing.SewedShape()
    if shape.IsNull():
        raise ValueError("五面缝合返回空形状")
    counts = {}
    for label, kind in (("faces", geometry["TopAbs_FACE"]), ("edges", geometry["TopAbs_EDGE"]), ("solids", geometry["TopAbs_SOLID"])):
        mapping = geometry["TopTools_IndexedMapOfShape"]()
        geometry["TopExp"].MapShapes_s(shape, kind, mapping)
        counts[label] = mapping.Extent()
    valid = geometry["BRepCheck_Analyzer"](shape).IsValid()
    if not valid or counts["faces"] != 5 or counts["solids"] != 0:
        raise ValueError("五面未构成有效开放壳，未把异常结果作为齿轮实体返回")
    free = sewing.NbFreeEdges()
    shared = counts["edges"]-free
    if shared != 4 or sewing.NbMultipleEdges() != 0:
        raise ValueError("五面未形成预期的四条共享边，请检查接边一致性或缝合容差")
    return GeneratedCadShell(shape, surfaces, {**counts, "valid": valid, "free_edges": free, "shared_edges": shared,
        "units": "m", "cad_internal_length_unit": "mm", "closed_solid": False, "sewing_tolerance_m": sewing_tolerance_m,
        "surface_interpolation": "统一均匀参数的三次C2插值",
        "continuity_scope": "各面内部C2；共享边已检查位置和拓扑，未将其宣称为跨面C1/C2认证", "metadata": dict(grid.get("metadata", {})),
        "limitations": ["未执行齿坯裁切、周期加工或实体封闭", "网格插值不替代实际机床设置、制造或接触检查"]})


def _shell(value):
    if not isinstance(value, GeneratedCadShell):
        raise TypeError("需要build_generated_shell返回的GeneratedCadShell")
    return value


def export_step(shell: GeneratedCadShell, path) -> Path:
    """导出毫米单位STEP；公开输入网格及返回三角网格仍为m。"""
    shell = _shell(shell)
    geometry = _ocp()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    settings = geometry["Interface_Static"]
    previous = settings.CVal_s("write.step.unit")
    settings.SetCVal_s("write.step.unit", "MM")
    try:
        writer = geometry["STEPControl_Writer"]()
        if writer.Transfer(shell.shape, geometry["STEPControl_AsIs"]) != geometry["IFSelect_RetDone"]:
            raise OSError("STEP形状传递失败")
        if writer.Write(str(target)) != geometry["IFSelect_RetDone"]:
            raise OSError("STEP文件写出失败")
    finally:
        if previous:
            settings.SetCVal_s("write.step.unit", previous)
    return target


def export_brep(shell: GeneratedCadShell, path) -> Path:
    """导出BREP；格式本身无单位，文件坐标为mm，须保留调用诊断。"""
    shell = _shell(shell)
    geometry = _ocp()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not geometry["BRepTools"].Write_s(shell.shape, str(target)):
        raise OSError("BREP文件写出失败")
    return target


def measure_surface_accuracy(shell: GeneratedCadShell, accuracy_grid: dict) -> dict:
    """检查均匀插值节点之外的SI点面距离；复合根区取底部三面的最小值。"""
    shell = _shell(shell)
    if not isinstance(accuracy_grid, dict) or accuracy_grid.get("units") != "m":
        raise ValueError("精度点须以m为单位")
    geometry = _ocp()
    groups = {}
    for name, group in accuracy_grid["groups"].items():
        if name == "root_transition":
            surfaces = [shell.surfaces[key] for key in ("left_fillet", "root_bottom", "right_fillet")]
        elif name in shell.surfaces:
            surfaces = [shell.surfaces[name]]
        else:
            raise ValueError("未知精度点区域")
        distances = [min(geometry["GeomAPI_ProjectPointOnSurf"](geometry["gp_Pnt"](*_point(point)), surface).LowerDistance()
                         for surface in surfaces)/1000 for point in group["points_m"]]
        if not distances:
            raise ValueError("精度点组不能为空")
        groups[name] = {"count": len(distances), "distances_m": distances, "maximum_error_m": max(distances),
                        "mean_error_m": sum(distances)/len(distances)}
    if not groups:
        raise ValueError("精度点组不能为空")
    return {"units": "m", "groups": groups, "count": sum(group["count"] for group in groups.values()),
            "maximum_error_m": max(group["maximum_error_m"] for group in groups.values()),
            "scope": "所给离散精度点到未裁切曲面距离，不是全域、制造或接触误差保证"}


def triangulate(shell: GeneratedCadShell, *, linear_deflection_m=5e-5, angular_deflection_rad=.1) -> dict:
    """三角化已有真实曲面壳，返回m及0起始索引；不增补实体边界。"""
    shell = _shell(shell)
    _positive(linear_deflection_m, "linear_deflection_m")
    _positive(angular_deflection_rad, "angular_deflection_rad")
    geometry = _ocp()
    geometry["BRepMesh_IncrementalMesh"](shell.shape, linear_deflection_m*1000, False, angular_deflection_rad, False).Perform()
    vertices, triangles, face_ids = [], [], []
    explorer = geometry["TopExp_Explorer"](shell.shape, geometry["TopAbs_FACE"])
    face_index = 0
    while explorer.More():
        face = geometry["TopoDS"].Face_s(explorer.Current())
        location = geometry["TopLoc_Location"]()
        mesh = geometry["BRep_Tool"].Triangulation_s(face, location)
        if mesh is not None:
            offset = len(vertices)
            transform = location.Transformation()
            for index in range(1, mesh.NbNodes()+1):
                point = mesh.Node(index).Transformed(transform)
                vertices.append((point.X()/1000, point.Y()/1000, point.Z()/1000))
            for index in range(1, mesh.NbTriangles()+1):
                ids = list(mesh.Triangle(index).Get())
                if face.Orientation() == geometry["TopAbs_REVERSED"]:
                    ids[1], ids[2] = ids[2], ids[1]
                triangles.append(tuple(offset+i-1 for i in ids))
                face_ids.append(face_index)
        face_index += 1
        explorer.Next()
    if not triangles:
        raise ValueError("CAD三角化没有返回有效三角形")
    return {"units": "m", "vertices_m": vertices, "triangles": triangles, "face_indices": face_ids,
            "metadata": {"closed_solid": False, "linear_deflection_m": linear_deflection_m,
                         "angular_deflection_rad": angular_deflection_rad,
                         "face_index_scope": "当前CAD面的遍历序号，不作为物理左右侧编号",
                         "orientation": "沿CAD曲面壳方向，不保证材料外法向"}}
