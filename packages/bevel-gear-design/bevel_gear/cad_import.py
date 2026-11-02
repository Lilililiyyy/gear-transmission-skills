"""中性STEP模型导入：按需使用OCP，形状内部统一为mm，公开度量为SI。

保留输入实际拓扑，不修补成实体；文件路径、哈希与原始头信息由调用方
私有记录。此模块只供几何交换和给定位姿查询，不求齿面或TCA/LTCA。
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path


@dataclass
class ImportedCadModel:
    """shape的内部长度单位为mm；诊断不保存来源身份或原始STEP头。"""
    shape: object
    diagnostics: dict


def _ocp():
    """导入本模块无需OCP，只有实际CAD调用才加载可选依赖。"""
    try:
        from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.Interface import Interface_Static
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRep import BRep_Tool
        from OCP.TopExp import TopExp, TopExp_Explorer
        from OCP.TopAbs import (TopAbs_COMPOUND, TopAbs_COMPSOLID, TopAbs_SOLID, TopAbs_SHELL,
                                TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX,
                                TopAbs_REVERSED)
        from OCP.TopTools import TopTools_IndexedMapOfShape
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        from OCP.BRepTools import BRepTools
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.TopLoc import TopLoc_Location
        from OCP.TopoDS import TopoDS
    except ImportError as error:
        raise ImportError("STEP导入和三角化需要可选OCP依赖") from error
    return locals()


def _describe(shape, ocp):
    """按真实拓扑分类；有面不等于实体，开放面不计算伪体积。"""
    maps = {}
    kinds = {"solids": "SOLID", "shells": "SHELL", "faces": "FACE", "edges": "EDGE", "vertices": "VERTEX"}
    for label, kind in kinds.items():
        mapping = ocp["TopTools_IndexedMapOfShape"]()
        ocp["TopExp"].MapShapes_s(shape, ocp["TopAbs_" + kind], mapping)
        maps[label] = mapping
    valid = ocp["BRepCheck_Analyzer"](shape).IsValid()
    types = {ocp["TopAbs_" + key]: value for key, value in
             (("COMPOUND", "compound"), ("COMPSOLID", "compsolid"), ("SOLID", "solid"),
              ("SHELL", "shell"), ("FACE", "face"), ("WIRE", "wire"), ("EDGE", "edge"), ("VERTEX", "vertex"))}
    covered = {key: ocp["TopTools_IndexedMapOfShape"]() for key in ("faces", "edges", "shells", "vertices")}
    volumes = []
    for index in range(1, maps["solids"].Extent() + 1):
        solid = maps["solids"].FindKey(index)
        for label in covered:
            ocp["TopExp"].MapShapes_s(solid, ocp["TopAbs_" + kinds[label]], covered[label])
        if valid:
            properties = ocp["GProp_GProps"]()
            ocp["BRepGProp"].VolumeProperties_s(solid, properties)
            volumes.append(abs(properties.Mass()) * 1e-9)
    closed_shells = sum(ocp["BRep_Tool"].IsClosed_s(maps["shells"].FindKey(index))
                        for index in range(1, maps["shells"].Extent() + 1))
    closed_solid = bool(valid and maps["solids"].Extent() > 0 and maps["shells"].Extent() > 0
                        and closed_shells == maps["shells"].Extent()
                        and all(covered[key].Extent() == maps[key].Extent() for key in covered))
    area = None
    if valid:
        areas = []
        for index in range(1, maps["faces"].Extent() + 1):
            properties = ocp["GProp_GProps"]()
            ocp["BRepGProp"].SurfaceProperties_s(maps["faces"].FindKey(index), properties)
            areas.append(abs(properties.Mass()) * 1e-6)
        area = math.fsum(areas)
    return {**{key: mapping.Extent() for key, mapping in maps.items()},
            "shape_type": types.get(shape.ShapeType(), "shape"), "valid": valid,
            "closed_shells": closed_shells, "closed_solid": closed_solid,
            "volume_m3": math.fsum(volumes) if volumes else None, "area_m2": area,
            "units": "m", "kernel_length_unit_m": .001, "cad_internal_length_unit": "mm",
            "volume_definition": "有效实体包围体积绝对值逐体相加；未求联合体，重叠实体会重复累计；无实体或拓扑无效时为null",
            "area_definition": "按含位置变换的唯一拓扑面逐面累计面积，不求几何重叠面并集",
            "limitations": ["导入不修补、封闭或重新生成齿面", "几何有效性和距离不代表制造、共轭或承载接触验证"]}


def load_step_model(path) -> ImportedCadModel:
    """读取中性STEP并显式归一到内部mm；不返回原始文件身份或头字段。"""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError("STEP输入文件不存在或不是普通文件")
    ocp = _ocp()
    reader = ocp["STEPControl_Reader"]()
    if reader.ReadFile(str(source)) != ocp["IFSelect_RetDone"]:
        raise ValueError("STEP读取失败，未返回模型")
    # OCP长度因子以mm为基准；读取后、转换前设置1，避免继承进程全局单位。
    reader.SetSystemLengthUnit(1.)
    roots = reader.NbRootsForTransfer()
    transferred = reader.TransferRoots()
    if roots <= 0 or transferred != roots:
        raise ValueError("STEP没有可转换根或部分根转换失败，未返回不完整模型")
    shape = reader.OneShape()
    if shape.IsNull():
        raise ValueError("STEP转换得到空形状")
    diagnostics = _describe(shape, ocp)
    diagnostics.update(roots_available=roots, roots_transferred=transferred,
                       reader_system_length_unit_mm=reader.SystemLengthUnit())
    return ImportedCadModel(shape, diagnostics)


def _model(model, ocp):
    if not isinstance(model, ImportedCadModel):
        raise TypeError("需要load_step_model返回的ImportedCadModel")
    if model.shape.IsNull() or not ocp["BRepCheck_Analyzer"](model.shape).IsValid():
        raise ValueError("CAD导出或三角化要求非空且拓扑有效的导入模型")
    return model


def triangulate_imported_model(model: ImportedCadModel, *, linear_deflection_m=5e-5,
                                angular_deflection_rad=.1) -> dict:
    """三角化实际BRep，应用位置和面方向；返回SI点及0起始索引。"""
    for name, value in (("linear_deflection_m", linear_deflection_m),
                        ("angular_deflection_rad", angular_deflection_rad)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name}须为正有限实数")
    ocp = _ocp()
    model = _model(model, ocp)
    mesher = ocp["BRepMesh_IncrementalMesh"](model.shape, linear_deflection_m * 1000,
                                           False, angular_deflection_rad, False)
    mesher.Perform()
    if not mesher.IsDone():
        raise ArithmeticError("导入CAD模型的三角化未完成")
    vertices, triangles, face_ids = [], [], []
    explorer = ocp["TopExp_Explorer"](model.shape, ocp["TopAbs_FACE"])
    face_index = 0
    while explorer.More():
        face = ocp["TopoDS"].Face_s(explorer.Current())
        location = ocp["TopLoc_Location"]()
        mesh = ocp["BRep_Tool"].Triangulation_s(face, location)
        if mesh is None or mesh.NbTriangles() == 0:
            raise ArithmeticError("至少一个CAD面未返回有效三角化，未导出不完整网格")
        offset = len(vertices)
        transform = location.Transformation()
        for index in range(1, mesh.NbNodes() + 1):
            point = mesh.Node(index).Transformed(transform)
            vertices.append((point.X() / 1000, point.Y() / 1000, point.Z() / 1000))
        for index in range(1, mesh.NbTriangles() + 1):
            ids = list(mesh.Triangle(index).Get())
            if face.Orientation() == ocp["TopAbs_REVERSED"]:
                ids[1], ids[2] = ids[2], ids[1]
            triangles.append(tuple(offset + value - 1 for value in ids))
            face_ids.append(face_index)
        face_index += 1
        explorer.Next()
    if not triangles:
        raise ValueError("导入模型没有可三角化曲面")
    topology = _describe(model.shape, ocp)
    return {"units": "m", "vertices_m": vertices, "triangles": triangles, "face_indices": face_ids,
            "metadata": {"closed_solid": topology["closed_solid"], "shape_type": topology["shape_type"],
                         "solid_count": topology["solids"], "kernel_length_unit_m": .001,
                         "topology_scope": "闭实体标志来自BRep拓扑；三角网格按面保留接边节点，未执行节点焊接",
                         "linear_deflection_m": linear_deflection_m, "angular_deflection_rad": angular_deflection_rad,
                         "face_index_scope": "当前CAD面遍历序号，不是齿轮工作面或物理区域编号",
                         "orientation": "保留CAD面拓扑方向，不重新推断材料外法向"}}


def export_model_step(model: ImportedCadModel, path) -> Path:
    """以中性毫米STEP重导模型，不复制原始STEP头。"""
    ocp = _ocp()
    model = _model(model, ocp)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    settings = ocp["Interface_Static"]
    previous = settings.CVal_s("write.step.unit")
    settings.SetCVal_s("write.step.unit", "MM")
    try:
        writer = ocp["STEPControl_Writer"]()
        if writer.Transfer(model.shape, ocp["STEPControl_AsIs"]) != ocp["IFSelect_RetDone"]:
            raise OSError("STEP形状转换失败")
        if writer.Write(str(target)) != ocp["IFSelect_RetDone"]:
            raise OSError("STEP模型写出失败")
    finally:
        if previous:
            settings.SetCVal_s("write.step.unit", previous)
    return target


def export_model_brep(model: ImportedCadModel, path) -> Path:
    """BREP无内置单位；输出坐标为mm，读取时须保留单位约定。"""
    ocp = _ocp()
    model = _model(model, ocp)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not ocp["BRepTools"].Write_s(model.shape, str(target)):
        raise OSError("BREP模型写出失败")
    return target
