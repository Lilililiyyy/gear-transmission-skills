"""可选VTK专业几何预览与VTP导出；不生成压力场或补造齿形。

导入模块仅用标准库，调用时才导入VTK。渲染始终离屏，不创建交互器、
不打开可见窗口。输入坐标为m，PNG文字中的尺寸换算为mm。
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path


def _vtk():
    try:
        import vtk
    except ImportError as error:
        raise ImportError("专业几何预览需要可选VTK依赖；核心网格和OCP功能可独立使用") from error
    return vtk


def _metadata(mesh):
    return {"length_unit": "m", "closed_solid": bool(mesh.get("metadata", {}).get("closed_solid", False)),
            "geometry_metadata": dict(mesh.get("metadata", {})),
            "preview_scope": "真实输入三角网格的专业几何预览，不生成接触压力或装饰齿形",
            "topology_note": "壳/实体标志沿用输入元数据，不由预览器补面或推定制造有效性"}


def _polydata(mesh):
    if not isinstance(mesh, dict) or mesh.get("units") != "m":
        raise ValueError("预览输入须为单位m的CAD三角网格")
    vertices, triangles = mesh.get("vertices_m", ()), mesh.get("triangles", ())
    if not vertices or not triangles:
        raise ValueError("预览网格须包含顶点及三角形")
    vtk = _vtk()
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    for point in vertices:
        if len(point) != 3 or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in point):
            raise ValueError("网格顶点须为三个有限SI坐标")
        points.InsertNextPoint(*point)
    cells = vtk.vtkCellArray()
    for triangle in triangles:
        if len(triangle) != 3 or any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(vertices) for i in triangle):
            raise ValueError("三角形须使用有效的0起始顶点索引")
        cells.InsertNextCell(3)
        for index in triangle:
            cells.InsertCellPoint(index)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)
    face_indices = mesh.get("face_indices")
    if face_indices is not None:
        if len(face_indices) != len(triangles) or any(isinstance(i, bool) or not isinstance(i, int) or i < 0 for i in face_indices):
            raise ValueError("face_indices须逐三角形给出非负整数面序号")
        values = vtk.vtkIntArray()
        values.SetName("face_indices")
        for index in face_indices:
            values.InsertNextValue(index)
        polydata.GetCellData().AddArray(values)
    metadata = vtk.vtkStringArray()
    metadata.SetName("geometry_metadata_json")
    metadata.InsertNextValue(json.dumps(_metadata(mesh), ensure_ascii=True))
    polydata.GetFieldData().AddArray(metadata)
    return vtk, polydata


def export_vtp(mesh: dict, path) -> Path:
    """保存真实三角拓扑、面序号、SI单位及开放壳说明，不依赖渲染窗口。"""
    vtk, polydata = _polydata(mesh)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(target))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    if writer.Write() != 1 or writer.GetErrorCode() != 0:
        raise OSError("VTP网格文件写出失败")
    return target


def _text_font(properties, vtk):
    # VTK内置西文字体不保证中文字形；存在系统中文字体时直接使用。
    fonts = Path(os.environ.get("WINDIR") or os.environ.get("SystemRoot") or ".")/"Fonts"
    for name in ("msyh.ttc", "simhei.ttf", "simsun.ttc"):
        candidate = fonts/name
        if candidate.is_file():
            properties.SetFontFamily(vtk.VTK_FONT_FILE)
            properties.SetFontFile(str(candidate))
            return
    properties.SetFontFamilyToArial()


def _fit_camera(camera, vertices, bounds, size):
    """按实际投影范围适配平行相机，并为顶部说明保留位置。"""
    center = [(bounds[2*i]+bounds[2*i+1])/2 for i in range(3)]
    extent = [bounds[2*i+1]-bounds[2*i] for i in range(3)]
    span = max(extent)
    if span <= 0:
        raise ValueError("网格包围尺寸为零，不能设置预览相机")
    direction = (1.3, -1.5, 1.1)
    length = math.sqrt(sum(value*value for value in direction))
    forward = tuple(-value/length for value in direction)
    right = (forward[1], -forward[0], 0.)
    right_length = math.hypot(right[0], right[1])
    right = tuple(value/right_length for value in right)
    up = (right[1]*forward[2]-right[2]*forward[1],
          right[2]*forward[0]-right[0]*forward[2],
          right[0]*forward[1]-right[1]*forward[0])
    projected = [(sum((point[i]-center[i])*right[i] for i in range(3)),
                  sum((point[i]-center[i])*up[i] for i in range(3))) for point in vertices]
    low = [min(point[i] for point in projected) for i in (0, 1)]
    high = [max(point[i] for point in projected) for i in (0, 1)]
    top = min(size[1]*.25, max(72., size[1]*.105))
    bottom = size[1]*.06
    available_width, available_height = size[0]*.84, size[1]-top-bottom
    scale = 1.08*.5*max((high[0]-low[0])*size[1]/available_width,
                       (high[1]-low[1])*size[1]/available_height)
    target_y = (size[1]-top+bottom)/2
    dx = (low[0]+high[0])/2
    dy = (low[1]+high[1])/2-2*scale*(target_y-size[1]/2)/size[1]
    focal = [center[i]+right[i]*dx+up[i]*dy for i in range(3)]
    camera.SetFocalPoint(*focal)
    camera.SetPosition(*(focal[i]+span*direction[i] for i in range(3)))
    camera.SetViewUp(*up)
    camera.ParallelProjectionOn()
    camera.SetParallelScale(scale)
    return extent


def render_cad_mesh(mesh: dict, path, *, size=(1400, 1100), show_edges=True,
                    lang="zh", annotate=True, cell_colors_rgb=None) -> Path:
    """离屏渲染白底灰蓝几何，边线表示边界/特征，不代表压力或接触斑。"""
    if len(size) != 2 or any(isinstance(v, bool) or not isinstance(v, int) or not 64 <= v <= 8192 for v in size):
        raise ValueError("图像宽高须为64至8192之间的整数")
    if not isinstance(show_edges, bool):
        raise ValueError("show_edges须为布尔值")
    if lang not in ("zh", "en"):
        raise ValueError("lang须为zh或en")
    if not isinstance(annotate, bool):
        raise ValueError("annotate须为布尔值；论文组合图可关闭内嵌文字")
    vtk, polydata = _polydata(mesh)
    if cell_colors_rgb is not None:
        if len(cell_colors_rgb) != len(mesh["triangles"]) or any(
            len(color) != 3 or any(isinstance(channel, bool) or not isinstance(channel, int) or not 0 <= channel <= 255
                                   for channel in color) for color in cell_colors_rgb):
            raise ValueError("cell_colors_rgb须逐三角片给出3个0至255整数颜色通道")
        colors = vtk.vtkUnsignedCharArray()
        colors.SetName("geometry_region_rgb")
        colors.SetNumberOfComponents(3)
        for color in cell_colors_rgb:
            colors.InsertNextTuple3(*color)
        polydata.GetCellData().AddArray(colors)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOff()
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.SplittingOn()
    normals.SetFeatureAngle(60.)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.ScalarVisibilityOff()
    if cell_colors_rgb is not None:
        mapper.SetScalarModeToUseCellFieldData()
        mapper.SelectColorArray("geometry_region_rgb")
        mapper.SetColorModeToDirectScalars()
        mapper.ScalarVisibilityOn()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    material = actor.GetProperty()
    material.SetColor(.58, .70, .77)
    material.SetAmbient(.30)
    material.SetDiffuse(.65)
    material.SetSpecular(.20)
    material.SetSpecularPower(28.)
    material.SetInterpolationToPhong()
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1., 1., 1.)
    renderer.TwoSidedLightingOn()
    renderer.AddActor(actor)
    if show_edges:
        edges = vtk.vtkFeatureEdges()
        edges.SetInputData(polydata)
        edges.BoundaryEdgesOn()
        edges.FeatureEdgesOn()
        edges.NonManifoldEdgesOn()
        edges.ManifoldEdgesOff()
        edges.SetFeatureAngle(30.)
        edge_mapper = vtk.vtkPolyDataMapper()
        edge_mapper.SetInputConnection(edges.GetOutputPort())
        edge_mapper.ScalarVisibilityOff()
        edge_actor = vtk.vtkActor()
        edge_actor.SetMapper(edge_mapper)
        edge_actor.GetProperty().SetColor(.16, .23, .27)
        edge_actor.GetProperty().SetLineWidth(1.3)
        edge_actor.GetProperty().LightingOff()
        renderer.AddActor(edge_actor)
    bounds = polydata.GetBounds()
    camera = renderer.GetActiveCamera()
    extent = _fit_camera(camera, mesh["vertices_m"], bounds, size)
    renderer.ResetCameraClippingRange()
    lights = vtk.vtkLightKit()
    lights.SetKeyLightIntensity(.9)
    lights.SetKeyToFillRatio(3.)
    lights.SetKeyToBackRatio(4.)
    lights.AddLightsToRenderer(renderer)
    metadata = _metadata(mesh)
    metadata["geometry_region_colors"] = cell_colors_rgb is not None
    caption = vtk.vtkTextActor()
    state = (("实体三角网格" if metadata["closed_solid"] else "开放曲面壳") if lang == "zh"
             else ("Solid triangular mesh" if metadata["closed_solid"] else "Open surface shell"))
    dimensions = " × ".join(f"{value*1000:.3f}" for value in extent)
    caption.SetInput((f"CAD几何预览\n{state}  ·  坐标单位 m\n包围尺寸 {dimensions} mm") if lang == "zh"
                     else f"CAD geometry\n{state}  ·  Coordinates in m\nBounding dimensions {dimensions} mm")
    caption.SetDisplayPosition(22, size[1]-18)
    properties = caption.GetTextProperty()
    properties.SetFontSize(max(12, round(size[1]*.018)))
    properties.SetColor(.16, .20, .23)
    properties.SetVerticalJustificationToTop()
    _text_font(properties, vtk)
    if annotate:
        renderer.AddActor2D(caption)
    window = vtk.vtkRenderWindow()
    window.OffScreenRenderingOn()
    if hasattr(window, "SetShowWindow"):
        window.SetShowWindow(False)
    window.SetSize(*size)
    window.SetMultiSamples(8)
    window.AddRenderer(renderer)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        if window.GetOffScreenRendering() != 1:
            raise RuntimeError("VTK未进入离屏模式，停止渲染")
        window.Render()
        show_window = bool(window.GetShowWindow()) if hasattr(window, "GetShowWindow") else None
        if window.GetOffScreenRendering() != 1 or show_window is True:
            raise RuntimeError("VTK未保持隐藏的离屏模式，停止导出")
        capabilities = window.ReportCapabilities()
        details = [line for line in capabilities.splitlines() if any(key in line.lower() for key in ("vendor string", "renderer string", "version string"))]
        metadata.update(language=lang, annotated=annotate, vtk_version=vtk.vtkVersion.GetVTKVersion(), renderer_class=renderer.GetClassName(),
                        render_window_class=window.GetClassName(), offscreen=bool(window.GetOffScreenRendering()),
                        show_window=show_window, renderer_details=details,
                        size=list(size), show_edges=show_edges, preview_material="灰蓝统一材质，无压力标量着色")
        capture = vtk.vtkWindowToImageFilter()
        capture.SetInput(window)
        capture.SetInputBufferTypeToRGB()
        capture.ReadFrontBufferOff()
        capture.Update()
        writer = vtk.vtkPNGWriter()
        writer.SetFileName(str(target))
        writer.SetInputConnection(capture.GetOutputPort())
        writer.AddText("cad_preview_metadata", json.dumps(metadata, ensure_ascii=True))
        writer.Write()
        if writer.GetErrorCode() != 0 or not target.is_file() or target.stat().st_size == 0:
            raise OSError("VTK离屏PNG写出失败")
    finally:
        window.Finalize()
    return target
