"""把真实接触多边形压力贴到实际CAD网格，使用Python/VTK离屏渲染。

压力按单元直接着色，不投影求新压力、不平滑、不补造接触椭圆。物理
坐标保持m；深度冲突只用VTK屏幕深度偏置，绝不挪动原压力坐标。
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tempfile


def _numeric():
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError('三维压力显示需要可选NumPy、Matplotlib与VTK') from error
    return np


def _rotation(value,np):
    value=np.asarray(value,dtype=float)
    if value.shape!=(3,3) or not np.isfinite(value).all() or not np.allclose(value.T@value,np.eye(3),rtol=0,atol=1e-9) or abs(np.linalg.det(value)-1)>1e-9:
        raise ValueError('压力世界坐标到网格坐标的姿态须为右手正交3×3旋转')
    return value


def prepare_pressure_overlay(polygons_m,pressures_pa,*,polygon_coordinates='mesh_local',
        mesh_rotation_rows=None,mesh_translation_m=(0.,0.,0.),forces_n=None,integration_areas_m2=None,
        application_points_m=None,
        area_definition='调用方显式提供的原接触积分面积'):
    """保存真实三维多边形及压力，可将world点按mesh姿态逆变换到局部。

    mesh_rotation_rows将网格局部点转到世界；本函数使用R.T@(x−t)。
    多边形可为三角、四角或更多顶点，不能以积分点替代有界区域。
    接触面积不会从显示多边形重新推算，半空间参考面积须由调用方给定。
    application_points_m为原求解力作用点，与polygon_coordinates使用相同坐标系；
    显示多边形质心只用于取景，不能替代该作用点作扭矩或平衡计算。
    """
    np=_numeric()
    if polygon_coordinates not in ('mesh_local','world'):raise ValueError('polygon_coordinates须为mesh_local或world')
    pressure=np.asarray(pressures_pa,dtype=float)
    if pressure.ndim!=1 or not len(pressure) or not np.isfinite(pressure).all() or np.any(pressure<0):
        raise ValueError('多边形压力须为逐单元有限非负Pa数组')
    if len(polygons_m)!=len(pressure):raise ValueError('每个真实多边形必须对应一个压力值')
    translation=np.asarray(mesh_translation_m,dtype=float)
    if translation.shape!=(3,) or not np.isfinite(translation).all():raise ValueError('mesh_translation_m须为有限三维m坐标')
    if polygon_coordinates=='world' and mesh_rotation_rows is None:raise ValueError('世界压力多边形须显式给整轮网格姿态，不能猜测坐标系')
    rotation=_rotation(np.eye(3) if mesh_rotation_rows is None else mesh_rotation_rows,np)
    polygons=[];normals=[]
    for value in polygons_m:
        points=np.asarray(value,dtype=float)
        if points.ndim!=2 or points.shape[1]!=3 or len(points)<3 or not np.isfinite(points).all():
            raise ValueError('压力区域须提供至少三个有限真实三维顶点，单位m')
        if polygon_coordinates=='world':points=(points-translation)@rotation
        normal=np.cross(points[1:]-points[0],np.roll(points,-1,axis=0)[1:]-points[0]).sum(axis=0)
        length=np.linalg.norm(normal)
        if length<=1e-24:raise ValueError('压力多边形退化，不能把点或线涂成接触区')
        polygons.append(points.tolist());normals.append((normal/length).tolist())
    if (forces_n is None)!=(integration_areas_m2 is None):raise ValueError('核对pA=f时须同时给原单元力与原积分面积')
    force=area=None
    if forces_n is not None:
        force=np.asarray(forces_n,dtype=float);area=np.asarray(integration_areas_m2,dtype=float)
        if force.shape!=pressure.shape or area.shape!=pressure.shape or not np.isfinite(force).all() or not np.isfinite(area).all() or np.any(force<0) or np.any(area<=0):
            raise ValueError('原接触力/积分面积须为逐单元有限非负力和正面积')
        if not np.allclose(pressure*area,force,rtol=2e-10,atol=1e-9):raise ValueError('原压力乘原积分面积与单元力不一致，停止三维显示')
    applications=None
    if application_points_m is not None:
        applications=np.asarray(application_points_m,dtype=float)
        if applications.shape!=(len(pressure),3) or not np.isfinite(applications).all():raise ValueError('原力作用点须为逐单元有限三维m数组')
        if polygon_coordinates=='world':applications=(applications-translation)@rotation
    active=pressure>0
    if not active.any():raise ValueError('没有正压力单元，未生成接触斑点占位图')
    weights=pressure[active] if force is None else force[active]
    centers=np.asarray([np.mean(p,axis=0) for p in polygons])[active]
    direction=np.asarray(normals)[active].copy();reference=direction[0]
    direction[np.einsum('ij,j->i',direction,reference)<0]*=-1
    average=np.average(direction,axis=0,weights=weights);average/=np.linalg.norm(average)
    identity=hashlib.sha256(json.dumps(dict(polygons_m=polygons,pressure_pa=pressure.tolist(),
        forces_n=None if force is None else force.tolist(),integration_areas_m2=None if area is None else area.tolist(),
        application_points_local_m=None if applications is None else applications.tolist()),
        sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')).hexdigest()
    return dict(schema='actual_pressure_overlay_v1',units='m',polygons_local_m=polygons,
        pressure_pa=pressure.tolist(),forces_n=None if force is None else force.tolist(),
        integration_areas_m2=None if area is None else area.tolist(),area_definition=area_definition,
        application_points_local_m=None if applications is None else applications.tolist(),
        active_mask=active.tolist(),mean_normal_local=average.tolist(),
        contact_center_local_m=np.average(centers,axis=0,weights=weights).tolist(),
        original_polygon_coordinates=polygon_coordinates,mesh_rotation_rows=rotation.tolist(),
        mesh_translation_m=translation.tolist(),source_identity_sha256=identity,
        pressure_values_modified=False,coordinates_displaced=False,
        pressure_interpolation='none',force_area_consistency_checked=force is not None)


def halfspace_pressure_overlay(result,*,tooth_offset=0,mesh_rotation_rows=None,mesh_translation_m=(0.,0.,0.)):
    """提取实际曲面四角，不把邻齿折叠到基齿；非零齿偏移须显式给整轮姿态。"""
    np=_numeric()
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('物理齿偏移tooth_offset须为整数')
    if not result.get('converged') or not result.get('fields'):raise ValueError('未接受的半解析候选不能贴到三维齿轮')
    offsets=result['selected_tooth_offsets']
    if tooth_offset not in offsets:raise ValueError('所选物理齿没有已接受压力场')
    patch=next(p for p in result['geometry']['patches'] if p['tooth_offset']==tooth_offset)
    if mesh_rotation_rows is None:
        if tooth_offset!=0:raise ValueError('非零齿偏移必须给整轮mesh_rotation_rows，patch姿态会把它折回基齿')
        mesh_rotation_rows=patch['frames']['pinion_rotation_rows']
    field=result['fields'][offsets.index(tooth_offset)]
    active=np.asarray(field['contact_mask'],dtype=bool)
    corners=np.asarray(patch['first_corner_points_world_m'],dtype=float)
    polygons=[corners[[j,j,j+1,j+1],[i,i+1,i+1,i]].tolist() for j,i in np.argwhere(active)]
    result_overlay=prepare_pressure_overlay(polygons,np.asarray(field['pressure_pa'])[active],polygon_coordinates='world',
        mesh_rotation_rows=mesh_rotation_rows,mesh_translation_m=mesh_translation_m,
        forces_n=np.asarray(field['forces_n'])[active],integration_areas_m2=np.full(int(active.sum()),field['cell_area_m2']),
        area_definition='半解析原参考接触切平面矩形面积，未以三维显示面积替代')
    if 'first_outward_normals_world' in patch:
        outward=np.asarray(patch['first_outward_normals_world'])[active]@np.asarray(mesh_rotation_rows)
        normal=np.average(outward,axis=0,weights=np.asarray(field['forces_n'])[active]);normal/=np.linalg.norm(normal)
        result_overlay['mean_normal_local']=normal.tolist()
    result_overlay.update(method='local_halfspace',tooth_offset=tooth_offset,
        driving_angle_rad=result['geometry']['driving_angle_rad'],
        pinion_torque_nm=result['geometry'].get('pinion_torque_nm'),
        model_scope='局部半空间压力，不含有限齿体/轮体/跨齿结构柔度')
    boundaries={entry['tooth_offset']:entry for entry in result.get('contact_boundary_checks',[])}
    result_overlay['all_tooth_context']=[dict(tooth_offset=offset,
        loaded_cell_count=int(np.count_nonzero(value['contact_mask'])),
        reference_plane_normal_force_n=float(np.sum(value['forces_n'])),
        boundary_checks=boundaries.get(offset,{})) for offset,value in zip(offsets,result['fields'])]
    return result_overlay


def finite_element_pressure_overlay(result,*,tooth_offset=0):
    """提取FE实际有界接触区域，支持任意顶点数；不将积分点扩成原粗面片。"""
    np=_numeric()
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('物理齿偏移tooth_offset须为整数')
    if not result.get('converged') or not result.get('contacts'):raise ValueError('未接受的FE结果不能贴到三维齿轮')
    contacts=[p for p in result['contacts'] if p.get('active') and p['force_n']>0 and p['tooth_offset']==tooth_offset]
    if not contacts:raise ValueError('所选FE物理齿没有承载多边形')
    if any('first_vertices_world_m' not in p for p in contacts):
        raise ValueError('FE压力须给实际积分区域的三维边界，不能把积分点涂满旧粗三角形')
    frames=result['mesh_reference_frames']
    overlay=prepare_pressure_overlay([p['first_vertices_world_m'] for p in contacts],[p['mean_pressure_pa'] for p in contacts],
        polygon_coordinates='world',mesh_rotation_rows=frames['pinion_rotation_rows'],
        mesh_translation_m=frames.get('translation_m',(0.,0.,0.)),
        forces_n=[p['force_n'] for p in contacts],integration_areas_m2=[p['area_m2'] for p in contacts],
        application_points_m=[p['application_point_world_m'] for p in contacts] if all('application_point_world_m' in p for p in contacts) else None,
        area_definition='FE结果显式保存的原接触积分区域面积')
    if all('normal_world' in p for p in contacts):
        normals=np.asarray([p['normal_world'] for p in contacts])@np.asarray(frames['pinion_rotation_rows'])
        normal=np.average(normals,axis=0,weights=[p['force_n'] for p in contacts]);normal/=np.linalg.norm(normal)
        overlay['mean_normal_local']=normal.tolist()
    overlay.update(method='finite_element',tooth_offset=tooth_offset,
        driving_angle_rad=result['driving_angle_rad'],pinion_torque_nm=result.get('pinion_torque_nm'),model_scope=result.get('model','三维线弹性FE接触'),
        contact_representation=result.get('contact_representation','explicit_contact_polygons'))
    if all('first_pressure_node_index' in p for p in contacts):
        indices=[p['first_pressure_node_index'] for p in contacts]
        overlay.update(pressure_node_indices=indices,active_pressure_node_count=len(set(indices)))
    all_contacts=[p for p in result['contacts'] if p.get('active') and p['force_n']>0]
    overlay['all_tooth_context']=[dict(tooth_offset=offset,loaded_dual_cell_count=len(part),
        active_pressure_node_count=len({p['first_pressure_node_index'] for p in part}) if all('first_pressure_node_index' in p for p in part) else None,
        normal_force_n=sum(p['force_n'] for p in part)) for offset in sorted({p['tooth_offset'] for p in all_contacts})
        for part in [[p for p in all_contacts if p['tooth_offset']==offset]]]
    return overlay


def _pressure_polydata(overlay,pressure_limit_mpa=None):
    from .cad_preview import _vtk
    from .halfspace_plotting import pressure_colour_scale
    np=_numeric();vtk=_vtk()
    if overlay.get('schema')!='actual_pressure_overlay_v1' or overlay.get('units')!='m':raise ValueError('需要显式单位m的实际压力多边形数据')
    pressure=np.asarray(overlay['pressure_pa']);active=pressure>0
    norm,cmap,scale=pressure_colour_scale(pressure,pressure_limit_mpa)
    points=vtk.vtkPoints();points.SetDataTypeToDouble();cells=vtk.vtkCellArray()
    p_values=vtk.vtkDoubleArray();p_values.SetName('pressure_pa')
    colour=vtk.vtkUnsignedCharArray();colour.SetName('pressure_rgb');colour.SetNumberOfComponents(3)
    area_values=vtk.vtkDoubleArray();area_values.SetName('original_integration_area_m2')
    force_values=vtk.vtkDoubleArray();force_values.SetName('force_n')
    applications=vtk.vtkDoubleArray();applications.SetName('application_point_local_m');applications.SetNumberOfComponents(3)
    pressure_nodes=vtk.vtkIdTypeArray();pressure_nodes.SetName('original_pressure_node_index')
    for index in np.flatnonzero(active):
        polygon=overlay['polygons_local_m'][index];cells.InsertNextCell(len(polygon))
        for point in polygon:cells.InsertCellPoint(points.InsertNextPoint(*point))
        p_values.InsertNextValue(float(pressure[index]))
        rgb=np.rint(np.asarray(cmap(norm(pressure[index]*1e-6)))[:3]*255).astype(int)
        colour.InsertNextTuple3(*map(int,rgb))
        if overlay['forces_n'] is not None:
            area_values.InsertNextValue(overlay['integration_areas_m2'][index]);force_values.InsertNextValue(overlay['forces_n'][index])
        if overlay.get('application_points_local_m') is not None:applications.InsertNextTuple3(*overlay['application_points_local_m'][index])
        if overlay.get('pressure_node_indices') is not None:pressure_nodes.InsertNextValue(overlay['pressure_node_indices'][index])
    poly=vtk.vtkPolyData();poly.SetPoints(points);poly.SetPolys(cells)
    poly.GetCellData().SetScalars(colour);poly.GetCellData().AddArray(p_values)
    if overlay['forces_n'] is not None:poly.GetCellData().AddArray(area_values);poly.GetCellData().AddArray(force_values)
    if overlay.get('application_points_local_m') is not None:poly.GetCellData().AddArray(applications)
    if overlay.get('pressure_node_indices') is not None:poly.GetCellData().AddArray(pressure_nodes)
    metadata=dict(**scale,units='m/Pa/N/m²',source_identity_sha256=overlay['source_identity_sha256'],
        source_polygon_count=len(pressure),rendered_positive_pressure_polygons=int(active.sum()),
        pressure_interpolation='none',colour_quantization='unsigned_byte_rgb',area_definition=overlay['area_definition'])
    field=vtk.vtkStringArray();field.SetName('pressure_metadata_json');field.InsertNextValue(json.dumps(metadata,ensure_ascii=True))
    poly.GetFieldData().AddArray(field)
    return vtk,poly,metadata


def export_pressure_overlay_vtp(overlay,path,*,pressure_limit_mpa=None):
    """保存真实多边形、原压力/力/积分面积和共用色标RGB，供专业后处理重读。"""
    vtk,poly,metadata=_pressure_polydata(overlay,pressure_limit_mpa)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(str(path));writer.SetInputData(poly);writer.SetDataModeToBinary()
    if writer.Write()!=1:raise OSError('三维压力VTP写出失败')
    return path


def render_pressure_overlay(mesh,overlay,path,*,pressure_limit_mpa=None,view='gear',
        view_direction_local=None,size=(2200,1700),body_opacity=1.,depth_bias_units=None):
    """在原CAD上渲染原压力层；gear为整轮，contact为同一几何的近景。"""
    from .cad_preview import _polydata
    np=_numeric()
    if view not in ('gear','contact'):raise ValueError('三维视图须为gear或contact')
    if not isinstance(body_opacity,(int,float)) or isinstance(body_opacity,bool) or not math.isfinite(body_opacity) or not 0<=body_opacity<=1:raise ValueError('实体显示透明度须为0–1有限值')
    if depth_bias_units is not None and (not isinstance(depth_bias_units,(int,float)) or isinstance(depth_bias_units,bool) or not math.isfinite(depth_bias_units)):raise ValueError('屏幕深度偏置须为有限数值或None自动值')
    if len(size)!=2 or any(isinstance(v,bool) or not isinstance(v,int) or not 64<=v<=8192 for v in size):raise ValueError('三维渲染尺寸须为64–8192整数像素')
    vtk,body=_polydata(mesh);_,pressure,metadata=_pressure_polydata(overlay,pressure_limit_mpa)
    # 仅为灰色CAD实体生成光照法线，不改顶点或三角面，也不插值压力颜色。
    body_normals=vtk.vtkPolyDataNormals();body_normals.SetInputData(body)
    body_normals.ConsistencyOn();body_normals.AutoOrientNormalsOff();body_normals.SplittingOn();body_normals.SetFeatureAngle(50.)
    body_normals.ComputePointNormalsOn();body_normals.ComputeCellNormalsOn()
    base_mapper=vtk.vtkPolyDataMapper();base_mapper.SetInputConnection(body_normals.GetOutputPort());base_mapper.ScalarVisibilityOff()
    base=vtk.vtkActor();base.SetMapper(base_mapper);base.GetProperty().SetColor(.77,.79,.81)
    base.GetProperty().SetAmbient(.27);base.GetProperty().SetDiffuse(.70);base.GetProperty().SetSpecular(.12)
    base.GetProperty().SetOpacity(body_opacity)
    base.GetProperty().SetInterpolationToPhong()
    pressure_mapper=vtk.vtkPolyDataMapper();pressure_mapper.SetInputData(pressure)
    pressure_mapper.SetScalarModeToUseCellData();pressure_mapper.SetColorModeToDirectScalars();pressure_mapper.ScalarVisibilityOn()
    pressure_mapper.InterpolateScalarsBeforeMappingOff()
    pressure_mapper.SetResolveCoincidentTopologyToPolygonOffset()
    pressure_mapper.SetRelativeCoincidentTopologyPolygonOffsetParameters(-2.,-2. if depth_bias_units is None else depth_bias_units)
    actor=vtk.vtkActor();actor.SetMapper(pressure_mapper);actor.GetProperty().LightingOff();actor.GetProperty().SetInterpolationToFlat()
    renderer=vtk.vtkRenderer();renderer.SetBackground(1.,1.,1.);renderer.AddActor(base);renderer.AddActor(actor)
    active=np.asarray(overlay['active_mask'],dtype=bool)
    pressure_points=np.concatenate([np.asarray(p) for p,keep in zip(overlay['polygons_local_m'],active) if keep])
    direction=np.array(overlay['mean_normal_local'] if view_direction_local is None else view_direction_local,dtype=float,copy=True)
    if direction.shape!=(3,) or not np.isfinite(direction).all() or np.linalg.norm(direction)<1e-12:raise ValueError('视线须为有限非零局部三维方向')
    direction/=np.linalg.norm(direction)
    up=np.array([0.,0.,1.])
    if abs(up@direction)>.92:up=np.array([0.,1.,0.])
    up-=direction*(up@direction);up/=np.linalg.norm(up);right=np.cross(up,direction)
    if view=='contact':
        projected=np.column_stack((pressure_points@right,pressure_points@up))
        _,vectors=np.linalg.eigh(np.cov(projected,rowvar=False));long=vectors[:,-1]
        if long[0]<0:long=-long
        right=right*long[0]+up*long[1];right/=np.linalg.norm(right);up=np.cross(direction,right)
    selected=np.asarray(mesh['vertices_m']) if view=='gear' else pressure_points
    center=(selected.min(axis=0)+selected.max(axis=0))/2
    coordinates=np.column_stack(((selected-center)@right,(selected-center)@up))
    limits=np.array([coordinates.min(axis=0),coordinates.max(axis=0)])
    center+=right*limits[:,0].mean()+up*limits[:,1].mean()
    extent=limits[1]-limits[0]
    scale=.57*max(float(extent[1]),float(extent[0])*size[1]/size[0],1e-6)
    distance=max(float(np.ptp(np.asarray(mesh['vertices_m']),axis=0).max()),scale)*4
    camera=renderer.GetActiveCamera();camera.SetFocalPoint(*center);camera.SetPosition(*(center+distance*direction));camera.SetViewUp(*up)
    camera.ParallelProjectionOn();camera.SetParallelScale(scale);renderer.ResetCameraClippingRange()
    lights=vtk.vtkLightKit();lights.SetKeyLightIntensity(.85);lights.SetKeyToFillRatio(3.);lights.AddLightsToRenderer(renderer)
    window=vtk.vtkRenderWindow();window.OffScreenRenderingOn()
    if hasattr(window,'SetShowWindow'):window.SetShowWindow(False)
    window.SetSize(*size);window.SetMultiSamples(0);window.AddRenderer(renderer)
    target=Path(path);target.parent.mkdir(parents=True,exist_ok=True)
    metadata.update(view=view,size=list(size),vtk_version=vtk.vtkVersion.GetVTKVersion(),
        geometry='调用方实际CAD三角网格，未造齿形',zero_pressure_geometry='grey',
        pressure_lighting=False,coordinate_displacement_m=0.,depth_bias='仅VTK屏幕深度偏置',
        body_opacity=body_opacity,
        camera=dict(position_local_m=list(camera.GetPosition()),focal_point_local_m=list(camera.GetFocalPoint()),
            view_up=list(camera.GetViewUp()),parallel_scale_m=camera.GetParallelScale()),
        source_mesh_vertices=len(mesh['vertices_m']),source_mesh_triangles=len(mesh['triangles']),
        contact_view_orientation='仅相机绕视线旋转，使实际压力区主轴水平' if view=='contact' else '网格局部z轴参考')
    try:
        window.Render()
        # 原CAD三角化可有明确弦差，真实NURBS压力多边形可能被弦面遮住。
        # 只按该精度在平行相机Z缓冲作显示偏置，绝不修改世界坐标或压力。
        clipping=camera.GetClippingRange();bits=max(16,int(window.GetDepthBufferSize()))
        quantum=(clipping[1]-clipping[0])/(2**bits-1)
        h=mesh.get('metadata',{}).get('linear_deflection_m',0.)
        allowance=0.
        if isinstance(h,(int,float)) and math.isfinite(h) and h>0:
            allowance=min(4*h,float(np.ptp(np.asarray(mesh['vertices_m']),axis=0).max())*.005)
        automatic=depth_bias_units is None
        effective_bias=-max(2.,math.ceil(allowance/max(quantum,1e-30))+2.) if automatic else depth_bias_units
        pressure_mapper.SetRelativeCoincidentTopologyPolygonOffsetParameters(-2.,effective_bias)
        metadata.update(depth_bias_units=effective_bias,depth_bias_automatic=automatic,
            depth_bias_equivalent_allowance_m=allowance if automatic else abs(effective_bias)*quantum,
            depth_buffer_bits=bits,camera_clipping_range_m=list(clipping),
            depth_bias_scope='补偿输入CAD网格弦差的屏幕显示；不用于读真实体间隙或变形')
        window.Render()
        if window.GetOffScreenRendering()!=1 or hasattr(window,'GetShowWindow') and window.GetShowWindow():
            raise RuntimeError('三维压力渲染未保持隐藏离屏模式')
        capture=vtk.vtkWindowToImageFilter();capture.SetInput(window);capture.SetInputBufferTypeToRGB();capture.ReadFrontBufferOff();capture.Update()
        writer=vtk.vtkPNGWriter();writer.SetFileName(str(target));writer.SetInputConnection(capture.GetOutputPort())
        writer.AddText('pressure_render_metadata',json.dumps(metadata,ensure_ascii=True));writer.Write()
        if writer.GetErrorCode()!=0 or not target.exists():raise OSError('三维压力PNG写出失败')
    finally:window.Finalize()
    return metadata


def plot_pressure_on_gear(mesh,overlay,output,*,lang='zh',pressure_limit_mpa=None,
        view_direction_local=None,stem='13_pressure_on_gear'):
    """183mm双栏实际齿轮总览/接触区近景，PNG/SVG/PDF及原压力三维VTP。"""
    from .workflow_plotting import _paper,_save,_text
    from .halfspace_plotting import pressure_colour_scale
    np=_numeric();norm,cmap,scale=pressure_colour_scale(overlay['pressure_pa'],pressure_limit_mpa)
    output=Path(output);rows=[]
    vtp=export_pressure_overlay_vtp(overlay,output/'data'/(stem+'.vtp'),pressure_limit_mpa=pressure_limit_mpa)
    for i,(polygon,p) in enumerate(zip(overlay['polygons_local_m'],overlay['pressure_pa'])):
        rows.append(dict(row_kind='actual_pressure_polygon',cell_index=i,pressure_pa=p,
            original_pressure_node_index=None if overlay.get('pressure_node_indices') is None else overlay['pressure_node_indices'][i],
            polygon_mesh_local_m_json=json.dumps(polygon),
            application_point_mesh_local_m_json=None if overlay.get('application_points_local_m') is None else json.dumps(overlay['application_points_local_m'][i]),
            force_n=None if overlay['forces_n'] is None else overlay['forces_n'][i],
            original_integration_area_m2=None if overlay['integration_areas_m2'] is None else overlay['integration_areas_m2'][i]))
    with tempfile.TemporaryDirectory(prefix='gear_pressure_render_') as directory:
        panels=[]
        for view in ('gear','contact'):
            path=Path(directory)/(view+'.png')
            record=render_pressure_overlay(mesh,overlay,path,pressure_limit_mpa=pressure_limit_mpa,
                view=view,view_direction_local=view_direction_local)
            panels.append((path,record))
        with _paper(lang) as plt:
            from matplotlib.cm import ScalarMappable
            fig,axes=plt.subplots(1,2,layout='constrained')
            method=overlay.get('method')
            nodal=overlay.get('contact_representation')=='nodal_dual_cells'
            name=(_text(lang,'整轮FE节点集总接触压力','Full-gear FE vertex-lumped nodal pressure') if nodal else
                _text(lang,'整轮FE接触压力','Full-gear FE contact pressure') if method=='finite_element' else
                _text(lang,'局部半解析接触压力','Local semi-analytical contact pressure') if method=='local_halfspace' else
                _text(lang,'实际齿面接触压力','Actual flank contact pressure'))
            title=name+'\n'+rf'$p_{{\max}}={scale["raw_maximum_pressure_mpa"]:.0f}$ MPa'
            if overlay.get('pinion_torque_nm') is not None:
                title+='; '+rf'$T_1={overlay["pinion_torque_nm"]:g}$ N·m, k={overlay.get("tooth_offset")}'
            fig.suptitle(title,fontsize=9)
            for ax,(path,record),letter in zip(axes,panels,'ab'):
                ax.imshow(plt.imread(path));ax.axis('off')
                ax.set_title(_text(lang,'实际齿轮与压力层','Actual gear with pressure') if record['view']=='gear' else
                    _text(lang,'真实齿面接触区近景','Actual flank contact detail'),loc='left',fontsize=8)
                ax.text(-.02,1.08,letter,transform=ax.transAxes,fontsize=9,weight='bold')
            bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=list(axes),fraction=.035,pad=.02,
                extend='max' if scale['over_limit_cell_count'] else 'neither')
            bar.set_label(_text(lang,'接触压力 p (MPa)','Contact pressure, p (MPa)'))
            caption=_text(lang,
                '真实CAD网格上的原接触多边形单元压力。蓝色低压、红色高压，未承载几何为灰色；压力层关闭光照且直接使用二维图相同的线性色表RGB，不做点压力插值或假想椭圆。右图只改变相机，不放大变形或加宽接触区。原积分面积和力不由渲染几何重新定义。三维几何栅格由VTK离屏计算后嵌入科研图，SVG/PDF不冒充CAD矢量；显示正确不等于压力网格或物理精度已验证。',
                'Original contact-polygon cell pressures on the actual CAD mesh. Blue denotes low and red high pressure; unloaded geometry remains grey. Pressure lighting is disabled and cell RGB uses the same linear map as 2D views, without point-pressure interpolation or an assumed ellipse. The right panel changes only the camera, not deformation scale or patch width. Rendering does not redefine integration areas or forces. VTK offscreen geometry rasters are embedded in the figure; SVG/PDF are not claimed as CAD vectors. Correct display does not establish pressure-mesh or physical accuracy.')
            if nodal:
                caption+=_text(lang,
                    f' 此图为vertex-lumped nodal pressure：{overlay.get("active_pressure_node_count")}个原承载节点拥有{len(overlay["pressure_pa"])}个原三角形中位对偶格，原节点为力作用点；格数不等于节点数。',
                    f' This is vertex-lumped nodal pressure: {overlay.get("active_pressure_node_count")} original loaded nodes own {len(overlay["pressure_pa"])} original triangular median-dual cells. Forces act at the original nodes; cell and node counts differ.')
            if overlay.get('all_tooth_context'):
                caption+=_text(lang,
                    f' 仅映射物理齿k={overlay.get("tooth_offset")}；T1为有符号小轮扭矩，其他齿载荷和边界状态单独保存在源JSON，不能由该齿截图推断整副已收敛。',
                    f' Only physical tooth k={overlay.get("tooth_offset")} is mapped. T1 is signed pinion torque. Other-tooth loads and boundary states are retained in source JSON; this tooth image cannot establish full-pair convergence.')
                edges=[entry for entry in overlay['all_tooth_context'] if entry.get('loaded_cell_count',0)>0 and
                    entry.get('boundary_checks',{}).get('touches_trimmed_or_model_validity_edge')]
                for entry in edges:
                    caption+=_text(lang,
                        f' 邻齿k={entry["tooth_offset"]}仍有{entry["loaded_cell_count"]}个承载格且承载域触及实际裁切/有效域边缘（参考平面合力{entry["reference_plane_normal_force_n"]:.2f} N），该全副限制保留。',
                        f' Tooth k={entry["tooth_offset"]} still has {entry["loaded_cell_count"]} loaded cells and its support reaches an actual trimmed/validity edge (reference-plane resultant {entry["reference_plane_normal_force_n"]:.2f} N); this full-pair limitation remains.')
            record=_save(fig,output,stem,lang,rows,caption,height_mm=95,extra=dict(**scale,
                pressure_source_identity_sha256=overlay['source_identity_sha256'],pressure_interpolation='none',
                coordinate_displacement_m=0.,render_panels=[p[1] for p in panels],
                model_scope=overlay.get('model_scope'),tooth_offset=overlay.get('tooth_offset'),
                contact_representation=overlay.get('contact_representation'),active_pressure_node_count=overlay.get('active_pressure_node_count'),
                all_tooth_context=overlay.get('all_tooth_context'),pressure_vtp_file=vtp.relative_to(output).as_posix()))
    return record
