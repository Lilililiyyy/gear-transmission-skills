"""接触压力的连续标量显示场；保留原求解值、原面积和原节点力矩定义。

FE沿实际原三角形作C0线性插值；半解析沿实际中心网格作P1重构。
本模块只后处理，不重新求接触，不作RGB混色、滤波或假想斑点拟合。
"""
from __future__ import annotations

import math


def _positive(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
        raise ValueError(name+'须为正有限数')
    return float(value)


def _node_key(point,tolerance_m):
    return tuple(int(round(float(x)/tolerance_m)) for x in point)


def reconstruct_fe_pressure_field(result,*,tooth_offset=0,coordinate_tolerance_m=1e-12):
    """从真实中位对偶格恢复原三角形，按原节点压力构成连续P1标量场。

    p_i=f_i/a_i且a_i=ΣA_triangle/3，故∫p_h dA=Σf_i。原节点仍是
    原数值模型的力矩积分点；显示场的面积力矩不重新定义原FE力矩。
    """
    import numpy as np
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('物理齿偏移须为整数')
    tolerance=_positive(coordinate_tolerance_m,'坐标合并容差m')
    if not result.get('converged') or result.get('contact_representation')!='nodal_dual_cells':
        raise ValueError('需要已接受的真实FE节点接触结果')
    nodes={}
    for node in result.get('nodal_contacts',[]):
        if node['tooth_offset']!=tooth_offset or not node.get('active') or node['force_n']<=0:continue
        point=np.asarray(node['first_point_world_m'],dtype=float)
        if point.shape!=(3,) or not np.isfinite(point).all():raise ValueError('原压力节点坐标须有限三维')
        key=_node_key(point,tolerance)
        if key in nodes:raise ValueError('原压力节点坐标重复，不能合并独立未知量')
        pressure=_positive(node['pressure_pa'],'原节点压力Pa');area=_positive(node['area_m2'],'原节点集总面积m²')
        force=_positive(node['force_n'],'原节点力N')
        if not math.isclose(pressure*area,force,rel_tol=2e-10,abs_tol=1e-9):raise ValueError('原节点pA=f不成立')
        nodes[key]=dict(point=point,pressure=pressure,force=force,area=area,node_id=node['first_pressure_node_index'])
    if not nodes:raise ValueError('所选物理齿没有承载压力节点')
    triangles={};dual_force=0.;dual_area=0.;dual_count=0
    for cell in result.get('contacts',[]):
        if cell['tooth_offset']!=tooth_offset or not cell.get('active') or cell['force_n']<=0:continue
        quad=np.asarray(cell['first_vertices_world_m'],dtype=float)
        if quad.shape!=(4,3) or not np.isfinite(quad).all():raise ValueError('连续FE显示需要实际中位对偶四边形')
        vertex=quad[0];original=np.asarray([vertex,2*quad[1]-vertex,2*quad[3]-vertex])
        if np.linalg.norm(original.mean(axis=0)-quad[2])>10*tolerance:raise ValueError('四边形不是原FE三角形的中位对偶格')
        area=float(np.linalg.norm(np.cross(original[1]-original[0],original[2]-original[0]))*.5)
        if not math.isfinite(area) or area<=0 or not math.isclose(area/3,cell['area_m2'],rel_tol=2e-8,abs_tol=1e-18):
            raise ValueError('恢复原三角形面积与中位对偶积分面积不一致')
        key=_node_key(vertex,tolerance)
        if key not in nodes or nodes[key]['node_id']!=cell['first_pressure_node_index']:
            raise ValueError('对偶格原顶点与原压力节点身份不一致')
        if not math.isclose(cell['mean_pressure_pa'],nodes[key]['pressure'],rel_tol=2e-10,abs_tol=1e-5):
            raise ValueError('原对偶格压力与原节点压力不一致')
        if not math.isclose(cell['mean_pressure_pa']*cell['area_m2'],cell['force_n'],rel_tol=2e-10,abs_tol=1e-9):
            raise ValueError('原对偶格pA=f不成立')
        keys=tuple(_node_key(p,tolerance) for p in original);triangle_key=tuple(sorted(keys))
        if len(set(keys))!=3:raise ValueError('坐标容差合并了退化三角形')
        if triangle_key not in triangles:triangles[triangle_key]=(original,area)
        dual_force+=float(cell['force_n']);dual_area+=float(cell['area_m2']);dual_count+=1
    if not triangles:raise ValueError('没有实际原三角形可用于连续显示')
    points=[];lookup={};indices=[];areas=[]
    for original,area in triangles.values():
        tri=[]
        for point in original:
            key=_node_key(point,tolerance)
            if key not in lookup:lookup[key]=len(points);points.append(point.tolist())
            elif np.linalg.norm(np.asarray(points[lookup[key]])-point)>5*tolerance:raise ValueError('合并坐标不一致')
            tri.append(lookup[key])
        indices.append(tri);areas.append(area)
    pressures=np.zeros(len(points));node_ids=np.full(len(points),-1,dtype=int);tributary=np.zeros(len(points))
    for key,index in lookup.items():
        if key in nodes:pressures[index]=nodes[key]['pressure'];node_ids[index]=nodes[key]['node_id']
    indices=np.asarray(indices,dtype=int);areas=np.asarray(areas)
    for j in range(3):np.add.at(tributary,indices[:,j],areas/3)
    for key,node in nodes.items():
        if key not in lookup or not math.isclose(tributary[lookup[key]],node['area'],rel_tol=3e-8,abs_tol=1e-17):
            raise ValueError('恢复原三角形未完整覆盖承载节点的原集总面积')
    total_force=sum(node['force'] for node in nodes.values())
    integral=float(np.sum(areas*pressures[indices].mean(axis=1)))
    if not math.isclose(dual_force,total_force,rel_tol=3e-9,abs_tol=1e-8) or not math.isclose(integral,total_force,rel_tol=3e-8,abs_tol=1e-8):
        raise ValueError('连续P1压力场与原节点力积分不守恒')
    return dict(points_world_m=np.asarray(points),triangles=indices,pressure_pa=pressures,
        triangle_integration_area_m2=areas,original_pressure_node_ids=node_ids,
        reconstruction='C0-P1原FE三角形标量插值；不滤波、不改原节点压力',
        diagnostics=dict(original_loaded_node_count=len(nodes),original_dual_cell_count=dual_count,
            recovered_triangle_count=len(indices),original_contact_area_m2=dual_area,
            display_positive_support_area_m2=float(areas.sum()),original_normal_force_n=total_force,
            reconstructed_scalar_integral_n=integral,relative_integral_error=(integral-total_force)/total_force,
            scalar_force_conserved=True,source_peak_pressure_pa=float(pressures.max()),
            reconstructed_peak_pressure_pa=float(pressures.max()),pressure_values_modified=False,
            original_moment_application='保留原节点；不以连续显示场的面积力矩替代',
            display_area_definition='原三角形内P1场正压支撑面积，含零压顶点过渡，不是原节点对偶承载面积'))


def _reconstruct_center_pressure_field(result,*,tooth_offset=0):
    """按实际曲面中心规则网格构成连续P1场，保留真实有效域和原参考面积。

    内部压力节点的P1面积权重等于dxdy；真实边界承载时不外推缺失域，
    因而显示积分可能与原常压力格不同，差值必须单列而非重新归一化。
    """
    import numpy as np
    if not result.get('converged') or tooth_offset not in result.get('selected_tooth_offsets',[]):
        raise ValueError('需要所选物理齿已接受的规则面元压力结果')
    offsets=result['selected_tooth_offsets']
    if len(set(offsets))!=len(offsets) or len(result.get('fields',[]))!=len(offsets):raise ValueError('物理齿偏移与原压力场布局不一致')
    index=offsets.index(tooth_offset);field=result['fields'][index]
    patches=[p for p in result.get('geometry',{}).get('patches',[]) if p['tooth_offset']==tooth_offset]
    if len(patches)!=1:raise ValueError('所选物理齿的实际几何patch须唯一存在')
    patch=patches[0]
    pressure=np.asarray(field['pressure_pa'],dtype=float);valid=np.asarray(field['valid_mask'])
    force=np.asarray(field['forces_n'],dtype=float);points=np.asarray(patch['first_points_world_m'],dtype=float)
    if pressure.ndim!=2 or pressure.shape!=valid.shape or force.shape!=pressure.shape or valid.dtype!=np.bool_:
        raise ValueError('规则面元压力、力及有效域布局不一致')
    if points.shape!=pressure.shape+(3,) or np.any(valid & ~np.isfinite(points).all(axis=-1)):
        raise ValueError('有效中心须有实际三维曲面点')
    if not np.array_equal(valid,np.asarray(patch['valid_mask'])) or not np.isfinite(pressure[valid]).all() or np.any(pressure[valid]<0):
        raise ValueError('压力场和实际几何有效域不一致或有非法压力')
    if not np.isfinite(force).all() or np.any(force<0) or np.any(force[~valid]!=0):raise ValueError('实际有效域外不得承载，原单元力须有限非负')
    cell_area=_positive(float(field['cell_area_m2']),'规则参考单元面积m²')
    if not np.allclose(pressure[valid]*cell_area,force[valid],rtol=2e-10,atol=1e-9):raise ValueError('原面元压力乘参考面积不等于原力')
    ny,nx=pressure.shape;all_triangles=[]
    for j in range(ny-1):
        for i in range(nx-1):
            for corners in (((j,i),(j,i+1),(j+1,i+1)),((j,i),(j+1,i+1),(j+1,i))):
                if all(valid[q] for q in corners) and any(pressure[q]>0 for q in corners):
                    all_triangles.append(tuple(y*nx+x for y,x in corners))
    if not all_triangles:raise ValueError('所选物理齿没有可在真实有效域内重构的压力场')
    full_indices=np.asarray(all_triangles,dtype=int);used=np.unique(full_indices)
    inverse=np.full(ny*nx,-1,dtype=int);inverse[used]=np.arange(len(used));triangles=inverse[full_indices]
    values=pressure.ravel()[used];tri_areas=np.full(len(triangles),cell_area*.5)
    integral=float(np.sum(tri_areas*values[triangles].mean(axis=1)));total=float(force.sum())
    if total<=0:raise ValueError('原面元法向力须为正')
    conserved=math.isclose(integral,total,rel_tol=2e-9,abs_tol=1e-8)
    return dict(points_world_m=points.reshape(-1,3)[used],triangles=triangles,pressure_pa=values,
        triangle_integration_area_m2=tri_areas,grid_node_indices=np.column_stack((used//nx,used%nx)),
        reconstruction='C0-P1实际规则中心网格标量插值；有效域外不补面',
        diagnostics=dict(original_active_cell_count=int(np.count_nonzero(pressure[valid]>0)),
            recovered_triangle_count=len(triangles),original_contact_area_m2=float(np.count_nonzero(pressure[valid]>0)*cell_area),
            display_positive_support_area_m2=float(tri_areas.sum()),original_normal_force_n=total,
            reconstructed_scalar_integral_n=integral,relative_integral_error=(integral-total)/total,
            scalar_force_conserved=conserved,source_peak_pressure_pa=float(pressure[valid].max()),
            reconstructed_peak_pressure_pa=float(values.max()),pressure_values_modified=False,
            original_moment_application='保留原常压力格及原力臂；显示重构不重新定义力矩',
            display_area_definition='参考切平面中心网格P1场支撑面积，不是原常压力格承载面积',
            boundary_integration_note='内部承载时原参考面积权重守恒；真实边界承载可能缺少外侧中心支承，不外推或重标定'))


def reconstruct_halfspace_pressure_field(result,*,tooth_offset=0):
    """半解析P0压力的实际中心网格P1显示；原积分与显示积分分列报告。"""
    return _reconstruct_center_pressure_field(result,tooth_offset=tooth_offset)


def reconstruct_distributed_fe_pressure_field(result,*,geometry=None,tooth_offset=0):
    """将整轮FE的原P0分布牵引压力作守恒中心P1显示，不重求接触。

    只共用规则网格的标量插值，不调用半空间柔度。每个原中心压力
    保持不变；积分权重须与原投影面元面积相等。缺少非承载邻域时
    拒绝连续显示，不补面、不重标定；实际接触面积仍由原P0定义。
    """
    import numpy as np
    if (result.get('contact_representation')!='distributed_P0_surface_cells' or
            result.get('halfspace_compliance_used') is not False):
        raise ValueError('需要明确的纯整轮FE分布P0压力，不能用节点压力或混合柔度冒充')
    original_geometry=result.get('geometry')
    if geometry is None:geometry=original_geometry
    if not isinstance(geometry,dict):raise ValueError('须提供同一压力小片的实际几何')
    source=dict(result,geometry=geometry)
    field=_reconstruct_center_pressure_field(source,tooth_offset=tooth_offset)
    diagnostics=field['diagnostics']
    if not diagnostics['scalar_force_conserved']:
        raise ValueError('承载中心缺少完整邻域；连续显示积分不守恒，须扩展或加密实际求解域')
    # 积分合计偶然相符仍不足以证明每个原压力自由度的面积权重正确。
    offsets=result['selected_tooth_offsets'];index=offsets.index(tooth_offset)
    original=result['fields'][index];pressure=np.asarray(original['pressure_pa'],dtype=float)
    area=float(original['cell_area_m2']);weights=np.zeros_like(pressure)
    node_weights=np.zeros(len(field['pressure_pa']))
    for column in range(3):
        np.add.at(node_weights,field['triangles'][:,column],field['triangle_integration_area_m2']/3)
    indices=field['grid_node_indices'];weights[indices[:,0],indices[:,1]]=node_weights
    loaded=pressure>0
    if not np.allclose(weights[loaded],area,rtol=2e-10,atol=area*1e-12):
        raise ValueError('原承载面元的中心P1面积权重不完整，不能外推成连续压力')
    field['reconstruction']='整轮FE原P0压力的守恒中心P1显示；不改变压力、不调用半空间柔度'
    diagnostics.update(source_contact_representation='distributed_P0_surface_cells',
        display_pressure_basis='center_P1',pressure_solver_basis='P0',
        all_loaded_center_weights_conserved=True,halfspace_compliance_used=False,
        display_area_definition='中心P1插值的正压支撑面积；原FE承载面积仍为正压P0面元面积之和',
        original_moment_application='保留原Gauss分布牵引及FE节点力矩；显示场不替代原力矩积分')
    return field


def project_pressure_field(field,pinion_request,frames):
    """只把真实三维显示场投影到齿宽/齿高坐标，不改变压力及原积分面积。"""
    import numpy as np
    from .contact_pattern import pinion_flank_coordinates
    from .generated_surface import _request
    request=_request(pinion_request)
    projected=pinion_flank_coordinates(field['points_world_m'],request.machine.pitch_cone_angle_rad,
        request.machine.reference_cone_distance_m,rotation_rows=frames['pinion_rotation_rows'],
        translation_m=frames.get('translation_m',(0.,0.,0.)))
    xy=np.column_stack((projected['width_position_m'],projected['normal_height_m']))
    triangles=np.asarray(field['triangles']).copy()
    first=xy[triangles[:,1]]-xy[triangles[:,0]];second=xy[triangles[:,2]]-xy[triangles[:,0]]
    signed=first[:,0]*second[:,1]-first[:,1]*second[:,0]
    if np.any(np.abs(signed)<1e-24):raise ValueError('实际场在子午投影中出现退化三角形')
    triangles[signed<0]=triangles[signed<0][:,[0,2,1]]
    return dict(field,points_width_height_m=xy,triangles=triangles,projection_definition=projected['projection_definition'])


def draw_pressure_field(ax,field,*,norm,cmap,level_count=256):
    """按原三角形线性标量作连续等值填色；不按顶点RGB作Gouraud混色。"""
    import numpy as np
    from matplotlib.tri import Triangulation
    if isinstance(level_count,bool) or not isinstance(level_count,int) or level_count<16:raise ValueError('云图色阶须为不小于16的整数')
    points=np.asarray(field['points_width_height_m'])*1000;values=np.asarray(field['pressure_pa'])*1e-6
    mesh=Triangulation(points[:,0],points[:,1],np.asarray(field['triangles']))
    maximum=float(values.max());limit=float(norm.vmax)
    levels=np.linspace(0.,limit,level_count+1)
    return ax.tricontourf(mesh,values,levels=levels,norm=norm,cmap=cmap,
        extend='max' if maximum>limit else 'neither',antialiased=False,zorder=2)


def pressure_field_rows(field,*,label):
    """导出连续显示场的原标量与连接关系，便于独立重绘；不替代求解源表。"""
    import json
    rows=[]
    for i,(point,xy,p) in enumerate(zip(field['points_world_m'],field['points_width_height_m'],field['pressure_pa'])):
        rows.append(dict(row_kind='continuous_scalar_vertex',result_label=label,display_vertex_index=i,
            world_x_m=float(point[0]),world_y_m=float(point[1]),world_z_m=float(point[2]),
            width_position_m=float(xy[0]),normal_height_m=float(xy[1]),pressure_pa=float(p)))
    for i,(triangle,area) in enumerate(zip(field['triangles'],field['triangle_integration_area_m2'])):
        rows.append(dict(row_kind='continuous_scalar_triangle',result_label=label,display_triangle_index=i,
            display_vertex_indices_json=json.dumps(list(map(int,triangle))),display_triangle_integration_area_m2=float(area)))
    return rows
