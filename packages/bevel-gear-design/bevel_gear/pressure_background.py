"""压力云图的零值与未知域显示；完整面蓝底必须有绑定原状态的覆盖证书。"""
from __future__ import annotations

from dataclasses import asdict,is_dataclass
import hashlib
import json
import math


DOMAIN_CERTIFICATE_SCHEMA='contact_pressure_domain_certificate_v1'
_METHODS={'full_gear_fe','local_halfspace'}
_CHECKS={
    'full_gear_fe':('input_geometry_bound','included_candidates_kkt_passed','omitted_candidates_excluded',
        'unresolved_regions_excluded','far_field_motion_bound_passed'),
    'local_halfspace':('input_geometry_bound','full_working_flank_covered','omitted_regions_excluded',
        'unresolved_regions_excluded','loaded_boundary_checks_passed')}


def _canonical(value):
    """只用于身份摘要；未知浮点保留为null，不写回物理场。"""
    if isinstance(value,dict):return {str(k):_canonical(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [_canonical(v) for v in value]
    if value is None or isinstance(value,(str,bool,int)):return value
    if isinstance(value,float):return value if math.isfinite(value) else None
    if callable(getattr(value,'as_dict',None)):return _canonical(value.as_dict())
    if callable(getattr(value,'tolist',None)):return _canonical(value.tolist())
    if callable(getattr(value,'item',None)):return _canonical(value.item())
    if is_dataclass(value):return _canonical(asdict(value))
    raise TypeError('压力覆盖身份包含不能序列化的对象：'+type(value).__name__)


def _digest(value):
    return hashlib.sha256(json.dumps(_canonical(value),sort_keys=True,separators=(',',':'),
        ensure_ascii=False,allow_nan=False).encode('utf-8')).hexdigest()


def pressure_domain_fingerprint(result,*,method,tooth_offset=0):
    """绑定同一物理齿的原压力、间隙、实际曲面位置与位姿，不修改原数据。"""
    if method not in _METHODS:raise ValueError('压力域方法须为full_gear_fe或local_halfspace')
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('物理齿偏移须为整数')
    if not result.get('converged'):raise ValueError('覆盖证书只能绑定已接受的压力结果')
    if method=='full_gear_fe' and result.get('contact_representation')=='distributed_P0_surface_cells':
        from .distributed_pressure_domain import distributed_pressure_fingerprint
        return distributed_pressure_fingerprint(result,tooth_offset=tooth_offset)
    if method=='full_gear_fe':
        rows=[row for row in result.get('nodal_contacts',[]) if row.get('tooth_offset')==tooth_offset]
        if not rows:
            shares=result.get('pair_force_shares',{})
            if not isinstance(shares,dict) or not any(str(key)==str(tooth_offset) for key in shares):
                raise ValueError('空候选物理齿必须属于原实际齿对请求窗口')
            if not result.get('nodal_contacts'):
                raise ValueError('空候选齿身份仍须绑定全轮实际原节点，不能使用空压力结果')
        names=('first_pressure_node_index','first_point_world_m','second_point_world_m','normal_world',
            'signed_gap_m','remaining_gap_m','force_n','pressure_pa','area_m2','active')
        values=[{name:row.get(name) for name in names} for row in sorted(rows,key=lambda r:r['first_pressure_node_index'])]
        cell_names=('first_pressure_node_index','first_vertices_world_m','application_point_world_m',
            'mean_pressure_pa','force_n','area_m2','active')
        cells=[{name:row.get(name) for name in cell_names} for row in result.get('contacts',[])
            if row.get('tooth_offset')==tooth_offset]
        payload=dict(method=method,tooth_offset=tooth_offset,driving_angle_rad=result['driving_angle_rad'],
            pinion_torque_nm=result['pinion_torque_nm'],frames=result['mesh_reference_frames'],
            elastic_te_rad=result.get('elastic_te_rad'),nodal_contacts=values,dual_cells=cells)
        if not rows:
            # 空局部候选不等于零压证明；远域排除仍须使用这些全轮原力重算响应。
            payload['empty_candidate_global_binding']=dict(pair_force_shares=result['pair_force_shares'],
                nodal_contacts=result['nodal_contacts'],contacts=result.get('contacts',[]))
    else:
        offsets=result.get('selected_tooth_offsets',[])
        if tooth_offset not in offsets:
            geometry=result.get('geometry',{})
            requested=geometry.get('build_parameters',{}).get('tooth_offsets',[])
            patches=[p for p in geometry.get('patches',[]) if p.get('tooth_offset')==tooth_offset]
            if tooth_offset not in requested or len(patches)!=1:
                raise ValueError('无有效格的齿必须属于原几何请求窗口且有唯一真实patch')
            if not offsets or len(set(offsets))!=len(offsets) or len(offsets)!=len(result.get('fields',[])):
                raise ValueError('未选齿身份仍须绑定完整非空的全局真实半解析解')
            payload=dict(method=method,tooth_offset=tooth_offset,driving_angle_rad=geometry['driving_angle_rad'],
                pinion_torque_nm=geometry['pinion_torque_nm'],theta_rad=result.get('theta_rad'),
                absent_selected_field_global_binding=dict(selected_tooth_offsets=offsets,
                    fields=result['fields'],geometry=geometry))
            return _digest(payload)
        if offsets.count(tooth_offset)!=1 or len(offsets)!=len(result.get('fields',[])):
            raise ValueError('覆盖证书找不到唯一的半解析物理齿压力场')
        field=result['fields'][offsets.index(tooth_offset)];geometry=result['geometry']
        patches=[p for p in geometry['patches'] if p['tooth_offset']==tooth_offset]
        if len(patches)!=1:raise ValueError('覆盖证书须有唯一实际曲面几何patch')
        patch=patches[0]
        payload=dict(method=method,tooth_offset=tooth_offset,driving_angle_rad=geometry['driving_angle_rad'],
            pinion_torque_nm=geometry['pinion_torque_nm'],theta_rad=result.get('theta_rad'),
            frames=patch['frames'],grid=patch.get('grid'),diagnostics=patch.get('diagnostics'),
            geometry={name:patch.get(name) for name in ('gaps_m','valid_mask','first_points_world_m',
                'first_corner_points_world_m','corner_valid_mask')},
            field={name:field.get(name) for name in ('pressure_pa','forces_n','valid_mask','contact_mask','cell_area_m2')})
    return _digest(payload)


def flank_background_fingerprint(background_polygons_m):
    """绑定本次显示使用的完整实际裁切齿面三角形，避免证书套到其他背景。"""
    if not isinstance(background_polygons_m,(list,tuple)) or not background_polygons_m:
        raise ValueError('覆盖证书须绑定完整实际齿面三角形')
    for polygon in background_polygons_m:
        if len(polygon)!=3 or any(len(point)!=2 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in point) for point in polygon):
            raise ValueError('实际齿面背景须为有限二维m坐标三角形')
    return _digest(background_polygons_m)


def global_fe_response_fingerprint(result):
    """绑定整轮响应所用全部节点载荷、位姿、TE和扭矩；不改变局部历史摘要。"""
    if result.get('contact_representation')=='distributed_P0_surface_cells':
        from .distributed_pressure_domain import distributed_response_fingerprint
        return distributed_response_fingerprint(result)
    return _digest(dict(nodal_contacts=result['nodal_contacts'],elastic_te_rad=result.get('elastic_te_rad'),
        frames=result['mesh_reference_frames'],pinion_torque_nm=result['pinion_torque_nm']))


def _partial_audit_identity(certificate,state_sha,method):
    """部分证书只接受正式审核运行中保持不变的几何、KKT及源码证据。"""
    checks=certificate.get('checks',{})
    required=('input_geometry_bound','audit_inputs_and_source_unchanged',
        'included_candidates_kkt_passed' if method=='full_gear_fe' else 'loaded_boundary_checks_passed')
    if not isinstance(checks,dict) or any(checks.get(name) is not True for name in required):
        raise ValueError('部分覆盖证书缺少几何、原KKT或源码不变检查')
    identity=certificate.get('audit_identity',{})
    if not isinstance(identity,dict):raise ValueError('部分覆盖证书缺少审计身份')
    source=identity.get('source_sha256_before')
    if (identity.get('unchanged') is not True or not isinstance(source,str) or len(source)!=64 or
            any(c not in '0123456789abcdef' for c in source) or identity.get('source_sha256_after')!=source or
            identity.get('state_fingerprint_before')!=state_sha or identity.get('state_fingerprint_after')!=state_sha):
        raise ValueError('部分覆盖证书的审计输入／源码身份不完整或不一致')
    response=certificate.get('fe_response_recomputed' if method=='full_gear_fe' else 'kkt_recomputed',{})
    failures='violations' if method=='full_gear_fe' else 'failures'
    if not isinstance(response,dict) or response.get('passed') is not True or response.get(failures)!=[]:
        raise ValueError('部分覆盖证书没有通过原KKT重算')
    if method=='full_gear_fe' and response.get('input_node_geometry_bound') is not True:
        raise ValueError('部分FE覆盖证书没有绑定原节点几何')
    omitted=('omitted_candidates_excluded','unresolved_regions_excluded','far_field_motion_bound_passed') if method=='full_gear_fe' else (
        'full_working_flank_covered','omitted_regions_excluded','unresolved_regions_excluded')
    if any(checks.get(name) is not False for name in omitted):
        raise ValueError('部分覆盖证书的完整面检查状态与未决域不一致')


def _partitioned_domain(certificate,background_polygons_m,result,pinion_request,method,tooth_offset,*,allow_complete=False):
    """验证中点四分树叶片完整分割，投影缓存与轻量汇总分开。

    每个原CAD三角形的叶片须正好覆盖一个无重复、无缺口的四叉树。
    仅比较总面积会漏掉等面积的重复／缺损，因此另外检查叶片路径区间。
    """
    import numpy as np
    from .contact_pattern import pinion_flank_coordinates
    from .generated_surface import _request
    region_data=certificate.get('regions')
    if not isinstance(region_data,dict) or not isinstance(region_data.get('regions'),list) or not region_data['regions']:
        raise ValueError('部分覆盖证书须保留完整原三角形分片')
    rows=region_data['regions'];count=len(rows);source_count=len(background_polygons_m)
    maximum=region_data.get('maximum_depth')
    if isinstance(maximum,bool) or not isinstance(maximum,int) or not 0<=maximum<=30:
        raise ValueError('覆盖分片最大深度须为0至30的整数')
    if region_data.get('input_triangle_count')!=source_count or region_data.get('terminal_region_count')!=count:
        raise ValueError('覆盖分片的原三角形或终端区域数量不一致')
    states=('original_solved_domain','excluded_by_distance_bound','unresolved')
    if any(not isinstance(row,dict) or row.get('state') not in states for row in rows):
        raise ValueError('覆盖分片含未定义的域状态')
    for row in rows:
        source=row.get('source_triangle_index');depth=row.get('depth')
        if (isinstance(source,bool) or not isinstance(source,int) or not 0<=source<source_count or
                isinstance(depth,bool) or not isinstance(depth,int) or not 0<=depth<=maximum):
            raise ValueError('覆盖分片的原三角形归属或深度不合法')
    try:world=np.asarray([row['vertices_world_m'] for row in rows],dtype=float)
    except (KeyError,TypeError,ValueError) as error:raise ValueError('覆盖分片须保留真实三维顶点') from error
    if world.shape!=(count,3,3) or not np.isfinite(world).all():raise ValueError('覆盖分片须为有限三维三角形')
    areas=np.linalg.norm(np.cross(world[:,1]-world[:,0],world[:,2]-world[:,0]),axis=1)*.5
    declared=np.asarray([row.get('area_m2') for row in rows],dtype=float)
    if not np.isfinite(declared).all() or np.any(areas<=0) or not np.allclose(areas,declared,rtol=1e-8,atol=1e-23):
        raise ValueError('覆盖分片的真实面积与声明不一致')
    if pinion_request is None:raise ValueError('部分覆盖证书投影需要同一小轮制造请求')
    request=_request(pinion_request)
    if method=='full_gear_fe' and result.get('contact_representation')=='distributed_P0_surface_cells':
        frames=next(p for p in result['geometry']['patches'] if p['tooth_offset']==tooth_offset)['frames']
    elif method=='full_gear_fe':frames=result['mesh_reference_frames']
    else:frames=next(p for p in result['geometry']['patches'] if p['tooth_offset']==tooth_offset)['frames']
    converted=pinion_flank_coordinates(world.reshape(-1,3),request.machine.pitch_cone_angle_rad,
        request.machine.reference_cone_distance_m,rotation_rows=frames['pinion_rotation_rows'],
        translation_m=frames.get('translation_m',(0.,0.,0.)))
    projected=np.column_stack((converted['width_position_m'],converted['normal_height_m'])).reshape(count,3,2)
    sources=np.asarray([row['source_triangle_index'] for row in rows],dtype=int)
    depths=np.asarray([row['depth'] for row in rows],dtype=int)
    if not np.array_equal(np.unique(sources),np.arange(source_count)):
        raise ValueError('覆盖分片缺少原CAD三角形')
    # 世界点与背景投影绑定至0.1nm；三维叶片面积仍用原世界坐标计算。
    coordinate_tolerance=1e-10
    order=np.argsort(sources,kind='stable');splits=np.searchsorted(sources[order],np.arange(source_count+1))
    original_area=0.
    for source in range(source_count):
        indices=order[splits[source]:splits[source+1]]
        candidates=projected[indices].reshape(-1,2);points=world[indices].reshape(-1,3)
        background=np.asarray(background_polygons_m[source]);distance=((candidates[:,None,:]-background[None,:,:])**2).sum(axis=2)
        nearest=np.argmin(distance,axis=0)
        if np.any(distance[nearest,np.arange(3)]>coordinate_tolerance**2):
            raise ValueError('覆盖分片没有保留原CAD三角形顶点')
        root=points[nearest];basis=(root[1:]-root[0]).T;gram=basis.T@basis
        root_area=float(np.linalg.norm(np.cross(root[1]-root[0],root[2]-root[0]))*.5)
        if root_area<=0:raise ValueError('覆盖分片的原CAD三角形退化')
        if len(indices)==1 and depths[indices[0]]==0:
            # 整个原三角形本身是唯一深0叶片：顶点一一绑定及面积已足以证明分割。
            # 保留原投影数组与求和顺序，不再对恒等坐标做逆矩阵和四分树恢复。
            if len(set(map(int,nearest)))!=3:
                raise ValueError('原CAD顶点匹配不是一一对应')
            if not math.isclose(float(areas[indices[0]]),root_area,rel_tol=1e-8,abs_tol=1e-22):
                raise ValueError('覆盖分片没有守恒原三角形面积')
            original_area+=root_area
            continue
        try:xy=(world[indices]-root[0])@basis@np.linalg.inv(gram)
        except np.linalg.LinAlgError as error:raise ValueError('覆盖分片的原三角形坐标不可逆') from error
        residual=world[indices]-(root[0]+xy@basis.T)
        if np.any(np.linalg.norm(residual,axis=2)>coordinate_tolerance):raise ValueError('覆盖分片离开原CAD三角形平面')
        bary=np.concatenate((1.-xy.sum(axis=2,keepdims=True),xy),axis=2)
        leaf_depth=depths[indices];codes=np.zeros(len(indices),dtype=np.int64)
        for level in range(int(leaf_depth.max())):
            active=np.flatnonzero(leaf_depth>level);value=bary[active];centre=value.mean(axis=1)
            vertex=np.argmax(centre,axis=1);child=np.where(centre.max(axis=1)>.5,vertex,3)
            codes[active]=codes[active]*4+child
            for corner in range(3):
                ids=np.flatnonzero(child==corner)
                if len(ids):value[ids]*=2;value[ids,:,corner]-=1
            ids=np.flatnonzero(child==3)
            if len(ids):value[ids]=1.-2*value[ids][:,:,[2,0,1]]
            bary[active]=value
        # 叶片必须是该四分树的整个三角形，而非同重心的任意小片。
        identity=np.eye(3);vertex_error=np.linalg.norm(bary[:,:,None,:]-identity[None,None,:,:],axis=3)
        nearest_vertex=np.argmin(vertex_error,axis=2)
        if np.any(vertex_error.min(axis=2)>2e-7) or np.any(np.sort(nearest_vertex,axis=1)!=np.arange(3)):
            raise ValueError('覆盖分片不是原CAD三角形的完整中点四分叶片')
        end_depth=int(leaf_depth.max());weights=4**(end_depth-leaf_depth)
        starts=codes*weights;ends=(codes+1)*weights;sequence=np.argsort(starts)
        if starts[sequence[0]]!=0 or ends[sequence[-1]]!=4**end_depth or np.any(starts[sequence[1:]]!=ends[sequence[:-1]]):
            raise ValueError('覆盖分片有缺损、重复或重叠')
        if not math.isclose(float(areas[indices].sum()),root_area,rel_tol=1e-8,abs_tol=1e-22):
            raise ValueError('覆盖分片没有守恒原三角形面积')
        original_area+=root_area
    actual_counts={state:sum(row['state']==state for row in rows) for state in states}
    actual_areas={state:math.fsum(float(area) for row,area in zip(rows,areas) if row['state']==state) for state in states}
    unknown_count=actual_counts['unresolved']
    expected_passed=unknown_count==0
    if (isinstance(certificate.get('unresolved_region_count'),bool) or isinstance(region_data.get('unresolved_region_count'),bool) or
            certificate.get('unresolved_region_count')!=unknown_count or region_data.get('unresolved_region_count')!=unknown_count or
            (unknown_count<=0 and not allow_complete) or region_data.get('passed') is not expected_passed):
        raise ValueError('部分覆盖证书的未决区域计数或状态不一致')
    declared_areas=region_data.get('area_m2',{})
    if not isinstance(declared_areas,dict) or set(declared_areas)!=set(states) or any(
            isinstance(declared_areas.get(state),bool) or not isinstance(declared_areas.get(state),(int,float)) or not math.isfinite(declared_areas[state]) or
            not math.isclose(declared_areas[state],actual_areas[state],rel_tol=1e-8,abs_tol=1e-22) for state in states):
        raise ValueError('覆盖分片面积分类汇总不一致')
    expected_sources=sorted({row['source_triangle_index'] for row in rows if row['state']=='unresolved'})
    if region_data.get('unresolved_source_triangles')!=expected_sources:
        raise ValueError('覆盖分片未决三角形归属汇总不一致')
    selected=np.asarray([row['state']=='unresolved' for row in rows])
    return dict(partition_summary=dict(input_triangle_count=source_count,terminal_region_count=count,
        state_region_counts=actual_counts,area_m2=actual_areas,total_area_m2=original_area,
        unresolved_area_fraction=actual_areas['unresolved']/original_area,
        coordinate_binding_tolerance_m=coordinate_tolerance,
        chord_error_allowance_m=region_data.get('chord_error_allowance_m'),maximum_depth=maximum,
        evidence_scope=region_data.get('evidence_scope'),
        partition_validation='原CAD顶点绑定；中点四分叶片路径完整且无重复；三维面积逐片及分类核对',
        partition_fingerprint=_digest([dict(source_triangle_index=row['source_triangle_index'],depth=row['depth'],
            state=row['state'],area_m2=row['area_m2'],vertices_world_m=row['vertices_world_m']) for row in rows])),
        _unresolved_polygons_m=projected[selected],_unresolved_source_indices=sources[selected])


def pressure_domain_metadata(info):
    """只导出身份、区域数量与面积；绘图缓存不进入JSON或源表。"""
    return {key:value for key,value in info.items() if not key.startswith('_')}


def _unresolved_display_geometry(info):
    """可选GEOS精确并集合并内部边；不简化、缓冲、吸附或扩大未知域。"""
    import numpy as np
    if '_unresolved_render_cache' in info:return info['_unresolved_render_cache']
    polygons=np.asarray(info['_unresolved_polygons_m']);sources=np.asarray(info['_unresolved_source_indices'])
    try:
        from shapely.geometry import Polygon
        from shapely.geometry.polygon import orient
        from shapely.ops import unary_union
    except ImportError:unary_union=None
    cache=[];merged_count=0
    for source in np.unique(sources):
        pieces=polygons[sources==source];entry=dict(source=int(source),triangles=pieces,rings=None)
        if unary_union is not None:
            union=unary_union([Polygon(triangle) for triangle in pieces])
            edge1=pieces[:,1]-pieces[:,0];edge2=pieces[:,2]-pieces[:,0]
            projected_area=float(np.abs(edge1[:,0]*edge2[:,1]-edge1[:,1]*edge2[:,0]).sum()*.5)
            if (union.is_valid and union.geom_type in ('Polygon','MultiPolygon') and
                    math.isclose(union.area,projected_area,rel_tol=1e-8,abs_tol=1e-20)):
                parts=[union] if union.geom_type=='Polygon' else list(union.geoms)
                entry['rings']=[]
                for part in parts:
                    shape=orient(part,sign=1.)
                    entry['rings'].append([np.asarray(shape.exterior.coords)]+
                        [np.asarray(interior.coords) for interior in shape.interiors])
                merged_count+=1
        cache.append(entry)
    info['_unresolved_render_cache']=cache
    info['partition_summary']['unresolved_display_geometry']=dict(
        engine='GEOS unary_union' if merged_count else 'original_triangles',
        merged_source_triangle_count=merged_count,source_triangle_count=len(cache),
        coordinate_snapping=False,buffer_applied=False,simplification_applied=False,
        union_area_verified_source_count=merged_count,
        rendering='vector_union_with_raster_fallback_for_original_triangles')
    return cache


def pressure_domain_background(result,background_polygons_m,*,method,tooth_offset=0,certificate=None,pinion_request=None):
    """检查外部覆盖证据的身份和必要检查项；绘图本身不重新证明非接触。"""
    if method=='full_gear_fe' and result.get('contact_representation')=='distributed_P0_surface_cells':
        from .distributed_pressure_domain import distributed_pressure_domain_background
        return distributed_pressure_domain_background(result,background_polygons_m,tooth_offset=tooth_offset,
            certificate=certificate,pinion_request=pinion_request)
    background_sha=flank_background_fingerprint(background_polygons_m)
    state_sha=pressure_domain_fingerprint(result,method=method,tooth_offset=tooth_offset)
    if method=='full_gear_fe':angle=result['driving_angle_rad'];torque=result['pinion_torque_nm']
    else:angle=result['geometry']['driving_angle_rad'];torque=result['geometry']['pinion_torque_nm']
    if certificate is None:
        return dict(full_flank_zero_certified=False,unknown_domain_present=True,method=method,tooth_offset=tooth_offset,
            complete_working_flank=False,certificate_passed=None,
            driving_angle_rad=angle,pinion_torque_nm=torque,
            state_fingerprint=state_sha,flank_background_fingerprint=background_sha,
            domain_state='partial_or_uncertified',zero_definition='仅实际已求解域中的零压力；完整面外域不假设为零')
    if not isinstance(certificate,dict) or certificate.get('schema')!=DOMAIN_CERTIFICATE_SCHEMA:
        raise ValueError('压力覆盖证书schema不受支持')
    for name,expected in (('method',method),('tooth_offset',tooth_offset),('state_fingerprint',state_sha),
            ('flank_background_fingerprint',background_sha)):
        if certificate.get(name)!=expected:raise ValueError('压力覆盖证书与当前状态／背景不匹配：'+name)
    for name,expected in (('driving_angle_rad',angle),('pinion_torque_nm',torque)):
        value=certificate.get(name)
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not math.isclose(value,expected,rel_tol=1e-12,abs_tol=1e-12):
            raise ValueError('压力覆盖证书相位或载荷不一致：'+name)
    if method=='full_gear_fe' and 'fe_response_state_fingerprint' in certificate:
        if certificate['fe_response_state_fingerprint']!=global_fe_response_fingerprint(result):
            raise ValueError('压力证书的完整FE响应输入与当前全轮原状态不同')
    partial=certificate.get('passed') is False and certificate.get('complete_working_flank') is False
    if certificate.get('passed') is not True or certificate.get('complete_working_flank') is not True:
        if partial:
            _partial_audit_identity(certificate,state_sha,method)
            partition=_partitioned_domain(certificate,background_polygons_m,result,pinion_request,method,tooth_offset)
            evidence=certificate.get('evidence_kind')
            if not isinstance(evidence,str) or not evidence.strip():raise ValueError('部分压力覆盖证书须说明证据类型')
            return dict(full_flank_zero_certified=False,unknown_domain_present=True,method=method,tooth_offset=tooth_offset,
                complete_working_flank=False,certificate_passed=False,
                driving_angle_rad=angle,pinion_torque_nm=torque,state_fingerprint=state_sha,
                flank_background_fingerprint=background_sha,evidence_kind=evidence,checks=dict(certificate['checks']),
                domain_state='partitioned_pressure_domain',
                zero_definition='仅已求解或距离界排除的真实分片显示零值蓝；未决分片保持灰色网纹',**partition)
        raise ValueError('完整齿面蓝色零值底需要已通过的全工作面覆盖证书')
    count=certificate.get('unresolved_region_count')
    if isinstance(count,bool) or not isinstance(count,int) or count!=0:raise ValueError('存在未认证区域，不能把完整齿面标成零压')
    if not isinstance(certificate.get('evidence_kind'),str) or not certificate['evidence_kind'].strip():
        raise ValueError('完整压力覆盖证书须说明证据类型')
    checks=certificate.get('checks',{})
    if not isinstance(checks,dict) or any(checks.get(name) is not True for name in _CHECKS[method]):
        raise ValueError('压力覆盖证书缺少本方法必要的几何、遗漏域或平衡检查')
    return dict(full_flank_zero_certified=True,unknown_domain_present=False,method=method,tooth_offset=tooth_offset,
        complete_working_flank=True,certificate_passed=True,
        driving_angle_rad=angle,pinion_torque_nm=torque,
        state_fingerprint=state_sha,flank_background_fingerprint=background_sha,
        evidence_kind=certificate['evidence_kind'],checks={name:checks[name] for name in _CHECKS[method]},
        domain_state='complete_pressure_domain',zero_definition='完整实际工作面中，重构正压区外为经认证的非接触零压力')


def pressure_background_rows(info,*,label):
    """背景状态单独导出；不把用于绘图的蓝色基底当作原求解压力表。"""
    row=dict(row_kind='pressure_domain_state',result_label=label,domain_state=info['domain_state'],
        method=info['method'],tooth_offset=info['tooth_offset'],driving_angle_rad=info['driving_angle_rad'],
        pinion_torque_nm=info['pinion_torque_nm'],
        full_flank_zero_certified=info['full_flank_zero_certified'],unknown_domain_present=info['unknown_domain_present'],
        state_fingerprint=info['state_fingerprint'],flank_background_fingerprint=info['flank_background_fingerprint'],
        zero_definition=info['zero_definition'])
    summary=info.get('partition_summary')
    if summary is not None:
        row.update(partition_fingerprint=summary['partition_fingerprint'],input_triangle_count=summary['input_triangle_count'],
            terminal_region_count=summary['terminal_region_count'],unresolved_area_fraction=summary['unresolved_area_fraction'])
        for state,area in summary['area_m2'].items():
            row[state+'_area_m2']=area;row[state+'_region_count']=summary['state_region_counts'][state]
    return [row]


def _zero_background_artist(polygons_m,colour,*,zorder):
    """零值域作一次复合路径填充，消除PDF逐三角形抗锯齿白缝。

    GEOS只求原三角形的精确并集，保留孔洞且核对投影面积；不用缓冲、
    吸附或简化。缺少GEOS或并集检查未通过时，仅将背景按导出DPI栅格化。
    """
    import numpy as np
    from matplotlib.collections import PathCollection,PolyCollection
    from matplotlib.path import Path
    pieces=np.asarray(polygons_m,dtype=float)
    try:
        from shapely.geometry import Polygon
        from shapely.geometry.polygon import orient
        from shapely.ops import unary_union
        from shapely.errors import GEOSException
    except ImportError:unary_union=None
    if unary_union is not None:
        try:union=unary_union([Polygon(triangle) for triangle in pieces])
        except (ValueError,GEOSException):union=None
        relative=pieces-pieces[:,:1];following=np.roll(relative,-1,axis=1)
        projected_area=float(np.abs((relative[:,:,0]*following[:,:,1]-relative[:,:,1]*following[:,:,0]).sum(axis=1)).sum()*.5)
        if (union is not None and union.is_valid and union.geom_type in ('Polygon','MultiPolygon') and
                math.isclose(union.area,projected_area,rel_tol=1e-8,abs_tol=1e-20)):
            paths=[]
            for part in ([union] if union.geom_type=='Polygon' else union.geoms):
                shape=orient(part,sign=1.)
                for ring in [shape.exterior,*shape.interiors]:paths.append(Path(np.asarray(ring.coords)*1000,closed=True))
            return PathCollection([Path.make_compound_path(*paths)],facecolors=[colour],edgecolors='none',
                antialiaseds=False,zorder=zorder)
    return PolyCollection(pieces*1000,facecolors=[colour],edgecolors='none',antialiaseds=False,
        rasterized=True,zorder=zorder)


def draw_pressure_background(ax,background_polygons_m,*,norm,cmap,info,known_polygons_m=()):
    """在实际裁切边界内绘制零压蓝底；未认证外域保持灰色网纹。"""
    import numpy as np
    from matplotlib.collections import PolyCollection
    zero=tuple(cmap(norm(0.)))
    partial=info['domain_state']=='partitioned_pressure_domain'
    if info['full_flank_zero_certified'] or partial:
        base=_zero_background_artist(background_polygons_m,zero,zorder=0)
    else:
        base=PolyCollection(np.asarray(background_polygons_m)*1000,facecolors='#D5DADD',edgecolors='#8C979E',
            linewidths=0.,hatch='///',antialiaseds=False,zorder=0)
    ax.add_collection(base)
    if partial:
        # 已验证分片完整分割全部CAD背景，蓝底与逐片涂蓝等价；灰片单独盖回。
        from matplotlib.path import Path
        from matplotlib.patches import PathPatch
        # 每组只裁切至自己的原CAD三角形，避免数万微片逐一遍历整面裁切路径。
        for entry in _unresolved_display_geometry(info):
            triangle=np.asarray(background_polygons_m[entry['source']])*1000
            clip=Path(np.vstack((triangle,triangle[0])),closed=True)
            if entry['rings'] is None:
                unresolved=PolyCollection(entry['triangles']*1000,facecolors='#D5DADD',
                    edgecolors='#8C979E',linewidths=0.,hatch='///',antialiaseds=False,zorder=4,rasterized=True)
                unresolved.set_clip_path(clip,ax.transData)
                # 自定义CAD路径不会自动与Axes默认边框相交，局部放大须同时裁窗。
                unresolved.set_clip_box(ax.bbox);ax.add_collection(unresolved)
            else:
                for rings in entry['rings']:
                    paths=[Path(ring*1000,closed=True) for ring in rings]
                    unresolved=PathPatch(Path.make_compound_path(*paths),facecolor='#D5DADD',edgecolor='#8C979E',
                        linewidth=0.,hatch='///',antialiased=False,zorder=4)
                    unresolved.set_clip_path(clip,ax.transData)
                    unresolved.set_clip_box(ax.bbox);ax.add_patch(unresolved)
    elif not info['full_flank_zero_certified'] and len(known_polygons_m):
        from matplotlib.path import Path
        paths=[Path(np.vstack((polygon,polygon[0])),closed=True) for polygon in np.asarray(background_polygons_m)*1000]
        known=_zero_background_artist(known_polygons_m,zero,zorder=.5)
        # 实际有效单元也只在原CAD裁切面内显示，不以矩形窗口补出齿面。
        known.set_clip_path(Path.make_compound_path(*paths),ax.transData)
        known.set_clip_box(ax.bbox)
        ax.add_collection(known)
    return base


def add_pressure_domain_legend(ax,info,*,norm,cmap,lang):
    """图级共用图例置于留白，避免窄小的全齿面定位图被图例遮住。"""
    from matplotlib.patches import Patch
    if lang not in ('zh','en'):raise ValueError('图例语言须为zh或en')
    fig=ax.figure
    previous=getattr(fig,'_pressure_domain_legend',None)
    unknown=info['unknown_domain_present'] or getattr(fig,'_pressure_unknown_domain',False)
    if previous is not None:previous.remove()
    fig._pressure_unknown_domain=unknown
    handles=[Patch(facecolor=cmap(norm(0.)),edgecolor='none',label='p = 0：非接触' if lang=='zh' else 'p = 0: no contact')]
    if unknown:
        handles.append(Patch(facecolor='#D5DADD',edgecolor='#8C979E',hatch='///',
            label='未决／未认证' if lang=='zh' else 'Unresolved / uncertified'))
    fig._pressure_domain_legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.47,.018),
        ncol=len(handles),fontsize=6,frameon=False,handlelength=1.2,columnspacing=1.5)
    fig._pressure_domain_legend.set_in_layout(False)
