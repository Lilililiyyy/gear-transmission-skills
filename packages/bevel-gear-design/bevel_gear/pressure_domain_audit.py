"""完整压力背景的几何排除审核；未决域保留，不修改原压力和掩码。

证书只覆盖输入实际工作面三角离散及显式弦差余量，不能替代连续CAD
严格区间证明、压力网格收敛或物理验证。远域由有限BRep最短距离的
1-Lipschitz界排除，近域必须被原求解域覆盖；不接受调用方手填通过。
"""
from __future__ import annotations

import math


def _progress(callback,stage,**values):
    if callback is not None:callback(dict(stage=stage,**values))


def _loaded_wheel_rotation(frames,theta_rad):
    """绕实际世界wheel_axis、共同原点转动，再作用原局部→世界位姿。"""
    import numpy as np
    theta=float(theta_rad);axis=np.asarray(frames['wheel_axis_world'],dtype=float)
    if axis.shape!=(3,) or not np.isfinite(axis).all() or not math.isfinite(theta) or abs(float(np.linalg.norm(axis))-1)>1e-10:
        raise ValueError('加载刚体旋转须有单位世界轴与有限角度')
    x,y,z=axis;skew=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
    rotation=np.eye(3)+math.sin(theta)*skew+(1-math.cos(theta))*(skew@skew)
    return rotation@np.asarray(frames['wheel_rotation_rows'],dtype=float)


def _kernel_sha256():
    from pathlib import Path
    import hashlib
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _finalize_certificate(certificate,result,*,state_before,source_before):
    """审核不许改变原状态，运行中源码变化不能拿旧结论冒充当前内核。"""
    from .pressure_background import pressure_domain_fingerprint
    state_after=pressure_domain_fingerprint(result,method=certificate['method'],tooth_offset=certificate['tooth_offset'])
    source_after=_kernel_sha256()
    unchanged=state_before==state_after and source_before==source_after
    certificate['audit_identity']=dict(source_sha256_before=source_before,source_sha256_after=source_after,
        state_fingerprint_before=state_before,state_fingerprint_after=state_after,unchanged=unchanged)
    certificate['checks']['audit_inputs_and_source_unchanged']=unchanged
    certificate['passed']=bool(certificate['passed'] and unchanged)
    return certificate


def _positive(value, name, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError(name + ('须为非负有限数' if zero else '须为正有限数'))
    return float(value)


def _selected_triangles(flank_mesh, triangle_indices):
    """统一OCP原生tuple与JSON回读list索引，禁止tuple被NumPy解释成多轴。"""
    import numpy as np
    vertices=np.asarray(flank_mesh['vertices_m'],dtype=float)
    indices=np.asarray(flank_mesh['triangles'],dtype=int)
    if indices.ndim!=2 or indices.shape[1]!=3:
        raise ValueError('实际齿面三角形索引须为N×3')
    return vertices[indices[np.asarray(triangle_indices,dtype=int)]]


class FiniteFaceDistanceOracle:
    """全部指定有限BRep面的点距，世界输入m，OCP内部mm；失败不返回零。"""
    def __init__(self, faces, rotation_rows, *, tolerance_m=1e-10):
        import numpy as np
        from OCP.BRep import BRep_Builder
        from OCP.TopoDS import TopoDS_Compound
        from OCP.Precision import Precision
        from .rigid_contact import RigidTransform
        self.rotation = np.asarray(RigidTransform(rotation_rows=rotation_rows).rotation_rows)
        self.tolerance_m = max(_positive(tolerance_m, '距离查询容差m'),float(Precision.Confusion_s())*.001)
        faces = tuple(faces)
        if not faces:
            raise ValueError('有限面距离查询须给非空真实BRep面')
        self.shape = TopoDS_Compound()
        builder = BRep_Builder(); builder.MakeCompound(self.shape)
        for face in faces:
            builder.Add(self.shape, face)
        self.query_count = 0

    def __call__(self, point_world_m):
        import numpy as np
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.gp import gp_Pnt
        point = np.asarray(point_world_m, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError('距离查询点须为有限三维m坐标')
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*(self.rotation.T @ point * 1000))).Vertex()
        query = BRepExtrema_DistShapeShape(vertex, self.shape)
        query.SetDeflection(self.tolerance_m*1000); query.Perform()
        self.query_count += 1
        if not query.IsDone():
            return dict(done=False, distance_m=None, error_bound_m=self.tolerance_m)
        distance = float(query.Value()) * .001
        return dict(done=math.isfinite(distance) and distance >= 0, distance_m=distance,
            error_bound_m=self.tolerance_m)


class TriangleDistanceOracle:
    """实际FE三角并集的欧氏最近距离；可选VTK加速，不筛法向/候选带。"""
    def __init__(self,triangles_world_m,*,tolerance_m=1e-10,backend='auto'):
        import numpy as np
        triangles=np.asarray(triangles_world_m,dtype=float)
        if triangles.ndim!=3 or triangles.shape[1:]!=(3,3) or not len(triangles) or not np.isfinite(triangles).all():
            raise ValueError('真实FE距离查询须给非空有限三角形N×3×3')
        if backend not in ('auto','exact','vtk'):raise ValueError('三角距离后端须为auto/exact/vtk')
        self.triangles=triangles;self.tolerance_m=_positive(tolerance_m,'三角距离误差m');self.query_count=0;self.locator=None
        if backend!='exact':
            try:
                import vtk
                from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray
                points=vtk.vtkPoints();points.SetData(numpy_to_vtk(triangles.reshape(-1,3),deep=True))
                ids=np.arange(len(triangles)*3,dtype=np.int64).reshape(-1,3)
                cells=vtk.vtkCellArray();cells.SetCells(len(ids),numpy_to_vtkIdTypeArray(np.column_stack((np.full(len(ids),3),ids)).ravel(),deep=True))
                poly=vtk.vtkPolyData();poly.SetPoints(points);poly.SetPolys(cells)
                self.locator=vtk.vtkStaticCellLocator();self.locator.SetDataSet(poly);self.locator.BuildLocator()
            except ImportError:
                if backend=='vtk':raise

    def __call__(self,point_world_m):
        import numpy as np
        from .contact_mesh import _closest_on_triangles
        point=np.asarray(point_world_m,dtype=float)
        if point.shape!=(3,) or not np.isfinite(point).all():raise ValueError('三角距离查询点须为有限三维m坐标')
        self.query_count+=1
        if self.locator is None:
            _,_,squared=_closest_on_triangles(point,self.triangles,np);distance=math.sqrt(max(0.,float(squared.min())))
        else:
            import vtk
            closest=[0.,0.,0.];index=vtk.mutable(0);sub=vtk.mutable(0);squared=vtk.mutable(0.)
            self.locator.FindClosestPoint(point.tolist(),closest,index,sub,squared)
            _,_,exact=_closest_on_triangles(point,self.triangles[int(index):int(index)+1],np)
            distance=math.sqrt(max(0.,float(exact[0])))
        return dict(done=math.isfinite(distance),distance_m=distance,error_bound_m=self.tolerance_m)


def _periodic_pinion_triangles(triangles,tooth_offset,pinion_teeth,wheel_teeth):
    """canonical槽0转到实际齿(-k)%z1；大轮对应k%z2，不改原相位。"""
    import numpy as np
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('物理齿偏移须为整数')
    for z in (pinion_teeth,wheel_teeth):
        if isinstance(z,bool) or not isinstance(z,int) or z<3:raise ValueError('实际齿数须为不少于3的整数')
    angle=-tooth_offset*2*math.pi/pinion_teeth;c,s=math.cos(angle),math.sin(angle)
    rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
    return np.asarray(triangles)@rotation.T,(-tooth_offset)%pinion_teeth,tooth_offset%wheel_teeth


def audit_triangle_pressure_domain(triangles_world_m, *, distance_oracle, solved_region_oracle,
        relative_motion_bound, chord_error_m, maximum_depth=8, minimum_triangle_area_m2=1e-15,
        maximum_regions=300000, tolerance_m=1e-10,progress=None):
    """逐片证明“原求解域覆盖”或“无接触”，返回未决实数与原片归属。

    距离函数须返回done/distance_m/error_bound_m；真实几何距离在整个
    三角形上的下界为d(重心)−最大顶点半径−查询误差−弦差。运动界为
    同一参考位姿到加载位姿的最大相对位移。回调是低层几何接口，高层
    入口会从真实接触状态/材料核/FE位移构造，不接收通过标志。
    """
    import numpy as np
    triangles = np.asarray(triangles_world_m, dtype=float)
    if triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not len(triangles) or not np.isfinite(triangles).all():
        raise ValueError('须给非空有限N×3×3实际工作面三角形m坐标')
    chord = _positive(chord_error_m, '弦差余量m', zero=True)
    tol = _positive(tolerance_m, '排除容差m')
    floor = _positive(minimum_triangle_area_m2, '最小细分面积m²')
    if isinstance(maximum_depth, bool) or not isinstance(maximum_depth, int) or maximum_depth < 0:
        raise ValueError('细分深度须为非负整数')
    if isinstance(maximum_regions, bool) or not isinstance(maximum_regions, int) or maximum_regions < len(triangles):
        raise ValueError('区域上限须不少于输入三角形数')
    rows = []; processed = 0
    stack = [(i, 0, triangle) for i, triangle in enumerate(triangles)]
    while stack:
        source, depth, triangle = stack.pop(); processed += 1
        if processed==1 or processed%1000==0:
            _progress(progress,'domain_regions',processed_regions=processed,pending_regions=len(stack),terminal_regions=len(rows))
        centre = triangle.mean(axis=0)
        radius = float(np.linalg.norm(triangle-centre, axis=1).max())
        area = float(np.linalg.norm(np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0]))*.5)
        if area <= 0:
            raise ValueError('实际工作面含退化三角形')
        covered = bool(solved_region_oracle(triangle))
        query = None; lower = None; motion = None
        if not covered:
            query = distance_oracle(centre)
            motion = _positive(relative_motion_bound(triangle), '相对运动界m', zero=True)
            if query.get('done') and query.get('distance_m') is not None:
                distance = _positive(query['distance_m'], '实际有限面距离m', zero=True)
                error = _positive(query.get('error_bound_m', 0.), '距离误差m', zero=True)
                lower = distance-radius-error-chord
        far = not covered and lower is not None and lower > motion+tol
        if covered or far:
            state = 'original_solved_domain' if covered else 'excluded_by_distance_bound'
        elif depth < maximum_depth and area > floor and processed+len(stack)+4 < maximum_regions:
            a,b,c = triangle; ab=(a+b)/2; bc=(b+c)/2; ca=(c+a)/2
            stack.extend((source,depth+1,np.array(piece)) for piece in ((a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)))
            continue
        else:
            state = 'unresolved'
        rows.append(dict(source_triangle_index=source, depth=depth, state=state, area_m2=area, vertices_world_m=triangle.tolist(),
            centre_world_m=centre.tolist(), radius_m=radius, distance_lower_bound_m=lower,
            relative_motion_upper_bound_m=motion, distance_query=query))
    areas = {name:math.fsum(r['area_m2'] for r in rows if r['state']==name)
        for name in ('original_solved_domain','excluded_by_distance_bound','unresolved')}
    unresolved = [r for r in rows if r['state']=='unresolved']
    _progress(progress,'domain_regions_completed',processed_regions=processed,terminal_regions=len(rows),unresolved_regions=len(unresolved))
    return dict(passed=not unresolved, unresolved_region_count=len(unresolved),
        input_triangle_count=len(triangles), terminal_region_count=len(rows), evaluated_region_count=processed,
        area_m2=areas, unresolved_source_triangles=sorted({r['source_triangle_index'] for r in unresolved}),
        minimum_excluded_clearance_margin_m=min((r['distance_lower_bound_m']-r['relative_motion_upper_bound_m']
            for r in rows if r['state']=='excluded_by_distance_bound'), default=None),
        chord_error_allowance_m=chord, maximum_depth=maximum_depth, regions=rows,
        evidence_scope='输入实际工作面三角离散；弦差是显式几何余量，非连续CAD严格区间证明')


def _grid_cover(patch, field, *, padding_m=0.):
    """三角形投影包围矩形全部位于原有效格并集才认可，宁可保留未决。"""
    import numpy as np
    grid = patch['grid'].as_dict() if hasattr(patch['grid'],'as_dict') else patch['grid']
    valid = np.asarray(patch['valid_mask'], dtype=bool)
    if not np.array_equal(valid, np.asarray(field['valid_mask'],dtype=bool)):
        raise ValueError('原压力场与原几何有效掩码不同')
    origin = np.asarray(patch['plane_origin_world_m']); basis=np.asarray(patch['plane_tangent_rows'])
    spacing=np.array([grid['spacing_x_m'],grid['spacing_y_m']]); start=np.asarray(grid['origin_xy_m'])
    def covered(triangle):
        xy=(triangle-origin)@basis.T
        lo=np.floor((xy.min(axis=0)-padding_m-start)/spacing+.5-1e-10).astype(int)
        hi=np.floor((xy.max(axis=0)+padding_m-start)/spacing+.5+1e-10).astype(int)
        if np.any(lo<0) or hi[0]>=valid.shape[1] or hi[1]>=valid.shape[0]: return False
        return bool(valid[lo[1]:hi[1]+1,lo[0]:hi[0]+1].all())
    return covered


def _halfspace_motion(result, offset, *, include_rigid_motion=True):
    """从原矩形力源推导任意片域位移上界，不从局部采样最大值猜全域界。"""
    import numpy as np
    from .halfspace_contact import rectangular_potential_integral
    index=result['selected_tooth_offsets'].index(offset); field=result['fields'][index]
    patch=next(p for p in result['geometry']['patches'] if p['tooth_offset']==offset)
    operator=result['operators'][index]; grid=field['grid']; forces=np.asarray(field['forces_n'])
    j,i=np.nonzero(forces>0); source_forces=forces[j,i]
    dx=grid['spacing_x_m']; dy=grid['spacing_y_m']; x0,y0=grid['origin_xy_m']
    sources=np.column_stack((x0+i*dx,y0+j*dy)); coefficient=operator['inverse_effective_modulus_per_pa']/(math.pi*dx*dy)
    origin=np.asarray(patch['plane_origin_world_m']); basis=np.asarray(patch['plane_tangent_rows'])
    # patch.frames属于最后几何线性化参考；只计从该参考到加载相位的剩余转角。
    theta=float(result.get('linearization_step_theta_rad',result['theta_rad']))
    def bound(triangle):
        xy=(triangle-origin)@basis.T
        nearest=np.maximum(xy.min(axis=0),np.minimum(sources,xy.max(axis=0)))
        delta=nearest-sources
        elastic=float(np.sum(source_forces*rectangular_potential_integral(delta[:,0],delta[:,1],dx/2,dy/2))*coefficient)
        # 距原点半径控制刚体转动；此上界允许方向最不利，不利用“弹性只张开”。
        radius=float(np.linalg.norm(triangle,axis=1).max())
        return elastic+(2*radius*abs(math.sin(theta/2)) if include_rigid_motion else 0.)
    return bound, dict(rectangular_source_count=len(source_forces), source_force_n=float(source_forces.sum()),
        inverse_effective_modulus_per_pa=operator['inverse_effective_modulus_per_pa'], remaining_rigid_angle_rad=theta,
        rigid_motion_already_applied=not include_rigid_motion,
        derivation='每个正矩形力源的核在片域投影包围盒最接近源中心处取最大值，逐源正和；已查询加载刚体位姿时不重复扣旋转位移')


def _kkt_halfspace(result, gap_tolerance_m):
    import numpy as np
    failures=[]
    for offset,field in zip(result['selected_tooth_offsets'],result['fields']):
        valid=np.asarray(field['valid_mask'],dtype=bool); f=np.asarray(field['forces_n'],dtype=float); p=np.asarray(field['pressure_pa'],dtype=float)
        gap=np.asarray(field['remaining_gap_m'],dtype=float); active=f>1e-12
        if not valid.any() or not np.isfinite(f).all() or not np.isfinite(p[valid]).all() or not np.isfinite(gap[valid]).all():
            failures.append(dict(tooth_offset=offset,reason='非有限或空有效域'));continue
        if np.any(f < -1e-12) or np.any(f[~valid] != 0) or np.any(gap[valid] < -gap_tolerance_m) or np.any(abs(gap[active]) > gap_tolerance_m):
            failures.append(dict(tooth_offset=offset,reason='原压力/间隙不满足KKT容差'))
        if not np.allclose(p[valid]*field['cell_area_m2'],f[valid],rtol=1e-10,atol=1e-10):
            failures.append(dict(tooth_offset=offset,reason='原pA与f不一致'))
    torque=None;target=None
    if 'geometry' in result:
        patches={p['tooth_offset']:p for p in result['geometry']['patches']}
        torque=sum(float(np.sum(np.asarray(field['forces_n'])[np.asarray(field['valid_mask'],dtype=bool)]*
            np.asarray(patches[offset]['wheel_lever_arms_m'])[np.asarray(field['valid_mask'],dtype=bool)]))
            for offset,field in zip(result['selected_tooth_offsets'],result['fields']))
        ratio=result['geometry']['pinion_to_wheel_ratio']
        target=result['geometry']['pinion_torque_nm']/ratio
        if not math.isclose(torque,target,rel_tol=1e-7,abs_tol=1e-7):
            failures.append(dict(reason='重构从动轴力矩不满足全局平衡',torque_nm=torque,target_nm=target))
    return dict(passed=not failures,gap_tolerance_m=gap_tolerance_m,failures=failures,
        wheel_torque_reconstructed_nm=torque,wheel_target_torque_nm=target)


def audit_generated_halfspace_pressure_domain(prepared, result, flank_mesh, pinion_request, *,
        tooth_offset=0, maximum_depth=8, maximum_regions=300000, tolerance_m=1e-10,progress=None):
    """对原真实半解析结果生成显示helper可读取的全背景覆盖证书。

    prepared须为同请求实际有限工作面。flank_mesh须为同CAD三角化，
    其metadata.linear_deflection_m用于显式弦差余量。任何未决使证书
    passed=False；不将原斜率拒绝、Newton未解或ROI外域直接设成零压。
    """
    import numpy as np
    from .generated_halfspace import PreparedGeneratedHalfspace
    from .generated_surface import _request
    from .contact_pattern import select_pinion_left_flank
    from .pressure_background import DOMAIN_CERTIFICATE_SCHEMA, pressure_domain_fingerprint, flank_background_fingerprint
    if not isinstance(prepared,PreparedGeneratedHalfspace) or not result.get('converged'):
        raise ValueError('须给实际有限面准备对象及已收敛原半解析压力')
    source_before=_kernel_sha256();state_before=pressure_domain_fingerprint(result,method='local_halfspace',tooth_offset=tooth_offset)
    _progress(progress,'halfspace_geometry_binding')
    request=_request(pinion_request)
    if prepared.tca.requests['pinion'] != request:raise ValueError('完整齿面与原压力的制造请求不同')
    index=result['selected_tooth_offsets'].index(tooth_offset); field=result['fields'][index]
    patch=next(p for p in result['geometry']['patches'] if p['tooth_offset']==tooth_offset)
    selected=select_pinion_left_flank(flank_mesh,request); coordinates=selected['coordinates']
    ids=selected['triangle_indices'];triangles=_selected_triangles(flank_mesh,ids)
    background=[[[coordinates['width_position_m'][i],coordinates['normal_height_m'][i]]
        for i in flank_mesh['triangles'][j]] for j in ids]
    rotation=np.asarray(patch['frames']['pinion_rotation_rows']); world=triangles@rotation.T
    deflection=_positive(flank_mesh.get('metadata',{}).get('linear_deflection_m'),'实际CAD三角化弦差m')
    first=FiniteFaceDistanceOracle(prepared.trimmed['pinion'].faces,rotation.tolist(),tolerance_m=tolerance_m)
    theta=float(result.get('linearization_step_theta_rad',result['theta_rad']))
    loaded_rotation=_loaded_wheel_rotation(patch['frames'],theta)
    second=FiniteFaceDistanceOracle(prepared.trimmed['wheel'].faces,loaded_rotation.tolist(),tolerance_m=tolerance_m)
    vertices=np.unique(world.reshape(-1,3),axis=0)
    distances=[first(v) for v in vertices]
    binding=all(d['done'] and d['distance_m']<=deflection+tolerance_m for d in distances)
    # 原有效格实际点也须确实在同一有限BRep，不能仅比较请求名称。
    valid=np.asarray(patch['valid_mask'],dtype=bool); points=np.asarray(patch['first_points_world_m'])[valid]
    probes=points[np.linspace(0,len(points)-1,min(25,len(points)),dtype=int)]
    pressure_distances=[first(v) for v in probes]
    binding=binding and bool(len(probes)) and all(d['done'] and d['distance_m']<=max(tolerance_m,1e-9) for d in pressure_distances)
    motion,motion_evidence=_halfspace_motion(result,tooth_offset,include_rigid_motion=False)
    _progress(progress,'halfspace_region_audit',input_triangles=len(world),geometry_bound=bool(binding))
    regions=audit_triangle_pressure_domain(world,distance_oracle=second,solved_region_oracle=_grid_cover(patch,field,padding_m=deflection),
        relative_motion_bound=motion,chord_error_m=deflection,maximum_depth=maximum_depth,
        maximum_regions=maximum_regions,tolerance_m=tolerance_m,progress=progress)
    gap_tol=max(2e-10,float(result.get('certificate',{}).get('gap_acceptance_m',tolerance_m))*1.01)
    kkt=_kkt_halfspace(result,gap_tol)
    boundaries=result.get('contact_boundary_checks',[])
    boundary_ok=bool(boundaries) and all(not b.get(name,True) for b in boundaries for name in
        ('touches_window_edge','touches_model_validity_edge','touches_unresolved_query_edge','touches_torque_direction_edge'))
    complete=bool(regions['passed'])
    checks=dict(input_geometry_bound=bool(binding),full_working_flank_covered=complete,
        omitted_regions_excluded=complete,unresolved_regions_excluded=complete,
        loaded_boundary_checks_passed=bool(boundary_ok and kkt['passed']))
    certificate=dict(schema=DOMAIN_CERTIFICATE_SCHEMA,method='local_halfspace',tooth_offset=tooth_offset,
        driving_angle_rad=result['geometry']['driving_angle_rad'],pinion_torque_nm=result['geometry']['pinion_torque_nm'],
        state_fingerprint=pressure_domain_fingerprint(result,method='local_halfspace',tooth_offset=tooth_offset),
        flank_background_fingerprint=flank_background_fingerprint(background),passed=all(checks.values()),
        complete_working_flank=complete,unresolved_region_count=regions['unresolved_region_count'],
        evidence_kind='实际有限BRep点距＋1-Lipschitz三角片下界＋原有效格覆盖；三角离散及显式弦差层级',
        checks=checks,regions=regions,motion_bound=motion_evidence,kkt_recomputed=kkt,
        geometry_binding=dict(checked_vertices=len(vertices),maximum_vertex_distance_m=max((d['distance_m'] for d in distances if d['done']),default=None),
            pressure_points_checked=len(probes),triangulation_deflection_allowance_m=deflection),
        loaded_opposing_rotation_rows=loaded_rotation.tolist(),
        limitations=['仅所选物理齿与当前相位，不认证其他齿或整个周期','不认证连续CAD严格覆盖、局部峰压精度、材料寿命或真实物理正确性'])
    return _finalize_certificate(certificate,result,state_before=state_before,source_before=source_before)


def _triangle_union_cover(triangles_world_m, origin, basis, *, padding_m=0.,backend='geos'):
    """GEOS普通三角并集覆盖；可显式选Python同数学差集作回归基准。

    不buffer、不snap、不取凸包，不使用假设无重叠的coverage_union。
    padding只将查询片扩为原有轴对齐包围矩形，不扩张已知求解域。
    """
    import numpy as np
    covers=(np.asarray(triangles_world_m)-origin)@basis.T
    lower=covers.min(axis=1); upper=covers.max(axis=1)
    if backend not in ('geos','python'):raise ValueError('三角并集后端须为geos或python')
    if backend=='geos':
        try:
            import shapely
            from shapely.geometry import Polygon,box
            from shapely.ops import unary_union
            from shapely.prepared import prep
        except ImportError as error:raise ImportError('完整FE压力域审核需要可选Shapely/GEOS专业几何引擎') from error
        polygons=[Polygon(triangle) for triangle in covers]
        if any(not polygon.is_valid for polygon in polygons):raise ValueError('投影FE三角形无效，不能修补后冒充原求解域')
        union=unary_union(polygons)
        if union.is_empty or not union.is_valid:raise ArithmeticError('GEOS原三角并集为空或无效，未进行buffer修补')
        prepared=prep(union)
        def geos_covered(triangle):
            xy=(triangle-origin)@basis.T
            if padding_m:
                lo=xy.min(axis=0)-padding_m;hi=xy.max(axis=0)+padding_m
                footprint=box(float(lo[0]),float(lo[1]),float(hi[0]),float(hi[1]))
            else:footprint=Polygon(xy)
            return bool(footprint.is_valid and footprint.area>1e-24 and prepared.covers(footprint))
        geos_covered.geometry_engine=dict(backend='GEOS',shapely_version=shapely.__version__,
            geos_version=getattr(shapely,'geos_version_string',None),operation='unary_union + prepared.covers',
            source_triangle_count=len(covers),union_area_m2=float(union.area),coordinate_snapping=False,buffer_applied=False)
        return geos_covered
    def cross(a,b):return float(a[0]*b[1]-a[1]*b[0])
    def area(poly):
        if len(poly)<3:return 0.
        return abs(sum(cross(poly[(i+1)%len(poly)]-poly[0],poly[i]-poly[0]) for i in range(1,len(poly)-1)))*.5
    def clip(poly,a,b,inside):
        out=[]
        for p,q in zip(poly,poly[1:]+poly[:1]):
            dp=cross(b-a,p-a);dq=cross(b-a,q-a)
            pin=dp>=0 if inside else dp<=0; qin=dq>=0 if inside else dq<=0
            if pin:out.append(p)
            if pin!=qin and abs(dp-dq)>1e-30:out.append(p+(q-p)*dp/(dp-dq))
        return out
    def covered(triangle):
        xy=(triangle-origin)@basis.T
        if padding_m:
            lo=xy.min(axis=0)-padding_m;hi=xy.max(axis=0)+padding_m
            xy=np.array([[lo[0],lo[1]],[hi[0],lo[1]],[hi[0],hi[1]],[lo[0],hi[1]]])
        total=area(list(xy))
        if total<=1e-24:return False
        candidates=np.flatnonzero(np.all(upper>=xy.min(axis=0),axis=1)&np.all(lower<=xy.max(axis=0),axis=1))
        remaining=[list(xy)]
        for index in candidates:
            cover=covers[index]
            if cross(cover[1]-cover[0],cover[2]-cover[0])<0:cover=cover[::-1]
            updated=[]
            for poly in remaining:
                carry=poly
                for a,b in zip(cover,np.roll(cover,-1,axis=0)):
                    outside=clip(carry,a,b,False)
                    if area(outside)>max(1e-24,total*1e-12):updated.append(outside)
                    carry=clip(carry,a,b,True)
                    if len(carry)<3:break
            remaining=updated
            if not remaining:return True
            if len(remaining)>1000:return False
        return sum(area(p) for p in remaining)<=max(1e-24,total*1e-10)
    covered.geometry_engine=dict(backend='python',operation='凸多边形逐片差集',source_triangle_count=len(covers))
    return covered


def _fe_displacement_evidence(prepared,result,*,gap_tolerance_m=1e-10):
    """按原节点力重新求FE响应并重构全部已有约束KKT，不信任旧passed字段。"""
    import numpy as np
    from .elasticity import ContactLoadChannel,contact_load_matrix
    rows=result['nodal_contacts'];forces=np.asarray([r['force_n'] for r in rows]);displacements={}
    errors=[];failures=[];frames=result['mesh_reference_frames'];theta=result['elastic_te_rad']
    for member,prefix in (('pinion','first'),('wheel','second')):
        channels=[ContactLoadChannel(tuple(r[prefix+'_node_indices']),tuple(r[prefix+'_weights']),tuple(r[prefix+'_force_normal_local'])) for r in rows]
        matrix=contact_load_matrix(prepared.systems[member],channels,matrix_format='csc',maximum_sparse_channels=100000)
        response=prepared.systems[member].solve(matrix@forces)
        displacements[member]=np.asarray(response['displacements_m'])
        nodes=np.asarray(prepared.systems[member].mesh.nodes_m); rotation=np.asarray(frames[member+'_rotation_rows'])
        for r in rows:
            point=np.asarray(r[prefix+'_weights'])@nodes[r[prefix+'_node_indices'],:]@rotation.T
            errors.append(float(np.linalg.norm(point-np.asarray(r[prefix+'_point_world_m']))))
    reconstruction=[]
    for r in rows:
        elastic=0.
        for member,prefix in (('pinion','first'),('wheel','second')):
            elastic+=float(np.asarray(r[prefix+'_weights'])@displacements[member][r[prefix+'_node_indices'],:]@np.asarray(r[prefix+'_force_normal_local']))
        gap=r['signed_gap_m']+elastic-r['wheel_lever_m']*theta
        reconstruction.append(abs(gap-r['remaining_gap_m']))
        if r['force_n'] < -1e-12 or gap < -gap_tolerance_m or (r['force_n']>1e-12 and abs(gap)>gap_tolerance_m):
            failures.append(dict(first_pressure_node_index=r['first_pressure_node_index'],remaining_gap_m=gap,force_n=r['force_n']))
        if not math.isclose(r['pressure_pa']*r['area_m2'],r['force_n'],rel_tol=1e-10,abs_tol=1e-9):
            failures.append(dict(first_pressure_node_index=r['first_pressure_node_index'],reason='原节点pA不等于f'))
    maxima={m:float(np.linalg.norm(u,axis=1).max()) for m,u in displacements.items()}
    radius=float(np.linalg.norm(np.asarray(prepared.systems['wheel'].mesh.nodes_m),axis=1).max())
    motion=sum(maxima.values())+2*radius*abs(math.sin(theta/2))
    ratio=prepared.models['pinion'].request.tooth_count/prepared.models['wheel'].request.tooth_count
    target=result['pinion_torque_nm']/ratio
    torque=sum(r['force_n']*r['wheel_lever_m'] for r in rows)
    torque_ok=math.isclose(torque,target,rel_tol=1e-7,abs_tol=1e-7)
    if not torque_ok:failures.append(dict(reason='重构从动轴力矩不满足全局平衡',torque_nm=torque,target_nm=target))
    return dict(passed=bool(not failures and max(reconstruction,default=math.inf)<=gap_tolerance_m),
        input_node_geometry_bound=bool(errors and max(errors)<=1e-10), maximum_node_position_error_m=max(errors,default=None),
        included_constraint_count=len(rows),maximum_kkt_reconstruction_error_m=max(reconstruction,default=None),
        violations=failures,maximum_displacement_m=maxima,wheel_radius_bound_m=radius,relative_motion_bound_m=motion,
        wheel_torque_reconstructed_nm=torque,wheel_target_torque_nm=target,
        definition='原节点力重算全轮FE位移；三角线性位移的范数上界取全部节点最大值，加从动轮刚体转角弦长')


def audit_generated_fe_pressure_domain(prepared,result,flank_mesh,pinion_request,*,tooth_offset=0,
        matching_tolerance_m=8e-7,maximum_depth=8,maximum_regions=300000,tolerance_m=1e-10,progress=None):
    """完整指定实际齿FE压力域证书；远域独立距界，近域须原候选节点三角并集覆盖。

    实际周期齿取小轮(-k)%z1、大轮k%z2；canonical背景先作对应周期
    旋转再进入原FE世界位姿。不改相位、齿偏移、压力或任何FE未知量。
    """
    import numpy as np
    from .generated_ltca import PreparedGeneratedLtca,select_generated_flank_facets
    from .generated_halfspace import _faces
    from .generated_surface import _request
    from .contact_pattern import select_pinion_left_flank
    from .pressure_background import DOMAIN_CERTIFICATE_SCHEMA,pressure_domain_fingerprint,flank_background_fingerprint
    if not isinstance(prepared,PreparedGeneratedLtca) or not result.get('converged'):
        raise ValueError('须给实际全轮FE准备对象及已收敛原节点压力')
    source_before=_kernel_sha256();state_before=pressure_domain_fingerprint(result,method='full_gear_fe',tooth_offset=tooth_offset)
    _progress(progress,'fe_geometry_binding')
    request=_request(pinion_request)
    if prepared.models['pinion'].request!=request:raise ValueError('完整齿面与FE制造请求不同')
    selected=select_pinion_left_flank(flank_mesh,request); coordinates=selected['coordinates']; ids=selected['triangle_indices']
    triangles=_selected_triangles(flank_mesh,ids)
    background=[[[coordinates['width_position_m'][i],coordinates['normal_height_m'][i]] for i in flank_mesh['triangles'][j]] for j in ids]
    frames=result['mesh_reference_frames'];rotation=np.asarray(frames['pinion_rotation_rows'])
    physical,pinion_index,wheel_index=_periodic_pinion_triangles(triangles,tooth_offset,request.tooth_count,prepared.models['wheel'].request.tooth_count)
    world=physical@rotation.T
    deflection=_positive(flank_mesh.get('metadata',{}).get('linear_deflection_m'),'实际CAD三角化弦差m')
    oracles={}
    loaded_rotation=_loaded_wheel_rotation(frames,result['elastic_te_rad'])
    for member,side in (('pinion','left_flank'),('wheel','right_flank')):
        model=prepared.models[member]
        tooth=pinion_index if member=='pinion' else wheel_index
        frame=frames[member+'_rotation_rows'] if member=='pinion' else loaded_rotation.tolist()
        oracles[member]=FiniteFaceDistanceOracle(_faces(model.shape,model.slot_face_groups[tooth][side]),frame,tolerance_m=tolerance_m)
    binding_queries=[oracles['pinion'](p) for p in np.unique(world.reshape(-1,3),axis=0)]
    binding=all(q['done'] and q['distance_m']<=deflection+tolerance_m for q in binding_queries)
    _progress(progress,'fe_response_reconstruction',geometry_bound=bool(binding))
    response=_fe_displacement_evidence(prepared,result,gap_tolerance_m=max(2e-10,tolerance_m*1.01))
    _progress(progress,'fe_original_domain',kkt_passed=response['passed'],input_node_geometry_bound=response['input_node_geometry_bound'])
    facets=select_generated_flank_facets(prepared.systems['pinion'].mesh,prepared.models['pinion'].shell.surfaces['left_flank'],
        tooth_count=request.tooth_count,tooth_index=pinion_index,matching_tolerance_m=matching_tolerance_m)
    candidates={r['first_pressure_node_index'] for r in result['nodal_contacts'] if r['tooth_offset']==tooth_offset}
    facet_ids=[tri for tri in facets['triangles'] if all(int(i) in candidates for i in tri)]
    if not facet_ids:raise ValueError('所选实际齿没有三顶点均属于原节点候选域的真实FE三角形')
    covered_triangles=np.asarray(facets['nodes_m'])[facet_ids]@rotation.T
    normals=np.asarray([r['normal_world'] for r in result['nodal_contacts'] if r['tooth_offset']==tooth_offset]);normal=normals.mean(axis=0);normal/=np.linalg.norm(normal)
    tangent=np.cross(normal,[1.,0.,0.] if abs(normal[0])<.9 else [0.,1.,0.]);tangent/=np.linalg.norm(tangent)
    basis=np.array([tangent,np.cross(normal,tangent)]);origin=covered_triangles[0,0]
    cover=_triangle_union_cover(covered_triangles,origin,basis,padding_m=deflection)
    elastic_bound=sum(response['maximum_displacement_m'].values())
    # 原接触约束采用转角一阶线性化。查询精确加载刚体位姿后，保留
    # |R(theta)q−q−theta*a×q|≤0.5*|q|*theta²，不能把它偷换成非线性接触通过。
    linearization_remainder=.5*response['wheel_radius_bound_m']*result['elastic_te_rad']**2
    loaded_motion_bound=elastic_bound+linearization_remainder
    _progress(progress,'fe_opposing_geometry',elastic_bound_m=elastic_bound,rotation_linearization_remainder_m=linearization_remainder)
    wheel_facets=select_generated_flank_facets(prepared.systems['wheel'].mesh,prepared.models['wheel'].shell.surfaces['right_flank'],
        tooth_count=prepared.models['wheel'].request.tooth_count,tooth_index=wheel_index,matching_tolerance_m=matching_tolerance_m)
    second_triangles=np.asarray(wheel_facets['nodes_m'])[wheel_facets['triangles']]@loaded_rotation.T
    fe_distance=TriangleDistanceOracle(second_triangles,tolerance_m=tolerance_m)
    def combined_distance(point):
        cad=oracles['wheel'](point);fe=fe_distance(point)
        done=bool(cad['done'] and fe['done'])
        return dict(done=done,distance_m=min(cad['distance_m'],fe['distance_m']) if done else None,
            error_bound_m=max(cad['error_bound_m'],fe['error_bound_m']),cad_distance_m=cad['distance_m'],fe_triangle_distance_m=fe['distance_m'])
    _progress(progress,'fe_region_audit',input_triangles=len(world),candidate_triangles=len(facet_ids),opposing_triangles=len(second_triangles),loaded_elastic_motion_bound_m=loaded_motion_bound)
    regions=audit_triangle_pressure_domain(world,distance_oracle=combined_distance,solved_region_oracle=cover,
        relative_motion_bound=lambda t:loaded_motion_bound,chord_error_m=deflection,
        maximum_depth=maximum_depth,maximum_regions=maximum_regions,tolerance_m=tolerance_m,progress=progress)
    complete=bool(regions['passed'])
    checks=dict(input_geometry_bound=bool(binding and response['input_node_geometry_bound']),included_candidates_kkt_passed=response['passed'],
        omitted_candidates_excluded=complete,unresolved_regions_excluded=complete,far_field_motion_bound_passed=complete)
    certificate=dict(schema=DOMAIN_CERTIFICATE_SCHEMA,method='full_gear_fe',tooth_offset=tooth_offset,
        driving_angle_rad=result['driving_angle_rad'],pinion_torque_nm=result['pinion_torque_nm'],
        state_fingerprint=pressure_domain_fingerprint(result,method='full_gear_fe',tooth_offset=tooth_offset),
        flank_background_fingerprint=flank_background_fingerprint(background),passed=all(checks.values()),complete_working_flank=complete,
        unresolved_region_count=regions['unresolved_region_count'],evidence_kind='实际有限BRep点距＋1-Lipschitz三角片下界＋原节点候选FE三角并集覆盖；三角离散层级',
        checks=checks,regions=regions,fe_response_recomputed=response,
        original_candidate_triangle_count=len(facet_ids),triangulation_deflection_allowance_m=deflection,
        near_domain_geometry_engine=cover.geometry_engine,
        physical_tooth_indices=dict(pinion=pinion_index,wheel=wheel_index),opposing_fe_triangle_count=len(second_triangles),
        opposing_distance_definition='逐查询取实际有限CAD面与真实FE对面三角并集距离之小；不忽略两种几何的差异',
        loaded_opposing_rotation_rows=loaded_rotation.tolist(),loaded_elastic_motion_bound_m=loaded_motion_bound,
        rotation_linearization_remainder_m=linearization_remainder,
        limitations=['只认证所选实际物理齿/当前相位的输入三角背景，不认证其他齿或周期','近域采用原节点线性三角面离散，不是连续接触严格证明','峰压网格收敛和压力物理验证另行报告'])
    return _finalize_certificate(certificate,result,state_before=state_before,source_before=source_before)
