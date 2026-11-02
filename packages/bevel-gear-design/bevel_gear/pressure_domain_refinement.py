"""只细分压力覆盖证书的未决片；保留原几何、KKT和完整分割证据。

本模块不重求接触、不改原压力、不扩张已知域。父证书先绑定同一原状态、
实际背景及审核源码；原已确认片保持原样，只有新距离界或原求解域覆盖
检查通过的真实子片可以从未决转为已知。
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import math
from pathlib import Path


def _source_sha256():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def refine_triangle_pressure_domain(parent_regions, *, distance_oracle,
        solved_region_oracle, relative_motion_bound, extra_depth=2,
        maximum_regions=300000, tolerance_m=1e-10, progress=None,
        targeted_region_indices=None):
    """仅在原未决片内继续四分；累计深度和原CAD三角归属保持不变。

    原证书的完整原三角分割须由高层先验证。本低层接口用于数值回归和
    显式几何适配；不会把调用者提供的布尔标志当作完整压力域证据。
    """
    import numpy as np
    from .pressure_domain_audit import audit_triangle_pressure_domain
    if isinstance(extra_depth, bool) or not isinstance(extra_depth, int) or not 1 <= extra_depth <= 12:
        raise ValueError('附加细分深度须为1至12的整数')
    if not isinstance(parent_regions, dict) or not isinstance(parent_regions.get('regions'), list):
        raise ValueError('父覆盖证书须保留真实分片')
    original = parent_regions['regions']
    if targeted_region_indices is None:
        selected = {i for i,r in enumerate(original) if r.get('state') == 'unresolved'}
    else:
        selected = set(targeted_region_indices)
        if (not selected or len(selected) != len(targeted_region_indices) or
                any(isinstance(i,bool) or not isinstance(i,int) or not 0 <= i < len(original) or original[i].get('state') != 'unresolved' for i in selected)):
            raise ValueError('定向索引须为父证书中不重复的真实未决片行号')
    selected_indices = sorted(selected)
    undecided = [original[i] for i in selected_indices]
    if not undecided:
        raise ValueError('父证书没有未决片，无需细分')
    allowed = ('original_solved_domain', 'excluded_by_distance_bound', 'unresolved')
    if any(r.get('state') not in allowed for r in original):
        raise ValueError('父证书含未定义的覆盖片状态')
    depth = parent_regions.get('maximum_depth')
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 0 or depth+extra_depth > 30:
        raise ValueError('累计细分深度须不大于30')
    source = np.asarray([r['vertices_world_m'] for r in undecided], dtype=float)
    result = audit_triangle_pressure_domain(source, distance_oracle=distance_oracle,
        solved_region_oracle=solved_region_oracle, relative_motion_bound=relative_motion_bound,
        chord_error_m=parent_regions['chord_error_allowance_m'], maximum_depth=extra_depth,
        maximum_regions=maximum_regions, tolerance_m=tolerance_m, progress=progress)
    rows = [deepcopy(r) for i,r in enumerate(original) if i not in selected]
    for piece in result['regions']:
        parent = undecided[piece['source_triangle_index']]
        row = deepcopy(piece)
        row['source_triangle_index'] = parent['source_triangle_index']
        row['depth'] = parent['depth']+piece['depth']
        rows.append(row)
    areas = {state: math.fsum(r['area_m2'] for r in rows if r['state'] == state) for state in allowed}
    before = math.fsum(r['area_m2'] for r in original)
    if not math.isclose(math.fsum(areas.values()), before, rel_tol=1e-10, abs_tol=1e-23):
        raise ArithmeticError('细分没有守恒原覆盖面积')
    unresolved = [r for r in rows if r['state'] == 'unresolved']
    margins = [r['distance_lower_bound_m']-r['relative_motion_upper_bound_m']
        for r in rows if r['state'] == 'excluded_by_distance_bound']
    return dict(passed=not unresolved, unresolved_region_count=len(unresolved),
        input_triangle_count=parent_regions['input_triangle_count'], terminal_region_count=len(rows),
        evaluated_region_count=parent_regions['evaluated_region_count']+result['evaluated_region_count'],
        area_m2=areas, unresolved_source_triangles=sorted({r['source_triangle_index'] for r in unresolved}),
        minimum_excluded_clearance_margin_m=min(margins, default=None),
        chord_error_allowance_m=parent_regions['chord_error_allowance_m'], maximum_depth=depth+extra_depth,
        regions=rows, evidence_scope=parent_regions.get('evidence_scope'),
        refinement_evaluated_region_count=result['evaluated_region_count'],
        refined_parent_region_count=len(undecided), targeted_parent_region_indices=selected_indices)


def _fe_oracles(prepared, result, certificate, offset, matching_tolerance_m, tolerance_m):
    """重建同一实际FE候选并集和加载对面，复用已绑定原响应的位移界。"""
    import numpy as np
    from .generated_ltca import select_generated_flank_facets
    from .generated_halfspace import _faces
    from .pressure_domain_audit import (FiniteFaceDistanceOracle, TriangleDistanceOracle,
        _loaded_wheel_rotation, _triangle_union_cover)
    frames = result['mesh_reference_frames']
    rotation = np.asarray(frames['pinion_rotation_rows'])
    z1 = prepared.models['pinion'].request.tooth_count
    z2 = prepared.models['wheel'].request.tooth_count
    # 原背景片已按实际齿变换；这里只恢复对应有限面及原候选三角形。
    pinion_index = (-offset) % z1; wheel_index = offset % z2
    kwargs = dict(matching_tolerance_m=matching_tolerance_m)
    first = select_generated_flank_facets(prepared.systems['pinion'].mesh,
        prepared.models['pinion'].shell.surfaces['left_flank'],
        tooth_count=z1, tooth_index=pinion_index, **kwargs)
    candidates = {r['first_pressure_node_index'] for r in result['nodal_contacts'] if r['tooth_offset'] == offset}
    facet_ids = [tri for tri in first['triangles'] if all(int(i) in candidates for i in tri)]
    if not facet_ids:
        raise ValueError('所选齿没有可复用的完整原FE候选三角域')
    triangles = np.asarray(first['nodes_m'])[facet_ids] @ rotation.T
    normals = np.asarray([r['normal_world'] for r in result['nodal_contacts'] if r['tooth_offset'] == offset])
    normal = normals.mean(axis=0); normal /= np.linalg.norm(normal)
    tangent = np.cross(normal, [1., 0., 0.] if abs(normal[0]) < .9 else [0., 1., 0.])
    tangent /= np.linalg.norm(tangent)
    basis = np.array([tangent, np.cross(normal, tangent)])
    cover = _triangle_union_cover(triangles, triangles[0, 0], basis,
        padding_m=certificate['triangulation_deflection_allowance_m'])
    loaded = _loaded_wheel_rotation(frames, result['elastic_te_rad'])
    if not np.allclose(loaded, certificate['loaded_opposing_rotation_rows'], rtol=0., atol=1e-14):
        raise ValueError('父证书加载位姿与当前原状态不同')
    wheel = prepared.models['wheel']
    cad = FiniteFaceDistanceOracle(_faces(wheel.shape, wheel.slot_face_groups[wheel_index]['right_flank']),
        loaded.tolist(), tolerance_m=tolerance_m)
    second = select_generated_flank_facets(prepared.systems['wheel'].mesh, wheel.shell.surfaces['right_flank'],
        tooth_count=z2, tooth_index=wheel_index, **kwargs)
    fe = TriangleDistanceOracle(np.asarray(second['nodes_m'])[second['triangles']] @ loaded.T,
        tolerance_m=tolerance_m)
    def distance(point):
        a, b = cad(point), fe(point)
        done = bool(a['done'] and b['done'])
        return dict(done=done, distance_m=min(a['distance_m'], b['distance_m']) if done else None,
            error_bound_m=max(a['error_bound_m'], b['error_bound_m']))
    response = certificate['fe_response_recomputed']
    remainder = .5*response['wheel_radius_bound_m']*result['elastic_te_rad']**2
    motion = sum(response['maximum_displacement_m'].values())+remainder
    if not math.isclose(motion, certificate['loaded_elastic_motion_bound_m'], rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError('父证书FE运动界与原响应不一致')
    return distance, cover, lambda triangle: motion


def _halfspace_oracles(prepared, result, certificate, offset, tolerance_m, *,
        geometry_engine='brep', opposing_deflection_m=1e-6):
    """从同一原矩形力源及有限对面重建半解析的距离/弹性上界。"""
    import numpy as np
    from .pressure_domain_audit import FiniteFaceDistanceOracle, TriangleDistanceOracle, _loaded_wheel_rotation, _grid_cover, _halfspace_motion
    patch = next(p for p in result['geometry']['patches'] if p['tooth_offset'] == offset)
    theta = float(result.get('linearization_step_theta_rad', result['theta_rad']))
    loaded = _loaded_wheel_rotation(patch['frames'], theta)
    if not np.allclose(loaded, certificate['loaded_opposing_rotation_rows'], rtol=0., atol=1e-14):
        raise ValueError('父证书加载位姿与原半解析状态不同')
    if geometry_engine == 'brep':
        distance = FiniteFaceDistanceOracle(prepared.trimmed['wheel'].faces, loaded.tolist(), tolerance_m=tolerance_m)
        oracle_evidence = dict(engine='OCP_finite_BRep',opposing_deflection_m=0.)
    elif geometry_engine == 'triangles':
        from OCP.BRep import BRep_Builder
        from OCP.TopoDS import TopoDS_Compound
        from .cad_import import ImportedCadModel, triangulate_imported_model
        from .pressure_background import _digest
        cache = getattr(prepared,'_pressure_refinement_distance_cache',None)
        if cache is None:
            cache = {}; setattr(prepared,'_pressure_refinement_distance_cache',cache)
        key = (opposing_deflection_m,tolerance_m)
        if key not in cache:
            compound = TopoDS_Compound(); builder = BRep_Builder(); builder.MakeCompound(compound)
            for face in prepared.trimmed['wheel'].faces:builder.Add(compound,face)
            mesh = triangulate_imported_model(ImportedCadModel(compound,{}),
                linear_deflection_m=opposing_deflection_m,angular_deflection_rad=.08)
            local = TriangleDistanceOracle(np.asarray(mesh['vertices_m'])[np.asarray(mesh['triangles'],dtype=int)],
                tolerance_m=tolerance_m+opposing_deflection_m,backend='vtk')
            cache[key]=(local,dict(engine='actual_finite_CAD_triangles_VTK',opposing_deflection_m=opposing_deflection_m,
                mesh_fingerprint=_digest(mesh),scope='原有限面对面三角距显式扣弦差；不作连续CAD严格区间证明'))
        local,oracle_evidence = cache[key]
        distance = lambda point: local(np.asarray(point) @ loaded)
    else:
        raise ValueError('细分几何引擎须为brep或triangles')
    if offset in result['selected_tooth_offsets']:
        field = result['fields'][result['selected_tooth_offsets'].index(offset)]
        motion, _ = _halfspace_motion(result, offset, include_rigid_motion=False)
        cover = _grid_cover(patch,field,padding_m=certificate['geometry_binding']['triangulation_deflection_allowance_m'])
    else:
        motion = lambda triangle:0.
        cover = lambda triangle:False
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.TopAbs import TopAbs_OUT
        from OCP.gp import gp_Pnt
        classifier = BRepClass3d_SolidClassifier(prepared.tca.cad_models['wheel'].shape)
        unsigned = distance
        def distance(point):
            classifier.Perform(gp_Pnt(*(np.asarray(point) @ loaded*1000.)),tolerance_m*1000.)
            if classifier.State() != TopAbs_OUT:return dict(done=False,distance_m=None,error_bound_m=tolerance_m)
            return unsigned(point)
    return distance, cover, motion, oracle_evidence


def refine_generated_pressure_domain(prepared, result, flank_mesh, pinion_request,
        parent_certificate, *, extra_depth=2, maximum_regions=300000,
        matching_tolerance_m=8e-7, tolerance_m=1e-10, progress=None,
        targeted_region_indices=None, geometry_engine='brep', opposing_deflection_m=1e-6):
    """深化已绑定原状态的部分覆盖证书；原压力及父证书均保持原样。

    prepared必须来自同制造输入和同FE网格。FE响应/KKT不重复求解，
    只复用父证书中已重算且绑定原状态的上界；每个新排除仍重新执行同一
    几何距离和运动界公式。输出保存父证书摘要及本模块前后源码摘要。
    """
    from .pressure_background import (_digest, pressure_domain_background, pressure_domain_fingerprint)
    from .pressure_domain_audit import _kernel_sha256, _finalize_certificate
    from .refined_fe_plotting import _background
    from .generated_surface import _request
    from .generated_ltca import PreparedGeneratedLtca
    from .generated_halfspace import PreparedGeneratedHalfspace
    if not isinstance(parent_certificate, dict) or parent_certificate.get('passed') is not False:
        raise ValueError('需要真实未决的部分覆盖父证书')
    method = parent_certificate.get('method'); offset = parent_certificate.get('tooth_offset')
    if isinstance(opposing_deflection_m,bool) or not isinstance(opposing_deflection_m,(int,float)) or not math.isfinite(opposing_deflection_m) or opposing_deflection_m<=0:
        raise ValueError('对面弦差须为正有限m')
    request = _request(pinion_request)
    if method == 'full_gear_fe':
        if not isinstance(prepared, PreparedGeneratedLtca) or prepared.models['pinion'].request != request:
            raise ValueError('FE细分准备对象与实际制造输入不一致')
    elif method == 'local_halfspace':
        if not isinstance(prepared, PreparedGeneratedHalfspace) or prepared.tca.requests['pinion'] != request:
            raise ValueError('半解析细分准备对象与实际制造输入不一致')
    else:
        raise ValueError('父覆盖证书方法不支持')
    source = _source_sha256(); audit_source = _kernel_sha256()
    parent_before = _digest(parent_certificate)
    state_before = pressure_domain_fingerprint(result, method=method, tooth_offset=offset)
    if parent_certificate.get('audit_identity', {}).get('source_sha256_after') != audit_source:
        raise ValueError('父审核源码身份与当前内核不一致，不能复用父证明')
    _, background, _ = _background(flank_mesh, pinion_request)
    info = pressure_domain_background(result, background, method=method, tooth_offset=offset,
        certificate=parent_certificate, pinion_request=pinion_request)
    if info.get('domain_state') != 'partitioned_pressure_domain':
        raise ValueError('父覆盖证书没有通过完整原三角分割审核')
    if method == 'full_gear_fe':
        if geometry_engine!='brep':raise ValueError('FE定向细分当前使用原BRep对照路线')
        distance,cover,motion=_fe_oracles(prepared,result,parent_certificate,offset,matching_tolerance_m,tolerance_m)
        oracle_evidence=dict(engine='OCP_finite_BRep_plus_actual_FE_triangles')
    else:
        distance,cover,motion,oracle_evidence=_halfspace_oracles(prepared,result,parent_certificate,offset,tolerance_m,
            geometry_engine=geometry_engine,opposing_deflection_m=opposing_deflection_m)
    regions = refine_triangle_pressure_domain(parent_certificate['regions'], distance_oracle=distance,
        solved_region_oracle=cover, relative_motion_bound=motion, extra_depth=extra_depth,
        maximum_regions=maximum_regions, tolerance_m=tolerance_m, progress=progress,
        targeted_region_indices=targeted_region_indices)
    certificate = deepcopy(parent_certificate)
    certificate['regions'] = regions
    certificate['unresolved_region_count'] = regions['unresolved_region_count']
    certificate['complete_working_flank'] = regions['passed']
    keys = ('omitted_candidates_excluded', 'unresolved_regions_excluded', 'far_field_motion_bound_passed') if method == 'full_gear_fe' else (
        'full_working_flank_covered', 'omitted_regions_excluded', 'unresolved_regions_excluded')
    for key in keys:
        certificate['checks'][key] = regions['passed']
    certificate['passed'] = all(certificate['checks'].values())
    parent_after = _digest(parent_certificate)
    certificate['refinement_identity'] = dict(parent_certificate_fingerprint=parent_before,
        parent_unchanged=parent_before == parent_after, source_sha256_before=source,
        source_sha256_after=_source_sha256(), state_fingerprint_before=state_before,
        state_fingerprint_after=pressure_domain_fingerprint(result, method=method, tooth_offset=offset),
        reused_parent_geometry_binding=True, reused_parent_kkt_evidence=True,
        parent_source_identity_bound=True, additional_depth=extra_depth,
        narrowed_only_parent_unresolved_regions=True,
        targeted_parent_region_indices=regions['targeted_parent_region_indices'], geometry_oracle=oracle_evidence)
    identity = certificate['refinement_identity']
    if (not identity['parent_unchanged'] or identity['source_sha256_before'] != identity['source_sha256_after'] or
            identity['state_fingerprint_before'] != identity['state_fingerprint_after']):
        raise ArithmeticError('细分期间父证书、原状态或细分源码改变')
    certificate = _finalize_certificate(certificate, result, state_before=state_before, source_before=audit_source)
    # 再审核新累计四分路径，防止丢片、重复、深度映射或面积分类错误。
    pressure_domain_background(result, background, method=method, tooth_offset=offset,
        certificate=certificate, pinion_request=pinion_request)
    return certificate
