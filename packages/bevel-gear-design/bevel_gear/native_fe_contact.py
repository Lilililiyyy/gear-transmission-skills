"""真实有限元边界的单相位受载接触外迭代，SI 单位。

给定已有整轮 FE 准备对象及明确的投影格/平面/整体姿态，每步重新构造
双边界三角交片、P0 牵引 W、平均 g0 和力臂，调用原互补 QP。不会以
半空间柔度替代体柔度，不接管调用方因子，不把本次接受等同于整域或收敛。
"""
from __future__ import annotations

from copy import deepcopy
import math

MEMBERS = ('pinion', 'wheel')


def _positive(value, name, *, maximum=None, integer=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value <= 0 or integer and not isinstance(value, int) or
            maximum is not None and value > maximum):
        raise ValueError(name+'须为有效正值，并满足该入口的上限')
    return value


def _mesh_identity(mesh):
    from .distributed_pressure_domain import numeric_array_identity
    return dict(nodes_m=numeric_array_identity(mesh.nodes_m),
        tetrahedra=numeric_array_identity(mesh.tetrahedra, integer=True),
        boundary_triangles=numeric_array_identity(mesh.boundary_triangles, integer=True))


def _validate_compacts(prepared, compact, offsets, *, boundary_verified, progress):
    """核对所有 compact 原节点、完整真实边界及整数索引，不建立大邻接字典。"""
    import numpy as np
    from types import SimpleNamespace
    from .native_fe_preparation import verified_tetra_boundary
    from .native_pressure_domain import _body_boundary_keys
    from .distributed_pressure_domain import numeric_array_identity
    if set(compact) != set(MEMBERS): raise ValueError('compact_facets须完整包含pinion/wheel')
    normalized = {m: {0: compact[m]} if 'nodes_m' in compact[m] else compact[m] for m in MEMBERS}
    identities = {}; audit = {}
    for member in MEMBERS:
        mesh = prepared.systems[member].mesh; nodes = np.asarray(mesh.nodes_m)
        if not boundary_verified:
            source = mesh if hasattr(mesh, 'boundary_tags') else SimpleNamespace(nodes_m=mesh.nodes_m,
                tetrahedra=mesh.tetrahedra, boundary_triangles=mesh.boundary_triangles,
                boundary_tags=np.zeros(len(mesh.boundary_triangles), dtype=np.int64))
            checked = verified_tetra_boundary(source, progress=progress)
            del checked
        _, boundary, packed = _body_boundary_keys(mesh)
        selections = []
        if any(k not in normalized[member] for k in offsets):
            raise ValueError('compact_facets缺少真实物理齿：'+member)
        for offset in offsets:
            part = normalized[member][offset]
            ids = np.asarray(part['global_node_ids']); triangles = np.asarray(part['triangles'])
            points = np.asarray(part['nodes_m'])
            if (part.get('name') != member or part.get('total_node_count') != len(nodes) or
                    ids.ndim != 1 or ids.dtype.kind not in 'iu' or not len(ids) or
                    ids.min() < 0 or ids.max() >= len(nodes) or len(np.unique(ids)) != len(ids) or
                    points.shape != (len(ids), 3) or not np.array_equal(points, nodes[ids])):
                raise ValueError('compact节点/全局编号没有绑定当前完整FE网格：'+member+' k='+str(offset))
            if (triangles.ndim != 2 or triangles.shape[1:] != (3,) or triangles.dtype.kind not in 'iu' or
                    not len(triangles) or triangles.min() < 0 or triangles.max() >= len(ids)):
                raise ValueError('compact三角形须为当前局部节点内的真实整数索引')
            keys = packed(ids[triangles]); where = np.searchsorted(boundary, keys)
            if np.any(where == len(boundary)) or not np.array_equal(boundary[where], keys):
                raise ValueError('compact含有非原四面体真实外边界三角形')
            selections.append(dict(tooth_offset=offset, global_node_ids=numeric_array_identity(ids, integer=True),
                triangles=numeric_array_identity(triangles, integer=True)))
        identities[member] = _mesh_identity(mesh)
        audit[member] = dict(actual_mesh=identities[member], compact_selections=selections,
            original_global_nodes_bound=True, complete_tetrahedral_boundary_verified=True)
    return normalized, identities, audit


def _coverage(system):
    import numpy as np
    rows = []
    for patch in system['geometry']['patches']:
        valid = np.asarray(patch['valid_mask'], dtype=bool)
        source = np.asarray(patch['source_valid_mask'], dtype=bool)
        partial = np.asarray(patch['native_partial_cell_mask'], dtype=bool)
        if source.shape != valid.shape or partial.shape != valid.shape:
            raise ValueError('源请求、完整和partial掩码布局不同')
        rows.append(dict(tooth_offset=patch['tooth_offset'], requested=int(source.sum()),
            complete=int((valid & source).sum()), missing=int((source & ~valid).sum()),
            partial=int(partial.sum()), unexpected=int((valid & ~source).sum())))
    passed = all(not row['missing'] and not row['partial'] and not row['unexpected'] for row in rows)
    passed = passed and all(a is None or a['partial_cell_count'] == 0 for a in system['typed_patch_arrays'])
    return dict(passed=passed, patches=rows, requested_cell_set_not_shrunk=True)


def _attach_angles(core, system, initial_reference):
    header = system['geometry']; step = float(core['equilibrium']['theta_rad'])
    reference = float(header['reference_follower_angle_rad']); driver = float(header['driving_angle_rad'])
    ratio = float(header['pinion_to_wheel_ratio'])
    if not all(math.isfinite(value) for value in (step, reference, driver, ratio)):
        raise ArithmeticError('原QP增量或当前参考姿态不是有限rad')
    if 'elastic_te_rad' in core: core['legacy_reference_increment_te_rad'] = core.pop('elastic_te_rad')
    loaded = reference+step
    core.update(linearization_step_theta_rad=step, reference_follower_angle_rad=reference,
        driving_angle_rad=driver, loaded_follower_angle_rad=loaded, total_loaded_te_rad=loaded+ratio*driver,
        cumulative_reference_follower_shift_rad=reference-initial_reference,
        cumulative_loaded_follower_shift_rad=loaded-initial_reference,
        native_unloaded_mesh_tangency_solved=False,
        native_TE_definition='总受载TE=实际从动角+齿数比×主动角；末步delta仅相对当次参考姿态，未求卸载FE网格相切')


def solve_native_fe_contact(prepared, geometry_template, *, compact_facets=None,
        maximum_outer_iterations=6, outer_step_tolerance_rad=1e-8, maximum_working_set=2048,
        maximum_fe_ray_closure_m=1e-10, mean_geometry_tolerance_m=1e-10,
        maximum_piece_penetration_m=2e-7, maximum_pinion_torque_relative_residual=0.02,
        batch_cells=128, include_body_fields=False, progress=None):
    """求单相位直齿/弧齿锥齿轮真实FE边界的P0受载接触，不写文件。

    geometry_template为SI字典：frames给两体整体R0/世界轴，另含主动角、
    参考从动角、主动有符号扭矩、z1/z2及patches。每patch明确物理齿k、
    grid、plane和requested_mask；source_material_reference_frames可保留
    材料显示Rk。旧间隙/压力/柔度不作为输入；每次调用均冷启动首个QP。

    compact_facets可显式给真实表面节点/全局编号/三角形；None时由完整
    四面体边界及有限CAD工作面提取。内部始终校验compact确属于prepared。
    先检查全部原片顶点，再在实际加载phi重建全部请求格核验平均间隙。
    返回accepted/result/system/geometry_outer_history/geometry_check等具名
    字段。accepted不含完整材料域、网格/相位收敛或纯弹性TE承诺。
    """
    import numpy as np
    from . import native_fe_operator as native
    from .generated_distributed_fe import solve_distributed_contact_system
    from .distributed_pressure_domain import distributed_body_identity
    from .pressure_background import _digest
    from .native_fe_preparation import extract_native_flank_compacts
    _positive(maximum_outer_iterations, 'maximum_outer_iterations', maximum=6, integer=True)
    _positive(outer_step_tolerance_rad, 'outer_step_tolerance_rad', maximum=1e-8)
    _positive(maximum_working_set, 'maximum_working_set', integer=True)
    _positive(maximum_fe_ray_closure_m, 'maximum_fe_ray_closure_m', maximum=1e-10)
    _positive(mean_geometry_tolerance_m, 'mean_geometry_tolerance_m', maximum=1e-10)
    # 默认仍0.2μm；工程P0平均接触案例可显式用至1μm，仅扩展点间隙验收。
    # 平均约束、柔度/力/外迭代更新公式和其它几何/平衡门限均不改变。
    _positive(maximum_piece_penetration_m, 'maximum_piece_penetration_m', maximum=1e-6)
    _positive(maximum_pinion_torque_relative_residual, 'maximum_pinion_torque_relative_residual')
    _positive(batch_cells, 'batch_cells', integer=True)
    if not isinstance(include_body_fields, bool): raise ValueError('include_body_fields须为布尔值')
    if set(prepared.systems) != set(MEMBERS): raise ValueError('prepared须为pinion/wheel两体原FE系统')
    names = ('frames', 'patches', 'driving_angle_rad', 'reference_follower_angle_rad', 'pinion_torque_nm', 'pinion_to_wheel_ratio')
    for name in names:
        if name not in geometry_template: raise ValueError('geometry_template缺少显式字段：'+name)
    template = {name: deepcopy(geometry_template[name]) for name in names}
    # 明确只取几何定义，旧pressure/gaps/柔度等额外载荷字段不参与新算子。
    allowed = ('tooth_offset', 'grid', 'plane', 'requested_mask', 'source_material_reference_frames')
    template['patches'] = [{name: deepcopy(p[name]) for name in allowed if name in p} for p in template['patches']]
    if any(any(name not in p for name in ('tooth_offset', 'grid', 'plane', 'requested_mask')) for p in template['patches']):
        raise ValueError('每个patch须明确tooth_offset/grid/plane/requested_mask')
    for name in ('driving_angle_rad', 'reference_follower_angle_rad', 'pinion_torque_nm', 'pinion_to_wheel_ratio'):
        value = template[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(name+'须为有限SI实数')
    if template['pinion_torque_nm'] == 0 or template['pinion_to_wheel_ratio'] <= 0:
        raise ValueError('主动有符号扭矩不能为零，z1/z2须为正')
    try:
        z1=prepared.models['pinion'].request.tooth_count;z2=prepared.models['wheel'].request.tooth_count
    except (AttributeError,KeyError,TypeError) as error:
        raise ValueError('prepared.models须声明两轮实际tooth_count，供原QP与TE使用同一齿数比') from error
    if any(isinstance(z,bool) or not isinstance(z,int) or z<3 for z in (z1,z2)):
        raise ValueError('两轮实际齿数须为不少于3的整数')
    actual_ratio=z1/z2
    if not math.isclose(template['pinion_to_wheel_ratio'],actual_ratio,rel_tol=2e-14,abs_tol=0.):
        raise ValueError('geometry_template齿数比与prepared.models实际z1/z2不一致，不能计算TE')
    # 该量由实际齿数派生，原QP扭矩控制和TE簿记使用同一个浮点比值。
    template['pinion_to_wheel_ratio']=actual_ratio
    offsets = [patch['tooth_offset'] for patch in template['patches']]
    if not offsets or len(set(offsets)) != len(offsets) or any(isinstance(k, bool) or not isinstance(k, int) for k in offsets):
        raise ValueError('patch须明确唯一整数物理齿偏移')
    extraction = None
    if compact_facets is None:
        compact_facets, extraction = extract_native_flank_compacts(prepared, offsets, progress=progress)
    compact, mesh_identity, compact_audit = _validate_compacts(prepared, compact_facets, offsets,
        boundary_verified=extraction is not None, progress=progress)
    input_fingerprint = _digest(dict(geometry_template=template, actual_native_mesh_numeric_identity=mesh_identity,
        compact_identity=compact_audit, maximum_outer_iterations=maximum_outer_iterations,
        outer_step_tolerance_rad=outer_step_tolerance_rad, maximum_working_set=maximum_working_set,
        maximum_fe_ray_closure_m=maximum_fe_ray_closure_m, mean_geometry_tolerance_m=mean_geometry_tolerance_m,
        maximum_piece_penetration_m=maximum_piece_penetration_m,
        maximum_pinion_torque_relative_residual=maximum_pinion_torque_relative_residual))
    system = native.build_native_fe_contact_operator(compact, template['patches'], template['frames'],
        template['reference_follower_angle_rad'], driving_angle_rad=template['driving_angle_rad'],
        pinion_torque_nm=template['pinion_torque_nm'], pinion_to_wheel_ratio=template['pinion_to_wheel_ratio'], batch_cells=batch_cells)
    base = system; history = []; previous_keys = None; checks = None; core = None; reasons = []
    for outer in range(maximum_outer_iterations):
        core = None
        coverage = _coverage(system)
        if not coverage['passed']:
            reasons.append('原请求压力格出现partial/缺失，未缩小请求域或补零'); break
        keys = {(r['tooth_offset'], r['row'], r['column']): i for i, r in enumerate(system['records'])}
        if len(keys) != len(system['records']): raise ValueError('当前operator物理格列键重复')
        warm = None if previous_keys is None else [keys[key] for key in previous_keys if key in keys]
        core = solve_distributed_contact_system(prepared, system, template['pinion_torque_nm'],
            maximum_working_set=maximum_working_set, maximum_active_mapping_error_m=maximum_fe_ray_closure_m,
            allow_loaded_rigid_overlap=True, initial_indices=warm, include_body_fields=True, progress=progress)
        for old, new in (('maximum_cad_fe_mapping_error_m', 'maximum_native_FE_ray_shape_closure_error_m'),
                ('maximum_active_cad_fe_mapping_error_m', 'maximum_active_native_FE_ray_shape_closure_error_m')):
            if old in core: core[new] = core.pop(old)
        if core.get('status') == 'active_CAD_FE_mapping_requires_refinement':
            core['status'] = 'native_FE_ray_shape_closure_failed'
        core.update(geometry_source='FE_boundary', integration_method=system['integration_method'], geometry=system['geometry'],
            operator_input_fingerprint=system['operator_input_fingerprint'], source_input_fingerprint=input_fingerprint,
            actual_native_mesh_numeric_identity=mesh_identity, native_closure_is_not_CAD_error=True,
            mesh_convergence_verified=False, external_tooth_coverage_verified=False)
        if core.get('converged') is not True or core.get('equilibrium', {}).get('certificate', {}).get('passed') is not True:
            reasons.append('当前原P0互补QP未接受'); break
        _attach_angles(core, system, template['reference_follower_angle_rad'])
        step = core['linearization_step_theta_rad']
        history.append(dict(outer_iteration=outer, reference_follower_angle_rad=core['reference_follower_angle_rad'],
            delta_theta_rad=step, loaded_follower_angle_rad=core['loaded_follower_angle_rad'],
            total_loaded_te_rad=core['total_loaded_te_rad'], cumulative_reference_follower_shift_rad=core['cumulative_reference_follower_shift_rad'],
            cumulative_loaded_follower_shift_rad=core['cumulative_loaded_follower_shift_rad'],
            operator_input_fingerprint=system['operator_input_fingerprint'],
            warm_physical_key_count=0 if warm is None else len(warm), fresh_compliance_cache=True,
            old_force_reused=False, source_candidate_coverage=coverage))
        if progress: progress(dict(stage='native_FE_contact_outer', **history[-1]))
        if abs(step) <= outer_step_tolerance_rad:
            checks = native.check_native_operator_pose(system, compact,
                {m: core['body_fields'][m]['displacements_m'] for m in MEMBERS}, step,
                core['equilibrium']['remaining_gaps_m'], batch_cells=batch_cells,
                mean_tolerance_m=mean_geometry_tolerance_m, maximum_penetration_m=maximum_piece_penetration_m, progress=progress)
            if checks.get('passed') is not True: reasons.append('真实原片顶点或有限φ平均几何审核未通过')
            break
        if outer+1 == maximum_outer_iterations:
            reasons.append('达到最多外迭代次数，原增量仍超限'); break
        previous_keys = [tuple(key) for key in core['active_cell_keys']]
        if len(set(previous_keys)) != len(previous_keys): raise ValueError('原QP活动物理格重复')
        system = native.build_native_operator_at_pose(base, compact,
            native.advance_native_follower_frame(system['frames'], step), core['loaded_follower_angle_rad'],
            batch_cells=batch_cells, progress=progress)
    if core is not None:
        core['geometry_outer_history'] = deepcopy(history)
        core['geometry_check'] = checks if checks is not None else dict(passed=False, reason='尚未进入完整有限φ审核')
        core['geometry_linearization_checked'] = checks is not None and checks.get('passed') is True
        residual = core.get('pinion_torque_relative_residual')
        if isinstance(residual, bool) or not isinstance(residual, (int, float)) or not math.isfinite(residual) or abs(residual) > maximum_pinion_torque_relative_residual:
            reasons.append('重构主动轴力矩相对名义扭矩超出门限或非有限')
        boundary = core.get('contact_boundary_checks')
        if not isinstance(boundary, list) or not boundary or any(row.get('touches_window_edge') is not False or
                row.get('touches_invalid_candidate_neighbor') is not False for row in boundary):
            reasons.append('承载区触及候选窗或无效邻格，需扩展真实候选域')
        if 'body_fields' in core:
            core['distributed_body_identity'] = distributed_body_identity(core)
            if not include_body_fields: core.pop('body_fields')
    accepted = bool(core is not None and not reasons and checks is not None and checks.get('passed') is True)
    if core is not None: core['native_contact_accepted'] = accepted
    return dict(accepted=accepted, result=core, system=system, geometry_outer_history=history, geometry_check=checks,
        input_fingerprint=input_fingerprint, preparation_audit=dict(extraction=extraction, compact_binding=compact_audit),
        rejection_reasons=reasons, limitations=['当前有限FE边界/小变形P1迹模型，不是大滑动或连续CAD点态接触定理',
            '原P0平均KKT与独立片顶点工程门限分开；未证明完整材料域和外齿零压',
            '未证明压力空间/相位收敛，未求卸载FE网格相切，末步delta不是纯弹性TE',
            '参考phi为调用方明确的几何初值；冷启动指未继承压力/C，不代表任意CAD无载窗口的大步路径已验证',
            'frames与声明主动/从动角的运动学一致性由调用方明确负责；本入口不反求CAD运动学关系'])
