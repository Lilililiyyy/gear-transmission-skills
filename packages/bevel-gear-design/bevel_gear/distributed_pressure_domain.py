"""新写P0分布牵引FE域审核：弱平均压力自由度与远场材料证据分开。

本模块绑定真实CAD位姿、全部P0载荷、同次FE位移/节点力及体网格。
已解小片的零P0值不等于逐点非穿透；Gauss残差和网格收敛单独保存。
不创建nodal_contacts，也不重新分解弹性矩阵。
"""
from __future__ import annotations

import hashlib
import json
import math


P0_REPRESENTATION = 'distributed_P0_surface_cells'
P0_CHECKS = ('input_geometry_bound', 'included_P0_average_kkt_passed',
    'body_response_bound', 'actual_loaded_pose_bound', 'full_working_flank_covered',
    'omitted_regions_excluded', 'unresolved_regions_excluded', 'far_field_motion_bound_passed')


def is_distributed_pressure_result(result):
    return isinstance(result, dict) and result.get('contact_representation') == P0_REPRESENTATION


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name+'须为有限实数')
    return float(value)


def numeric_array_identity(values, *, integer=False):
    """用固定小端dtype和shape摘要绑定数组；NumPy缺席时小型身份测试仍可运行。"""
    dtype = '<i8' if integer else '<f8'
    try:
        import numpy as np
    except ImportError:
        import struct
        def flatten(value):
            if isinstance(value, (list, tuple)):
                children = [flatten(item) for item in value]
                if children and any(shape != children[0][0] for shape, _ in children):
                    raise ValueError('身份数组不能为不规则嵌套')
                shape = (len(value),)+(children[0][0] if children else ())
                return shape, b''.join(raw for _, raw in children)
            number = _finite(value, '身份数组值')
            if integer and number != int(number):
                raise ValueError('整数身份数组不能截断浮点索引')
            return (), struct.pack('<q' if integer else '<d', int(number) if integer else number)
        shape, raw = flatten(values)
    else:
        original = np.asarray(values)
        if original.dtype.kind not in ('i', 'u', 'f') or not np.isfinite(original).all():
            raise ValueError('身份数组须为有限数值')
        if integer and np.any(original != original.astype(np.int64)):
            raise ValueError('整数身份数组不能截断浮点索引')
        array = np.asarray(original, dtype=dtype, order='C')
        shape, raw = array.shape, array.tobytes(order='C')
    digest = hashlib.sha256(json.dumps(dict(dtype=dtype, shape=list(shape)), separators=(',', ':')).encode('ascii')+b'\0'+raw).hexdigest()
    return dict(dtype=dtype, shape=list(shape), sha256=digest)


def distributed_body_identity(result):
    """原数组与外存NPZ身份使用同协议；两者同时提供时必须精确一致。"""
    from .pressure_background import _digest
    bodies = result.get('body_fields')
    supplied = result.get('distributed_body_identity')
    if bodies is None:
        if not isinstance(supplied, dict) or set(supplied) != {'pinion', 'wheel'}:
            raise ValueError('P0域审核须绑定同次u/F原数组或其明确外存身份')
        identity = supplied
    else:
        if not isinstance(bodies, dict) or set(bodies) != {'pinion', 'wheel'}:
            raise ValueError('P0全轮响应必须包含pinion/wheel')
        identity = {}
        for member, response in bodies.items():
            displacement = numeric_array_identity(response['displacements_m'])
            force = numeric_array_identity(response['nodal_forces_n'])
            if displacement['shape'] != force['shape'] or len(force['shape']) != 2 or force['shape'][1] != 3 or not force['shape'][0]:
                raise ValueError('P0体响应须为同一完整非空N×3的u(m)/F(N)')
            try:
                import numpy as np
            except ImportError:
                maximum_u = max(math.sqrt(sum(float(value)**2 for value in row)) for row in response['displacements_m'])
            else:
                maximum_u = float(np.linalg.norm(np.asarray(response['displacements_m'], dtype=float), axis=1).max())
            identity[member] = dict(displacement_m=displacement, nodal_force_n=force, maximum_displacement_norm_m=maximum_u,
                response_diagnostics={name: response[name] for name in ('strain_energy_j', 'external_work_j',
                    'maximum_free_residual_n', 'relative_free_residual')})
        if supplied is not None and _digest(identity) != _digest(supplied):
            raise ValueError('外存P0响应身份与当前实际u/F数组不一致')
    for response in identity.values():
        if _finite(response.get('maximum_displacement_norm_m'), '实际位移范数上界m') < 0:
            raise ValueError('实际位移范数上界不能为负')
        for name in ('displacement_m', 'nodal_force_n'):
            item = response.get(name, {})
            if (item.get('dtype') != '<f8' or len(item.get('shape', [])) != 2 or item['shape'][1] != 3
                    or not isinstance(item.get('sha256'), str) or len(item['sha256']) != 64):
                raise ValueError('P0体响应外存身份缺少明确N×3小端数组和SHA')
    return identity


def distributed_loaded_pose(result):
    """最终phi只取最后参考姿态加最后增量，total TE只作为独立历史量绑定。"""
    from .native_pressure_domain import native_geometry_mode
    geometry = result['geometry']
    reference = _finite(geometry['reference_follower_angle_rad'], '最后参考从动角rad')
    step = _finite(result['linearization_step_theta_rad'], '最后从动角增量rad')
    loaded = _finite(result['loaded_follower_angle_rad'], '实际受载从动角rad')
    if not math.isclose(loaded, reference+step, rel_tol=0., abs_tol=2e-14):
        raise ValueError('P0实际从动角必须为最后reference+step，不能重复叠加total TE')
    if native_geometry_mode(result):
        # 原生分支尚可未求卸载FE网格相切；增量角不能冒充纯弹性或总TE。
        equilibrium_step = _finite(result['equilibrium']['theta_rad'], '原QP从动角增量rad')
        if not math.isclose(step, equilibrium_step, rel_tol=0., abs_tol=2e-14):
            raise ValueError('原生FE位姿增量须与原QP theta逐项对应')
        driving = _finite(geometry['driving_angle_rad'], '主动角rad')
        ratio = _finite(geometry['pinion_to_wheel_ratio'], '小轮/大轮传动比')
        total = _finite(result['total_loaded_te_rad'], '实际总受载TE rad')
        if not math.isclose(total, loaded+ratio*driving, rel_tol=0., abs_tol=2e-14):
            raise ValueError('原生总受载TE须由实际从动角与理想传动角计算')
        return dict(driving_angle_rad=driving, pinion_torque_nm=_finite(geometry['pinion_torque_nm'], '主动扭矩Nm'),
            last_reference_follower_angle_rad=reference, last_step_theta_rad=step,
            loaded_follower_angle_rad=loaded, total_loaded_te_rad=total,
            increment_definition='相对于当次参考刚体角的QP增量，不是纯弹性TE',
            native_unloaded_mesh_tangency_solved=result.get('native_unloaded_mesh_tangency_solved') is True)
    return dict(driving_angle_rad=_finite(geometry['driving_angle_rad'], '主动角rad'),
        pinion_torque_nm=_finite(geometry['pinion_torque_nm'], '主动扭矩Nm'),
        last_reference_follower_angle_rad=reference, last_step_theta_rad=step,
        loaded_follower_angle_rad=loaded, total_elastic_te_rad=_finite(result['elastic_te_rad'], '总TE rad'))


def distributed_response_fingerprint(result):
    """绑定全部真实P0压力、几何、全轮u/F与平均约束；不借用节点接触协议。"""
    from .pressure_background import _digest
    from .native_pressure_domain import native_result_identity
    if not is_distributed_pressure_result(result) or result.get('converged') is not True:
        raise ValueError('需要已接受的真实P0分布牵引结果')
    if result.get('pressure_filter_applied') is not False or result.get('halfspace_compliance_used') is not False:
        raise ValueError('P0 FE覆盖身份只接受未滤波的真实FE结构柔度结果')
    requested = result.get('requested_tooth_offsets')
    patches = result['geometry'].get('patches', [])
    if (not isinstance(requested, (list, tuple)) or not requested or len(set(requested)) != len(requested)
            or any(isinstance(k, bool) or not isinstance(k, int) for k in requested)
            or sorted(requested) != sorted(p['tooth_offset'] for p in patches)):
        raise ValueError('P0请求窗口必须逐齿绑定唯一真实几何patch')
    selected, fields = result.get('selected_tooth_offsets', []), result.get('fields', [])
    if (not selected or len(set(selected)) != len(selected) or len(selected) != len(fields)
            or any(k not in requested or f.get('tooth_offset') != k for k, f in zip(selected, fields))):
        raise ValueError('P0有效压力场须与selected齿索引逐项对应')
    payload = dict(method='full_gear_fe', representation=P0_REPRESENTATION,
        pose=distributed_loaded_pose(result), requested_tooth_offsets=requested,
        selected_tooth_offsets=selected, geometry=result['geometry'], fields=fields,
        equilibrium=result['equilibrium'], cells=result.get('cells', []),
        quadrature_order=result['quadrature_order'], active_cell_keys=result.get('active_cell_keys', []),
        body_identity=distributed_body_identity(result), mesh_identity=result.get('distributed_mesh_identity'),
        body_response_checks=result.get('body_response_checks'), geometry_check=result.get('geometry_check'),
        geometry_linearization_checked=result.get('geometry_linearization_checked'),
        geometry_outer_history=result.get('geometry_outer_history'))
    # 只在显式协议字段存在时扩展身份，缺省旧NURBS payload/SHA保持原样。
    explicit = native_result_identity(result)
    if explicit:payload['explicit_geometry_protocol'] = explicit
    return _digest(payload)


def distributed_pressure_fingerprint(result, *, tooth_offset=0):
    """空selected齿仍绑定完整同次P0响应；请求窗口之外的齿明确拒绝。"""
    from .pressure_background import _digest
    if isinstance(tooth_offset, bool) or not isinstance(tooth_offset, int) or tooth_offset not in result.get('requested_tooth_offsets', []):
        raise ValueError('P0物理齿必须属于实际请求窗口，缺场不构成零压证明')
    return _digest(dict(method='full_gear_fe', representation=P0_REPRESENTATION,
        tooth_offset=tooth_offset, full_response_fingerprint=distributed_response_fingerprint(result)))


def distributed_average_kkt_evidence(result):
    """重核原每片合力/平均间隙标量KKT；逐Gauss和连续点约束分别报告。"""
    equilibrium = result['equilibrium']; certificate = equilibrium['certificate']
    forces, gaps = equilibrium['forces_n'], equilibrium['remaining_gaps_m']
    if not len(forces) or len(forces) != len(gaps):
        raise ValueError('P0原标量合力和平均间隙长度必须一致且非空')
    force_tolerance = _finite(certificate['force_acceptance_n'], '力容差N')
    gap_tolerance = _finite(certificate['gap_acceptance_m'], '平均间隙容差m')
    complementarity = sum(abs(float(f)*float(g)) for f, g in zip(forces, gaps))
    active_stationarity = max((abs(float(g)) for f, g in zip(forces, gaps) if f > force_tolerance), default=math.inf)
    failures = []
    if min(forces) < -force_tolerance: failures.append('原P0合力违反非负约束')
    if min(gaps) < -gap_tolerance: failures.append('原P0积分平均间隙违反非负约束')
    if active_stationarity > gap_tolerance: failures.append('原P0承载小片平均驻点误差超限')
    if complementarity > certificate['complementarity_acceptance_nm']: failures.append('原P0平均互补误差超限')
    if certificate.get('passed') is not True: failures.append('原P0离散证书未通过')
    return dict(passed=not failures, failures=failures, force_count=len(forces),
        minimum_force_n=float(min(forces)), minimum_average_gap_m=float(min(gaps)),
        active_average_stationarity_m=active_stationarity, complementarity_nm=complementarity,
        original_gauss_minimum_gap_m=result.get('minimum_gauss_remaining_gap_m'),
        actual_loaded_gauss_minimum_gap_m=result.get('geometry_check', {}).get('actual_loaded_minimum_gauss_gap_m'),
        pointwise_nonpenetration_certified=False, pressure_mesh_convergence_certified=False,
        definition='原P0小片合力和积分平均间隙的标量KKT；不冒充逐节点或逐Gauss非穿透')


def matrix_free_body_response(mesh, material, displacements_m, nodal_forces_n, fixed_nodes, *, chunk_size=20000):
    """分块回装BᵀDBu验证Ku−F，不构造全局K、不因子分解、不重新求u。"""
    import numpy as np
    from .elasticity import _constitutive, _element_kinematics
    nodes, tetra = np.asarray(mesh.nodes_m, dtype=float), np.asarray(mesh.tetrahedra, dtype=np.int64)
    u, force = np.asarray(displacements_m, dtype=float), np.asarray(nodal_forces_n, dtype=float)
    if u.shape != nodes.shape or force.shape != nodes.shape or not np.isfinite(u).all() or not np.isfinite(force).all():
        raise ValueError('独立Ku−F重建必须使用完整同网格u/F')
    _, constitutive = _constitutive(material, np)
    internal = np.zeros_like(nodes)
    for start in range(0, len(tetra), chunk_size):
        indices = tetra[start:start+chunk_size]
        b, volumes = _element_kinematics(nodes[indices], np)
        element_u = u[indices].reshape((-1, 12))
        strain = np.einsum('eai,ei->ea', b, element_u)
        stress = strain@constitutive.T
        element_force = np.einsum('eai,ea,e->ei', b, stress, volumes).reshape((-1, 4, 3))
        np.add.at(internal, indices.ravel(), element_force.reshape((-1, 3)))
    fixed = np.asarray(fixed_nodes, dtype=np.int64)
    free = np.ones(len(nodes), dtype=bool); free[fixed] = False
    residual = internal-force
    relative = float(np.linalg.norm(residual[free])/max(float(np.linalg.norm(force[free])), 1.))
    fixed_max = float(np.abs(u[fixed]).max()) if len(fixed) else math.inf
    return dict(passed=bool(relative <= 1e-8 and fixed_max <= 1e-14),
        relative_free_residual=relative, maximum_free_residual_n=float(np.abs(residual[free]).max()),
        maximum_fixed_displacement_m=fixed_max, maximum_displacement_norm_m=float(np.linalg.norm(u, axis=1).max()),
        strain_energy_j=.5*float(np.sum(u*internal)), external_half_work_j=.5*float(np.sum(u*force)),
        displacement_identity=numeric_array_identity(u), nodal_force_identity=numeric_array_identity(force),
        mesh_nodes_identity=numeric_array_identity(nodes), mesh_tetrahedra_identity=numeric_array_identity(tetra, integer=True),
        matrix_assembled=False, factorization_performed=False,
        definition='原CST四面体分块应变/应力/节点内力回装，固定孔壁自由度单独排除')


def restore_distributed_geometry(geometry):
    """将JSON几何恢复为原系统重建所需grid/角点数组，不修改原输入对象。"""
    import numpy as np
    from .halfspace_contact import RectangularContactGrid
    restored = dict(geometry)
    restored['patches'] = []
    for original in geometry['patches']:
        patch = dict(original)
        if isinstance(patch['grid'], dict):
            patch['grid'] = RectangularContactGrid(**patch['grid'])
        for name in ('first_corner_points_world_m', 'second_corner_points_world_m'):
            patch[name] = np.asarray(patch[name], dtype=float)
        restored['patches'].append(patch)
    return restored


def audit_distributed_response(prepared, system, result, *, body_response_evidence=None):
    """用真实W重算Wf、平均弹性开口和合力矩，并绑定同次u/网格/Ku残差。

生产回调直接提供原prepared/system。离线重建只有几何/网格时，须明确
提供分块Ku−F证据，且u/F和网格数组身份逐项吻合；两条路径均不求新u。
"""
    import numpy as np
    from .pressure_background import _digest
    from .native_pressure_domain import native_geometry_mode, audit_native_partition_geometry
    records, matrices = system['records'], system['matrices']
    native_proof = None; native_levers = None; native_gaps = None
    if native_geometry_mode(result):
        native_gaps, native_keys, native_proof = audit_native_partition_geometry(system,
            {m: prepared.systems[m].mesh for m in ('pinion', 'wheel')}, result=result, include_recomputed_values=True)
        if native_proof.get('passed') is not True or native_proof.get('W_reconstructed_from_original_pieces') is not True:
            raise ArithmeticError('原生FE分片平均g0/W/刚体力臂独立重建未通过')
        native_levers = native_proof.pop('recomputed_lever_arms_m')
        if not np.array_equal(native_keys, np.asarray([(r['tooth_offset'], r['row'], r['column']) for r in records])):
            raise ValueError('原生FE独立分片格序与原W列不一致')
    forces = np.asarray(result['equilibrium']['forces_n'], dtype=float)
    if len(records) != len(forces) or set(matrices) != {'pinion', 'wheel'}:
        raise ValueError('实际W列、P0力向量和两体必须完整对应')
    force_from_fields = np.concatenate([np.asarray(field['forces_n'])[np.asarray(field['valid_mask'], dtype=bool)] for field in result['fields']])
    if not np.array_equal(force_from_fields, forces):
        raise ValueError('原P0二维合力与实际W列排序不一致')
    response = distributed_body_identity(result)
    checks, openings, matrix_identity = {}, [], {}; native_energies = []
    for member, matrix in matrices.items():
        body, saved = prepared.systems[member], result['body_fields'][member]
        u = np.asarray(saved['displacements_m'], dtype=float).reshape(-1)
        stored_force = np.asarray(saved['nodal_forces_n'], dtype=float).reshape(-1)
        reconstructed = np.asarray(matrix@forces).reshape(-1)
        error = float(np.linalg.norm(reconstructed-stored_force)/max(float(np.linalg.norm(stored_force)), 1.))
        openings.append(np.asarray(matrix.T@u).reshape(-1))
        matrix = matrix.copy(); matrix.sum_duplicates(); matrix.sort_indices()
        matrix_identity[member] = dict(format='CSC', shape=list(matrix.shape),
            data=numeric_array_identity(matrix.data), indices=numeric_array_identity(matrix.indices, integer=True),
            indptr=numeric_array_identity(matrix.indptr, integer=True))
        if body_response_evidence is not None:
            mechanical = body_response_evidence[member]
            if (mechanical['displacement_identity'] != response[member]['displacement_m'] or
                    mechanical['nodal_force_identity'] != response[member]['nodal_force_n'] or
                    mechanical['mesh_nodes_identity'] != numeric_array_identity(body.mesh.nodes_m) or
                    mechanical['mesh_tetrahedra_identity'] != numeric_array_identity(body.mesh.tetrahedra, integer=True)):
                raise ValueError('离线Ku−F证据没有绑定相同完整u/F和体网格')
            residual = mechanical['relative_free_residual']
            fixed_max = mechanical['maximum_fixed_displacement_m']
        else:
            if not hasattr(body, 'stiffness') or not hasattr(body, 'free_dofs'):
                raise ValueError('没有原K时必须提供真实分块Ku−F证据，不能仅复制原残差字段')
            internal = np.asarray(body.stiffness@u).reshape(-1)
            free = body.free_dofs
            residual = float(np.linalg.norm((internal-stored_force)[free])/max(float(np.linalg.norm(stored_force[free])), 1.))
            fixed_max = float(np.abs(u[body.fixed_dofs]).max())
        checks[member] = dict(passed=bool(error <= 1e-10 and residual <= 1e-8 and fixed_max <= 1e-14),
            Wf_relative_error=error, Ku_minus_F_relative_free_residual=residual,
            maximum_fixed_displacement_m=fixed_max,
            maximum_displacement_norm_m=float(np.linalg.norm(u.reshape((-1, 3)), axis=1).max()),
            response_identity=response[member], mesh_nodes_identity=numeric_array_identity(body.mesh.nodes_m),
            mesh_tetrahedra_identity=numeric_array_identity(body.mesh.tetrahedra, integer=True))
        if native_proof is not None:
            if body_response_evidence is not None:
                if 'strain_energy_j' not in mechanical or 'external_half_work_j' not in mechanical:
                    raise ValueError('原生离线Ku证据还须提供同次应变能和外功，不能仅复制残差')
                strain = _finite(mechanical['strain_energy_j'], '独立FE应变能J')
                half_work = _finite(mechanical['external_half_work_j'], '独立FE半外功J')
            else:
                strain = .5*float(u@internal); half_work = .5*float(u@stored_force)
            energy_scale = max(abs(strain), abs(half_work), 1e-30)
            error_energy = max(abs(strain-half_work), abs(strain-saved['strain_energy_j']),
                abs(2*half_work-saved['external_work_j']))/energy_scale
            checks[member].update(strain_energy_j=strain, external_half_work_j=half_work,
                original_energy_relative_error=error_energy)
            checks[member]['passed'] = checks[member]['passed'] and error_energy <= 1e-8
            native_energies.append(strain)
    opening = sum(openings)
    saved_opening = np.concatenate([np.asarray(field['elastic_deflection_m'])[np.asarray(field['valid_mask'], dtype=bool)] for field in result['fields']])
    gaps = native_gaps if native_proof is not None else np.asarray([record['gap_m'] for record in records])
    levers = native_levers['wheel'] if native_proof is not None else np.asarray([record['levers_m']['wheel'] for record in records])
    remaining = gaps+opening-levers*result['equilibrium']['theta_rad']
    remaining_error = float(np.abs(remaining-np.asarray(result['equilibrium']['remaining_gaps_m'])).max())
    opening_error = float(np.abs(opening-saved_opening).max())
    pinion_levers = native_levers['pinion'] if native_proof is not None else np.asarray([record['levers_m']['pinion'] for record in records])
    pinion_torque = float(pinion_levers@forces)
    wheel_torque = float(levers@forces)
    torque = result['geometry']['pinion_torque_nm']
    ratio = result['geometry']['pinion_to_wheel_ratio']
    torque_error = max(abs(pinion_torque-torque)/abs(torque), abs(wheel_torque-torque/ratio)/abs(torque/ratio))
    wheel_error = abs(wheel_torque-torque/ratio)/abs(torque/ratio)
    pinion_residual = (pinion_torque-torque)/abs(torque)
    pinion_reconstruction_error = abs(pinion_residual-result['pinion_torque_relative_residual'])
    weak = distributed_average_kkt_evidence(result)
    passed = bool(all(check['passed'] for check in checks.values()) and remaining_error <= 2e-10
        and opening_error <= 2e-10 and wheel_error <= 1e-8 and pinion_reconstruction_error <= 1e-10 and weak['passed'])
    evidence = dict(schema='distributed_P0_response_binding_v1', passed=passed, body_checks=checks,
        W_matrix_identity=matrix_identity, full_response_fingerprint=distributed_response_fingerprint(result),
        maximum_average_gap_reconstruction_error_m=remaining_error, maximum_WTu_opening_error_m=opening_error,
        maximum_torque_relative_error=torque_error, reconstructed_pinion_torque_nm=pinion_torque,
        reconstructed_wheel_torque_nm=wheel_torque, weak_average_kkt=weak,
        wheel_control_balance_relative_error=wheel_error, nominal_pinion_torque_relative_residual=pinion_residual,
        pinion_reconstruction_error=pinion_reconstruction_error, pinion_torque_exactly_controlled=False,
        actual_loaded_pose=distributed_loaded_pose(result),
        quadrature_record_fingerprint=native_proof['record_geometry_fingerprint'] if native_proof is not None else _digest(records),
        factorization_performed=False, displacement_resolved=False,
        definition='真实同列Wf、Wᵀu、原P0平均间隙与两体Ku−F；逐Gauss非穿透另列')
    if native_proof is not None:
        summed = math.fsum(native_energies)
        condensed = _finite(result['equilibrium']['strain_energy_j'], '原凝聚FE应变能J')
        energy_error = max(abs(summed-condensed), abs(summed-_finite(result['fe_strain_energy_j'], '原两体FE应变能J')))/max(abs(summed),abs(condensed),1e-30)
        evidence.update(native_geometry_recomputed=native_proof,
            recomputed_mean_gap_identity=numeric_array_identity(gaps),
            recomputed_lever_identity={m: numeric_array_identity(v) for m, v in native_levers.items()},
            total_recomputed_FE_strain_energy_j=summed, condensed_energy_relative_error=energy_error,
            definition='真实FE交片平均g0/W/力臂独立重建，原Wf/Wᵀu与两体Ku−F/能量；非NURBS采样')
        evidence['passed'] = evidence['passed'] and energy_error <= 1e-8
    return evidence


def distributed_pressure_domain_background(result, background_polygons_m, *, tooth_offset=0,
        certificate=None, pinion_request=None):
    """正式读取P0专用域证据，保持其离散平均语义，不套节点KKT要求。"""
    from .pressure_background import (DOMAIN_CERTIFICATE_SCHEMA, flank_background_fingerprint,
        _partitioned_domain)
    state = distributed_pressure_fingerprint(result, tooth_offset=tooth_offset)
    background = flank_background_fingerprint(background_polygons_m)
    pose = distributed_loaded_pose(result)
    base = dict(method='full_gear_fe', contact_representation=P0_REPRESENTATION, tooth_offset=tooth_offset,
        driving_angle_rad=pose['driving_angle_rad'], pinion_torque_nm=pose['pinion_torque_nm'],
        state_fingerprint=state, flank_background_fingerprint=background,
        pointwise_nonpenetration_certified=False, pressure_mesh_convergence_certified=False,
        actual_loaded_gauss_minimum_gap_m=result.get('geometry_check', {}).get('actual_loaded_minimum_gauss_gap_m'),
        zero_definition='原有效P0压力自由度及守恒中心P1显示域的离散零值；远域另经整轮距离/位移/材料排除，逐Gauss非穿透独立验证')
    if certificate is None:
        return dict(base, full_flank_zero_certified=False, unknown_domain_present=True,
            complete_working_flank=False, certificate_passed=None, domain_state='partial_or_uncertified')
    if certificate.get('schema') != DOMAIN_CERTIFICATE_SCHEMA or certificate.get('contact_representation') != P0_REPRESENTATION:
        raise ValueError('P0压力必须使用明确的分布牵引域证书，不能复用节点接触证书')
    for key, expected in (('method', 'full_gear_fe'), ('tooth_offset', tooth_offset),
            ('state_fingerprint', state), ('flank_background_fingerprint', background)):
        if certificate.get(key) != expected:
            raise ValueError('P0域证书与实际状态/背景不一致：'+key)
    for key in ('driving_angle_rad', 'pinion_torque_nm'):
        if not math.isclose(_finite(certificate.get(key), key), pose[key], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('P0域证书相位或载荷不匹配')
    checks = certificate.get('checks', {})
    required = ('input_geometry_bound', 'included_P0_average_kkt_passed', 'body_response_bound',
        'actual_loaded_pose_bound', 'audit_inputs_and_source_unchanged')
    if any(checks.get(key) is not True for key in required):
        raise ValueError('P0证书缺少几何、弱平均KKT、真实u/F绑定或源码不变证据')
    identity = certificate.get('audit_identity', {})
    if (identity.get('unchanged') is not True or identity.get('state_fingerprint_before') != state
            or identity.get('state_fingerprint_after') != state or
            identity.get('source_sha256_before') != identity.get('source_sha256_after')):
        raise ValueError('P0证书审核前后输入或执行源码不一致')
    response = certificate.get('distributed_response_recomputed', {})
    if response.get('passed') is not True or response.get('full_response_fingerprint') != distributed_response_fingerprint(result):
        raise ValueError('P0证书未绑定全部同次u/F、实际W或弱平均约束')
    count = certificate.get('unresolved_region_count')
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError('P0证书未决计数不合法')
    complete = certificate.get('passed') is True and certificate.get('complete_working_flank') is True
    if complete:
        if count or any(checks.get(key) is not True for key in P0_CHECKS):
            raise ValueError('P0完整齿面离散零底仍有未决区域或缺少远域证明')
        # 完整证书同样逐片检查原CAD四分树，不能只信总面积和passed标签。
        partition = _partitioned_domain(certificate, background_polygons_m, result, pinion_request,
            'full_gear_fe', tooth_offset, allow_complete=True)
        return dict(base, full_flank_zero_certified=True, unknown_domain_present=False,
            complete_working_flank=True, certificate_passed=True, domain_state='complete_pressure_domain',
            evidence_kind=certificate['evidence_kind'], checks=dict(checks), **partition)
    if certificate.get('passed') is not False or certificate.get('complete_working_flank') is not False or count <= 0:
        raise ValueError('P0部分覆盖的通过状态与未决区域不一致')
    for key in ('full_working_flank_covered', 'omitted_regions_excluded', 'unresolved_regions_excluded', 'far_field_motion_bound_passed'):
        if checks.get(key) is not False:
            raise ValueError('P0部分覆盖不能保留完整远域通过标志')
    partition = _partitioned_domain(certificate, background_polygons_m, result, pinion_request,
        'full_gear_fe', tooth_offset)
    return dict(base, full_flank_zero_certified=False, unknown_domain_present=True,
        complete_working_flank=False, certificate_passed=False, domain_state='partitioned_pressure_domain',
        evidence_kind=certificate['evidence_kind'], checks=dict(checks), **partition)


class DistributedPressureDomainAudit:
    """实际完整CAD背景与同次P0 FE响应的轻量域上下文；不持有新FE分解。

prepared_surface提供models/tca。内部从同模型重建完整参考工作面；若传
flank_mesh则其有序背景必须与该完整面一致。目标可用完整整轮三角边界，
或附全原单元AABB覆盖证明的候选子网；每查询检查裁剪盒完整性，失败
自动回退实际整轮BRep点距。所有弦差均明确为工程设置。
"""
    def __init__(self, prepared_surface, *, flank_mesh=None, source_deflection_m=1e-6,
            opposing_triangles_m=None, opposing_complete_box_m=None, opposing_identity=None,
            target_deflection_m=1e-6, cad_identity=None, tolerance_m=1e-10):
        import numpy as np
        from .surface_workflow import _flank_mesh
        from .refined_fe_plotting import _background
        from .pressure_background import _digest, flank_background_fingerprint
        self.prepared = prepared_surface
        self.models = getattr(prepared_surface, 'models', None) or prepared_surface.tca.cad_models
        self.source_deflection = source_deflection_m; self.target_deflection = target_deflection_m
        self.tolerance = tolerance_m; self.cad_identity = {} if cad_identity is None else dict(cad_identity)
        for model in self.models.values():
            if model.diagnostics.get('valid') is not True or model.diagnostics.get('solids') != 1:
                raise ValueError('P0材料域需要两轮各一个实际有效闭实体')
        canonical = _flank_mesh(prepared_surface, linear_deflection_m=source_deflection_m, angular_deflection_rad=.03)
        actual_sampling = self.models['pinion'].grid['patches']['left_flank']
        canonical.setdefault('metadata', {})['source_surface_sampling'] = {
            key: actual_sampling[key] for key in ('width_count', 'profile_count')}
        selected, background, _ = _background(canonical, self.models['pinion'].request)
        if flank_mesh is not None:
            _, supplied, _ = _background(flank_mesh, self.models['pinion'].request)
            if flank_background_fingerprint(supplied) != flank_background_fingerprint(background):
                raise ValueError('外部背景未精确绑定同模型完整实际齿面；请使用context.flank_mesh')
        self.flank_mesh = canonical; self.background = background
        self.background_fingerprint = flank_background_fingerprint(background)
        ids = np.asarray(canonical['triangles'], dtype=np.int64)[np.asarray(selected['triangle_indices'], dtype=np.int64)]
        self.triangles = np.asarray(canonical['vertices_m'])[ids]
        self.source_geometry_identity = _digest(dict(request=self.models['pinion'].request,
            blank=self.models['pinion'].settings, grid=self.models['pinion'].grid,
            full_background=self.background_fingerprint, cad_identity=self.cad_identity))
        self.target = None
        self.target_identity = opposing_identity
        self.complete_box = None if opposing_complete_box_m is None else np.asarray(opposing_complete_box_m, dtype=float)
        if self.complete_box is not None and self.complete_box.shape != (2, 3):
            raise ValueError('完整候选盒须为2×3的实际m坐标')
        if opposing_triangles_m is not None:
            if (not isinstance(opposing_identity, dict) or opposing_identity.get('target_cad_sha256') != self.cad_identity.get('wheel') or
                    not opposing_identity.get('original_whole_mesh_sha256')):
                raise ValueError('目标三角必须绑定同一整轮CAD及原完整网格身份')
            if self.complete_box is not None and opposing_identity.get('all_original_cell_aabbs_tested') is not True:
                raise ValueError('裁剪目标须有全原单元AABB覆盖证据，不能用近邻面子集')
            if self.complete_box is None and opposing_identity.get('covers_complete_closed_boundary') is not True:
                raise ValueError('没有裁剪盒时须提供整个闭体边界')
            from .triangle_surface_distance import TriangleSurfaceDistance
            self.target = TriangleSurfaceDistance(opposing_triangles_m, numerical_tolerance_m=tolerance_m)

    def _loaded_frames(self, result):
        from .generated_tca import generated_tca_frames
        pose = distributed_loaded_pose(result)
        return generated_tca_frames(self.prepared.tca, pose['driving_angle_rad'], pose['loaded_follower_angle_rad'])

    def _point_distance(self, point, wheel_rotation):
        import numpy as np
        local = np.asarray(point)@wheel_rotation
        if self.target is not None:
            target = self.target
            target.locator.FindClosestPoint(local, target.closest, target.cellid, target.subid, target.distance2)
            distance = math.sqrt(float(target.distance2))
            complete = self.complete_box is None or bool(np.all(local-distance-self.tolerance >= self.complete_box[0]-1e-14)
                and np.all(local+distance+self.tolerance <= self.complete_box[1]+1e-14))
            if complete:
                return dict(done=True, distance_m=distance, error_bound_m=self.target_deflection+self.tolerance,
                    query_kind='完整整轮三角边界点距；逐查询完整候选盒检查', complete_target_search=True,
                    closest_target_cell=int(target.cellid), target_chord_setting_m=self.target_deflection)
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.gp import gp_Pnt
        query = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*(local*1000.))).Vertex(), self.models['wheel'].shape)
        query.SetDeflection(self.tolerance*1000.); query.Perform()
        return dict(done=bool(query.IsDone()), distance_m=float(query.Value())*.001 if query.IsDone() else None,
            error_bound_m=self.tolerance, query_kind='实际完整整轮BRep点距', complete_target_search=True,
            target_chord_setting_m=0.)

    def geometry_partition(self, result, *, tooth_offset=0, maximum_depth=2, maximum_regions=300000,
            response_evidence=None, progress=None):
        """几何/位移先行结果；没有W/Ku绑定时不能把本结果当正式压力域证书。"""
        import numpy as np
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.TopAbs import TopAbs_OUT
        from OCP.gp import gp_Pnt
        from .pressure_domain_audit import _grid_cover, audit_triangle_pressure_domain, _periodic_pinion_triangles
        from .pressure_material_domain import exact_vertex_components
        distributed_pressure_fingerprint(result, tooth_offset=tooth_offset)
        frames = self._loaded_frames(result)
        pinion_rotation, wheel_rotation = (np.asarray(frames[m+'_rotation_rows']) for m in ('pinion', 'wheel'))
        local, _, _ = _periodic_pinion_triangles(self.triangles, tooth_offset,
            self.models['pinion'].request.tooth_count, self.models['wheel'].request.tooth_count)
        world = local@pinion_rotation.T
        patch = next(p for p in result['geometry']['patches'] if p['tooth_offset'] == tooth_offset)
        from .native_pressure_domain import native_geometry_mode
        if native_geometry_mode(result) and 'source_material_reference_frames' not in patch:
            raise ValueError('原生FE完整齿面审核须提供参考材料齿Rk，不能拿整体R0代替')
        material_frames = patch['source_material_reference_frames'] if native_geometry_mode(result) else patch['frames']
        if not np.allclose(self.triangles@np.asarray(material_frames['pinion_rotation_rows']).T, world, atol=1e-12, rtol=0):
            raise ValueError('P0实际物理齿与参考齿面的周期位姿不一致')
        if tooth_offset in result['selected_tooth_offsets']:
            field = result['fields'][result['selected_tooth_offsets'].index(tooth_offset)]
            cover = _grid_cover(patch, field, padding_m=self.source_deflection)
        else:
            cover = lambda triangle: False
        if 'body_fields' in result:
            maxima = {m: float(np.linalg.norm(result['body_fields'][m]['displacements_m'], axis=1).max()) for m in self.models}
        else:
            if response_evidence is None or response_evidence.get('full_response_fingerprint') != distributed_response_fingerprint(result):
                raise ValueError('外存响应的几何域计算仍须绑定真实同次位移上界')
            maxima = {m: value['maximum_displacement_norm_m'] for m, value in distributed_body_identity(result).items()}
        motion = sum(maxima.values())
        regions = audit_triangle_pressure_domain(world,
            distance_oracle=lambda point: self._point_distance(point, wheel_rotation), solved_region_oracle=cover,
            relative_motion_bound=lambda triangle: motion, chord_error_m=self.source_deflection,
            maximum_depth=maximum_depth, maximum_regions=maximum_regions, tolerance_m=self.tolerance, progress=progress)
        rows = regions['regions']
        far = [i for i, row in enumerate(rows) if row['state'] == 'excluded_by_distance_bound']
        classifier = BRepClass3d_SolidClassifier(self.models['wheel'].shape)
        components = []
        for component, members in enumerate(exact_vertex_components([row['vertices_world_m'] for row in rows], far)):
            representative = min(members, key=lambda i: rows[i]['distance_lower_bound_m']-motion)
            point = np.asarray(rows[representative]['centre_world_m'])@wheel_rotation
            classifier.Perform(gp_Pnt(*(point*1000.)), self.tolerance*1000.)
            outside = classifier.State() == TopAbs_OUT
            components.append(dict(component_index=component, region_indices=members,
                material_state=str(classifier.State()), outside=bool(outside), exact_shared_vertices_only=True,
                representative_region=representative, representative_target_local_m=point.tolist()))
            for index in members:
                rows[index]['whole_material_component'] = component
                if not outside:
                    rows[index]['state'] = 'unresolved'
                    rows[index]['rigid_material_inside_or_on_not_assumed_contact'] = True
        unknown = [row for row in rows if row['state'] == 'unresolved']
        regions.update(passed=not unknown, unresolved_region_count=len(unknown),
            unresolved_source_triangles=sorted({row['source_triangle_index'] for row in unknown}),
            area_m2={name: math.fsum(row['area_m2'] for row in rows if row['state'] == name)
                for name in ('original_solved_domain', 'excluded_by_distance_bound', 'unresolved')})
        return dict(regions=regions, loaded_opposing_rotation_rows=wheel_rotation.tolist(),
            maximum_displacements_m=maxima, relative_elastic_norm_bound_m=motion, material_components=components,
            actual_loaded_pose=distributed_loaded_pose(result), geometry_only=True,
            absent_selected_field=tooth_offset not in result['selected_tooth_offsets'],
            near_domain_definition='原有效P0压力自由度的已求解显示域；平均约束不等于逐点非穿透',
            strict_continuous_cad_proof=False)

    def audit(self, result, *, tooth_offset=0, response_evidence, maximum_depth=2, maximum_regions=300000, progress=None):
        """正式P0域证书；真实响应绑定缺失或不符时明确拒绝。"""
        from pathlib import Path
        from .pressure_background import DOMAIN_CERTIFICATE_SCHEMA
        before = distributed_pressure_fingerprint(result, tooth_offset=tooth_offset)
        source_before = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        if response_evidence.get('passed') is not True or response_evidence.get('full_response_fingerprint') != distributed_response_fingerprint(result):
            raise ValueError('P0完整域入口需要同一次真实Wf/Ku−F及平均约束验收')
        geometric = self.geometry_partition(result, tooth_offset=tooth_offset, maximum_depth=maximum_depth,
            maximum_regions=maximum_regions, response_evidence=response_evidence, progress=progress)
        regions = geometric['regions']; complete = regions['passed']; pose = geometric['actual_loaded_pose']
        weak = distributed_average_kkt_evidence(result)
        checks = dict(input_geometry_bound=True, included_P0_average_kkt_passed=weak['passed'],
            body_response_bound=True, actual_loaded_pose_bound=bool(result.get('geometry_linearization_checked') and result.get('geometry_check', {}).get('passed')),
            full_working_flank_covered=complete, omitted_regions_excluded=complete,
            unresolved_regions_excluded=complete, far_field_motion_bound_passed=complete)
        after = distributed_pressure_fingerprint(result, tooth_offset=tooth_offset)
        source_after = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        unchanged = before == after and source_before == source_after
        checks['audit_inputs_and_source_unchanged'] = unchanged
        certificate = dict(schema=DOMAIN_CERTIFICATE_SCHEMA, method='full_gear_fe',
            contact_representation=P0_REPRESENTATION, tooth_offset=tooth_offset,
            driving_angle_rad=pose['driving_angle_rad'], pinion_torque_nm=pose['pinion_torque_nm'],
            state_fingerprint=before, flank_background_fingerprint=self.background_fingerprint,
            passed=all(checks.values()), complete_working_flank=complete, checks=checks,
            unresolved_region_count=regions['unresolved_region_count'], regions=regions,
            distributed_response_recomputed=response_evidence, weak_average_kkt_recomputed=weak,
            geometry_binding=dict(source_geometry_identity=self.source_geometry_identity, cad_identity=self.cad_identity,
                source_working_flank_generated_from_same_model=True, source_chord_setting_m=self.source_deflection,
                target_identity=self.target_identity),
            loaded_opposing_rotation_rows=geometric['loaded_opposing_rotation_rows'],
            maximum_body_displacements_m=geometric['maximum_displacements_m'],
            whole_material_components=geometric['material_components'],
            audit_identity=dict(state_fingerprint_before=before, state_fingerprint_after=after,
                source_sha256_before=source_before, source_sha256_after=source_after, unchanged=unchanged),
            evidence_kind='原P0自由度/中心P1显示已解域＋完整整轮距离、同次u范数界及真实材料OUT；平均KKT与逐Gauss分开',
            strict_continuous_cad_proof=False,
            limitations=['CAD弦差为明确工程设置，非独立连续曲面全域严格误差上界',
                'P0有效小片的已求解状态不证明每个Gauss点非穿透',
                '刚体CAD材料内且可能被弹性打开的区域保留未知，须进一步实际u映射审核'])
        distributed_pressure_domain_background(result, self.background, tooth_offset=tooth_offset,
            certificate=certificate, pinion_request=self.models['pinion'].request)
        return certificate


def audit_distributed_external_window(prepared_surface, result, *, response_evidence,
        external_offsets=None, tolerance_m=1e-10):
    """完整外齿工作面到真实最终整轮CAD的正距离，扣同次两体全局u范数。

全局P1形函数为凸权，因此任意材料迹位移范数不超过所有体节点最大值。
此函数使用实际最终phi，额外刚体位移不再重复扣入total TE。刚体材料
内/ON或未分离时保留失败，不能据此判接触，也不能补零。
"""
    import numpy as np
    from pathlib import Path
    from OCP.BRep import BRep_Builder, BRep_Tool
    from OCP.TopoDS import TopoDS_Compound, TopoDS
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopTools import TopTools_IndexedMapOfShape
    from OCP.TopExp import TopExp
    from OCP.TopAbs import TopAbs_VERTEX, TopAbs_OUT
    from OCP.gp import gp_Trsf
    from .generated_halfspace import _faces
    from .generated_tca import generated_tca_frames
    before = distributed_response_fingerprint(result)
    if response_evidence.get('passed') is not True or response_evidence.get('full_response_fingerprint') != before:
        raise ValueError('外齿排除必须绑定原Wf/Ku与同次真实P0载荷响应')
    requested = sorted(result['requested_tooth_offsets'])
    if requested != list(range(min(requested), max(requested)+1)):
        raise ValueError('外齿窗口不能跳过内部物理齿')
    offsets = [min(requested)-1, max(requested)+1] if external_offsets is None else list(external_offsets)
    if len(offsets) != len(set(offsets)) or any(isinstance(k, bool) or not isinstance(k, int) or k in requested for k in offsets):
        raise ValueError('外齿须为请求窗口之外不重复的真实整数齿偏移')
    models = getattr(prepared_surface, 'models', None) or prepared_surface.tca.cad_models
    pose = distributed_loaded_pose(result)
    frames = generated_tca_frames(prepared_surface.tca, pose['driving_angle_rad'], pose['loaded_follower_angle_rad'])
    identity = distributed_body_identity(result)
    maxima = {m: value['maximum_displacement_norm_m'] for m, value in identity.items()}
    motion = sum(maxima.values())
    def transform(shape, rotation):
        value = gp_Trsf(); value.SetValues(*(float(x) for row in rotation for x in (*row, 0.)))
        return BRepBuilderAPI_Transform(shape, value, True).Shape()
    target = transform(models['wheel'].shape, frames['wheel_rotation_rows'])
    classifier = BRepClass3d_SolidClassifier(target)
    rows = []
    source_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for offset in offsets:
        tooth = (-offset) % models['pinion'].request.tooth_count
        faces = _faces(models['pinion'].shape, models['pinion'].slot_face_groups[tooth]['left_flank'])
        compound = TopoDS_Compound(); builder = BRep_Builder(); builder.MakeCompound(compound)
        for face in faces: builder.Add(compound, face)
        source = transform(compound, frames['pinion_rotation_rows'])
        query = BRepExtrema_DistShapeShape(source, target)
        query.SetDeflection(tolerance_m*1000.); query.Perform()
        done = bool(query.IsDone()); raw = float(query.Value())*.001 if done else None
        inner = bool(query.InnerSolution()) if done else None
        probes = []
        for face in faces:
            shape = transform(face, frames['pinion_rotation_rows'])
            vertices = TopTools_IndexedMapOfShape(); TopExp.MapShapes_s(shape, TopAbs_VERTEX, vertices)
            for i in range(1, vertices.Extent()+1):
                point = BRep_Tool.Pnt_s(TopoDS.Vertex_s(vertices.FindKey(i)))
                classifier.Perform(point, tolerance_m*1000.)
                probes.append(dict(material_state=str(classifier.State()), outside=classifier.State() == TopAbs_OUT,
                    source_point_world_m=(np.asarray(point.Coord())*.001).tolist()))
        gap = raw-motion-tolerance_m if raw is not None else None
        passed = bool(done and inner is False and gap > 0 and probes and all(p['outside'] for p in probes))
        rows.append(dict(tooth_offset=offset, actual_pinion_tooth=tooth, raw_CAD_distance_m=raw,
            minimum_remaining_gap_m=gap, relative_elastic_norm_bound_m=motion,
            complete_actual_finite_faces=len(faces), inner_solution=inner, actual_material_probes=probes,
            passed=passed, evidence_kind='真实最终CAD整片到整轮距离−同次全局P1体位移范数界；实际材料OUT'))
    after = distributed_response_fingerprint(result)
    source_after = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    unchanged = before == after and source_sha == source_after
    return dict(schema='distributed_P0_external_window_v1', method='full_gear_fe',
        contact_representation=P0_REPRESENTATION, requested_tooth_offsets=requested, external_tooth_offsets=offsets,
        rows=rows, passed=bool(unchanged and all(row['passed'] for row in rows)),
        full_response_fingerprint=before, body_response_identity=identity,
        maximum_body_displacements_m=maxima, actual_loaded_pose=pose,
        response_W_Ku_bound=True, actual_geometry_checked=True, original_state_unchanged=before == after,
        source_sha256_before=source_sha, source_sha256_after=source_after,
        strict_continuous_cad_proof=False,
        interpretation='几何正距与同次FE迹位移界排除外齿；P0平均KKT、逐Gauss和网格收敛分别验收')


def distributed_gauss_mean_gaps(prepared_surface, result, *, batch_size=4096, progress=None):
    """同原参考位姿/任意plane basis重采样Gauss平均g0；无FE映射或分解。"""
    import numpy as np
    from .native_pressure_domain import native_geometry_mode
    if native_geometry_mode(result):
        raise ValueError('FE_boundary须重建真实FE交片平均g0，不能使用NURBS Gauss或中心间隙')
    from .generated_halfspace import _initial_uv
    from .nurbs_batch import sample_pair_batch
    order = result['quadrature_order']
    locations, weights = np.polynomial.legendre.leggauss(order)
    gauss = [(float(a), float(b), float(wa*wb/4)) for b, wb in zip(locations, weights) for a, wa in zip(locations, weights)]
    values, keys = [], []
    for patch in result['geometry']['patches']:
        valid = np.asarray(patch['valid_mask'], dtype=bool); indices = list(zip(*np.nonzero(valid)))
        if not indices: continue
        grid = patch['grid'].as_dict() if hasattr(patch['grid'], 'as_dict') else patch['grid']
        dx, dy = grid['spacing_x_m'], grid['spacing_y_m']; ox, oy = grid['origin_xy_m']
        origin, basis = np.asarray(patch['plane_origin_world_m']), np.asarray(patch['plane_tangent_rows'])
        normal, frames = np.asarray(patch['plane_outward_normal_world']), patch['frames']
        seeds = {m: _initial_uv(prepared_surface.trimmed[m].surface, np.asarray(frames[m+'_rotation_rows']),
            origin, basis, patch['seed_uv'][m]) for m in ('pinion', 'wheel')}
        cell_count = max(1, batch_size//len(gauss))
        for start in range(0, len(indices), cell_count):
            chosen = indices[start:start+cell_count]
            xy = np.asarray([(ox+(i+a*.5)*dx, oy+(j+b*.5)*dy) for j, i in chosen for a, b, _ in gauss])
            samples, failures = sample_pair_batch(prepared_surface, frames, origin, basis, normal,
                patch['seed_uv'], xy, 2e-11, .3, .3, seeds)
            if any(sample is None for sample in samples):
                raise ArithmeticError('原有效P0小片Gauss重新查询失败，未补零：'+str(next(f for sample, f in zip(samples, failures) if sample is None)))
            gaps = np.asarray([sample['gap_m'] for sample in samples]).reshape((-1, len(gauss)))
            values.extend((gaps@np.asarray([g[2] for g in gauss])).tolist())
            keys.extend((patch['tooth_offset'], int(j), int(i)) for j, i in chosen)
        if progress: progress(dict(stage='P0_Gauss_mean_gap_reconstructed', tooth_offset=patch['tooth_offset'], cells=len(indices)))
    return np.asarray(values), np.asarray(keys, dtype=np.int64)


def audit_distributed_response_archive(state_path, *, run_root, prepared_surface, meshes,
        materials, recheck_body=True, progress=None):
    """读取原工作流NPZ/W并复验同次响应，返回加回原u/F的派生视图和证据。

原JSON与NPZ均只读。外存文件必须匹配原状态receipt SHA、输入指纹与物理
格索引；Wf、Wᵀu、真实Gauss平均g0和两体力矩重新计算。默认还分块重算
Ku−F；不构造/分解全局K。返回视图可保留distributed_body_identity后移除
body_fields，身份不变，供流式周期后处理使用。
"""
    from pathlib import Path
    import numpy as np
    from scipy import sparse
    from .refined_fe_postprocessing import _relative_file
    from .generated_tca import generated_tca_frames
    from .pressure_background import _digest
    state_path, root = Path(state_path).resolve(), Path(run_root).resolve()
    def sha(path):
        value = hashlib.sha256()
        with Path(path).open('rb') as stream:
            while block := stream.read(4*1024*1024): value.update(block)
        return value.hexdigest()
    state_sha = sha(state_path)
    if state_path.suffix == '.gz':
        import gzip
        with gzip.open(state_path, 'rt', encoding='utf-8-sig') as stream:state = json.load(stream)
    else:state = json.loads(state_path.read_text(encoding='utf-8-sig'))
    if state.get('accepted') is not True or state.get('source_identity_unchanged') is not True:
        raise ValueError('P0响应档案只接受原业务已接受、源身份不变的状态')
    result = dict(state['result']); receipt = state['response_receipt']
    from .native_pressure_domain import native_geometry_mode
    if native_geometry_mode(result):
        raise ValueError('原生FE档案须使用audit_native_response_archive并提供真实分片系统，不走NURBS档案分支')
    if receipt.get('passed') is not True or receipt.get('input_fingerprint') != state['input_fingerprint']:
        raise ValueError('P0原响应receipt未通过或不属于本输入')
    first = _relative_file(root, receipt['bodies']['pinion']['response_file'])
    receipt_path = first.parent/'receipt.json'
    final_receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
    if final_receipt.get('source_state_sha256') != state_sha or final_receipt['bodies'] != receipt['bodies']:
        raise ValueError('外存响应receipt未绑定本原始JSON字节')
    models = getattr(prepared_surface, 'models', None) or prepared_surface.tca.cad_models
    body_fields, arrays, matrices, files = {}, {}, {}, {}
    force = np.asarray(result['equilibrium']['forces_n'], dtype=float)
    expected_keys = np.asarray([(field['tooth_offset'], int(j), int(i)) for field in result['fields']
        for j, i in zip(*np.nonzero(np.asarray(field['valid_mask'], dtype=bool)))], dtype=np.int64)
    for member in ('pinion', 'wheel'):
        item = receipt['bodies'][member]
        response_path = _relative_file(root, item['response_file']); W_path = _relative_file(root, item['load_matrix_file'])
        if sha(response_path) != item['response_sha256'] or sha(W_path) != item['load_matrix_sha256']:
            raise ValueError('原P0响应或W的外存SHA不一致')
        with np.load(response_path, allow_pickle=False) as saved:
            arrays[member] = {name: saved[name] for name in saved.files}
        values = arrays[member]
        if not np.array_equal(values['force_plane_n'], force) or not np.array_equal(values['contact_cell_keys'], expected_keys):
            raise ValueError('外存原W列/P0力/实际格索引不匹配，不能跨压力网格复用')
        matrix = sparse.load_npz(W_path)
        if matrix.shape != tuple(item['load_matrix_shape']) or matrix.shape != (3*len(meshes[member].nodes_m), len(force)):
            raise ValueError('外存W维度与原网格/P0自由度不一致')
        matrices[member] = matrix
        saved_checks = result['body_response_checks'][member]
        body_fields[member] = dict(displacements_m=values['displacements_m'], nodal_forces_n=values['nodal_forces_n'],
            **{key: saved_checks[key] for key in ('strain_energy_j', 'external_work_j', 'maximum_free_residual_n', 'relative_free_residual')})
        files[member] = dict(response_sha256=item['response_sha256'], W_sha256=item['load_matrix_sha256'])
    result['body_fields'] = body_fields
    result['distributed_body_identity'] = distributed_body_identity(result)
    reference = generated_tca_frames(prepared_surface.tca, result['geometry']['driving_angle_rad'],
        result['geometry']['reference_follower_angle_rad'])
    opening = np.zeros(len(force)); levers = {}; body_checks = {}; W_identity = {}
    for member in ('pinion', 'wheel'):
        matrix, saved = matrices[member], arrays[member]
        u, force_nodes = saved['displacements_m'].reshape(-1), saved['nodal_forces_n'].reshape(-1)
        Wf = np.asarray(matrix@force).reshape(-1)
        Wf_error = float(np.linalg.norm(Wf-force_nodes)/max(float(np.linalg.norm(force_nodes)), 1.))
        opening += np.asarray(matrix.T@u).reshape(-1)
        nodes = np.asarray(meshes[member].nodes_m, dtype=float)
        rotation = np.asarray(reference[member+'_rotation_rows'])
        world = nodes@rotation.T
        velocity_local = np.cross(reference[member+'_axis_world'], world)@rotation
        levers[member] = -np.asarray(matrix.T@velocity_local.reshape(-1)).reshape(-1)
        bore = models[member].settings.bore_diameter_m*.5
        tolerance = max(models[member].settings.boolean_tolerance_m*2, bore*1e-8)
        fixed = np.flatnonzero(np.abs(np.linalg.norm(nodes[:, :2], axis=1)-bore) <= tolerance)
        expected_free = np.setdiff1d(np.arange(3*len(nodes)), (3*fixed[:, None]+np.arange(3)).ravel())
        if not np.array_equal(saved['free_dofs'], expected_free):
            raise ValueError('外存Ku−F自由度与原孔壁约束不一致')
        stored_relative = float(np.linalg.norm(saved['ku_minus_f_n'].reshape(-1)[expected_free])/max(float(np.linalg.norm(force_nodes[expected_free])), 1.))
        if recheck_body:
            mechanical = matrix_free_body_response(meshes[member], materials[member], saved['displacements_m'], saved['nodal_forces_n'], fixed)
        else:
            mechanical = dict(passed=stored_relative <= 1e-8, relative_free_residual=stored_relative,
                maximum_fixed_displacement_m=float(np.abs(saved['displacements_m'][fixed]).max()),
                definition='原工作流在原K上保存的Ku−F向量；读取时重新计算其自由度范数')
        body_checks[member] = dict(passed=bool(mechanical['passed'] and Wf_error <= 1e-10 and stored_relative <= 1e-8),
            Wf_relative_error=Wf_error, Ku_minus_F_relative_free_residual=mechanical['relative_free_residual'],
            stored_Ku_minus_F_relative_free_residual=stored_relative,
            maximum_displacement_norm_m=result['distributed_body_identity'][member]['maximum_displacement_norm_m'],
            response_identity=result['distributed_body_identity'][member], mechanical_recheck=mechanical,
            mesh_nodes_identity=numeric_array_identity(nodes), mesh_tetrahedra_identity=numeric_array_identity(meshes[member].tetrahedra, integer=True),
            model_semantic_identity=_digest(dict(request=models[member].request, settings=models[member].settings, grid=models[member].grid)),
            original_archive_files=files[member])
        W_identity[member] = dict(shape=list(matrix.shape), data=numeric_array_identity(matrix.data),
            indices=numeric_array_identity(matrix.indices, integer=True), indptr=numeric_array_identity(matrix.indptr, integer=True))
        if progress: progress(dict(stage='archived_P0_body_bound', member=member, Wf_relative_error=Wf_error,
            Ku_minus_F_relative_free_residual=mechanical['relative_free_residual']))
    mean_gaps, gap_keys = distributed_gauss_mean_gaps(prepared_surface, result, progress=progress)
    if not np.array_equal(gap_keys, expected_keys):
        raise ValueError('Gauss平均间隙的原格顺序与原W列不一致')
    remaining = mean_gaps+opening-levers['wheel']*result['equilibrium']['theta_rad']
    gap_error = float(np.abs(remaining-np.asarray(result['equilibrium']['remaining_gaps_m'])).max())
    saved_opening = np.concatenate([np.asarray(field['elastic_deflection_m'])[np.asarray(field['valid_mask'], dtype=bool)] for field in result['fields']])
    opening_error = float(np.abs(opening-saved_opening).max())
    torque = result['geometry']['pinion_torque_nm']; ratio = result['geometry']['pinion_to_wheel_ratio']
    pinion_torque, wheel_torque = (float(levers[m]@force) for m in ('pinion', 'wheel'))
    torque_error = max(abs(pinion_torque-torque)/abs(torque), abs(wheel_torque-torque/ratio)/abs(torque/ratio))
    wheel_error = abs(wheel_torque-torque/ratio)/abs(torque/ratio)
    pinion_residual = (pinion_torque-torque)/abs(torque)
    pinion_reconstruction_error = abs(pinion_residual-result['pinion_torque_relative_residual'])
    weak = distributed_average_kkt_evidence(result)
    source = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    evidence = dict(schema='distributed_P0_response_binding_v1',
        passed=bool(all(b['passed'] for b in body_checks.values()) and gap_error <= 2e-10 and opening_error <= 2e-10
            and wheel_error <= 1e-8 and pinion_reconstruction_error <= 1e-10 and weak['passed']),
        full_response_fingerprint=distributed_response_fingerprint(result), body_checks=body_checks,
        W_matrix_identity=W_identity, maximum_average_gap_reconstruction_error_m=gap_error,
        maximum_WTu_opening_error_m=opening_error, maximum_torque_relative_error=torque_error,
        reconstructed_pinion_torque_nm=pinion_torque, reconstructed_wheel_torque_nm=wheel_torque,
        wheel_control_balance_relative_error=wheel_error, nominal_pinion_torque_relative_residual=pinion_residual,
        pinion_reconstruction_error=pinion_reconstruction_error, pinion_torque_exactly_controlled=False,
        recomputed_mean_gap_identity=numeric_array_identity(mean_gaps), weak_average_kkt=weak,
        actual_loaded_pose=distributed_loaded_pose(result), source_state_sha256=state_sha,
        response_receipt_sha256=sha(receipt_path), input_fingerprint=state['input_fingerprint'],
        archived_files_identity=files, original_state_unchanged=sha(state_path) == state_sha,
        audit_source_sha256=source, factorization_performed=False, displacement_resolved=False,
        independent_Ku_recheck=recheck_body,
        definition='原NPZ/W完整身份与格序；Wf/Wᵀu/真实Gauss平均g0/两体力矩；Ku−F原向量及可选独立分块复验')
    if not evidence['original_state_unchanged']:
        raise ArithmeticError('验证期间原P0状态字节改变')
    return result, evidence
