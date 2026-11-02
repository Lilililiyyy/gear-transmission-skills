"""整轮FE周期审核扩展：同一次重算位移复用于压力域和外齿KKT检查。

求解业务通过显式观察回调传入prepared，不修改接触方程或原状态。
外齿检查仍属于同一离散节点候选模型，不冒充连续CAD全域证明。
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

from .pressure_domain_batch import PressureDomainBatchAudit


def _source_sha256():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def recompute_fe_response_with_fields(prepared, result, *, gap_tolerance_m=2e-10):
    """用原完整节点力重算两轮FE；KKT公式与单次域审核相同，同时返回u。"""
    import numpy as np
    from .elasticity import ContactLoadChannel, contact_load_matrix
    rows = result['nodal_contacts']; forces = np.asarray([r['force_n'] for r in rows])
    displacements = {}; errors = []; failures = []
    frames = result['mesh_reference_frames']; theta = result['elastic_te_rad']
    for member, prefix in (('pinion', 'first'), ('wheel', 'second')):
        channels = [ContactLoadChannel(tuple(r[prefix+'_node_indices']), tuple(r[prefix+'_weights']),
            tuple(r[prefix+'_force_normal_local'])) for r in rows]
        matrix = contact_load_matrix(prepared.systems[member], channels, matrix_format='csc', maximum_sparse_channels=100000)
        response = prepared.systems[member].solve(matrix @ forces)
        displacements[member] = np.asarray(response['displacements_m'])
        nodes = np.asarray(prepared.systems[member].mesh.nodes_m)
        rotation = np.asarray(frames[member+'_rotation_rows'])
        for r in rows:
            point = np.asarray(r[prefix+'_weights']) @ nodes[r[prefix+'_node_indices'], :] @ rotation.T
            errors.append(float(np.linalg.norm(point-np.asarray(r[prefix+'_point_world_m']))))
    reconstruction = []
    for r in rows:
        elastic = 0.
        for member, prefix in (('pinion', 'first'), ('wheel', 'second')):
            elastic += float(np.asarray(r[prefix+'_weights']) @ displacements[member][r[prefix+'_node_indices'], :] @ np.asarray(r[prefix+'_force_normal_local']))
        gap = r['signed_gap_m']+elastic-r['wheel_lever_m']*theta
        reconstruction.append(abs(gap-r['remaining_gap_m']))
        if r['force_n'] < -1e-12 or gap < -gap_tolerance_m or (r['force_n'] > 1e-12 and abs(gap) > gap_tolerance_m):
            failures.append(dict(first_pressure_node_index=r['first_pressure_node_index'], remaining_gap_m=gap, force_n=r['force_n']))
        if not math.isclose(r['pressure_pa']*r['area_m2'], r['force_n'], rel_tol=1e-10, abs_tol=1e-9):
            failures.append(dict(first_pressure_node_index=r['first_pressure_node_index'], reason='原节点pA不等于f'))
    maxima = {m: float(np.linalg.norm(u, axis=1).max()) for m, u in displacements.items()}
    radius = float(np.linalg.norm(np.asarray(prepared.systems['wheel'].mesh.nodes_m), axis=1).max())
    motion = sum(maxima.values())+2*radius*abs(math.sin(theta/2))
    ratio = prepared.models['pinion'].request.tooth_count/prepared.models['wheel'].request.tooth_count
    target = result['pinion_torque_nm']/ratio
    torque = sum(r['force_n']*r['wheel_lever_m'] for r in rows)
    if not math.isclose(torque, target, rel_tol=1e-7, abs_tol=1e-7):
        failures.append(dict(reason='重构从动轴力矩不满足全局平衡', torque_nm=torque, target_nm=target))
    evidence = dict(passed=bool(not failures and max(reconstruction, default=math.inf) <= gap_tolerance_m),
        input_node_geometry_bound=bool(errors and max(errors) <= 1e-10), maximum_node_position_error_m=max(errors, default=None),
        included_constraint_count=len(rows), maximum_kkt_reconstruction_error_m=max(reconstruction, default=None),
        violations=failures, maximum_displacement_m=maxima, wheel_radius_bound_m=radius, relative_motion_bound_m=motion,
        wheel_torque_reconstructed_nm=torque, wheel_target_torque_nm=target,
        definition='原节点力重算全轮FE位移；三角线性位移范数界取全部节点最大值，加从动刚体转角弦长')
    return evidence, displacements


class FullGearPressureDomainBatchAudit(PressureDomainBatchAudit):
    """保留当前相位u数组，避免压力域和外候选检查重复分解或重复求响应。"""
    def _fe_evidence(self, result):
        from .pressure_background import _digest
        key = _digest(dict(nodal_contacts=result['nodal_contacts'], elastic_te_rad=result['elastic_te_rad'],
            frames=result['mesh_reference_frames'], pinion_torque_nm=result['pinion_torque_nm']))
        if key != self._last_response_key:
            self._last_response, self._full_displacements = recompute_fe_response_with_fields(self.prepared, result,
                gap_tolerance_m=max(2e-10, self.tolerance_m*1.01))
            self._last_response_key = key; self._response_count += 1
        return self._last_response, key

    def audit(self, result, **kwargs):
        before = _source_sha256()
        certificate = super().audit(result, **kwargs)
        after = _source_sha256()
        if before != after:
            raise ArithmeticError('整轮FE批量响应审核源码在计算中改变')
        certificate['fe_response_extension_identity'] = dict(source_sha256_before=before, source_sha256_after=after,
            unchanged=True, full_displacements_retained=True)
        return certificate

    def audit_external_tooth_constraints(self, result, *, tooth_offsets=(-3, 3), maximum_distance_m=.0002,
            maximum_normal_angle_rad=.3, gap_tolerance_m=1e-10):
        """将外齿试力置零，逐原离散候选检查全轮响应后的剩余间隙。"""
        import numpy as np
        from .nodal_contact import pair_nodal_contacts
        from .pressure_background import _digest
        evidence, response_key = self._fe_evidence(result)
        if evidence['passed'] is not True or evidence['input_node_geometry_bound'] is not True:
            raise ValueError('原完整FE响应/KKT未通过，不能审核外齿')
        offsets = tuple(tooth_offsets)
        if not offsets or len(set(offsets)) != len(offsets) or any(isinstance(k, bool) or not isinstance(k, int) or k in self.tooth_offsets for k in offsets):
            raise ValueError('外齿须为显式请求窗口以外的不重复整数')
        frames = result['mesh_reference_frames']; direction = 1 if result['pinion_torque_nm'] > 0 else -1
        records = []
        for offset in offsets:
            first = self._flank_facets('pinion', offset); second = self._flank_facets('wheel', offset)
            pairs = pair_nodal_contacts(first, second, frames['pinion_rotation_rows'], frames['wheel_rotation_rows'],
                maximum_distance_m=maximum_distance_m, maximum_normal_angle_rad=maximum_normal_angle_rad)
            actual = []
            for pair in pairs:
                normal = np.asarray(pair['normal_world'])
                ell1 = float(normal @ np.cross(frames['pinion_axis_world'], pair['first_point_world_m']))
                ell2 = float(-normal @ np.cross(frames['wheel_axis_world'], pair['second_point_world_m']))
                if direction*ell1 <= 1e-10 or direction*ell2 <= 1e-10:
                    continue
                elastic = 0.
                for member, prefix in (('pinion', 'first'), ('wheel', 'second')):
                    elastic += float(np.asarray(pair[prefix+'_weights']) @ self._full_displacements[member][pair[prefix+'_node_indices'], :] @ np.asarray(pair[prefix+'_force_normal_local']))
                gap = pair['signed_gap_m']+elastic-ell2*result['elastic_te_rad']
                actual.append(dict(first_pressure_node_index=pair['first_pressure_node_index'], signed_gap_m=pair['signed_gap_m'],
                    elastic_normal_displacement_m=elastic, remaining_gap_m=gap, pinion_lever_m=ell1, wheel_lever_m=ell2, force_n=0.))
            local_maxima = {member: float(np.linalg.norm(self._full_displacements[member][np.unique(facets['triangles'])], axis=1).max())
                for member, facets in (('pinion', first), ('wheel', second))}
            records.append(dict(tooth_offset=offset, available_geometric_pairs=len(pairs), accepted_signed_load_pairs=len(actual),
                minimum_remaining_gap_m=min((r['remaining_gap_m'] for r in actual), default=None),
                zero_outer_force_kkt_passed=all(r['remaining_gap_m'] >= -gap_tolerance_m for r in actual),
                local_facet_displacement_max_m=local_maxima, actual_candidate_checks=actual))
        return dict(schema='bevel_external_fe_constraints_v1', method='full_gear_fe',
            driving_angle_rad=result['driving_angle_rad'], pinion_torque_nm=result['pinion_torque_nm'],
            outer_offsets=list(offsets), response_state_fingerprint=response_key,
            whole_state_fingerprint=_digest(result), actual_full_fe_response_recomputed=True,
            candidate_distance_window_m=maximum_distance_m, gap_tolerance_m=gap_tolerance_m,
            passed=all(r['zero_outer_force_kkt_passed'] for r in records), rows=records,
            interpretation='与原FE同实际节点候选、法向及有符号力臂筛选；外齿零试力的离散KKT可行性，不是连续CAD全域接触完备证明')


class FePressureDomainObserver:
    """可直接传给FE业务的state_observer；原state只读，证书独立落盘。

    expected_phase_count显式给定完整相位数组长度；None时只累积证据，
    不宣称完整周期审核已完成。分片仍使用原全局相位索引。
    """
    def __init__(self, output, *, flank_mesh, pinion_request, tooth_offsets=(-2,-1,0,1,2),
            outer_offsets=(-3,3), maximum_depth=2, maximum_distance_m=.0002,
            maximum_normal_angle_rad=.3, engine='triangles', expected_phase_count=None):
        if expected_phase_count is not None and (isinstance(expected_phase_count,bool) or
                not isinstance(expected_phase_count,int) or expected_phase_count<1):
            raise ValueError('预期相位数量须为正整数或None')
        self.output = Path(output); self.flank_mesh = flank_mesh; self.pinion_request = pinion_request
        self.offsets = tuple(tooth_offsets); self.outer_offsets = tuple(outer_offsets)
        self.maximum_depth = maximum_depth; self.maximum_distance_m = maximum_distance_m
        self.maximum_normal_angle_rad = maximum_normal_angle_rad; self.engine = engine
        self.expected_phase_count = expected_phase_count
        self.context = None

    def __call__(self, prepared, phase_index, state_path):
        import json
        from .refined_fe_workflow import _write
        path = Path(state_path); before = hashlib.sha256(path.read_bytes()).hexdigest()
        state = json.loads(path.read_text(encoding='utf-8-sig'))
        if state.get('accepted') is not True or state.get('phase_index') != phase_index:
            raise ValueError('观察回调只接收同相位已接受的完整FE状态')
        if (isinstance(phase_index,bool) or not isinstance(phase_index,int) or phase_index<0 or
                self.expected_phase_count is not None and phase_index>=self.expected_phase_count):
            raise ValueError('相位索引超出显式预期相位数组')
        manifest_path = self.output/'manifest.json'
        manifest = json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.is_file() else dict(
            schema='bevel_cycle_domain_batch_v1', method='full_gear_fe', complete=False, certificates=[],
            failures=[], expected_phase_count=self.expected_phase_count, expected_tooth_offsets=list(self.offsets))
        if manifest.get('expected_phase_count')!=self.expected_phase_count or manifest.get('expected_tooth_offsets')!=list(self.offsets):
            raise ValueError('已有观察清单的相位数量或齿对窗口不同')
        if any(r['phase_index'] == phase_index for r in manifest['certificates']):
            raise ValueError('观察输出已有同相位证书，拒绝覆盖重复相位')
        if self.context is None:
            self.context = FullGearPressureDomainBatchAudit(prepared, self.flank_mesh, self.pinion_request,
                method='full_gear_fe', engine=self.engine, tooth_offsets=self.offsets)
        elif self.context.prepared is not prepared:
            raise ValueError('同一观察对象不能暗换FE准备对象')
        result = state['result']; records = []
        for offset in self.offsets:
            certificate = self.context.audit(result, tooth_offset=offset, maximum_depth=self.maximum_depth)
            folder = self.output/f'phase_{phase_index:04d}'/('k_m'+str(-offset) if offset < 0 else 'k_'+str(offset))
            target = folder/'certificate.json'; _write(certificate, target)
            records.append(dict(phase_index=phase_index, tooth_offset=offset, source_state_sha256=before,
                certificate_file=target.relative_to(self.output).as_posix(), certificate_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                passed=certificate['passed'], partial_domain_evidence_available=bool(not certificate['passed'] and certificate['checks']['input_geometry_bound']
                    and certificate['checks']['included_candidates_kkt_passed'] and certificate['checks']['audit_inputs_and_source_unchanged']),
                checks=certificate['checks'], unresolved_region_count=certificate['unresolved_region_count'],
                area_m2=certificate['regions']['area_m2'], original_state_unchanged=True))
        external = self.context.audit_external_tooth_constraints(result, tooth_offsets=self.outer_offsets,
            maximum_distance_m=self.maximum_distance_m, maximum_normal_angle_rad=self.maximum_normal_angle_rad)
        external.update(phase_index=phase_index, source_state_sha256=before)
        _write(external, self.output/f'phase_{phase_index:04d}'/'external_constraints.json')
        if hashlib.sha256(path.read_bytes()).hexdigest() != before:
            raise ArithmeticError('观察回调改变了原FE状态文件')
        manifest['certificates'].extend(records)
        manifest['complete'] = (self.expected_phase_count is not None and
            len(manifest['certificates']) == self.expected_phase_count*len(self.offsets))
        _write(manifest, manifest_path)
        return dict(phase_index=phase_index, original_state_unchanged=True, certificate_count=len(records),
            external_zero_force_kkt_passed=external['passed'])
