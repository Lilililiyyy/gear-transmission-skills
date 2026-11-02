"""外齿窗口的实际几何排除：完整CAD工作面和真实全轮FE边界的双重界。

不把“没有近域节点候选”当作非接触证明。全轮位移界来自已绑定原状态
且独立重算的域证书，外齿实际有限面及FE三角形分别逐域核查。
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path


def _source_sha256():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _verify_global_fe_response_binding(result, certificate):
    """外域位移界来自全轮原力，不能只比较一个局部齿的节点摘要。"""
    from .pressure_background import _digest
    expected = _digest(dict(nodal_contacts=result['nodal_contacts'],elastic_te_rad=result['elastic_te_rad'],
        frames=result['mesh_reference_frames'],pinion_torque_nm=result['pinion_torque_nm']))
    if certificate.get('fe_response_state_fingerprint') != expected:
        raise ValueError('证书完整FE响应输入与当前原状态不同或缺失')
    return expected


def fe_surface_clearance_bound(source_triangles_world_m, target_distance_oracle,
        relative_motion_bound_m, *, maximum_depth=6, maximum_regions=300000,
        tolerance_m=1e-10, progress=None):
    """对全部真实源FE三角形作1-Lipschitz距离下界，无CAD弦差替代。"""
    from .pressure_domain_audit import audit_triangle_pressure_domain
    if isinstance(relative_motion_bound_m, bool) or not isinstance(relative_motion_bound_m, (int, float)) or not math.isfinite(relative_motion_bound_m) or relative_motion_bound_m < 0:
        raise ValueError('相对弹性与旋转余项界须为非负有限m')
    return audit_triangle_pressure_domain(source_triangles_world_m,
        distance_oracle=target_distance_oracle, solved_region_oracle=lambda triangle: False,
        relative_motion_bound=lambda triangle: relative_motion_bound_m,
        chord_error_m=0., maximum_depth=maximum_depth, maximum_regions=maximum_regions,
        tolerance_m=tolerance_m, progress=progress)


class ExternalFeWindowAudit:
    """无LU副本的外齿几何审核；整轮实际FE边界定位器按局部框架复用。"""
    def __init__(self, models, meshes, *, requested_tooth_offsets,
            external_tooth_offsets=(-3, 3), matching_tolerance_m=8e-7,
            tolerance_m=1e-10, engine='triangles'):
        import numpy as np
        from .contact_mesh import boundary_facets
        from .pressure_domain_audit import TriangleDistanceOracle
        from .pressure_domain_batch import _mesh_sha256, _brep_sha256
        if not isinstance(models, dict) or set(models) != {'pinion', 'wheel'} or not isinstance(meshes, dict) or set(meshes) != set(models):
            raise ValueError('外齿审核须提供同一实际双轮CAD和双轮FE网格')
        requested = tuple(requested_tooth_offsets); external = tuple(external_tooth_offsets)
        if not requested or not external or any(isinstance(k, bool) or not isinstance(k, int) for k in requested+external) or len(set(external)) != len(external) or set(requested) & set(external):
            raise ValueError('外齿必须是实际请求窗口外的不重复整数')
        if engine not in ('triangles', 'exact'):
            raise ValueError('外FE距离引擎须为triangles或exact')
        self.models = models; self.meshes = meshes; self.requested = requested; self.external = external
        self.matching_tolerance_m = matching_tolerance_m; self.tolerance_m = tolerance_m
        self.boundaries = {m: boundary_facets(mesh) for m, mesh in meshes.items()}
        target = self.boundaries['wheel']
        self.target_distance = TriangleDistanceOracle(np.asarray(target['nodes_m'])[target['triangles']],
            tolerance_m=tolerance_m, backend='vtk' if engine == 'triangles' else 'exact')
        self.mesh_sha256 = {m: _mesh_sha256(mesh) for m, mesh in meshes.items()}
        self.cad_sha256 = {m: _brep_sha256(model.shape) for m, model in models.items()}
        self._facets = {}; self._source_before = _source_sha256()

    def audit(self, result, domain_certificate, *, maximum_depth=6,
            maximum_regions=300000, progress=None):
        """绑定原状态和FE响应，给每个外齿真正的非负剩余距离下界。"""
        import numpy as np
        from OCP.BRep import BRep_Builder
        from OCP.TopoDS import TopoDS_Compound
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from .generated_halfspace import _faces
        from .generated_ltca import select_generated_flank_facets
        from .pressure_domain_audit import _loaded_wheel_rotation
        from .pressure_domain_batch import _mesh_sha256
        from .pressure_background import pressure_domain_fingerprint, _digest
        from .rigid_contact import _transform, RigidTransform
        if self._source_before != _source_sha256():
            raise ValueError('外窗口审核源码改变')
        if not result.get('converged') or domain_certificate.get('method') != 'full_gear_fe':
            raise ValueError('需要原已收敛FE结果及同状态整轮FE域证书')
        if sorted(int(k) for k in result.get('pair_force_shares',{})) != sorted(self.requested):
            raise ValueError('实际已求解齿对窗口与外齿审核请求不一致')
        expected = pressure_domain_fingerprint(result, method='full_gear_fe', tooth_offset=domain_certificate['tooth_offset'])
        if domain_certificate.get('state_fingerprint') != expected:
            raise ValueError('外窗口审核证书与原压力状态不同')
        _verify_global_fe_response_binding(result,domain_certificate)
        response = domain_certificate.get('fe_response_recomputed', {})
        checks = domain_certificate.get('checks', {})
        if (response.get('passed') is not True or response.get('input_node_geometry_bound') is not True or
                checks.get('included_candidates_kkt_passed') is not True or checks.get('audit_inputs_and_source_unchanged') is not True):
            raise ValueError('外窗口审核需要原完整FE响应与KKT的有效重算证据')
        cached_meshes = domain_certificate.get('batch_identity', {}).get('geometry', {}).get('fe_mesh_sha256')
        if cached_meshes is not None and cached_meshes != self.mesh_sha256:
            raise ValueError('外窗口FE网格与原响应证书不同')
        if any(_mesh_sha256(mesh) != self.mesh_sha256[m] for m, mesh in self.meshes.items()):
            raise ValueError('外窗口审核期间实际FE网格改变')
        frames = result['mesh_reference_frames']; first_rotation = np.asarray(frames['pinion_rotation_rows'])
        loaded = _loaded_wheel_rotation(frames, result['elastic_te_rad'])
        remainder = .5*response['wheel_radius_bound_m']*result['elastic_te_rad']**2
        motion = sum(response['maximum_displacement_m'].values())+remainder
        wheel = _transform(self.models['wheel'].shape, RigidTransform(rotation_rows=loaded.tolist()), .001)
        def target_distance(point):
            return self.target_distance(np.asarray(point) @ loaded)
        rows = []
        for offset in self.external:
            model = self.models['pinion']; tooth = (-offset) % model.request.tooth_count
            compound = TopoDS_Compound(); builder = BRep_Builder(); builder.MakeCompound(compound)
            for face in _faces(model.shape, model.slot_face_groups[tooth]['left_flank']):
                builder.Add(compound, face)
            source = _transform(compound, RigidTransform(rotation_rows=first_rotation.tolist()), .001)
            query = BRepExtrema_DistShapeShape(source, wheel)
            query.SetDeflection(self.tolerance_m*1000.); query.Perform()
            done = query.IsDone()
            cad_distance = float(query.Value())*.001 if done else None
            inner = bool(query.InnerSolution()) if done else None
            cad_margin = cad_distance-motion-self.tolerance_m if done and not inner else None
            if offset not in self._facets:
                self._facets[offset] = select_generated_flank_facets(self.meshes['pinion'], model.shell.surfaces['left_flank'],
                    tooth_count=model.request.tooth_count, tooth_index=tooth, boundary=self.boundaries['pinion'],
                    matching_tolerance_m=self.matching_tolerance_m)
            facets = self._facets[offset]
            triangles = np.asarray(facets['nodes_m'])[facets['triangles']] @ first_rotation.T
            if len(triangles):
                region = fe_surface_clearance_bound(triangles, target_distance, motion,
                    maximum_depth=maximum_depth, maximum_regions=maximum_regions,
                    tolerance_m=self.tolerance_m, progress=progress)
                fe_margin = region['minimum_excluded_clearance_margin_m'] if region['passed'] else None
            else:
                region = None; fe_margin = None
            passed = bool(cad_margin is not None and cad_margin > self.tolerance_m and
                fe_margin is not None and fe_margin > self.tolerance_m)
            rows.append(dict(tooth_offset=offset, passed=passed, cad_query_done=bool(done),
                cad_source_partly_inside_wheel=inner, cad_working_face_to_loaded_whole_wheel_distance_m=cad_distance,
                cad_remaining_clearance_lower_bound_m=cad_margin, fe_remaining_clearance_lower_bound_m=fe_margin,
                minimum_remaining_gap_m=min(cad_margin, fe_margin) if passed else None,
                relative_motion_bound_m=motion, source_fe_working_triangle_count=len(triangles), fe_region_audit=region))
        after = _source_sha256()
        if after != self._source_before:
            raise ArithmeticError('外窗口审核源码在计算中改变')
        return dict(schema='bevel_external_fe_geometric_window_v1', method='full_gear_fe',
            driving_angle_rad=result['driving_angle_rad'], pinion_torque_nm=result['pinion_torque_nm'],
            requested_tooth_offsets=list(self.requested), outer_offsets=list(self.external),
            passed=all(r['passed'] for r in rows), rows=rows,
            response_state_fingerprint=expected, source_domain_certificate_fingerprint=_digest(domain_certificate),
            source_identity=dict(source_sha256_before=self._source_before, source_sha256_after=after, unchanged=True),
            geometry_identity=dict(cad_brep_sha256=self.cad_sha256, fe_mesh_sha256=self.mesh_sha256),
            loaded_elastic_motion_bound_m=motion, rotation_linearization_remainder_m=remainder,
            scope='原FE全边界的三角距离下界＋实际CAD整工作面到加载整轮实体距离，扣同状态全轮弹性运动界；不证明所选窗内峰压收敛')


def audit_external_fe_window(models, meshes, result, domain_certificate, *,
        requested_tooth_offsets, external_tooth_offsets=(-3, 3), maximum_depth=6,
        maximum_regions=300000, engine='triangles', progress=None):
    """单状态便捷接口；周期场景复用ExternalFeWindowAudit类节省几何准备。"""
    context = ExternalFeWindowAudit(models, meshes, requested_tooth_offsets=requested_tooth_offsets,
        external_tooth_offsets=external_tooth_offsets, engine=engine)
    return context.audit(result, domain_certificate, maximum_depth=maximum_depth,
        maximum_regions=maximum_regions, progress=progress)
