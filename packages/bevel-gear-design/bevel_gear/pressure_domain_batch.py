"""批量压力域审核：缓存刚体不变的几何绑定，使用有弦差余量的实际面网格。

三角距离不是连续CAD严格证明。对面采用同一有限工作面的真实三角化，
查询下界显式扣除对面弦差；源面弦差仍由逐片审核扣除。FE还取原实际
对面FE三角距离之小，避免细CAD表面掩盖粗FE面几何偏差。
"""
from __future__ import annotations

import hashlib
import io
import math
from pathlib import Path


def _source_sha256():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def bounded_triangle_distance(distance_m, chord_error_m, tolerance_m=1e-10):
    """由真实三角网格距离给CAD距离下界；弦差只减不加。"""
    values = (distance_m, chord_error_m, tolerance_m)
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('距离、弦差和查询余量须为非负有限m')
    return max(0., distance_m-chord_error_m-tolerance_m)


def _brep_sha256(shape):
    """哈希真实BRep几何，不含可变的三角化缓存。"""
    from OCP.BRepTools import BRepTools
    from OCP.TopTools import TopTools_FormatVersion
    stream = io.BytesIO()
    BRepTools.Write_s(shape, stream, False, False, TopTools_FormatVersion.TopTools_FormatVersion_VERSION_1)
    return hashlib.sha256(stream.getvalue()).hexdigest()


def _mesh_sha256(mesh):
    import numpy as np
    digest = hashlib.sha256()
    for name in ('nodes_m', 'tetrahedra', 'boundary_triangles', 'boundary_tags'):
        values = np.asarray(getattr(mesh, name))
        digest.update(str(values.dtype).encode()); digest.update(str(values.shape).encode()); digest.update(values.tobytes())
    return digest.hexdigest()


class PressureDomainBatchAudit:
    """同准备对象下复用几何绑定与单相位FE响应；audit不重求压力。

    engine='triangles' 按需需要VTK；engine='brep'使用原OCP有限面距离，
    保留较慢的对照路线。构造或导入本模块均不会自动导入VTK。
    """
    def __init__(self, prepared, flank_mesh, pinion_request, *, method,
            opposing_deflection_m=1e-6, engine='triangles', matching_tolerance_m=8e-7,
            tolerance_m=1e-10, tooth_offsets=(-2,-1,0,1,2)):
        import numpy as np
        from .generated_surface import _request
        from .generated_ltca import PreparedGeneratedLtca
        from .generated_halfspace import PreparedGeneratedHalfspace
        from .refined_fe_plotting import _background
        from .pressure_domain_audit import _selected_triangles, _kernel_sha256
        from .pressure_background import _digest
        if method not in ('full_gear_fe', 'local_halfspace') or engine not in ('triangles', 'brep'):
            raise ValueError('方法或压力域几何引擎不支持')
        for value in (opposing_deflection_m, matching_tolerance_m, tolerance_m):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError('压力域几何误差须为正有限m')
        request = _request(pinion_request)
        if method == 'full_gear_fe':
            if not isinstance(prepared, PreparedGeneratedLtca) or prepared.models['pinion'].request != request:
                raise ValueError('FE准备对象与同制造输入不一致')
            models = prepared.models
        else:
            if not isinstance(prepared, PreparedGeneratedHalfspace) or prepared.tca.requests['pinion'] != request:
                raise ValueError('半解析准备对象与同制造输入不一致')
            models = prepared.tca.cad_models
            if models is None:
                raise ValueError('批量覆盖审核须复用实际有限实体')
        selected, background, _ = _background(flank_mesh, pinion_request)
        self.prepared = prepared; self.models = models; self.flank_mesh = flank_mesh
        self.request = request; self.method = method; self.engine = engine
        self.opposing_deflection_m = opposing_deflection_m
        self.matching_tolerance_m = matching_tolerance_m; self.tolerance_m = tolerance_m
        self.background = background
        self.triangles = _selected_triangles(flank_mesh, selected['triangle_indices'])
        self.deflection = flank_mesh.get('metadata', {}).get('linear_deflection_m')
        if not isinstance(self.deflection, (int, float)) or self.deflection <= 0:
            raise ValueError('源CAD三角背景必须保存实际弦差')
        self.background_sha = _digest(background)
        self._own_source = _source_sha256(); self._core_source = _kernel_sha256()
        self._geometry = self._geometry_identity()
        self._before_cache_brep_sha256 = dict(self._geometry['brep'])
        self._bindings = {}; self._pinion_oracles = {}; self._facets = {}; self._cad_oracles = {}; self._fe_oracles = {}
        self._last_response_key = None; self._last_response = None
        self._response_count = 0; self._binding_count = 0
        self._point_count = len(np.unique(self.triangles.reshape(-1, 3), axis=0))
        offsets = tuple(tooth_offsets)
        if not offsets or len(set(offsets)) != len(offsets) or any(isinstance(k, bool) or not isinstance(k, int) for k in offsets):
            raise ValueError('批量缓存齿对窗口须为不重复整数')
        self.tooth_offsets = offsets
        # OCP网格缓存和面检查可改变BRep的辅助标志。先完成所有缓存准备，
        # 再冻结实际几何身份；后续相位仅作刚体映射，不能持续接受身份变化。
        for offset in offsets:
            self._cad_oracle(offset, np.eye(3))
            self._binding(offset)
        self._geometry = self._geometry_identity()
        for binding in self._bindings.values():
            binding['geometry_fingerprint'] = self._geometry['fingerprint']

    def _geometry_identity(self):
        from .pressure_background import _digest
        identity = dict(brep={member: _brep_sha256(model.shape) for member, model in self.models.items()},
            background_fingerprint=self.background_sha, source_deflection_m=self.deflection,
            opposing_deflection_m=self.opposing_deflection_m, engine=self.engine)
        if self.method == 'full_gear_fe':
            identity['fe_mesh_sha256'] = {m: _mesh_sha256(system.mesh) for m, system in self.prepared.systems.items()}
            identity['fe_material_and_support'] = {m: _digest(dict(material=system.material,
                fixed_dofs=system.fixed_dofs)) for m, system in self.prepared.systems.items()}
        identity['fingerprint'] = _digest(identity)
        return identity

    def _faces(self, member, offset):
        from .generated_halfspace import _faces
        if self.method == 'local_halfspace':
            return self.prepared.trimmed[member].faces
        model = self.models[member]
        tooth = (-offset if member == 'pinion' else offset) % model.request.tooth_count
        return _faces(model.shape, model.slot_face_groups[tooth]['left_flank' if member == 'pinion' else 'right_flank'])

    def _binding(self, offset):
        """只在材料局部框架检查一次；跨相位刚体旋转保持所有欧氏距离。"""
        import numpy as np
        from .pressure_domain_audit import FiniteFaceDistanceOracle, _periodic_pinion_triangles
        if self.method == 'local_halfspace' and self._bindings and offset not in self._bindings:
            reference = next(iter(self._bindings))
            self._bindings[offset] = dict(self._bindings[reference])
            self._pinion_oracles[offset] = self._pinion_oracles[reference]
            return self._bindings[offset]
        if offset not in self._bindings:
            local = self.triangles
            if self.method == 'full_gear_fe':
                local, _, _ = _periodic_pinion_triangles(local, offset,
                    self.request.tooth_count, self.models['wheel'].request.tooth_count)
            oracle = FiniteFaceDistanceOracle(self._faces('pinion', offset), np.eye(3).tolist(), tolerance_m=self.tolerance_m)
            self._pinion_oracles[offset] = oracle
            queries = [oracle(p) for p in np.unique(local.reshape(-1, 3), axis=0)]
            passed = all(q['done'] and q['distance_m'] <= self.deflection+self.tolerance_m for q in queries)
            self._bindings[offset] = dict(passed=bool(passed), checked_vertices=len(queries),
                maximum_vertex_distance_m=max((q['distance_m'] for q in queries if q['done']), default=None),
                invariance='同一材料齿面与背景的刚体等距映射，非跨齿几何猜测',
                geometry_fingerprint=self._geometry['fingerprint'], source_deflection_m=self.deflection)
            self._binding_count += 1
        return self._bindings[offset]

    def _cad_oracle(self, offset, rotation):
        """有限对面局部网格只建一次；查询点逆刚体变换，无插值距离场。"""
        import numpy as np
        from .pressure_domain_audit import FiniteFaceDistanceOracle, TriangleDistanceOracle
        key = offset if self.method == 'full_gear_fe' else 0
        if key not in self._cad_oracles:
            faces = self._faces('wheel', offset)
            if self.engine == 'brep':
                oracle = FiniteFaceDistanceOracle(faces, np.eye(3).tolist(), tolerance_m=self.tolerance_m)
                evidence = dict(engine='OCP_finite_BRep_distance', opposing_chord_allowance_m=0.)
            else:
                from OCP.BRep import BRep_Builder
                from OCP.TopoDS import TopoDS_Compound
                from .cad_import import ImportedCadModel, triangulate_imported_model
                from .pressure_background import _digest
                compound = TopoDS_Compound(); builder = BRep_Builder(); builder.MakeCompound(compound)
                for face in faces:
                    builder.Add(compound, face)
                mesh = triangulate_imported_model(ImportedCadModel(compound, {}),
                    linear_deflection_m=self.opposing_deflection_m, angular_deflection_rad=.08)
                triangles = np.asarray(mesh['vertices_m'])[np.asarray(mesh['triangles'], dtype=int)]
                oracle = TriangleDistanceOracle(triangles, tolerance_m=self.tolerance_m+self.opposing_deflection_m, backend='vtk')
                evidence = dict(engine='actual_finite_CAD_triangles_VTK', opposing_chord_allowance_m=self.opposing_deflection_m,
                    triangle_count=len(triangles), mesh_fingerprint=_digest(mesh),
                    scope='实际裁切有限面的三角离散；显式扣对面弦差，不宣称连续CAD严格区间界')
            self._cad_oracles[key] = (oracle, evidence)
        oracle, evidence = self._cad_oracles[key]
        frame = np.asarray(rotation)
        return lambda point: oracle(np.asarray(point) @ frame), evidence

    def _flank_facets(self, member, offset):
        from .generated_ltca import select_generated_flank_facets
        key = member, offset
        if key not in self._facets:
            model = self.models[member]
            self._facets[key] = select_generated_flank_facets(self.prepared.systems[member].mesh,
                model.shell.surfaces['left_flank' if member == 'pinion' else 'right_flank'],
                tooth_count=model.request.tooth_count,
                tooth_index=(-offset if member == 'pinion' else offset) % model.request.tooth_count,
                matching_tolerance_m=self.matching_tolerance_m)
        return self._facets[key]

    def _fe_evidence(self, result):
        from .pressure_background import _digest
        from .pressure_domain_audit import _fe_displacement_evidence
        key = _digest(dict(nodal_contacts=result['nodal_contacts'], elastic_te_rad=result['elastic_te_rad'],
            frames=result['mesh_reference_frames'], pinion_torque_nm=result['pinion_torque_nm']))
        if key != self._last_response_key:
            self._last_response = _fe_displacement_evidence(self.prepared, result,
                gap_tolerance_m=max(2e-10, self.tolerance_m*1.01))
            self._last_response_key = key; self._response_count += 1
        return self._last_response, key

    def audit(self, result, *, tooth_offset=0, maximum_depth=2, maximum_regions=300000, progress=None):
        """逐原片证明覆盖或排除；缓存只节省重复验证，不直接制造排除片。"""
        import numpy as np
        from .pressure_domain_audit import (_kernel_sha256, _finalize_certificate, _loaded_wheel_rotation,
            _periodic_pinion_triangles, _triangle_union_cover, TriangleDistanceOracle,
            _grid_cover, _halfspace_motion, _kkt_halfspace, audit_triangle_pressure_domain)
        from .pressure_background import DOMAIN_CERTIFICATE_SCHEMA, pressure_domain_fingerprint, pressure_domain_background
        if isinstance(tooth_offset, bool) or not isinstance(tooth_offset, int):
            raise ValueError('实际齿偏移须为整数')
        if tooth_offset not in self.tooth_offsets:
            raise ValueError('实际齿不属于本批量缓存的显式请求窗口')
        if not result.get('converged'):
            raise ValueError('只审核已收敛原压力结果')
        own_before = _source_sha256(); core_before = _kernel_sha256()
        if own_before != self._own_source or core_before != self._core_source:
            raise ValueError('批量审核的源码身份改变')
        actual_geometry = self._geometry_identity()
        if actual_geometry != self._geometry:
            changed = [key for key in actual_geometry if actual_geometry[key] != self._geometry.get(key)]
            raise ValueError('批量审核的实际CAD、材料或FE网格身份改变：'+','.join(changed))
        state_before = pressure_domain_fingerprint(result, method=self.method, tooth_offset=tooth_offset)
        binding = self._binding(tooth_offset)
        common = dict(schema=DOMAIN_CERTIFICATE_SCHEMA, method=self.method, tooth_offset=tooth_offset,
            state_fingerprint=state_before, flank_background_fingerprint=self.background_sha)
        if self.method == 'full_gear_fe':
            frames = result['mesh_reference_frames']; rotation = np.asarray(frames['pinion_rotation_rows'])
            physical, pinion_index, wheel_index = _periodic_pinion_triangles(self.triangles, tooth_offset,
                self.request.tooth_count, self.models['wheel'].request.tooth_count)
            world = physical @ rotation.T
            response, response_key = self._fe_evidence(result)
            first = self._flank_facets('pinion', tooth_offset)
            candidates = {r['first_pressure_node_index'] for r in result['nodal_contacts'] if r['tooth_offset'] == tooth_offset}
            facet_ids = [tri for tri in first['triangles'] if all(int(i) in candidates for i in tri)]
            if facet_ids:
                covered = np.asarray(first['nodes_m'])[facet_ids] @ rotation.T
                normals = np.asarray([r['normal_world'] for r in result['nodal_contacts'] if r['tooth_offset'] == tooth_offset])
                normal = normals.mean(axis=0); normal /= np.linalg.norm(normal)
                tangent = np.cross(normal, [1., 0., 0.] if abs(normal[0]) < .9 else [0., 1., 0.]); tangent /= np.linalg.norm(tangent)
                cover = _triangle_union_cover(covered, covered[0, 0], np.array([tangent, np.cross(normal, tangent)]), padding_m=self.deflection)
            else:
                # 没有完整候选三角形时没有任何原已知近域，必须逐片真实排除。
                cover = lambda triangle: False
            loaded = _loaded_wheel_rotation(frames, result['elastic_te_rad'])
            cad, opposing = self._cad_oracle(tooth_offset, loaded)
            if tooth_offset not in self._fe_oracles:
                second = self._flank_facets('wheel', tooth_offset)
                triangles = np.asarray(second['nodes_m'])[second['triangles']]
                self._fe_oracles[tooth_offset] = TriangleDistanceOracle(triangles,
                    tolerance_m=self.tolerance_m, backend='vtk' if self.engine == 'triangles' else 'auto')
            fe_oracle = self._fe_oracles[tooth_offset]
            def distance(point):
                a = cad(point); b = fe_oracle(np.asarray(point) @ loaded)
                done = bool(a['done'] and b['done'])
                return dict(done=done, distance_m=min(a['distance_m'], b['distance_m']) if done else None,
                    error_bound_m=max(a['error_bound_m'], b['error_bound_m']))
            remainder = .5*response['wheel_radius_bound_m']*result['elastic_te_rad']**2
            motion_bound = sum(response['maximum_displacement_m'].values())+remainder
            regions = audit_triangle_pressure_domain(world, distance_oracle=distance, solved_region_oracle=cover,
                relative_motion_bound=lambda triangle: motion_bound, chord_error_m=self.deflection,
                maximum_depth=maximum_depth, maximum_regions=maximum_regions, tolerance_m=self.tolerance_m, progress=progress)
            complete = regions['passed']
            checks = dict(input_geometry_bound=bool(binding['passed'] and response['input_node_geometry_bound']),
                included_candidates_kkt_passed=response['passed'], omitted_candidates_excluded=complete,
                unresolved_regions_excluded=complete, far_field_motion_bound_passed=complete)
            common.update(driving_angle_rad=result['driving_angle_rad'], pinion_torque_nm=result['pinion_torque_nm'],
                fe_response_recomputed=response, fe_response_state_fingerprint=response_key,
                original_candidate_triangle_count=len(facet_ids), triangulation_deflection_allowance_m=self.deflection,
                physical_tooth_indices=dict(pinion=pinion_index, wheel=wheel_index),
                loaded_elastic_motion_bound_m=motion_bound, rotation_linearization_remainder_m=remainder,
                opposing_distance_definition='同一加载位姿中，取有弦差余量的实际有限CAD三角距与原FE对面三角距之小')
        else:
            patch = next(p for p in result['geometry']['patches'] if p['tooth_offset'] == tooth_offset)
            rotation = np.asarray(patch['frames']['pinion_rotation_rows'])
            world = self.triangles @ rotation.T
            theta = float(result.get('linearization_step_theta_rad', result['theta_rad']))
            loaded = _loaded_wheel_rotation(patch['frames'], theta)
            distance, opposing = self._cad_oracle(tooth_offset, loaded)
            if tooth_offset in result['selected_tooth_offsets']:
                offset = result['selected_tooth_offsets'].index(tooth_offset)
                field = result['fields'][offset]
                points = np.asarray(patch['first_points_world_m'], dtype=float)[np.asarray(patch['valid_mask'], dtype=bool)]
                probes = points[np.linspace(0, len(points)-1, min(25, len(points)), dtype=int)]
                pressure_queries = [self._pinion_oracles[tooth_offset](p @ rotation) for p in probes]
                pressure_bound = bool(len(probes)) and all(q['done'] and q['distance_m'] <= max(self.tolerance_m, 1e-9) for q in pressure_queries)
                motion, motion_evidence = _halfspace_motion(result, tooth_offset, include_rigid_motion=False)
                cover = _grid_cover(patch, field, padding_m=self.deflection)
            else:
                # 本模型各齿半空间彼此独立，未建有效格的齿无自身力源；仍须
                # 对整个实际有限齿面逐片排除，缺失格不构成已知零压域。
                probes = []; pressure_bound = True
                motion = lambda triangle: 0.
                cover = lambda triangle: False
                motion_evidence = dict(source_force_n=0., rectangular_source_count=0,
                    rigid_motion_already_applied=True, absent_selected_field=True,
                    derivation='局部半空间没有跨齿弹性耦合；零自身力源产生零自身位移，实际共同加载刚体位姿另行查询')
                from OCP.BRepClass3d import BRepClass3d_SolidClassifier
                from OCP.TopAbs import TopAbs_OUT
                from OCP.gp import gp_Pnt
                classifier = BRepClass3d_SolidClassifier(self.models['wheel'].shape)
                unsigned = distance
                def distance(point):
                    classifier.Perform(gp_Pnt(*(np.asarray(point) @ loaded*1000.)), self.tolerance_m*1000.)
                    if classifier.State() != TopAbs_OUT:
                        return dict(done=False, distance_m=None, error_bound_m=self.tolerance_m)
                    return unsigned(point)
            regions = audit_triangle_pressure_domain(world, distance_oracle=distance,
                solved_region_oracle=cover,
                relative_motion_bound=motion, chord_error_m=self.deflection,
                maximum_depth=maximum_depth, maximum_regions=maximum_regions, tolerance_m=self.tolerance_m, progress=progress)
            kkt = _kkt_halfspace(result, max(2e-10, float(result.get('certificate', {}).get('gap_acceptance_m', self.tolerance_m))*1.01))
            boundaries = result.get('contact_boundary_checks', [])
            boundary_ok = bool(boundaries) and all(not b.get(name, True) for b in boundaries for name in
                ('touches_window_edge', 'touches_model_validity_edge', 'touches_unresolved_query_edge', 'touches_torque_direction_edge'))
            complete = regions['passed']
            checks = dict(input_geometry_bound=bool(binding['passed'] and pressure_bound), full_working_flank_covered=complete,
                omitted_regions_excluded=complete, unresolved_regions_excluded=complete,
                loaded_boundary_checks_passed=bool(boundary_ok and kkt['passed']))
            common.update(driving_angle_rad=result['geometry']['driving_angle_rad'], pinion_torque_nm=result['geometry']['pinion_torque_nm'],
                kkt_recomputed=kkt, motion_bound=motion_evidence,
                geometry_binding=dict(**binding, pressure_points_checked=len(probes), pressure_points_bound=pressure_bound,
                    triangulation_deflection_allowance_m=self.deflection))
        common.update(passed=all(checks.values()), complete_working_flank=complete, checks=checks,
            unresolved_region_count=regions['unresolved_region_count'], regions=regions,
            loaded_opposing_rotation_rows=loaded.tolist(),
            evidence_kind='同一实际有限工作面三角距离＋显式双面弦差余量＋原求解域/KKT；刚体不变几何绑定与同相位FE响应缓存',
            limitations=['只审核本相位和实际齿的离散背景，不认证时间连续极值或峰压网格收敛',
                'CAD弦差为显式离散误差余量，不是NURBS连续区间严格证明', '压力和物理验证分别报告'])
        own_after = _source_sha256()
        common['batch_identity'] = dict(source_sha256_before=own_before, source_sha256_after=own_after,
            unchanged=own_before == own_after, geometry=self._geometry, binding=binding,
            cache_preparation_brep_sha256_before=self._before_cache_brep_sha256,
            opposing_geometry=opposing, geometry_binding_count=self._binding_count,
            fe_response_reconstruction_count=self._response_count,
            cache_scope='同实际CAD/背景/FE网格；几何刚体等距复用，同原状态完整力源响应复用')
        if own_before != own_after:
            raise ArithmeticError('批量审核源码在计算中改变')
        common = _finalize_certificate(common, result, state_before=state_before, source_before=core_before)
        pressure_domain_background(result, self.background, method=self.method, tooth_offset=tooth_offset,
            certificate=common, pinion_request=self.request)
        return common
