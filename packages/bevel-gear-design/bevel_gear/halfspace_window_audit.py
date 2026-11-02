"""局部半空间的外邻齿排除；只适用于现有无跨齿弹性耦合模型。"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path


class ExternalHalfspaceWindowAudit:
    """复用同一实际双轮，以整片有限工作面到加载轮体的正分离排除外齿。"""
    def __init__(self, prepared, *, requested_tooth_offsets, external_tooth_offsets=None,
            tolerance_m=1e-10):
        from .generated_halfspace import PreparedGeneratedHalfspace
        from .pressure_domain_batch import _brep_sha256
        if not isinstance(prepared, PreparedGeneratedHalfspace) or prepared.tca.cad_models is None:
            raise ValueError('外齿排除需要同实际有限实体的半空间准备对象')
        requested = tuple(requested_tooth_offsets)
        if (not requested or len(set(requested)) != len(requested) or
                any(isinstance(k, bool) or not isinstance(k, int) for k in requested)):
            raise ValueError('请求齿对窗口须为不重复整数')
        outer = tuple(external_tooth_offsets) if external_tooth_offsets is not None else (min(requested)-1, max(requested)+1)
        if (not outer or len(set(outer)) != len(outer) or set(outer) & set(requested) or
                any(isinstance(k, bool) or not isinstance(k, int) for k in outer)):
            raise ValueError('外邻齿须为请求窗口之外的不重复整数')
        if isinstance(tolerance_m, bool) or not isinstance(tolerance_m, (int, float)) or not math.isfinite(tolerance_m) or tolerance_m <= 0:
            raise ValueError('实际有限面距离余量须为正有限m')
        self.prepared = prepared; self.models = prepared.tca.cad_models
        self.requested = requested; self.outer = outer; self.tolerance_m = tolerance_m
        self.geometry_identity = {m: _brep_sha256(model.shape) for m, model in self.models.items()}
        self.source_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def audit(self, result):
        """真实非负余隙才通过；失败、内部解及接触均不能当作零载证明。"""
        from OCP.BRep import BRep_Builder
        from OCP.TopoDS import TopoDS_Compound
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from .generated_halfspace import _faces
        from .generated_tca import generated_tca_frames
        from .pressure_domain_audit import _loaded_wheel_rotation
        from .pressure_domain_batch import _brep_sha256
        from .rigid_contact import _transform, RigidTransform
        if not result.get('converged') or result.get('active_contact_domain_accepted') is not True:
            raise ValueError('外邻齿排除只能使用已接受的半空间原结果')
        geometry = result['geometry']
        offsets = geometry.get('build_parameters', {}).get('tooth_offsets',
            [p['tooth_offset'] for p in geometry['patches']])
        if sorted(offsets) != sorted(self.requested):
            raise ValueError('实际已求解窗口与外齿排除请求不一致')
        if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != self.source_sha256:
            raise ValueError('外齿排除源码身份改变')
        if {m: _brep_sha256(model.shape) for m, model in self.models.items()} != self.geometry_identity:
            raise ValueError('外齿排除的实际CAD身份改变')
        frames = generated_tca_frames(self.prepared.tca, geometry['driving_angle_rad'],
            geometry['reference_follower_angle_rad'], tooth_offset=0)
        theta = result.get('linearization_step_theta_rad', result['theta_rad'])
        loaded = _loaded_wheel_rotation(frames, theta)
        wheel = _transform(self.models['wheel'].shape, RigidTransform(rotation_rows=loaded.tolist()), .001)
        rows = []
        for k in self.outer:
            model = self.models['pinion']; tooth = (-k) % model.request.tooth_count
            compound = TopoDS_Compound(); builder = BRep_Builder(); builder.MakeCompound(compound)
            for face in _faces(model.shape, model.slot_face_groups[tooth]['left_flank']):
                builder.Add(compound, face)
            first = _transform(compound, RigidTransform(rotation_rows=frames['pinion_rotation_rows']), .001)
            query = BRepExtrema_DistShapeShape(first, wheel)
            query.SetDeflection(self.tolerance_m/.001); query.Perform()
            done = bool(query.IsDone()); inner = bool(query.InnerSolution()) if done else None
            distance = float(query.Value())*.001 if done else None
            margin = distance-self.tolerance_m if done and not inner else None
            passed = margin is not None and margin > 0.
            rows.append(dict(tooth_offset=k, passed=bool(passed), minimum_remaining_gap_m=margin,
                evidence_kind='实际有限外齿工作面至加载整轮实体的数值距离，扣查询余量并拒绝内部解；局部半空间无跨齿耦合',
                cad_query_done=done, source_partly_inside_wheel=inner,
                rigid_distance_m=distance, distance_allowance_m=self.tolerance_m))
        return dict(schema='bevel_external_halfspace_window_v1', method='local_halfspace',
            passed=all(r['passed'] for r in rows), rows=rows,
            requested_tooth_offsets=list(self.requested), outer_offsets=list(self.outer),
            geometry_identity=self.geometry_identity, source_sha256=self.source_sha256,
            interpretation='外齿零试力的有限面几何排除，不代表有限齿体弹性或峰压精度验证')
