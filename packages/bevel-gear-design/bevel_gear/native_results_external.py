"""无原压力field的物理齿独立排除；不补造零压力数据。"""
import hashlib

from .distributed_pressure_domain import distributed_response_fingerprint,distributed_body_identity,distributed_loaded_pose


def native_field_partition(result):
    """区分原输入窗口、真正压力field和无field的空请求。"""
    import numpy as np
    nominal=list(result['requested_tooth_offsets']);selected=list(result['selected_tooth_offsets']);fields=result['fields'];represented=[f['tooth_offset'] for f in fields]
    if any(isinstance(k,bool) or not isinstance(k,int) for k in nominal+selected) or not nominal or not represented or len(set(nominal))!=len(nominal) or represented!=selected or len(set(represented))!=len(represented) or not set(represented)<=set(nominal):raise ValueError('原输入/selected/pressure field物理齿布局不一致')
    nominal=sorted(nominal)
    zero_fields=[]
    for field in fields:
        valid=np.asarray(field['valid_mask'])
        if valid.dtype!=np.bool_ or valid.ndim!=2 or not valid.any():raise ValueError('已表示field须含真实有效原格；不能以空field充当零压')
        force=np.asarray(field['forces_n'],dtype=float)
        if force.shape!=valid.shape or not np.isfinite(force).all() or np.any(force<0):raise ValueError('原field合力数组非法')
        if not np.any(force>0):zero_fields.append(field['tooth_offset'])
    missing=[k for k in nominal if k not in represented]
    for offset in missing:
        patches=[p for p in result['geometry']['patches'] if p['tooth_offset']==offset]
        if len(patches)!=1:raise ValueError('无field的原请求缺少唯一实际geometry patch')
        mask=np.asarray(patches[0].get('source_valid_mask',patches[0]['valid_mask']))
        if mask.dtype!=np.bool_ or mask.ndim!=2 or mask.any():raise ValueError('非空原请求缺pressure field，不能当作零压或外齿放行')
    return dict(nominal_requested_offsets=nominal,represented_field_offsets=represented,unrepresented_nominal_offsets=missing,zero_pressure_field_offsets=zero_fields)


def audit_native_unrepresented_teeth(prepared_surface, result, *, response_evidence,
        external_offsets=None, tolerance_m=1e-10):
    """完整外齿工作面到真实最终整轮CAD的正距离，扣同次两体全局u范数。

全局P1形函数为凸权，因此任意材料迹位移范数不超过所有体节点最大值。
此函数使用实际最终phi，额外刚体位移不再重复扣入total TE。刚体材料
内/ON或未分离时保留失败，不能据此判接触，也不能补零。
"""
    partition = native_field_partition(result)
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
    requested = partition['nominal_requested_offsets']
    represented = partition['represented_field_offsets']
    if external_offsets is None and requested != list(range(min(requested), max(requested)+1)):
        raise ValueError('外齿窗口不能跳过内部物理齿')
    offsets = [min(requested)-1, max(requested)+1] if external_offsets is None else list(external_offsets)
    if len(offsets) != len(set(offsets)) or any(isinstance(k, bool) or not isinstance(k, int) or k in represented for k in offsets):
        raise ValueError('独立排除对象须为无原pressure field的不重复整数齿偏移')
    models = getattr(prepared_surface, 'models', None) or prepared_surface.tca.cad_models
    teeth = models['pinion'].request.tooth_count
    if len({(-k)%teeth for k in offsets})!=len(offsets) or {(-k)%teeth for k in offsets}&{(-k)%teeth for k in represented}:raise ValueError('独立排除对象重复实际物理齿或已表示pressure field')
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
    return dict(schema='native_FE_unrepresented_physical_teeth_v1', method='full_gear_fe',
        contact_representation=P0_REPRESENTATION, requested_tooth_offsets=requested, external_tooth_offsets=offsets,
        represented_field_tooth_offsets=represented, unrepresented_nominal_offsets=partition['unrepresented_nominal_offsets'],
        rows=rows, passed=bool(unchanged and all(row['passed'] for row in rows)),
        full_response_fingerprint=before, body_response_identity=identity,
        maximum_body_displacements_m=maxima, actual_loaded_pose=pose,
        response_W_Ku_bound=True, actual_geometry_checked=True, original_state_unchanged=before == after,
        source_sha256_before=source_sha, source_sha256_after=source_after,
        strict_continuous_cad_proof=False,
        interpretation='几何正距与同次FE迹位移界排除外齿；P0平均KKT、逐Gauss和网格收敛分别验收')
