"""用专业有限BRep射线查询认证压力邻接排除域，不改变已求得的压力。"""
from __future__ import annotations

import hashlib
import math


class FiniteFaceRayOracle:
    """有限工作面集合的精确射线查询；公开长度m，OCP内部mm。"""
    def __init__(self,faces,*,tolerance_m=1e-10):
        import numpy as np
        from OCP.BRep import BRep_Builder
        from OCP.BRepBndLib import BRepBndLib
        from OCP.Bnd import Bnd_Box
        from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
        from OCP.TopoDS import TopoDS_Compound
        if isinstance(tolerance_m,bool) or not isinstance(tolerance_m,(int,float)) or not math.isfinite(tolerance_m) or tolerance_m<=0:
            raise ValueError('有限面射线容差须为正有限长度m')
        faces=tuple(faces)
        if not faces:raise ValueError('须给至少一个实际有限工作面')
        self.shape=TopoDS_Compound();builder=BRep_Builder();builder.MakeCompound(self.shape)
        for face in faces:builder.Add(self.shape,face)
        box=Bnd_Box();BRepBndLib.AddOptimal_s(self.shape,box,False,True)
        if box.IsVoid() or box.IsWhole():raise ValueError('实际工作面边界盒须为有限非空域')
        bounds=box.Get()
        if not all(math.isfinite(v) for v in bounds):raise ValueError('实际工作面边界盒含非有限坐标')
        self.corners_mm=np.array([[x,y,z] for x in (bounds[0],bounds[3]) for y in (bounds[1],bounds[4]) for z in (bounds[2],bounds[5])])
        self.engine=IntCurvesFace_ShapeIntersector();self.engine.Load(self.shape,tolerance_m*1000)
        self.tolerance_m=float(tolerance_m)

    def query(self,world_point_m,world_direction,*,rotation_rows=None):
        """查询一条整直线与全部选定有限工作面，零交点保留IsDone证据。

        rotation_rows把本体局部m坐标转至world。以所有实际面边界盒顶点
        到射线原点的最远距离选双向线段，故不会因猜测线长漏掉面交点。
        """
        import numpy as np
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.gp import gp_Lin,gp_Pnt,gp_Dir
        point=np.asarray(world_point_m,dtype=float);direction=np.asarray(world_direction,dtype=float)
        rotation=np.eye(3) if rotation_rows is None else np.asarray(rotation_rows,dtype=float)
        if point.shape!=(3,) or direction.shape!=(3,) or rotation.shape!=(3,3) or not np.isfinite(point).all() or not np.isfinite(direction).all() or not np.isfinite(rotation).all():
            raise ValueError('射线须给有限3D点/方向和3×3旋转矩阵')
        if not np.allclose(rotation.T@rotation,np.eye(3),rtol=0,atol=1e-10) or abs(float(np.linalg.det(rotation))-1)>1e-10:
            raise ValueError('射线换轴须为正交且行列式+1的旋转')
        length=float(np.linalg.norm(direction))
        if length<=1e-15:raise ValueError('射线方向不能为零')
        local=rotation.T@point*1000;normal=rotation.T@(direction/length)
        span=float(np.linalg.norm(self.corners_mm-local,axis=1).max())+max(1e-3,100*self.tolerance_m*1000)
        line=gp_Lin(gp_Pnt(*tuple(local)),gp_Dir(*tuple(normal)))
        self.engine.Perform(line,-span,span)
        done=bool(self.engine.IsDone());count=int(self.engine.NbPnt()) if done else None
        distance_m=None;distance_done=None
        if done and count==0:
            # 面内重合直线可能没有孤立交点；独立最短距离防止误判为空域。
            edge=BRepBuilderAPI_MakeEdge(line,-span,span).Edge()
            distance=BRepExtrema_DistShapeShape(edge,self.shape);distance.Perform();distance_done=bool(distance.IsDone())
            if distance_done:distance_m=float(distance.Value())*.001
        outside=bool(done and count==0 and distance_done and distance_m>self.tolerance_m)
        return dict(done=done,intersection_count=count,outside_selected_finite_faces=outside,
            finite_line_face_distance_m=distance_m,distance_query_done=distance_done,
            local_line_origin_m=(local*.001).tolist(),local_unit_direction=normal.tolist(),
            guaranteed_line_parameter_interval_m=[-span*.001,span*.001],tolerance_m=self.tolerance_m,
            intersection_points_local_m=[(np.asarray(self.engine.Pnt(i).Coord())*.001).tolist() for i in range(1,count+1)] if done else None,
            method='OCP IntCurvesFace_ShapeIntersector，全部实际有限工作面及真实裁切域')


def _field_identity(result):
    """绑定压力、力、有效/活跃掩码；后验分类不能改这些物理未知量。"""
    import numpy as np
    values={}
    for offset,field in zip(result['selected_tooth_offsets'],result['fields']):
        for name in ('pressure_pa','forces_n','valid_mask','contact_mask'):
            array=np.ascontiguousarray(field[name])
            values[str(offset)+':'+name]=hashlib.sha256(array.tobytes()).hexdigest()
    return values


def certify_halfspace_contact_boundaries(prepared,result,*,original_source_hashes=None,tolerance_m=1e-10):
    """对活跃邻接未解单元作有限BRep外域认证，保留原拒绝输入不修改。

    只接受每个被排除单元至少一个必需中心/四角射线对至少一体实际
    工作面IsDone且零交点的情形。任何已有交点而未解的单元仍拒绝。
    不改变间隙、压力、力臂、有效格或模型门限；认证报告保留原拒绝、
    原candidate场SHA及原源码身份。原始完整candidate由调用方归档。
    """
    import numpy as np
    from .generated_halfspace import _surface_sample
    if result.get('status')!='active_contact_reaches_unaccepted_domain_boundary' or result.get('converged'):
        return result
    candidate=result.get('candidate')
    if not candidate or not candidate.get('converged'):return result
    geometry=candidate['geometry']
    if not geometry['build_parameters']['check_cell_corners']:raise ValueError('有限域认证须采用中心及四角的原严格排除策略')
    before=_field_identity(candidate);patches={p['tooth_offset']:p for p in geometry['patches']}
    oracles={member:FiniteFaceRayOracle(trimmed.faces,tolerance_m=tolerance_m) for member,trimmed in prepared.trimmed.items()}
    checks=[];new_boundaries=[];binding_error=0.
    for offset,field,boundary in zip(candidate['selected_tooth_offsets'],candidate['fields'],candidate['contact_boundary_checks']):
        patch=patches[offset];origin=np.asarray(patch['plane_origin_world_m']);basis=np.asarray(patch['plane_tangent_rows']);normal=np.asarray(patch['plane_outward_normal_world'])
        valid_positions=np.argwhere(patch['valid_mask'])
        for number in np.linspace(0,len(valid_positions)-1,min(3,len(valid_positions)),dtype=int):
            j,i=valid_positions[number]
            for member,prefix in (('pinion','first'),('wheel','second')):
                point,*_=_surface_sample(prepared.trimmed[member].surface,patch[prefix+'_uv'][j,i])
                rotation=np.asarray(patch['frames'][member+'_rotation_rows'])
                binding_error=max(binding_error,float(np.linalg.norm(rotation@point-patch[prefix+'_points_world_m'][j,i])))
        if binding_error>1e-10:raise ValueError('有限域oracle的实际曲面与原候选压力格不一致')
        padded=np.pad(field['contact_mask'],1,constant_values=False)
        adjacent=padded[1:-1,:-2]|padded[1:-1,2:]|padded[:-2,1:-1]|padded[2:,1:-1]
        targets=np.argwhere(adjacent&patch['unresolved_query_rejected_mask']);grid=patch['grid'];patch_proven=True
        point_cache={}
        for j,i in targets:
            rays=[]
            for dy,dx in ((0.,0.),(-.5,-.5),(-.5,.5),(.5,-.5),(.5,.5)):
                key=(float(j+dy),float(i+dx))
                if key not in point_cache:
                    xy=np.array([grid.origin_xy_m[0]+(i+dx)*grid.spacing_x_m,grid.origin_xy_m[1]+(j+dy)*grid.spacing_y_m])
                    world=origin+basis.T@xy
                    point_cache[key]={member:oracle.query(world,normal,rotation_rows=patch['frames'][member+'_rotation_rows'])
                                      for member,oracle in oracles.items()}
                rays.append(dict(cell_offset=[dy,dx],members=point_cache[key]))
            proven=any(value['outside_selected_finite_faces'] for ray in rays for value in ray['members'].values())
            patch_proven=patch_proven and proven
            checks.append(dict(tooth_offset=offset,row=int(j),column=int(i),exclusion_proven=proven,required_point_queries=rays))
        updated=dict(boundary)
        if targets.size and patch_proven:
            updated['touches_unresolved_query_edge']=False;updated['touches_trimmed_domain_edge']=True
            updated['unresolved_edge_certified_as_finite_brep_outside']=True
        new_boundaries.append(updated)
    after=_field_identity(candidate)
    if before!=after:raise ArithmeticError('有限域后验认证意外改变物理场')
    accepted=not any(b['touches_window_edge'] or b['touches_model_validity_edge'] or b['touches_unresolved_query_edge'] or
                     b['touches_torque_direction_edge'] for b in new_boundaries)
    report=dict(accepted=bool(accepted),original_status=result['status'],original_candidate_converged=True,
        original_source_hashes=original_source_hashes,original_boundary_checks=candidate['contact_boundary_checks'],
        certified_boundary_checks=new_boundaries,field_sha256_before=before,field_sha256_after=after,
        original_candidate_fields_unchanged=True,original_input_not_mutated=True,source_surface_binding_error_m=binding_error,
        excluded_neighbor_cells_checked=len(checks),all_cell_exclusions_proven=all(c['exclusion_proven'] for c in checks),checks=checks,
        definition='仅证明原严格中心/四角规则下被排除单元属于实际有限BRep投影外域；不接受有交点的未知格，不改压力或有效格',
        raw_candidate_preservation='原拒绝输入不原地修改；调用方须保留原始candidate档案及其源码SHA')
    if accepted:
        restored=dict(candidate);restored.update(candidate=None,status='converged_with_finite_brep_boundary_certificate',
            contact_boundary_checks=new_boundaries,active_contact_domain_accepted=True,boundary_certification=report)
        return restored
    rejected=dict(result);rejected['boundary_certification']=report;return rejected
