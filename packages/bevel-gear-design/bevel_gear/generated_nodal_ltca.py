"""实际展成整轮FE的高分辨率节点接触：真实体柔度，不调用半空间核。"""
from __future__ import annotations
from dataclasses import asdict
import math


def generated_refinement_tracks(prepared_tca,*,driving_angle_rad=0.,tooth_offsets=(-1,0,1),
        half_widths_m=(.006,.0018),track_spacing_m=.00015,sample_spacing_m=.0005):
    """由实际NURBS参数导数生成面上加密带；它只控制网格，不定义接触域。"""
    import numpy as np
    from OCP.gp import gp_Pnt,gp_Vec
    from .generated_tca import solve_generated_tca
    from .generated_ltca import _rz
    for value in (*half_widths_m,track_spacing_m,sample_spacing_m):
        if isinstance(value,bool) or not math.isfinite(value) or value<=0:raise ValueError('加密带尺寸须为正m')
    reference=solve_generated_tca(prepared_tca,driving_angle_rad,tooth_offset=0)
    if not reference['converged']:raise ArithmeticError('加密参考齿须实际相切成功')
    tracks={m:[] for m in ('pinion','wheel')};records=[]
    for offset in tooth_offsets:
        sample=solve_generated_tca(prepared_tca,driving_angle_rad,tooth_offset=offset)
        seed=sample if sample['converged'] else reference
        for member,side in (('pinion','left_flank'),('wheel','right_flank')):
            surface=prepared_tca.shells[member].surfaces[side]
            # prepared_tca的surfaces为本轮工作面，参数含两轮独立u/v。
            uv=np.asarray(seed[member+'_uv'])
            p,du,dv=gp_Pnt(),gp_Vec(),gp_Vec();surface.D1(*uv,p,du,dv)
            jac=np.column_stack((du.Coord(),dv.Coord()))*.001
            e1=jac[:,0]/np.linalg.norm(jac[:,0]);e2=jac[:,1]-e1*(e1@jac[:,1]);e2/=np.linalg.norm(e2)
            inverse=np.linalg.pinv(jac)@np.column_stack((e1,e2));bounds=surface.Bounds()
            rotation=_rz(((-offset if member=='pinion' else offset)%prepared_tca.requests[member].tooth_count)*2*math.pi/prepared_tca.requests[member].tooth_count,np)
            xs=np.linspace(-half_widths_m[0],half_widths_m[0],math.ceil(2*half_widths_m[0]/sample_spacing_m)+1)
            ys=np.linspace(-half_widths_m[1],half_widths_m[1],math.ceil(2*half_widths_m[1]/track_spacing_m)+1)
            for y in ys:
                track=[]
                for x in xs:
                    u,v=uv+inverse@np.array([x,y])
                    if bounds[0]<=u<=bounds[1] and bounds[2]<=v<=bounds[3]:
                        point=rotation@(np.asarray(surface.Value(float(u),float(v)).Coord())*.001)
                        if not track or np.linalg.norm(point-track[-1])>1e-12:track.append(point.tolist())
                if len(track)>=2:tracks[member].append(track)
            records.append(dict(member=member,tooth_offset=offset,tca_seed_converged=sample['converged'],
                role='只用于真实面上的网格控制；正式接触节点须另行几何筛选'))
    return dict(tracks=tracks,diagnostics=records,half_widths_m=list(half_widths_m),track_spacing_m=track_spacing_m)


def solve_generated_nodal_fe_ltca(prepared,driving_angle_rad,pinion_torque_nm,*,tooth_offsets=(-1,0,1),
        maximum_distance_m=.0002,maximum_normal_angle_rad=.3,matching_tolerance_m=2e-7,
        mesh_phase_tolerance_m=1e-10,maximum_mesh_phase_iterations=15,maximum_candidates=30000,
        maximum_loaded_angle_rad=.005,maximum_pinion_torque_relative_residual=.02,
        maximum_working_set=2048,include_fields=False,progress=None):
    """独立节点变量、实际面配对、全轮FE列生成KKT；输出真实对偶压力格。"""
    import numpy as np
    from .generated_ltca import PreparedGeneratedLtca,select_generated_flank_facets
    from .generated_tca import solve_generated_tca,generated_tca_frames
    from .contact_mesh import boundary_facets
    from .nodal_contact import pair_nodal_contacts,nodal_pressure_cells
    from .elasticity import ContactLoadChannel,contact_load_matrix
    from .fe_contact_solver import solve_fe_contact_equilibrium
    if not isinstance(prepared,PreparedGeneratedLtca):raise TypeError('须传入真实整轮FE准备对象')
    if not math.isfinite(pinion_torque_nm) or pinion_torque_nm==0:raise ValueError('载荷须为非零有符号Nm')
    offsets=tuple(tooth_offsets);direction=1 if pinion_torque_nm>0 else -1
    key=('nodal',driving_angle_rad,direction,offsets,maximum_distance_m,maximum_normal_angle_rad,matching_tolerance_m,mesh_phase_tolerance_m)
    if key not in prepared._contact_cache:
        tca={k:solve_generated_tca(prepared.tca,driving_angle_rad,tooth_offset=k) for k in offsets}
        valid=[r for r in tca.values() if r['converged']]
        if not valid:raise ArithmeticError('没有有效无载相切姿态')
        reference=min(valid,key=lambda r:direction*r['follower_increment_angle_rad']);follower=reference['follower_increment_angle_rad']
        facets={}
        for member,side in (('pinion','left_flank'),('wheel','right_flank')):
            model=prepared.models[member];mesh=prepared.systems[member].mesh
            boundary=boundary_facets(mesh)
            for k in offsets:
                tooth=(-k if member=='pinion' else k)%model.request.tooth_count
                facets[member,k]=select_generated_flank_facets(mesh,model.shell.surfaces[side],
                    tooth_count=model.request.tooth_count,tooth_index=tooth,boundary=boundary,matching_tolerance_m=matching_tolerance_m)
        correction=0.;history=[]
        for _ in range(maximum_mesh_phase_iterations):
            frames=generated_tca_frames(prepared.tca,driving_angle_rad,follower+correction)
            records=[]
            for k in offsets:
                pairs=pair_nodal_contacts(facets['pinion',k],facets['wheel',k],frames['pinion_rotation_rows'],frames['wheel_rotation_rows'],
                    maximum_distance_m=maximum_distance_m,maximum_normal_angle_rad=maximum_normal_angle_rad)
                for pair in pairs:
                    normal=np.asarray(pair['normal_world'])
                    ell1=float(normal@np.cross(frames['pinion_axis_world'],pair['first_point_world_m']))
                    ell2=float(-normal@np.cross(frames['wheel_axis_world'],pair['second_point_world_m']))
                    if direction*ell1>1e-10 and direction*ell2>1e-10:
                        records.append(dict(pair,tooth_offset=k,pinion_lever_m=ell1,wheel_lever_m=ell2))
            if not records:raise ArithmeticError('没有可承受指定有符号载荷的真实节点配对')
            if len(records)>maximum_candidates:raise ValueError(f'真实节点候选{len(records)}超过显式上限{maximum_candidates}')
            if len({r['first_pressure_node_index'] for r in records})!=len(records):raise ValueError('第一侧节点重复，拒绝冗余接触变量')
            gaps=np.asarray([r['signed_gap_m'] for r in records]);lever=np.asarray([r['wheel_lever_m'] for r in records])
            history.append(dict(correction_rad=correction,minimum_gap_m=float(gaps.min()),candidates=len(records)))
            if abs(gaps.min())<=mesh_phase_tolerance_m:break
            correction+=direction*float(np.min(gaps/(direction*lever)))
            if abs(correction)>.02:raise ArithmeticError('网格初相位修正超过0.02rad')
        else:raise ArithmeticError('节点接触网格初相位未收敛')
        matrices={}
        for member,prefix in (('pinion','first'),('wheel','second')):
            channels=[ContactLoadChannel(tuple(r[prefix+'_node_indices']),tuple(r[prefix+'_weights']),tuple(r[prefix+'_force_normal_local'])) for r in records]
            matrices[member]=contact_load_matrix(prepared.systems[member],channels,matrix_format='csc',maximum_sparse_channels=maximum_candidates)
        prepared._contact_cache[key]=dict(records=records,matrices=matrices,gaps=gaps,lever=lever,frames=frames,
            correction=correction,follower=follower,tca=tca,reference=reference,history=history,facets=facets)
    cache=prepared._contact_cache[key];records=cache['records'];count=len(records)
    free={m:matrix[prepared.systems[m].free_dofs] for m,matrix in cache['matrices'].items()}
    saved_columns=cache.setdefault('compliance_columns',{})
    def columns(indices):
        indices=list(map(int,indices));missing=[i for i in indices if i not in saved_columns]
        if missing:
            value=np.zeros((count,len(missing)))
            for member,w in free.items():
                system=prepared.systems[member]
                for start in range(0,len(missing),24):
                    selected=missing[start:start+24]
                    displacement=system.factorization.solve(w[:,selected].toarray())
                    value[:,start:start+len(selected)]+=w.T@displacement
            for j,index in enumerate(missing):saved_columns[index]=value[:,j].copy()
        return np.column_stack([saved_columns[i] for i in indices])
    ratio=prepared.tca.requests['pinion'].tooth_count/prepared.tca.requests['wheel'].tooth_count
    target=pinion_torque_nm/ratio
    if progress:progress(dict(stage='nodal_qp',candidates=count))
    equilibrium=solve_fe_contact_equilibrium(cache['gaps'],lever_arms_m=cache['lever'],torque_nm=target,
        compliance_columns=columns,maximum_working_set=maximum_working_set,
        allow_initial_penetration=bool(cache['gaps'].min()<0),gap_tolerance_m=1e-10)
    result=dict(converged=False,status=equilibrium['status'],contact_representation='nodal_dual_cells',contacts=[],
        model='真实整轮三维线弹性FE+节点—面集总接触；不含半空间柔度',
        driving_angle_rad=driving_angle_rad,pinion_torque_nm=pinion_torque_nm,wheel_target_torque_nm=target,
        preparation=prepared.diagnostics,equilibrium=equilibrium,mesh_reference_frames=cache['frames'],
        mesh_unloaded_correction_rad=cache['correction'],mesh_phase_history=cache['history'],candidate_count=count,
        tca={str(k):r for k,r in cache['tca'].items()},pressure_peak_validated=False,discretization_converged=None)
    if not equilibrium['converged']:return result
    forces=np.asarray(equilibrium['forces_n']);theta=equilibrium['theta_rad']
    pinion_reconstructed=float(sum(r['pinion_lever_m']*f for r,f in zip(records,forces)))
    residual=(pinion_reconstructed-pinion_torque_nm)/abs(pinion_torque_nm)
    if abs(theta)>maximum_loaded_angle_rad or abs(residual)>maximum_pinion_torque_relative_residual:
        result.update(status='nodal_rotation_or_pinion_torque_limit',candidate=dict(equilibrium=equilibrium,pinion_relative_residual=residual));return result
    active=np.asarray(equilibrium['active_indices'],dtype=int);nodal_rows=[];maps={k:{} for k in offsets}
    active_set=set(map(int,active))
    for index,(r,f) in enumerate(zip(records,forces)):
        value=dict(r,force_n=float(f),pressure_pa=float(f/r['area_m2']),active=index in active_set,
            remaining_gap_m=float(equilibrium['remaining_gaps_m'][index]))
        nodal_rows.append(value);maps[r['tooth_offset']][r['first_pressure_node_index']]=value
    cells=[]
    for k in offsets:cells.extend(nodal_pressure_cells(cache['facets']['pinion',k],maps[k],cache['frames']['pinion_rotation_rows'],tooth_offset=k))
    if not math.isclose(sum(c['force_n'] for c in cells),forces.sum(),rel_tol=1e-10,abs_tol=1e-8):raise ArithmeticError('对偶压力格与节点总力不守恒')
    summaries={};fields={};energies=0.
    for member,matrix in cache['matrices'].items():
        response=prepared.systems[member].solve(matrix@forces);energies+=response['strain_energy_j']
        summaries[member]={k:response[k] for k in ('strain_energy_j','maximum_free_residual_n','total_force_balance_n')}
        if include_fields:fields[member]=response
    ca=columns(active)[active,:];la=cache['lever'][active]
    stiffness=float(la@np.linalg.solve(ca,la))
    result.update(converged=True,status='converged',contacts=cells,nodal_contacts=nodal_rows,
        active_nodal_count=len(active),active_contact_count=len(cells),total_normal_force_n=float(forces.sum()),
        active_area_m2=float(sum(c['area_m2'] for c in cells)),maximum_cell_mean_pressure_pa=max(c['mean_pressure_pa'] for c in cells),
        pair_force_shares={str(k):float(sum(f for r,f in zip(records,forces) if r['tooth_offset']==k)) for k in offsets},
        pinion_torque_reconstructed_nm=pinion_reconstructed,pinion_torque_balance_relative_residual=residual,
        unloaded_tca_te_rad=cache['reference']['unloaded_te_rad'],mesh_correction_rad=cache['correction'],elastic_te_rad=theta,
        total_relative_te_rad=cache['reference']['unloaded_te_rad']+cache['correction']+theta,
        incremental_mesh_stiffness_nm_per_rad=stiffness,elastic_fields=fields,stress_summary=summaries,
        fe_energy_j=energies,fe_to_condensed_energy_relative_error=(energies/equilibrium['strain_energy_j']-1),
        limitations=['节点压力及中位对偶格为明确的集总接触离散，不是旧面片压力平滑',
            '压力力矩作用点为节点积分点，不能以对偶四边形中心替代',
            '仍需体网格/接触网格收敛；固定孔壁、无摩擦、小变形，不含轴承箱体'])
    return result
