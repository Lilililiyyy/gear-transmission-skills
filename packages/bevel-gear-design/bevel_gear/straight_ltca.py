"""理想直齿对象的三维实体有限元承载接触，新写小变形离散模型。

Gmsh四面体→三维线弹性K→功共轭C→非负力/间隙/扭矩平衡。
工作面与根区按straight_cad的建模边界解释，工程根延伸不是切削
圆角。面片平均压力不作为网格收敛Hertz峰值或ISO齿根应力使用。
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from pathlib import Path

from .geometry import GeometryResult
from .elasticity import (ElasticMaterial, ElasticSystem, ContactLoadChannel,
    mesh_step_tetra, condense_contact_flexibility, contact_load_matrix)
from .straight_tca import solve_straight_tca, straight_pair_frames
from .contact_equilibrium import solve_contact_equilibrium


@dataclass
class PreparedStraightLtca:
    """两轮实际CAD、体网格和受约束弹性系统，供多个载荷/转角复用。"""
    geometry: GeometryResult
    models: dict
    systems: dict
    diagnostics: dict
    _contact_cache: dict = field(default_factory=dict,repr=False)


def prepare_straight_ltca(geometry, *, pinion_cad_settings, wheel_cad_settings,
                          pinion_material, wheel_material, directory,
                          element_size_m=.004, full_gears=True, maximum_nodes=50000,
                          contact_refined_size_m=None, refinement_radius_m=.002,
                          refinement_driving_angle_rad=0.):
    """生成可重用的两轮实体FE；孔壁全部平移固定，扇区侧面保持自由。

    默认整轮实体，凝聚矩阵保留同轮多齿的非对角耦合。显式设置False
    使用单齿扇区仅作边界敏感性试算；孔壁很远时自由侧面的扇区可能
    过软，须由小变形/小转角检查拒绝，不能据此代替整轮结果。
    """
    from .straight_cad import build_straight_tooth_sector, build_straight_gear, StraightCadSettings
    from .cad_import import ImportedCadModel, export_model_step
    from .elasticity import _numeric
    np, _, _ = _numeric()
    if not isinstance(geometry, GeometryResult):
        raise TypeError("有限元LTCA需要GeometryResult")
    if not isinstance(full_gears,bool):
        raise ValueError("full_gears须为布尔值")
    if geometry.inputs.spiral_angle_rad != 0:
        raise NotImplementedError("此对象的实体FE接触限理想直齿，弧齿须使用实际展成齿面独立接入")
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    models,systems={},{}
    refine_tca=solve_straight_tca(geometry,refinement_driving_angle_rad) if contact_refined_size_m is not None else None
    if refine_tca is not None and not refine_tca['converged']:
        raise ValueError("局部加密参考姿态没有有效解析接触线")
    for member,settings,material in (("pinion",pinion_cad_settings,pinion_material),("wheel",wheel_cad_settings,wheel_material)):
        if isinstance(settings,dict):settings=StraightCadSettings(**settings)
        if isinstance(material,dict):material=ElasticMaterial(**material)
        builder=build_straight_gear if full_gears else build_straight_tooth_sector
        model=builder(geometry,settings,member=member)
        path=directory/(member+"_elastic_body.step")
        export_model_step(ImportedCadModel(model.shape,model.diagnostics),path)
        refine_points=[]
        if refine_tca is not None:
            # 明确围绕给定姿态的解析接触线加密；不声称覆盖任意转角。
            delta=refine_tca[member+'_polar_angle_rad']
            lo=model.boundaries['end_inner']['z_m']
            hi=model.boundaries['end_outer']['z_m']
            for z in np.linspace(lo,hi,max(9,int((hi-lo)/refinement_radius_m)*2+1)):
                refine_points.append(model.surface.profile.evaluate(float(z)/math.cos(delta),delta,'right'))
        mesh=mesh_step_tetra(path,element_size_m=element_size_m,maximum_nodes=maximum_nodes,
            refinement_points_m=refine_points,refined_size_m=contact_refined_size_m,
            refinement_radius_m=refinement_radius_m if refine_points else None,
            refinement_polyline=bool(refine_points))
        nodes=np.asarray(mesh.nodes_m)
        radius=np.linalg.norm(nodes[:,:2],axis=1)
        bore=model.boundaries['bore']['radius_m']
        tolerance=max(settings.tolerance_m*2,bore*1e-8)
        fixed=np.flatnonzero(np.abs(radius-bore)<=tolerance)
        system=ElasticSystem(mesh,material,fixed,maximum_dofs=3*maximum_nodes)
        system.diagnostics.update(bore_radius_m=bore,bore_selection_tolerance_m=tolerance,
            cad_volume_m3=model.diagnostics['volume_m3'],
            relative_mesh_volume_error=(system.diagnostics['volume_m3']/model.diagnostics['volume_m3']-1))
        models[member],systems[member]=model,system
    diagnostics=dict(full_gears=full_gears,element_size_m=element_size_m,
        contact_refined_size_m=contact_refined_size_m,refinement_driving_angle_rad=refinement_driving_angle_rad,
        members={m:s.diagnostics for m,s in systems.items()},
        support_model="两轮孔壁完全固定；未包含轴、轴承与箱体柔度",
        root_model="显式工程径向根连接；不是实际切削根圆角",
        neighbour_coupling="整轮三维实体耦合" if full_gears else "各单齿扇区独立；省略邻齿之间的轮体耦合",
        limitations=["低阶四面体与离散接触面须独立加密；没有预设工业精度", "本模型是理想直齿有限元近似，非工艺齿面认证"])
    return PreparedStraightLtca(geometry,models,systems,diagnostics)


def solve_straight_fe_ltca(prepared, driving_angle_rad, pinion_torque_nm, *,
                          tooth_offsets=(-1,0,1), maximum_distance_m=2e-4,
                          maximum_normal_angle_rad=.25, maximum_channels=240,
                          mesh_phase_tolerance_m=1e-9, maximum_mesh_phase_iterations=8,
                          maximum_loaded_angle_rad=.005,
                          maximum_pinion_torque_relative_residual=.02, include_fields=True):
    """在给定主动角与扭矩下求离散承载接触；工作侧为正转right。

    先独立求解析相切TCA，再对有限元三角面的离散初相位作显式修正。
    该数值相位与弹性滞后分开输出，禁止把网格齿形误差称加工TE。
    柔度包含两体实体变形，不附加未经给定的经验啮合刚度。
    从动轴满足KKT力矩约束，主动轴重构残差另按显式门限核验。未接受
    的加载候选只保存在candidate，不发布压力、弹性TE或应力场为成功结果。
    """
    from .elasticity import _numeric
    from .contact_mesh import select_straight_flank_facets, pair_facet_contacts
    np, _, _ = _numeric()
    if not isinstance(prepared,PreparedStraightLtca):
        raise TypeError("需要已准备的PreparedStraightLtca")
    for name,value in (("主动角rad",driving_angle_rad),("小轮扭矩Nm",pinion_torque_nm)):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError(name+"须为有限实数")
    if pinion_torque_nm<=0:raise ValueError("当前承载分支要求正小轮扭矩")
    for name,value in (("候选距离",maximum_distance_m),("网格相位残差",mesh_phase_tolerance_m),
                       ("小转角上限",maximum_loaded_angle_rad),
                       ("主动轴重构扭矩相对残差上限",maximum_pinion_torque_relative_residual)):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
            raise ValueError(name+"须为正有限数")
    if not isinstance(include_fields,bool):
        raise ValueError("include_fields须为布尔值")
    if isinstance(maximum_channels,bool) or not isinstance(maximum_channels,int) or not 1<=maximum_channels<=512:
        raise ValueError("maximum_channels须为1至512的整数")
    if isinstance(maximum_mesh_phase_iterations,bool) or not isinstance(maximum_mesh_phase_iterations,int) or maximum_mesh_phase_iterations<1:
        raise ValueError("网格相位迭代上限须为正整数")
    offsets=list(tooth_offsets)
    if not offsets or len(set(offsets))!=len(offsets) or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets):
        raise ValueError("配对齿偏移须为非空且不重复的整数序列")
    cache_key=(driving_angle_rad,tuple(offsets),maximum_distance_m,maximum_normal_angle_rad,maximum_channels,mesh_phase_tolerance_m,maximum_mesh_phase_iterations)
    cache_hit=cache_key in prepared._contact_cache
    if not cache_hit:
        g=prepared.geometry
        ratio=g.pinion.teeth/g.wheel.teeth
        tca={}
        facets={}
        for k in offsets:
            result=solve_straight_tca(g,driving_angle_rad,tooth_offset=k,flank="right")
            tca[k]=result
            if not result['converged']:continue
            for member,count in (("pinion",g.pinion.teeth),("wheel",g.wheel.teeth)):
                model=prepared.models[member]
                tooth=(k if member=="pinion" else -k)%count if model.tooth_index is None else model.tooth_index
                facets[member,k]=select_straight_flank_facets(prepared.systems[member].mesh,model.surface,
                    flank="right",tooth_index=tooth,matching_tolerance_m=max(3e-6,model.settings.tolerance_m*5))
        usable=[k for k in offsets if tca[k]['converged']]
        if not usable:
            raise ArithmeticError("所选齿对在实际有效齿廓内均无相切解")
        nominal=-driving_angle_rad*ratio
        def geometry_at(correction):
            records=[]
            for k in usable:
                frames=straight_pair_frames(g,driving_angle_rad,nominal+correction,tooth_offset=k)
                base=straight_pair_frames(g,driving_angle_rad,nominal+correction,tooth_offset=0)
                a=frames['pinion_rotation_rows'] if prepared.models['pinion'].tooth_index is not None else base['pinion_rotation_rows']
                b=frames['wheel_rotation_rows'] if prepared.models['wheel'].tooth_index is not None else base['wheel_rotation_rows']
                pairs=pair_facet_contacts(facets['pinion',k],facets['wheel',k],a,b,
                    maximum_distance_m=maximum_distance_m,maximum_normal_angle_rad=maximum_normal_angle_rad)
                axis1=np.array(frames['pinion_axis_world'])
                axis2=np.array(frames['wheel_axis_world'])
                for pair in pairs:
                    normal=np.array(pair['normal_world'])
                    lever1=float(normal@np.cross(axis1,pair['first_point_world_m']))
                    lever2=float(-normal@np.cross(axis2,pair['second_point_world_m']))
                    if lever1<=1e-10 or lever2<=1e-10:continue
                    record=dict(pair,tooth_offset=k,pinion_lever_m=lever1,wheel_lever_m=lever2)
                    records.append(record)
            if not records:raise ArithmeticError("实际实体网格在候选窗口内没有法向相对的工作侧面")
            if len(records)>maximum_channels:
                raise ValueError(f"候选接触面片{len(records)}个超过显式上限{maximum_channels}；须调整网格或候选距离，不静默抽点")
            return records
        correction=0.
        phase_history=[]
        for _ in range(maximum_mesh_phase_iterations):
            records=geometry_at(correction)
            gaps=np.array([r['signed_gap_m'] for r in records])
            lever=np.array([r['wheel_lever_m'] for r in records])
            minimum=float(gaps.min())
            phase_history.append(dict(follower_correction_rad=correction,minimum_signed_gap_m=minimum,candidate_count=len(records)))
            if abs(minimum)<=mesh_phase_tolerance_m:break
            correction+=float(np.min(gaps/lever))
            if abs(correction)>.02:
                raise ArithmeticError("网格初相位修正超出0.02rad，不能把大范围齿形失配当离散误差消除")
        else:
            raise ArithmeticError("有限元面片的未加载相位未达到显式间隙残差门槛")
        if float(gaps.min()) < -mesh_phase_tolerance_m:
            raise ArithmeticError("网格未加载工作面仍有超限负间隙")
        n=len(records)
        compliance=np.zeros((n,n))
        member_channels={}
        condensation_diagnostics={}
        for member in ("pinion","wheel"):
            first=member=="pinion"
            prefix='first' if first else 'second'
            channels=[ContactLoadChannel(tuple(r[prefix+'_node_indices']),tuple(r[prefix+'_weights']),
                tuple(r[prefix+'_force_normal_local'])) for r in records]
            member_channels[member]=channels
            system=prepared.systems[member]
            groups=[list(range(n))] if prepared.models[member].tooth_index is None else [
                [i for i,r in enumerate(records) if r['tooth_offset']==k] for k in usable]
            condensation_diagnostics[member]=[]
            for indices in groups:
                if not indices:continue
                # 单体投影可半正定；真正的接触稳定性在两体C相加后由KKT层检查。
                result=condense_contact_flexibility(system,[channels[i] for i in indices],require_positive_definite=False)
                compliance[np.ix_(indices,indices)]+=np.array(result['compliance_m_per_n'])
                condensation_diagnostics[member].append({k:v for k,v in result.items() if k!='compliance_m_per_n'})
        prepared._contact_cache[cache_key]=(g,ratio,tca,facets,usable,nominal,correction,phase_history,records,gaps,lever,n,compliance,member_channels,condensation_diagnostics)
    g,ratio,tca,facets,usable,nominal,correction,phase_history,records,gaps,lever,n,compliance,member_channels,condensation_diagnostics=prepared._contact_cache[cache_key]
    target_wheel_torque=pinion_torque_nm/ratio
    equilibrium=solve_contact_equilibrium(gaps.tolist(),compliance.tolist(),lever_arms_m=lever.tolist(),
        torque_nm=target_wheel_torque,allow_initial_penetration=bool(float(gaps.min())<0),gap_tolerance_m=max(1e-10,mesh_phase_tolerance_m*2))
    f=np.array(equilibrium.forces_n)
    contacts=[]
    # 导出积分面片的实际FE顶点；整轮采用同一轴相位，扇区才使用各齿局部帧。
    base_frames=straight_pair_frames(g,driving_angle_rad,nominal+correction,tooth_offset=0)
    for i,r in enumerate(records):
        tooth_frames=straight_pair_frames(g,driving_angle_rad,nominal+correction,tooth_offset=r['tooth_offset'])
        vertices={}
        for member,prefix in (('pinion','first'),('wheel','second')):
            frames=base_frames if prepared.models[member].tooth_index is None else tooth_frames
            rotation=np.array(frames[member+'_rotation_rows'])
            local=prepared.systems[member].nodes[list(r[prefix+'_node_indices'])]
            vertices[prefix+'_vertices_world_m']=(local@rotation.T).tolist()
        contacts.append(dict(r,force_n=float(f[i]),mean_pressure_pa=float(f[i]/r['area_m2']),
            remaining_gap_m=equilibrium.remaining_gaps_m[i],active=i in equilibrium.active_indices,**vertices))
    elastic_fields={}
    stress_summary={}
    for member in ("pinion","wheel"):
        system=prepared.systems[member]
        channels=member_channels[member]
        groups=[('whole',list(range(n)))] if prepared.models[member].tooth_index is None else [
            (str(k),[i for i,r in enumerate(records) if r['tooth_offset']==k]) for k in usable]
        elastic_fields[member]={}
        stress_summary[member]={}
        for label,indices in groups:
            if not indices:continue
            loads=contact_load_matrix(system,[channels[i] for i in indices])@f[indices]
            response=system.solve(loads)
            stress_summary[member][label]=dict(
                maximum_element_von_mises_pa=max(response['element_von_mises_pa']),
                maximum_displacement_m=float(np.linalg.norm(response['displacements_m'],axis=1).max()),
                maximum_free_residual_n=response['maximum_free_residual_n'],
                total_force_balance_n=response['total_force_balance_n'],strain_energy_j=response['strain_energy_j'])
            if include_fields:elastic_fields[member][label]=response
    active=equilibrium.active_indices
    if active:
        ca=compliance[np.ix_(active,active)]
        la=lever[active]
        incremental_stiffness=float(la@np.linalg.solve(ca,la))
    else:incremental_stiffness=None
    theta=equilibrium.theta_rad
    model_valid=theta is not None and abs(theta)<=maximum_loaded_angle_rad
    pinion_reconstructed=float(sum(r['pinion_lever_m']*force for r,force in zip(records,f)))
    pinion_relative_residual=(pinion_reconstructed-pinion_torque_nm)/abs(pinion_torque_nm)
    pinion_balance_valid=abs(pinion_relative_residual)<=maximum_pinion_torque_relative_residual
    accepted=equilibrium.converged and model_valid and pinion_balance_valid
    candidate=None if accepted else dict(theta_rad=theta,incremental_mesh_stiffness_nm_per_rad=incremental_stiffness,
        contacts=contacts,stress_summary=stress_summary,elastic_fields=elastic_fields,equilibrium=asdict(equilibrium),
        pinion_torque_reconstructed_nm=pinion_reconstructed,wheel_torque_reconstructed_nm=equilibrium.resultant_load,
        interpretation='未通过的候选，不作为已接受加载物理结果')
    status=(equilibrium.status if not equilibrium.converged else 'small_rotation_limit' if not model_valid
            else 'pinion_torque_residual_limit' if not pinion_balance_valid else 'converged')
    loaded_theta=theta if accepted else None
    result=dict(converged=accepted,contact_representation='facets',
        status=status,
        model='理想直齿三维线弹性FE＋三角面片无摩擦小转角接触',
        geometry_flexibility_cache_hit=cache_hit,discretization_converged=None,pressure_peak_validated=False,
        driving_angle_rad=driving_angle_rad,pinion_torque_nm=pinion_torque_nm,wheel_target_torque_nm=target_wheel_torque,
        pinion_torque_reconstructed_nm=pinion_reconstructed if accepted else None,
        pinion_torque_balance_residual_nm=pinion_reconstructed-pinion_torque_nm if accepted else None,
        pinion_torque_balance_relative_residual=pinion_relative_residual if accepted else None,
        wheel_torque_reconstructed_nm=equilibrium.resultant_load if accepted else None,
        torque_constraint='从动真实轴力臂严格KKT；目标按名义齿数比换算，主动重构扭矩独立核验',
        analytic_unloaded_te_rad=tca[usable[0]]['unloaded_transmission_error_rad'],mesh_unloaded_correction_rad=correction,
        elastic_follower_lag_rad=loaded_theta,
        total_follower_te_rad=None if not accepted else correction+theta,
        incremental_mesh_stiffness_nm_per_rad=incremental_stiffness if accepted else None,
        total_normal_force_n=float(f.sum()) if accepted else None,
        active_contact_count=len(active) if accepted else None,candidate_count=n,
        active_area_m2=float(sum(records[i]['area_m2'] for i in active)) if accepted else None,
        maximum_cell_mean_pressure_pa=max((c['mean_pressure_pa'] for c in contacts),default=0.) if accepted else None,
        pair_force_shares={str(k):float(sum(force for r,force in zip(records,f) if r['tooth_offset']==k)) for k in usable} if accepted else None,
        contacts=contacts if accepted else [],equilibrium=asdict(equilibrium),tca={str(k):v for k,v in tca.items()},
        mesh_phase_history=phase_history,condensation=condensation_diagnostics,
        stress_summary=stress_summary if accepted else {},elastic_fields=elastic_fields if accepted else {},candidate=candidate,
        quality_checks={'small_rotation_passed':model_valid,'pinion_torque_balance_passed':pinion_balance_valid,
            'maximum_pinion_torque_relative_residual':maximum_pinion_torque_relative_residual,
            'candidate_pinion_torque_relative_residual':pinion_relative_residual},
        preparation=prepared.diagnostics,
        limitations=['面片平均压力不是已收敛Hertz峰值；常应力单元峰值不是ISO齿根应力',
            '无摩擦、小变形、小转角线性化；须检查网格、候选窗口和接触边界收敛',
            '工程根区与固定孔壁决定本算例柔度，未作工艺齿面或物理试验标定',
            '网格初相位与弹性TE分开报告，不能将离散齿形误差当加工TE',
            prepared.diagnostics['neighbour_coupling']])
    return result
