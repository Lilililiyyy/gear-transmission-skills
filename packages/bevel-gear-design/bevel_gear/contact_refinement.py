"""从真实有限曲面主轴窄带生成本体SI加密折线，不求压力或限定接触域。"""
from __future__ import annotations
from dataclasses import asdict,is_dataclass
import math


def _positive(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
        raise ValueError(name+'须为正有限数')
    return float(value)


def continuous_valid_tracks(points,valid):
    """行/列仅连接相邻有效采样点；遇到任何invalid立即断开，不跨空档。"""
    import numpy as np
    points=np.asarray(points,dtype=float);valid=np.asarray(valid)
    if valid.ndim!=2 or valid.dtype!=np.bool_ or points.shape!=valid.shape+(3,):raise ValueError('实际曲面点须对应二维布尔有效域')
    if not np.isfinite(points[valid]).all():raise ValueError('有效曲面点必须有限，不能外推或补零')
    lines=[];records=[];used=np.zeros_like(valid)
    def scan(indices,axis,fixed):
        run=[]
        def close():
            if len(run)>=2:
                xyz=points[tuple(np.asarray(run).T)]
                if np.any(np.linalg.norm(np.diff(xyz,axis=0),axis=1)<=1e-12):raise ArithmeticError('相邻有效网格点在实际曲面重合')
                lines.append(xyz.tolist())
                records.append(dict(axis=axis,fixed_index=fixed,start_index=list(run[0]),end_index=list(run[-1]),point_count=len(run)))
                for key in run:used[key]=True
        for key in indices:
            if valid[key]:run.append(key)
            else:close();run=[]
        close()
    ny,nx=valid.shape
    for j in range(ny):scan([(j,i) for i in range(nx)],'row',j)
    for i in range(nx):scan([(j,i) for j in range(ny)],'column',i)
    return lines,records,dict(valid_point_count=int(valid.sum()),used_valid_point_count=int(used.sum()),
        isolated_valid_point_count=int((valid&~used).sum()),invalid_sample_bridged=False,
        continuity_definition='相邻有效采样段；折线弦不是整段精确NURBS曲线或材料覆盖证书')


def _verify_surface_nodes(prepared,patch,member,prefix,world,valid):
    """独立重取每个有效UV的CAD原值与有限面分类，验证真实点后才输出控制。"""
    import numpy as np
    surface=prepared.trimmed[member].surface;bounds=surface.Bounds()
    uv=np.asarray(patch[prefix+'_uv'],dtype=float)[valid]
    R=np.asarray(patch['frames'][member+'_rotation_rows']);expected=[];classified=0
    for u,v in uv:
        if not math.isfinite(float(u)) or not math.isfinite(float(v)):raise ArithmeticError('有效面点缺少有限原UV')
        local=np.asarray(surface.Value(bounds[0]+u*(bounds[1]-bounds[0]),bounds[2]+v*(bounds[3]-bounds[2])).Coord())*.001
        if not prepared.trimmed[member].inside((float(u),float(v)),local):raise ArithmeticError('输出控制点不在实际有限BRep工作面')
        classified+=1;expected.append(R@local)
    error=float(np.linalg.norm(np.asarray(expected).reshape((-1,3))-world[valid],axis=1).max()) if len(uv) else 0.
    if error>1e-10:raise ArithmeticError('实际控制点与原NURBS UV重取值不一致')
    return dict(independent_CAD_point_checks=classified,finite_BRep_classification_checks=classified,
        maximum_original_surface_point_error_m=error)


def create_contact_refinement(prepared_surface,driving_angles_rad,*,tooth_offsets=(-1,0,1),
        half_widths_m=(.008,.0006),longitudinal_step_m=.00025,transverse_step_m=.0001,
        refined_size_m=75e-6,refinement_core_radius_m=.0003,refinement_radius_m=.0012,
        working_side_sign=-1,alignment_step_m=1e-4,alignment_fallback_to_surface_basis=True,
        ray_backend='batch',verify_surface_points=True,progress=None):
    """生成主轴窄带的纯SI区域文档与逐相位几何诊断。

    prepared_surface必须持有同制造请求的真实闭实体/有限工作面。工作侧
    符号只适配公共几何参考位姿选择，内部单位幅值不代表计算机械载荷。
    行/列的连续有效面点变为折线，每相位统一用整轮R0.T映回本体。
    输出可叠加原全轨迹h100背景；它不是压力求解域、加载覆盖或收敛证书。
    """
    import numpy as np
    from .contact_plane import build_principal_aligned_patches
    from .generated_tca import generated_tca_frames
    angles=list(driving_angles_rad);offsets=list(tooth_offsets)
    if not angles or any(isinstance(a,bool) or not isinstance(a,(int,float)) or not math.isfinite(a) for a in angles):raise ValueError('源相位须为非空有限rad列表')
    if any(b<=a for a,b in zip(angles,angles[1:])):raise ValueError('源相位须严格递增，不能遗漏失败相位后冒充完整请求')
    if not offsets or len(set(offsets))!=len(offsets) or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets):raise ValueError('控制齿对索引须是不重复整数')
    if len(half_widths_m)!=2:raise ValueError('窄带须给长短两个半宽m')
    a,b=(_positive(x,'加密带半宽m') for x in half_widths_m)
    dx=_positive(longitudinal_step_m,'沿带采样步长m');dy=_positive(transverse_step_m,'跨带采样步长m')
    fine=_positive(refined_size_m,'局部单元尺寸m');core=_positive(refinement_core_radius_m,'近场core半径m');reach=_positive(refinement_radius_m,'影响半径m')
    _positive(alignment_step_m,'主轴差分步长m')
    if not fine<=core<reach:raise ValueError('局部控制须满足fine≤core<radius')
    if isinstance(working_side_sign,bool) or working_side_sign not in (-1,1):raise ValueError('工作侧符号须为±1')
    if not isinstance(verify_surface_points,bool) or not isinstance(alignment_fallback_to_surface_basis,bool):raise ValueError('几何核验/fallback开关须为布尔值')
    if ray_backend not in ('scalar','batch'):raise ValueError('真实曲面射线后端须为scalar或batch')
    models=getattr(prepared_surface.tca,'cad_models',None)
    if not isinstance(models,dict) or set(models)!={'pinion','wheel'}:raise ValueError('需要持有真实双轮CAD及制造身份的有限面准备对象')
    manufacturing={}
    for member,model in models.items():
        if not is_dataclass(model.request) or not is_dataclass(model.settings):raise ValueError('实际制造request/blank须为标准化数据类')
        manufacturing[member]=dict(request=asdict(model.request),blank=asdict(model.settings))
    nx=2*math.ceil(a/dx)+1;ny=2*math.ceil(b/dy)+1;dx=2*a/(nx-1);dy=2*b/(ny-1)
    if nx*ny>500000:raise ValueError('控制采样单片超过显式五十万点上限')
    tracks={m:[] for m in models};phases=[]
    for phase_index,angle in enumerate(angles):
        geometry=build_principal_aligned_patches(prepared_surface,angle,float(working_side_sign),step_m=alignment_step_m,
            fallback_to_surface_basis=alignment_fallback_to_surface_basis,tooth_offsets=offsets,
            grid_shape=(nx,ny),spacing_m=(dx,dy),cover_complete_working_flank=False,check_cell_corners=False,
            follower_correction_rad=1e-10,ray_backend=ray_backend)
        frame0=generated_tca_frames(prepared_surface.tca,angle,geometry['reference_follower_angle_rad'],tooth_offset=0)
        translation=np.asarray(frame0.get('translation_m',[0.,0.,0.]));records=[];added={m:0 for m in models}
        for patch in geometry['patches']:
            valid=np.asarray(patch['valid_mask'],dtype=bool);row=dict(tooth_offset=patch['tooth_offset'],members={})
            for member,prefix in (('pinion','first'),('wheel','second')):
                R0=np.asarray(frame0[member+'_rotation_rows'])
                if not np.allclose(R0.T@R0,np.eye(3),rtol=0,atol=1e-12):raise ArithmeticError('整轮参考frame不是正交旋转')
                world=np.asarray(patch[prefix+'_points_world_m'],dtype=float)
                proof=_verify_surface_nodes(prepared_surface,patch,member,prefix,world,valid) if verify_surface_points else dict(independent_CAD_point_checks=0,verification_explicitly_disabled=True)
                local=(world-translation)@R0
                lines,segments,statistics=continuous_valid_tracks(local,valid)
                tracks[member].extend(lines);added[member]+=len(lines)
                error=float(np.linalg.norm(local[valid]@R0.T+translation-world[valid],axis=1).max()) if valid.any() else 0.
                wrong=(world[valid]-translation)@np.asarray(patch['frames'][member+'_rotation_rows'])
                wrong_delta=float(np.linalg.norm(local[valid]-wrong,axis=1).max()) if valid.any() else 0.
                row['members'][member]=dict(**statistics,track_count=len(lines),segments=segments,
                    body_frame_roundtrip_error_m=error,maximum_error_if_patch_frame_were_used_m=wrong_delta,**proof)
            records.append(row)
        if not all(added.values()):raise ArithmeticError('本相位没有连续有限控制带，不能跳过失败相位后输出完整计划')
        phase=dict(phase_index=phase_index,driving_angle_rad=angle,reference_follower_angle_rad=geometry['reference_follower_angle_rad'],
            whole_gear_frames=frame0,grid_alignment=geometry['grid_alignment'],patches=records,added_track_counts=added,
            source_geometry_valid=True,loaded_coverage_certified=False)
        phases.append(phase)
        if progress:progress(dict(stage='contact_refinement_phase_completed',phase_index=phase_index,driving_angle_rad=angle,
            added_track_counts=added,valid_point_count=sum(r['members']['pinion']['valid_point_count'] for r in records)))
    regions={m:[dict(refinement_polylines_m=tracks[m],refined_size_m=fine,
        refinement_core_radius_m=core,refinement_radius_m=reach)] for m in models}
    document=dict(schema='bevel_refinement_regions_v1',units='m',coordinate_system='gear_body_local',
        source_driver_angles_rad=angles,manufacturing=manufacturing,regions=regions,
        description='真实有限曲面主轴窄带的采样折线控制；叠加原全轨迹背景，不限定压力求解域，不构成实际载荷或收敛证明')
    diagnostics=dict(schema='bevel_contact_refinement_diagnostics_v1',requested_angles_rad=angles,requested_tooth_offsets=offsets,
        half_widths_m=[a,b],actual_grid_shape=[nx,ny],actual_sampling_steps_m=[dx,dy],working_side_sign=working_side_sign,
        working_side_identity='仅选择公共几何的工作侧和参考角；未计算机械载荷/压力',phases=phases,
        completed_phase_count=len(phases),all_requested_phases_generated=len(phases)==len(angles),
        independent_surface_verification_enabled=verify_surface_points,track_counts={m:len(lines) for m,lines in tracks.items()},
        pressure_solved=False,stiffness_assembled=False,mesh_generated=False,
        actual_loaded_region_coverage_certified=False,continuous_cycle_coverage_certified=False)
    return dict(document=document,diagnostics=diagnostics)
