"""从原已接受P0状态汇总整轮响应；不重算接触、不填补失败相位。"""
from __future__ import annotations

import math


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name+'须为有限数值')
    return float(value)


def load_distributed_fe_cycle(summary_file):
    """逐状态核对输入指纹及实际文件SHA，然后汇总，接受状态不补造。"""
    import hashlib
    import json
    from pathlib import Path
    path=Path(summary_file);raw=path.read_bytes();summary=json.loads(raw.decode('utf-8-sig'))
    if summary.get('schema')!='bevel_distributed_fe_workflow_v1':raise ValueError('需要分布FE工作流汇总文件')
    states=[];bindings=[]
    for row in summary['states']:
        state_path=(path.parent/row['state_file']).resolve()
        if not state_path.is_relative_to(path.parent.resolve()):raise ValueError('原状态路径不得越出计算结果目录')
        state_raw=state_path.read_bytes();digest=hashlib.sha256(state_raw).hexdigest()
        expected=row.get('source_state_sha256') or row.get('state_sha256')
        if expected is None or expected!=digest:raise ValueError('原逐相位状态缺少或不符合SHA绑定')
        state=json.loads(state_raw.decode('utf-8-sig'))
        if (state['phase_index']!=row['phase_index'] or state['driving_angle_rad']!=row['driving_angle_rad']
                or state.get('accepted') is not row.get('accepted')
                or state.get('input_fingerprint')!=summary.get('input_fingerprint')
                or state.get('source_identity_unchanged') is not True):
            raise ValueError('原状态的相位、输入、验收或计算源身份与汇总不一致')
        states.append(state);bindings.append(dict(phase_index=state['phase_index'],source_state_sha256=digest))
    result=summarize_distributed_fe_cycle(states)
    result.update(input_fingerprint=summary.get('input_fingerprint'),
        source_summary_sha256=hashlib.sha256(raw).hexdigest(),source_states=bindings)
    return result


def summarize_distributed_fe_cycle(states):
    """汇总从动轴TE、切线刚度、原面元面积与参考平面合力分担。

    输入为原 ``bevel_distributed_fe_state_v1`` 状态序列。失败状态保留
    断点；不以零力或零TE替代。各齿参考平面不同，因此F_plane分担
    仅是各原压力积分量的归一化，并非共同法向的矢量合力分担。
    一个齿距的端点覆盖也不等于相位极值收敛。
    """
    if not isinstance(states, (list, tuple)) or not states:
        raise ValueError('需要非空原始分布FE状态序列')
    rows=[]; seen=set(); ratio=None; geometry_sources=set()
    for state in states:
        if state.get('schema')!='bevel_distributed_fe_state_v1':
            raise ValueError('输入须采用原分布FE状态schema')
        index=state.get('phase_index')
        if isinstance(index,bool) or not isinstance(index,int) or index<0 or index in seen:
            raise ValueError('原相位索引须为互不重复的非负整数')
        seen.add(index); angle=_finite(state.get('driving_angle_rad'),'主动角rad')
        result=state.get('result',{})
        row=dict(phase_index=index,driving_angle_rad=angle,accepted=state.get('accepted') is True,
            status=result.get('status'))
        if not row['accepted']:
            rows.append(row);continue
        if (result.get('converged') is not True or result.get('geometry_linearization_checked') is not True
                or result.get('contact_representation')!='distributed_P0_surface_cells'
                or result.get('halfspace_compliance_used') is not False
                or result.get('equilibrium',{}).get('certificate',{}).get('passed') is not True):
            raise ValueError('标记已接受的状态缺少原FE、几何或KKT通过证据')
        geometry=result['geometry']; local_ratio=_finite(geometry['pinion_to_wheel_ratio'],'齿数比')
        if local_ratio<=0 or (ratio is not None and local_ratio!=ratio):
            raise ValueError('各相位必须使用相同的正齿数比')
        ratio=local_ratio
        if _finite(geometry['driving_angle_rad'],'原几何主动角')!=angle:
            raise ValueError('原状态角与原几何角不一致')
        nominal=-angle*ratio
        total=_finite(result['loaded_follower_angle_rad'],'实际受载从动角')-nominal
        native=result.get('geometry_source')=='FE_boundary'
        geometry_sources.add('FE_boundary' if native else 'CAD_surface_gauss')
        if native:
            reference=_finite(result['reference_follower_angle_rad'],'当前原生参考从动角')
            step=_finite(result['linearization_step_theta_rad'],'末次原生QP增量')
            reported=_finite(result['total_loaded_te_rad'],'原生总受载TE')
            if (not math.isclose(total,reported,abs_tol=1e-12,rel_tol=1e-9) or
                    not math.isclose(reference+step,result['loaded_follower_angle_rad'],abs_tol=1e-12,rel_tol=1e-9)):
                raise ArithmeticError('原生总TE或当前参考角+末次增量与实际加载角不一致')
            # 没有求得卸载FE网格相切时，不能将末次增量命名为纯弹性TE。
            row.update(geometry_source='FE_boundary',linearization_step_theta_rad=step,
                reference_follower_angle_rad=reference,native_unloaded_mesh_tangency_solved=False)
        else:
            tca=geometry['reference_tca']
            unloaded=_finite(tca['follower_increment_angle_rad'],'首接触从动角')-nominal
            reference=_finite(result['unloaded_reference_follower_angle_rad'],'原无载参考从动角')
            correction=reference-nominal-unloaded
            elastic=_finite(result['elastic_te_rad'],'从动轴弹性TE')
            if not math.isclose(total,unloaded+correction+elastic,abs_tol=1e-12,rel_tol=1e-9):
                raise ArithmeticError('直接受载角与无载/参考修正/弹性TE之和不一致')
            row.update(unloaded_te_rad=unloaded,reference_phase_correction_rad=correction,elastic_te_rad=elastic)
        stiffness=_finite(result['incremental_stiffness_nm_per_rad'],'从动轴切线刚度')
        if stiffness<=0:raise ValueError('已接受有载状态的切线刚度须为正')
        offsets=result['selected_tooth_offsets'];fields=result['fields']
        if len(offsets)!=len(fields) or len(set(offsets))!=len(offsets):
            raise ValueError('原齿对索引与压力场不能一一对应')
        pair=[]
        for offset,field in zip(offsets,fields):
            if field.get('tooth_offset')!=offset:raise ValueError('压力场齿对索引错位')
            area=_finite(field['cell_area_m2'],'原参考面元面积')
            if area<=0:raise ValueError('原参考面元面积须为正')
            pressure=field['pressure_pa'];force=field['forces_n'];mask=field['contact_mask']
            if not (len(pressure)==len(force)==len(mask)):raise ValueError('原二维压力/力/接触掩码形状不同')
            values=[];pressures=[];active=0
            for pp,ff,mm in zip(pressure,force,mask):
                if not (len(pp)==len(ff)==len(mm)):raise ValueError('原二维压力场行布局不一致')
                for p,f,m in zip(pp,ff,mm):
                    p=_finite(p,'原P0压力Pa');f=_finite(f,'原平面合力N')
                    if p<0 or f<0 or not math.isclose(f,p*area,rel_tol=1e-10,abs_tol=1e-10):
                        raise ArithmeticError('原P0压力与平面合力不满足非负pA关系')
                    if not isinstance(m,bool):raise ValueError('原接触掩码须为布尔值')
                    values.append(f);pressures.append(p);active+=int(m)
            pair.append(dict(tooth_offset=offset,force_plane_n=math.fsum(values),
                peak_pressure_pa=max(pressures,default=0.),active_P0_area_plane_m2=active*area))
        total_force=math.fsum(p['force_plane_n'] for p in pair)
        if total_force<=0 or not math.isclose(total_force,_finite(result['total_force_plane_n'],'原总平面合力'),rel_tol=1e-10,abs_tol=1e-9):
            raise ArithmeticError('原逐齿平面合力之和与总量不一致')
        total_area=math.fsum(p['active_P0_area_plane_m2'] for p in pair)
        if not math.isclose(total_area,_finite(result['contact_area_plane_m2'],'原P0总接触面积'),rel_tol=1e-10,abs_tol=1e-15):
            raise ArithmeticError('原逐齿P0承载面积之和与总量不一致')
        peak=max(p['peak_pressure_pa'] for p in pair)
        if not math.isclose(peak,_finite(result['maximum_pressure_pa'],'原峰压'),rel_tol=1e-10,abs_tol=1e-3):
            raise ArithmeticError('原二维压力最大值与峰压字段不一致')
        for p in pair:p['force_plane_fraction']=p['force_plane_n']/total_force
        row.update(total_te_rad=total,wheel_axis_stiffness_nm_per_rad=stiffness,maximum_pressure_pa=peak,
            active_P0_area_plane_m2=total_area,total_force_plane_n=total_force,pair_force_plane_shares=pair,
            pinion_torque_relative_residual=_finite(result['pinion_torque_relative_residual'],'主动轴力矩残差'))
        if native:
            row['minimum_piece_vertex_gap_m']=_finite(result['geometry_check']['loaded_FE_phi_geometry_check']['minimum_remaining_gap_m'],
                '有限phi下真实交片顶点最小间隙')
        else:
            row['minimum_gauss_remaining_gap_m']=_finite(result['minimum_gauss_remaining_gap_m'],'原Gauss最小间隙')
        rows.append(row)
    rows.sort(key=lambda r:r['driving_angle_rad'])
    if len({r['driving_angle_rad'] for r in rows})!=len(rows):raise ValueError('不同相位索引不能使用重复角度')
    if len(geometry_sources)>1:raise ValueError('同一周期不能混用CAD积分与实际FE边界几何')
    return dict(schema='bevel_distributed_fe_cycle_response_v1',rows=rows,
        selected_state_count=len(rows),accepted_state_count=sum(r['accepted'] for r in rows),
        selected_states_completed=all(r['accepted'] for r in rows),phase_resolution_converged=False,
        geometry_source=next(iter(geometry_sources),None),
        units=dict(angle='rad',pressure='Pa',area='m²',force_plane='N',stiffness='N·m/rad'),
        force_share_definition='各齿原参考平面∫p dA_plane的归一化；不冒称共同法向矢量分担',
        te_definition='总TE=实际从动角+z1/z2主动角；原生FE只报告总量与末次增量，未求卸载网格相切时不分解纯弹性TE')


def plot_distributed_fe_cycle(data,output,*,lang='zh'):
    """同一原P0周期的TE、从动轴刚度及平面力分担，保留失败断点。"""
    from .workflow_plotting import _paper,_save,_text,_panel
    from .publication_style import BLUE,RED,COLORS
    if data.get('schema')!='bevel_distributed_fe_cycle_response_v1':raise ValueError('需要分布FE周期响应汇总')
    rows=data['rows'];x=[r['driving_angle_rad']*180/math.pi for r in rows];arcsec=180/math.pi*3600
    native=data.get('geometry_source')=='FE_boundary'
    with _paper(lang) as plt:
        fig,grid=plt.subplots(2,2,sharex=True,layout='constrained');axes=grid.ravel()
        if native:
            axes[0].plot(x,[r.get('total_te_rad',math.nan)*arcsec for r in rows],'o-',ms=2.3,color=BLUE)
            axes[0].set_ylabel(_text(lang,'受载从动轴总传动误差\nΔθ₂ (arcsec)','Total loaded follower TE\nΔθ₂ (arcsec)'))
            axes[1].plot(x,[r.get('maximum_pressure_pa',math.nan)*1e-6 for r in rows],'o-',ms=2.3,color=RED)
            axes[1].set_ylabel(_text(lang,'原P0峰值接触压力\npₘₐₓ (MPa)','Original P0 peak pressure\npₘₐₓ (MPa)'))
        else:
            axes[0].plot(x,[r.get('unloaded_te_rad',math.nan)*arcsec for r in rows],'o-',ms=2.3,color=BLUE)
            axes[0].set_ylabel(_text(lang,'无载传动误差\nΔθ₂,₀ (arcsec)','Unloaded transmission error\nΔθ₂,₀ (arcsec)'))
            for key,color,zh,en in (('elastic_te_rad',RED,'弹性TE','Elastic TE'),('total_te_rad',BLUE,'总TE','Total TE')):
                axes[1].plot(x,[r.get(key,math.nan)*arcsec for r in rows],'o-',ms=2.3,color=color,label=_text(lang,zh,en))
            axes[1].set_ylabel(_text(lang,'受载从动轴传动误差\nΔθ₂ (arcsec)','Loaded follower transmission error\nΔθ₂ (arcsec)'));axes[1].legend(ncol=2,fontsize=6)
        axes[2].plot(x,[r.get('wheel_axis_stiffness_nm_per_rad',math.nan)*1e-3 for r in rows],'o-',ms=2.3,color=BLUE)
        axes[2].set_ylabel(_text(lang,'从动轴切线刚度\nKθ₂ (kN·m/rad)','Follower tangent stiffness\nKθ₂ (kN·m/rad)'))
        offsets=sorted({p['tooth_offset'] for r in rows for p in r.get('pair_force_plane_shares',[]) if p['force_plane_n']>0})
        for i,k in enumerate(offsets):
            fractions=[]
            for r in rows:
                if not r['accepted']:fractions.append(math.nan);continue
                fractions.append(next((p['force_plane_fraction']*100 for p in r['pair_force_plane_shares'] if p['tooth_offset']==k),0.))
            axes[3].plot(x,fractions,'o-',ms=2.3,color=COLORS[i%len(COLORS)],label=f'k={k}')
        axes[3].set_ylabel(_text(lang,'参考平面力分担\n','Reference-plane force share\n')+r'$F_{\mathrm{plane},k}/\sum_k F_{\mathrm{plane},k}$ (%)')
        axes[3].legend(ncol=max(1,len(offsets)),fontsize=6)
        for ax in axes[2:]:ax.set_xlabel(_text(lang,'主动轮增量角 θ₁ (deg)','Pinion increment angle θ₁ (deg)'))
        for ax,letter in zip(axes,'abcd'):_panel(ax,letter)
        source=[]
        for r in rows:
            base={k:v for k,v in r.items() if k!='pair_force_plane_shares'}
            for p in r.get('pair_force_plane_shares',[]):base.update({f'pair_{p["tooth_offset"]}_{key}':value for key,value in p.items() if key!='tooth_offset'})
            source.append(base)
        note=(_text(lang,'同一原生FE周期；总TE由实际加载角直接计算，未求卸载网格相切，不分解纯弹性项。曲线连接原计算点，失败相位保留断点；分担为各局部参考平面压力积分的归一化。',
                'One native FE cycle. Total TE is calculated directly from the loaded angle; no unloaded FE tangency or elastic-only split is claimed. Lines join computed samples, with gaps at rejected phases. Shares normalize pressure integrals on local reference planes.') if native else
            _text(lang,'同一整轮P0有限元接受状态；TE为从动轴量。曲线连接实际计算点，失败相位留断点。分担是各局部参考平面的压力积分，非共同法向矢量分担；完整SI源数据保留无载及参考相位修正。',
                'Accepted states of one full-gear P0 FE analysis. TE refers to the follower shaft. Lines join computed samples; rejected phases remain gaps. Shares normalize pressure integrals on each local reference plane, rather than a common normal-force vector. SI data retain unloaded TE and reference-phase corrections.'))
        return _save(fig,output,'18_loaded_cycle_response',lang,source,note,
            width_mm=190,height_mm=130,extra=dict(source_schema=data['schema'],phase_resolution_converged=False))
