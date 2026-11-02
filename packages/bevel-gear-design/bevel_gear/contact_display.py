"""四类接触压力的工程显示：半解析在前，每幅独立色标。

原始离散解始终单列；FE可选合力守恒压力恢复及随接触轨迹的相位插值。
插值包络是显示估计，不是新增有限元平衡解，不供强度验收。
"""
from pathlib import Path
import gzip,json,math,hashlib
import numpy as np

from .contact_recovery import raster_grid,raw_raster,recover_sequence


def prepare_pressure_comparison(rows, *, view='engineering', recovery_source_file=None,
        shape=(161,641), subdivisions=16, smoothing_length_m=None):
    """返回按半解析瞬时/周期、FE瞬时/周期排列的四个SI显示场。"""
    if view not in ('raw','engineering'):raise ValueError('显示方式只允许raw或engineering')
    by_method={row['method']:row for row in rows}
    if set(by_method)!={'full_gear_fe','local_halfspace'}:raise ValueError('须同时提供FE和局部半解析')
    panels=[];outline=by_method['local_halfspace']['instantaneous']['outline_geometry']
    grid=raster_grid(outline,shape)
    for method in ('local_halfspace','full_gear_fe'):
        for key,kind in (('instantaneous','instantaneous'),('cycle_maximum','sampled_cycle_maximum')):
            data=by_method[method][key]
            if data['method']!=method or data['kind']!=kind or not data['outline_geometry'].equals(outline):
                raise ValueError('四个输入须为同一材料齿面及明确的方法/图种类')
            panels.append(dict(data=data,method=method,kind=kind))
    report=dict(view=view,actual_original_pressure_modified=False,new_FE_solves=0,
        display_pressure_is_raw=(view=='raw'),physical_validation_added=False,
        axis_units='m',pressure_units='Pa',display_grid_shape=list(shape),colour_scale='independent')
    # 半解析仍使用实际P1源的瞬时场和129个真实相位最大值，不借给FE。
    for panel in panels[:2]:panel['pressure_pa']=raw_raster(panel['data'],grid)
    if view=='engineering':
        if recovery_source_file is None:raise ValueError('工程显示需要原P0网格源，不能由图片反推或拟合斑点')
        file=Path(recovery_source_file);raw=file.read_bytes();seq=json.loads(gzip.decompress(raw))
        if seq.get('schema')!='bevel_pressure_recovery_source_v1':raise ValueError('压力恢复源协议不符')
        source_cycle=panels[3]['data'];identity={(f['driver_angle_rad'],f['tooth_offset']):f['source_state_sha256'] for f in source_cycle['source_fields']}
        source_angles=sorted({f['driver_angle_rad'] for f in source_cycle['source_fields']})
        if seq['driver_angles_rad']!=source_angles:raise ValueError('P0恢复源与原始FE周期的实际相位不同')
        for source in seq['fields']:
            if np.any(np.asarray(source['pressure_pa'])>0):
                if identity.get((source['driver_angle_rad'],source['tooth_offset']))!=source['source_state_sha256']:
                    raise ValueError('恢复源不是同一真实FE状态/物理齿，不混用新旧网格')
        instant,cycle,diagnostics=recover_sequence(seq,grid,subdivisions=subdivisions,smoothing_length_m=smoothing_length_m)
        panels[2]['pressure_pa']=instant;panels[3]['pressure_pa']=cycle
        report.update(recovery=diagnostics,recovery_source_sha256=hashlib.sha256(raw).hexdigest(),
            FE_display_pressure_is_raw=False,FE_cycle_is_phase_interpolated=True)
    else:
        for panel in panels[2:]:panel['pressure_pa']=raw_raster(panel['data'],grid)
    summaries=[]
    for panel in panels:
        panel['pressure_pa']=np.where(grid['flank_mask'],panel['pressure_pa'],0.)
        data=panel['data'];peak=float(panel['pressure_pa'].max())
        if not math.isfinite(peak) or peak<=0:raise ValueError('显示须有真实非零压力源，不生成占位图')
        # 原始视图色标覆盖全部源节点峰值，不能因显示格漏过窄峰而截短色标。
        scale_peak=max(peak,float(data['source_peak_pressure_pa'])) if view=='raw' else peak
        panel['limit_pa']=math.ceil(scale_peak/100e6)*100e6
        summary=dict(method=panel['method'],kind=panel['kind'],raw_peak_pressure_pa=float(data['source_peak_pressure_pa']),
            display_peak_pressure_pa=peak,colour_limit_pa=panel['limit_pa'],
            solved_phase_count=len({f['driver_angle_rad'] for f in data['source_fields']}),
            postprocessed=(view=='engineering' and panel['method']=='full_gear_fe'))
        summaries.append(summary)
    report['panels']=summaries
    return dict(panels=panels,grid=grid,report=report)


def plot_pressure_comparison(rows, output, *, lang='zh', stem='17_contact_composite',
        display_scope='certified',comparison_evidence=None,source_compression='gzip',
        instantaneous_zoom=False,view='engineering',recovery_source_file=None,
        prepared=None,subdivisions=16,smoothing_length_m=None):
    """绘制四幅独立MPa色标云图；FE工程估计和原始求解峰值分别标注。"""
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.patches import Patch
    from .publication_style import contact_pressure_colormap
    from .workflow_plotting import _paper,_save,_text
    from .contact_composite import _draw_geometry
    bundle=prepare_pressure_comparison(rows,view=view,recovery_source_file=recovery_source_file,
        subdivisions=subdivisions,smoothing_length_m=smoothing_length_m) if prepared is None else prepared
    panels=bundle['panels'];grid=bundle['grid'];report=bundle['report'];cmap=contact_pressure_colormap()
    output=Path(output);output.mkdir(parents=True,exist_ok=True);(output/'data').mkdir(exist_ok=True)
    with _paper(lang) as plt:
        fig,axes=plt.subplots(2,2,figsize=(200/25.4,(146 if instantaneous_zoom else 112)/25.4))
        fig.subplots_adjust(left=.075,right=.94,top=.84,bottom=.29 if instantaneous_zoom else .23,
            hspace=1.65 if instantaneous_zoom else .95,wspace=.50)
        for index,(ax,panel) in enumerate(zip(axes.ravel(),panels)):
            data=panel['data'];pressure=panel['pressure_pa'];limit=panel['limit_pa']
            norm=Normalize(0.,limit*1e-6);levels=np.linspace(0.,limit*1e-6,97)
            _draw_geometry(ax,data['outline_geometry'],facecolor=cmap(0.),edgecolor='#555D65',linewidth=.5,zorder=0)
            outline_patch=ax.patches[-1]
            image=ax.contourf(grid['xx_m']*1e3,grid['yy_m']*1e3,pressure*1e-6,
                levels=levels,norm=norm,cmap=cmap,antialiased=False,zorder=1)
            # Matplotlib 3.6的等值集合尚无统一裁切接口，逐集合设置同一真实轮廓。
            if hasattr(image,'set_clip_path'):image.set_clip_path(outline_patch)
            else:
                for collection in image.collections:collection.set_clip_path(outline_patch)
            _draw_geometry(ax,data['unknown_geometry'],facecolor='#D5DADD' if display_scope=='certified' else 'none',
                edgecolor='#89959C',hatch='////',linewidth=0.,zorder=2)
            _draw_geometry(ax,data['outline_geometry'],facecolor='none',edgecolor='#555D65',linewidth=.5,zorder=3)
            x0,y0,x1,y1=data['outline_geometry'].bounds;ax.set(xlim=(x0*1e3,x1*1e3),ylim=(y0*1e3,y1*1e3));ax.set_aspect('equal')
            method=_text(lang,'局部半解析','Semi-analytical') if panel['method']=='local_halfspace' else _text(lang,'整轮有限元','Full-gear FE')
            if panel['kind']=='instantaneous':kind=_text(lang,'瞬时压力','Instantaneous pressure')
            else:kind=_text(lang,'周期最大压力','Cycle-maximum pressure')
            if view=='engineering' and panel['method']=='full_gear_fe':
                kind=_text(lang,'瞬时压力恢复','Recovered instantaneous pressure') if index==2 else _text(lang,'相位插值最大压力','Phase-interpolated maximum')
            ax.set_title(method+'\n'+kind,loc='left',fontsize=8,pad=7)
            summary=report['panels'][index]
            label_peak=summary['raw_peak_pressure_pa'] if view=='raw' else summary['display_peak_pressure_pa']
            label=rf'$p_\mathrm{{max}}={label_peak*1e-6:.0f}$ MPa'
            ax.text(.99,1.015,label,transform=ax.transAxes,ha='right',va='bottom',fontsize=6)
            ax.annotate('abcd'[index],(0,1),xycoords='axes fraction',xytext=(-9,10),textcoords='offset points',ha='right',weight='bold',fontsize=9,annotation_clip=False)
            ax.set_xlabel(_text(lang,r'齿宽位置 $w_\mathrm{f}$ (mm)',r'Face-width position $w_\mathrm{f}$ (mm)'),fontsize=7)
            ax.set_ylabel(_text(lang,'法向高度 h (mm)','Normal height h (mm)'),fontsize=7);ax.tick_params(labelsize=6)
            cax=make_axes_locatable(ax).append_axes('right',size='3%',pad=.07)
            bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,ticks=[0.,limit*5e-7,limit*1e-6])
            bar.ax.tick_params(labelsize=5.5,length=2);bar.set_label('p (MPa)',fontsize=6,labelpad=3)
            if instantaneous_zoom and index in (0,2):
                positive=pressure>summary['display_peak_pressure_pa']*1e-4
                ys,xs=np.where(positive);padx=.0005;pady=.00015
                box=(grid['x_m'][xs.min()]-padx,grid['x_m'][xs.max()]+padx,grid['y_m'][ys.min()]-pady,grid['y_m'][ys.max()]+pady)
                parent=ax.get_position();detail=fig.add_axes([parent.x0+.035,.51 if index==0 else .105,parent.width*.80,.08])
                detail.contourf(grid['xx_m']*1e3,grid['yy_m']*1e3,pressure*1e-6,levels=levels,norm=norm,cmap=cmap)
                detail.set(xlim=(box[0]*1e3,box[1]*1e3),ylim=(box[2]*1e3,box[3]*1e3));detail.set_aspect('equal');detail.tick_params(labelsize=5)
                detail.set_title(_text(lang,'压力区细节（相同色标）','Pressure detail (same scale)'),fontsize=6)
        title=_text(lang,'锥齿轮接触压力对照','Bevel-gear contact pressure comparison')
        actual_angle=float(panels[0]['data']['source_fields'][0]['driver_angle_rad'])
        conditions=rf'$\alpha={actual_angle:g}$ rad'
        if comparison_evidence is not None and 'signed_pinion_torque_nm' in comparison_evidence:
            conditions+=rf', $T_1={comparison_evidence["signed_pinion_torque_nm"]:g}$ N·m'
        fig.suptitle(title+'\n'+conditions,fontsize=10,y=.98)
        n_fe=report['panels'][3]['solved_phase_count'];n_sa=report['panels'][1]['solved_phase_count']
        raw_fe=report['panels'][3]['raw_peak_pressure_pa']*1e-6
        count=report.get('recovery',{}).get('temporal_recovery',{}).get('display_phase_count',n_fe)
        if view=='engineering':
            note=_text(lang,
                f'半解析：{n_sa}个求解相位。FE：{n_fe}个求解相位 → {count}个显示相位；下排为合力守恒恢复/相位插值，非新增FE解。\n每幅使用独立线性MPa色标；FE原始周期峰值{raw_fe:.1f} MPa另行保留。',
                f'Semi-analytical: {n_sa} solves. FE: {n_fe} solves → {count} display phases; lower row uses force-preserving recovery / phase interpolation.\nIndependent linear MPa scales. The original FE cycle peak, {raw_fe:.1f} MPa, is retained separately.')
        else:
            note=_text(lang,f'原始离散解：半解析{n_sa}相位，FE {n_fe}相位；不恢复压力、不进行相位插值。每幅独立线性色标。',
                f'Raw discrete results: semi-analytical {n_sa} phases, FE {n_fe} phases. No recovery or phase interpolation. Independent scales.')
        fig.text(.075,.08 if instantaneous_zoom else .105,note,fontsize=6,va='bottom',linespacing=1.5,color='#40474C')
        if any(p['data']['unknown_geometry'].area>0 for p in panels):
            fig.legend(handles=[Patch(facecolor='#D5DADD',edgecolor='#89959C',hatch='////',label=_text(lang,'原采样未解析区域','Unresolved in original samples'))],loc='lower center',bbox_to_anchor=(.5,.015),fontsize=6)
        caption=_text(lang,'半解析瞬时、周期最大；FE瞬时、周期最大。各图独立色标。工程显示不覆盖原始解，不用于峰压精度或强度验收。',
            'Semi-analytical instantaneous/cycle maximum, then FE instantaneous/cycle maximum. Independent scales; engineering displays do not replace raw solutions or validate peak pressure.')
        record=_save(fig,output,stem,lang,report['panels'],caption,width_mm=200,height_mm=146 if instantaneous_zoom else 112,
            extra=dict(schema='bevel_contact_display_v1',display=report,colour_scale_shared_by_all_panels=False,original_pressure_values_modified=False))
    np.savez_compressed(output/'data'/f'{stem}_SI.npz',x_m=grid['x_m'],y_m=grid['y_m'],flank_mask=grid['flank_mask'],pressure_pa=np.stack([p['pressure_pa'] for p in panels]))
    (output/'data'/f'{stem}_display.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return record
