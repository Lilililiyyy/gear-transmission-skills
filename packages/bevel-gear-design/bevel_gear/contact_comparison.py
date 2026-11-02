"""同一物理齿、相位和载荷下展示FE面片与局部半解析压力，保持模型语义独立。"""
from __future__ import annotations
import json
import math


def contact_method_comparison_data(bundle, surface_results, flank_mesh, pinion_request, *,
        fe_phase_index=None, surface_phase_index=None, tooth_offset=0):
    """校验相位/有符号载荷/制造请求，提取原压力多边形，不统一两种物理面积。"""
    import numpy as np
    from .contact_pattern import contact_pattern_data
    from .halfspace_plotting import halfspace_contact_pattern_data
    from .workflow_plotting import _loaded_results
    if bundle['case']['manufacturing']['pinion']['request']!=pinion_request:
        raise ValueError('两方法展示须使用同一显式小轮制造请求')
    fe_results=_loaded_results(bundle)
    if fe_phase_index is None:fe_phase_index=len(fe_results)//2
    if surface_phase_index is None:surface_phase_index=len(surface_results)//2
    for index,count in ((fe_phase_index,len(fe_results)),(surface_phase_index,len(surface_results))):
        if isinstance(index,bool) or not isinstance(index,int) or not 0<=index<count:
            raise ValueError('方法比较相位索引不存在')
    if len(fe_results)<3:raise ValueError('FE完整齿面数据入口要求至少三个已解姿态')
    # 旧FE提取器接收三个不同索引；这里只展示明确选择的一个，不复制姿态。
    indices=list(dict.fromkeys((fe_phase_index,0,len(fe_results)-1,len(fe_results)//2,1)))[:3]
    fe=contact_pattern_data(bundle,flank_mesh,phase_indices=indices,tooth_offset=tooth_offset)
    hs=halfspace_contact_pattern_data(surface_results,flank_mesh,pinion_request,
        phase_indices=[surface_phase_index],tooth_offset=tooth_offset)
    first=next(p for p in fe['phases'] if p['phase_index']==fe_phase_index);second=hs['phases'][0]
    if not math.isclose(first['driver_rad'],second['driver_rad'],abs_tol=1e-12,rel_tol=0.):
        raise ValueError('两方法必须展示同一实际小轮相位')
    if not math.isclose(first['pinion_torque_nm'],second['pinion_torque_nm'],abs_tol=1e-10,rel_tol=1e-12):
        raise ValueError('两方法必须使用同一有符号小轮载荷')
    if not first['facets'] or not second['polygons_m']:raise ValueError('所选物理齿须在两方法中均有承载结果')
    state=surface_results[surface_phase_index]
    materials=bundle['case'].get('materials',{})
    if materials and state.get('operators'):
        inverse=sum((1-materials[m]['poisson_ratio']**2)/materials[m]['young_modulus_pa'] for m in ('pinion','wheel'))
        if any(not math.isclose(operator['inverse_effective_modulus_per_pa'],inverse,rel_tol=1e-10,abs_tol=0.) for operator in state['operators']):
            raise ValueError('两方法的合成弹性材料不同，不能标作同工况')
    field=state['fields'][state['selected_tooth_offsets'].index(tooth_offset)]
    models=[dict(method='full_gear_fe',polygons_m=[json.loads(row['polygon_width_height_m_json']) for row in first['facets']],
        pressures_pa=[row['pressure_pa'] for row in first['facets']],forces_n=[row['force_n'] for row in first['facets']],
        areas_m2=[row['physical_area_m2'] for row in first['facets']],area_definition='实际三维FE面片面积',
        force_n=first['normal_force_n'],force_symbol='n',
        limitation='粗面片均压，压力网格未收敛'),
        dict(method='local_halfspace',polygons_m=second['polygons_m'],pressures_pa=second['pressure_pa'],
        forces_n=np.asarray(field['forces_n'])[field['contact_mask']].tolist(),
        areas_m2=[field['cell_area_m2']]*second['active_cell_count'],area_definition='参考接触平面矩形面积',
        force_n=second['plane_normal_force_n'],force_symbol='n,plane',
        limitation='不含齿根/轮体及跨齿结构柔度')]
    rows=[]
    for model in models:
        p=np.asarray(model['pressures_pa']);a=np.asarray(model['areas_m2']);f=np.asarray(model['forces_n'])
        if not np.isfinite(p).all() or np.any(p<0) or not np.allclose(p*a,f,rtol=1e-10,atol=1e-9):
            raise ValueError('模型压力/原积分面积/力不一致，不能作比较图')
        points=np.asarray(model['polygons_m']).reshape(-1,2)
        _,axes=np.linalg.eigh(np.cov(points,rowvar=False));sizes=np.ptp((points-points.mean(axis=0))@axes,axis=0)
        model.update(active_count=len(p),peak_pressure_mpa=float(p.max()*1e-6),
            integration_area_mm2=float(a.sum()*1e6),support_length_mm=float(sizes[1]*1000),
            support_width_mm=float(sizes[0]*1000),mean_pressure_mpa=float(f.sum()/a.sum()*1e-6))
        for polygon,pressure,area,force in zip(model['polygons_m'],p,a,f):
            rows.append(dict(row_kind='actual_loaded_element',method=model['method'],
                driver_rad=first['driver_rad'],pinion_torque_nm=first['pinion_torque_nm'],tooth_offset=tooth_offset,
                polygon_width_height_m_json=json.dumps(polygon),pressure_pa=float(pressure),
                original_integration_area_m2=float(area),force_n=float(force),area_definition=model['area_definition']))
        rows.append(dict(row_kind='method_summary',method=model['method'],active_count=model['active_count'],
            peak_pressure_pa=model['peak_pressure_mpa']*1e6,integration_area_m2=model['integration_area_mm2']*1e-6,
            support_length_m=model['support_length_mm']*.001,support_width_m=model['support_width_mm']*.001,
            force_n=model['force_n'],area_definition=model['area_definition']))
    for polygon in hs['background_polygons_m']:
        rows.append(dict(row_kind='actual_flank_triangle',polygon_width_height_m_json=json.dumps(polygon)))
    return dict(models=models,rows=rows,background_polygons_m=hs['background_polygons_m'],
        driving_angle_rad=first['driver_rad'],pinion_torque_nm=first['pinion_torque_nm'],tooth_offset=tooth_offset,
        interpretation='同工况下不同模型的结果展示，不是同柔度、同离散方法的精度对标')


def plot_contact_method_comparison(bundle,surface_results,output,*,flank_mesh,pinion_request,
        lang='zh',fe_phase_index=None,surface_phase_index=None,tooth_offset=0,pressure_limit_mpa=None):
    """两方法各显示完整工作面和等比例放大，统一线性色标、显式超量程。"""
    import numpy as np
    from .halfspace_plotting import pressure_colour_scale
    from .workflow_plotting import _paper,_panel,_save,_text
    data=contact_method_comparison_data(bundle,surface_results,flank_mesh,pinion_request,
        fe_phase_index=fe_phase_index,surface_phase_index=surface_phase_index,tooth_offset=tooth_offset)
    norm,cmap,scale=pressure_colour_scale([p for model in data['models'] for p in model['pressures_pa']],pressure_limit_mpa)
    background=np.asarray(data['background_polygons_m'])*1000;bounds=background.reshape(-1,2)
    with _paper(lang) as plt:
        from matplotlib.collections import PolyCollection
        from matplotlib.cm import ScalarMappable
        fig,axes=plt.subplots(2,2,layout='constrained',gridspec_kw=dict(width_ratios=[1.4,1.]))
        for i,model in enumerate(data['models']):
            polygons=np.asarray(model['polygons_m'])*1000;points=polygons.reshape(-1,2)
            method=_text(lang,'整轮FE：面片平均压力','Full-gear FE: facet-mean pressure') if i==0 else _text(lang,'局部半解析：单元压力','Local semi-analytical: cell pressure')
            for j,ax in enumerate(axes[i]):
                ax.add_collection(PolyCollection(background,facecolors='#E9EDF0',edgecolors='none',zorder=0))
                ax.add_collection(PolyCollection(polygons,array=np.asarray(model['pressures_pa'])*1e-6,
                    cmap=cmap,norm=norm,edgecolors='#34444C' if i==0 else 'none',linewidths=.4,antialiaseds=False,zorder=2))
                ax.axhline(0.,ls=(0,(3,3)),color='#93A1AA',lw=.45,zorder=1)
                if j:
                    span=np.ptp(points,axis=0);margin=np.maximum(span*.12,.12)
                    ax.set_xlim(points[:,0].min()-margin[0],points[:,0].max()+margin[0])
                    ax.set_ylim(points[:,1].min()-margin[1],points[:,1].max()+margin[1])
                    title=_text(lang,'局部放大（等比例）','Local view (equal aspect)')+'\n'+rf'$p_{{\max}}={model["peak_pressure_mpa"]:.0f}$ MPa'
                else:
                    ax.set_xlim(bounds[:,0].min()-.4,bounds[:,0].max()+.4)
                    ax.set_ylim(bounds[:,1].min()-.25,bounds[:,1].max()+.25)
                    title=method+'\n'+_text(lang,f'{model["active_count"]}个承载单元；',f'{model["active_count"]} loaded elements; ')+rf'$F_{{{model["force_symbol"]}}}={model["force_n"]:.0f}$ N'
                ax.set_title(title,loc='left',fontsize=7,pad=5)
                ax.set_ylabel(_text(lang,'节锥法向高度 h (mm)','Pitch-cone normal\nheight, h (mm)'));ax.set_aspect('equal',adjustable='box')
                ax.annotate('abcd'[i*2+j],(0,1),xycoords='axes fraction',xytext=(-6,18),textcoords='offset points',
                    ha='right',va='bottom',fontsize=9,weight='bold',annotation_clip=False)
        for ax in axes[-1]:ax.set_xlabel(_text(lang,r'齿宽位置 $w_\mathrm{f}$ (mm)',r'Face-width position, $w_\mathrm{f}$ (mm)'))
        fig.suptitle(_text(lang,'相同工况，不同柔度与接触离散模型','Same case, different compliance and contact models')+
            f'\n'+rf'$\theta_1={data["driving_angle_rad"]*180/math.pi:.2f}^\circ$, '+rf'$T_1={data["pinion_torque_nm"]:g}$ N·m, k={tooth_offset}',fontsize=9)
        bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),pad=.025,shrink=.88,
            extend='max' if scale['over_limit_cell_count'] else 'neither')
        bar.set_label(_text(lang,'接触压力 (MPa)','Contact pressure (MPa)'))
        if scale['over_limit_cell_count']:bar.ax.set_title(_text(lang,f'红色 >{scale["display_pressure_limit_mpa"]:g}',f'Red >{scale["display_pressure_limit_mpa"]:g}'),fontsize=6,pad=7)
        caption=_text(lang,
            '同一小轮物理齿、相位与有符号载荷下，分别展示整轮FE粗面片均压和局部半空间半解析单元压力。两行共用同一线性色标，红色和上端箭头明确表示超出显示量程，原始峰值仍标出。灰色为完整实际工作面，右列仅等比例放大，不插值/加宽压力区。FE使用原三维面片面积，半空间使用参考切平面矩形面积；两者力与压力定义分别保存。FE仅有少量承载面片，未解析连续压力斑点，不能据此称压力精度对标通过；半空间不含有限齿体/轮体柔度。此为瞬时压力，不是整周期累计着色印迹；两种方法的真实压力均未完成物理验证。',
            'The same physical pinion tooth, phase and signed load are shown for coarse full-gear FE facet-mean pressures and local half-space cell pressures. Both rows share one linear scale; red and the extended colorbar explicitly denote over-range cells, while raw peaks remain labelled. Grey shows the actual full flank. Right panels are equal-aspect local views without interpolation or widening. FE retains original 3D facet areas; the half-space retains reference-plane rectangle areas. The few loaded FE facets do not resolve a continuous patch, so this is not validated pressure agreement. The half-space excludes finite tooth/body compliance. These are instantaneous pressures, not cycle-accumulated marking patterns. Neither method has established actual-flank physical validation.')
        return _save(fig,output,'12_contact_method_comparison',lang,data['rows'],caption,height_mm=115,
            extra=dict(**scale,geometry_aspect='equal',interpolation='none',models=[{k:v for k,v in m.items()
                if k not in ('polygons_m','pressures_pa','forces_n','areas_m2')} for m in data['models']],
                driving_angle_rad=data['driving_angle_rad'],tooth_offset=tooth_offset,
                pinion_torque_nm=data['pinion_torque_nm'],comparison_scope=data['interpretation']))
