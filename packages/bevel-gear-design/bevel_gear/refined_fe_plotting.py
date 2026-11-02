"""真实FE节点压力的连续标量云图、完整裁切齿面与等比例局部图。

每个压力节点拥有若干原FE三角形的中位对偶四边形；一个节点不能按
对偶格数重复计数。默认按原FE三角形作C0-P1标量插值；原求解压力和
积分定义不改，不作RGB混色、滤波、假想椭圆拟合或坐标扩宽。
"""
from __future__ import annotations

import json
import math


def _finite(value,name,*,positive=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or positive and value<=0:
        raise ValueError(name+'须为有限'+('正' if positive else '')+'实数')
    return float(value)


def _background(flank_mesh,pinion_request):
    from .contact_pattern import select_pinion_left_flank
    selected=select_pinion_left_flank(flank_mesh,pinion_request)
    coordinates=selected['coordinates'];polygons=[];rows=[]
    for index in selected['triangle_indices']:
        triangle=flank_mesh['triangles'][index]
        polygon=[[coordinates['width_position_m'][i],coordinates['normal_height_m'][i]] for i in triangle]
        polygons.append(polygon)
        rows.append(dict(row_kind='actual_trimmed_flank_triangle',cad_triangle_index=index,
            cad_face_index=flank_mesh['face_indices'][index],polygon_width_height_m_json=json.dumps(polygon)))
    if not polygons:raise ValueError('缺少实际完整裁切齿面，不能以矩形补造背景')
    return selected,polygons,rows


def _connectivity(polygons,forces,areas,node_ids,pressures):
    """按真实三维共边统计承载域，不以凸包或插值把分离区域接通。"""
    import numpy as np
    parent=list(range(len(polygons)));owners={};jumps=[]
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,polygon in enumerate(polygons):
        keys=[tuple(int(round(float(v)/1e-12)) for v in point) for point in polygon]
        for a,b in zip(keys,keys[1:]+keys[:1]):
            edge=tuple(sorted((a,b)))
            if edge in owners:
                j=owners[edge];parent[find(i)]=find(j)
                if node_ids[i]!=node_ids[j]:jumps.append(abs(pressures[i]-pressures[j])/max((pressures[i]+pressures[j])*.5,1e-30))
            else:owners[edge]=i
    groups={}
    for i in range(len(polygons)):groups.setdefault(find(i),[]).append(i)
    components=sorted((dict(dual_cell_count=len(indices),active_pressure_node_count=len({node_ids[i] for i in indices}),
        normal_force_n=sum(forces[i] for i in indices),original_integration_area_m2=sum(areas[i] for i in indices))
        for indices in groups.values()),key=lambda v:v['normal_force_n'],reverse=True)
    return dict(edge_connected_component_count=len(components),components=components,
        shared_coordinate_tolerance_m=1e-12,
        adjacent_node_relative_pressure_jump_median=None if not jumps else float(np.median(jumps)),
        adjacent_node_relative_pressure_jump_p95=None if not jumps else float(np.quantile(jumps,.95)),
        pressure_jump_definition='共边不同承载节点的|p_i-p_j|/((p_i+p_j)/2)，含真实梯度，不直接等于误差')


def refined_fe_contact_data(results,flank_mesh,pinion_request,*,labels=None,tooth_offset=0):
    """保留真实对偶格/原力节点/积分面积，核对每节点与同齿载荷后提取SI数据。"""
    from .contact_pattern import pinion_flank_coordinates
    from .generated_surface import _request
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('tooth_offset须为整数物理齿偏移')
    if not isinstance(results,(list,tuple)) or not 1<=len(results)<=3:raise ValueError('需要一至三个实际FE结果')
    labels=list(labels) if labels is not None else [str(i+1) for i in range(len(results))]
    if len(labels)!=len(results) or len(set(labels))!=len(labels) or any(not isinstance(label,str) or not label.strip() for label in labels):
        raise ValueError('每个FE结果须有不同的非空文字标签')
    request=_request(pinion_request)
    selected,background,rows=_background(flank_mesh,pinion_request)
    models=[];state=None
    for result,label in zip(results,labels):
        if not result.get('converged') or result.get('contact_representation')!='nodal_dual_cells':
            raise ValueError('需要已接受的nodal_dual_cells真实FE结果，不能把旧面片或积分点冒充细网格')
        angle=_finite(result['driving_angle_rad'],'driving_angle_rad');torque=_finite(result['pinion_torque_nm'],'pinion_torque_nm')
        if state is None:state=(angle,torque)
        elif not math.isclose(angle,state[0],rel_tol=0.,abs_tol=1e-12) or not math.isclose(torque,state[1],rel_tol=1e-12,abs_tol=1e-10):
            raise ValueError('同图FE细化结果须为同一实际相位与有符号扭矩')
        frames=result['mesh_reference_frames'];cells=[];nodes={}
        for node in result.get('nodal_contacts',[]):
            if node['tooth_offset']==tooth_offset and node.get('active') and node['force_n']>0:
                index=node['first_pressure_node_index']
                if isinstance(index,bool) or not isinstance(index,int) or index<0 or index in nodes:raise ValueError('原承载压力节点索引须为唯一非负整数')
                pressure=_finite(node['pressure_pa'],'nodal_pressure_pa',positive=True)
                area=_finite(node['area_m2'],'nodal_area_m2',positive=True);force=_finite(node['force_n'],'nodal_force_n',positive=True)
                if not math.isclose(pressure*area,force,rel_tol=2e-10,abs_tol=1e-9):raise ValueError('原压力节点的pA=f不一致')
                application=node['first_point_world_m']
                if len(application)!=3:raise ValueError('原压力节点须保留三维力作用点')
                for value in application:_finite(value,'first_point_world_m')
                nodes[index]=node
        if not nodes:raise ValueError('所选物理齿没有原承载节点，不能只凭着色四边形声称FE压力')
        grouped={i:dict(force=0.,area=0.,count=0) for i in nodes}
        for index,contact in enumerate(result['contacts']):
            if not contact.get('active') or contact['tooth_offset']!=tooth_offset or contact['force_n']<=0:continue
            node_id=contact['first_pressure_node_index']
            if node_id not in nodes:raise ValueError('承载对偶格找不到其原压力节点')
            world=contact.get('first_vertices_world_m')
            if not isinstance(world,(list,tuple)) or len(world)<3:raise ValueError('每个FE对偶格须保留实际三维边界，不能以节点中心补造')
            pressure=_finite(contact['mean_pressure_pa'],'mean_pressure_pa',positive=True)
            area=_finite(contact['area_m2'],'area_m2',positive=True);force=_finite(contact['force_n'],'force_n',positive=True)
            if not math.isclose(pressure*area,force,rel_tol=2e-10,abs_tol=1e-9):raise ValueError('对偶格原压力乘原积分面积与力不一致')
            node=nodes[node_id]
            if not math.isclose(pressure,node['pressure_pa'],rel_tol=2e-10,abs_tol=1e-5):raise ValueError('对偶格压力与其原节点压力不一致')
            application=contact.get('application_point_world_m')
            if application is None or len(application)!=3 or any(abs(_finite(a,'application_point_world_m')-b)>1e-12 for a,b in zip(application,node['first_point_world_m'])):
                raise ValueError('须保留原节点力作用点，不能以对偶格质心替代')
            converted=pinion_flank_coordinates(world,request.machine.pitch_cone_angle_rad,request.machine.reference_cone_distance_m,
                rotation_rows=frames['pinion_rotation_rows'],translation_m=frames.get('translation_m',(0.,0.,0.)))
            polygon=list(map(list,zip(converted['width_position_m'],converted['normal_height_m'])))
            cells.append(dict(pressure_pa=pressure,force_n=force,area_m2=area,node_id=node_id,
                polygon_width_height_m=polygon,polygon_world_m=world,application_point_world_m=application))
            grouped[node_id]['force']+=force;grouped[node_id]['area']+=area;grouped[node_id]['count']+=1
            rows.append(dict(row_kind='actual_loaded_nodal_dual_cell',result_label=label,contact_index=index,
                tooth_offset=tooth_offset,first_pressure_node_index=node_id,pressure_pa=pressure,force_n=force,
                original_integration_area_m2=area,polygon_width_height_m_json=json.dumps(polygon),
                polygon_world_m_json=json.dumps(world),application_point_world_m_json=json.dumps(application)))
        for index,node in nodes.items():
            actual=grouped[index]
            if not math.isclose(actual['force'],node['force_n'],rel_tol=2e-9,abs_tol=1e-8) or not math.isclose(actual['area'],node['area_m2'],rel_tol=2e-9,abs_tol=1e-18):
                raise ValueError('某压力节点的对偶格原力/面积总和与节点集总量不一致')
        total=sum(c['force_n'] for c in cells)
        expected=result['pair_force_shares'].get(str(tooth_offset),result['pair_force_shares'].get(tooth_offset,0.))
        if not math.isclose(total,expected,rel_tol=2e-9,abs_tol=1e-8):raise ValueError('所选物理齿原接触力与齿对分担不一致')
        connectivity=_connectivity([c['polygon_world_m'] for c in cells],[c['force_n'] for c in cells],
            [c['area_m2'] for c in cells],[c['node_id'] for c in cells],[c['pressure_pa'] for c in cells])
        model=dict(label=label,cells=cells,active_pressure_node_count=len(nodes),loaded_dual_cell_count=len(cells),
            peak_pressure_mpa=max(c['pressure_pa'] for c in cells)*1e-6,normal_force_n=total,
            integration_area_m2=sum(c['area_m2'] for c in cells),connectivity=connectivity,
            pressure_peak_validated=result.get('pressure_peak_validated',False),
            discretization_converged=result.get('discretization_converged'),
            pressure_definition='vertex-lumped nodal pressure：同节点各中位对偶格采用原p_node，非平滑面片压力')
        rows.append(dict(row_kind='nodal_contact_summary',result_label=label,
            **{k:v for k,v in model.items() if k not in ('cells','connectivity','label')},connectivity_json=json.dumps(connectivity)))
        models.append(model)
    return dict(models=models,rows=rows,background_polygons_m=background,tooth_offset=tooth_offset,
        driving_angle_rad=state[0],pinion_torque_nm=state[1],
        projection_definition=selected['coordinates']['projection_definition'],
        physical_area_definition='原三维FE三角形的节点中位对偶积分面积；子午投影不重新定义面积',
        face_identification=selected['face_identification'])


def plot_refined_fe_contact_pattern(results,output,*,pinion_request,flank_mesh,labels=None,
        tooth_offset=0,lang='zh',pressure_limit_mpa=None,zoom_window_m=None,stem='14_refined_fe_contact',display='continuous',
        domain_certificates=None):
    """原压力连续显示；domain_certificates逐结果绑定完整齿面的非接触证据。"""
    import numpy as np
    from .halfspace_plotting import pressure_colour_scale
    from .workflow_plotting import _paper,_save,_text
    from .pressure_field import reconstruct_fe_pressure_field,project_pressure_field,draw_pressure_field,pressure_field_rows
    from .pressure_background import pressure_domain_background,pressure_background_rows,draw_pressure_background,add_pressure_domain_legend,pressure_domain_metadata
    if display not in ('continuous','cells'):raise ValueError('display须为continuous或cells')
    data=refined_fe_contact_data(results,flank_mesh,pinion_request,labels=labels,tooth_offset=tooth_offset)
    certificates=[None]*len(results) if domain_certificates is None else list(domain_certificates)
    if len(certificates)!=len(results):raise ValueError('覆盖证书须与FE结果逐一对应')
    domains=[pressure_domain_background(result,data['background_polygons_m'],method='full_gear_fe',
        tooth_offset=tooth_offset,certificate=certificate,pinion_request=pinion_request) for result,certificate in zip(results,certificates)]
    for model,domain in zip(data['models'],domains):data['rows'].extend(pressure_background_rows(domain,label=model['label']))
    fields=[]
    if display=='continuous':
        for result,model in zip(results,data['models']):
            field=project_pressure_field(reconstruct_fe_pressure_field(result,tooth_offset=tooth_offset),pinion_request,result['mesh_reference_frames'])
            fields.append(field);data['rows'].extend(pressure_field_rows(field,label=model['label']))
    if zoom_window_m is not None and (len(zoom_window_m)!=2 or any(_finite(v,'zoom_window_m',positive=True)<=0 for v in zoom_window_m)):
        raise ValueError('局部裁窗须为两个正的物理m尺寸')
    norm,cmap,scale=pressure_colour_scale([c['pressure_pa'] for m in data['models'] for c in m['cells']],pressure_limit_mpa)
    background=np.asarray(data['background_polygons_m'])*1000;bounds=background.reshape(-1,2)
    zoom_records=[]
    with _paper(lang) as plt:
        from matplotlib.collections import PolyCollection
        from matplotlib.cm import ScalarMappable
        from matplotlib.patches import Rectangle
        fig,axes=plt.subplots(len(data['models']),2,squeeze=False,gridspec_kw=dict(width_ratios=[.8,1.5]))
        # 等比例的扁长齿面不参与自动压缩布局；预留轴名、图例及色标位置。
        fig.subplots_adjust(left=.085,right=.855,bottom=.19,top=.81,wspace=.45,hspace=1.15)
        for i,model in enumerate(data['models']):
            polygons=[np.asarray(c['polygon_width_height_m'])*1000 for c in model['cells']]
            values=np.asarray([c['pressure_pa']*1e-6 for c in model['cells']]);points=np.concatenate(polygons)
            if display=='continuous':points=np.asarray(fields[i]['points_width_height_m'])*1000
            if zoom_window_m is None:
                span=np.ptp(points,axis=0);margin=np.maximum(span*.1,.08)
                low=points.min(axis=0)-margin;high=points.max(axis=0)+margin
            else:
                peak=model['cells'][int(np.argmax(values))]
                position=np.asarray(peak['polygon_width_height_m']).mean(axis=0)*1000
                half=np.asarray(zoom_window_m)*500;low=position-half;high=position+half
            zoom_records.append(dict(result_label=model['label'],width_height_bounds_m=[(low*.001).tolist(),(high*.001).tolist()]))
            known=[c['polygon_width_height_m'] for c in model['cells']]
            if display=='continuous':known=np.asarray(fields[i]['points_width_height_m'])[fields[i]['triangles']]
            for j,ax in enumerate(axes[i]):
                draw_pressure_background(ax,data['background_polygons_m'],norm=norm,cmap=cmap,info=domains[i],known_polygons_m=known)
                if display=='continuous':draw_pressure_field(ax,fields[i],norm=norm,cmap=cmap)
                else:ax.add_collection(PolyCollection(polygons,array=values,cmap=cmap,norm=norm,
                    edgecolors='none',antialiaseds=False,zorder=2))
                ax.axhline(0.,ls=(0,(3,3)),color='#8C989F',lw=.5,zorder=1)
                if j:
                    ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1])
                    title=_text(lang,'局部等比例放大','Local view, equal aspect')
                    if zoom_window_m is not None:title=_text(lang,'峰压附近物理裁窗','Physical crop near peak')+f'  {zoom_window_m[0]*1000:g} × {zoom_window_m[1]*1000:g} mm'
                    title+='\n'+rf'$p_{{\max}}={model["peak_pressure_mpa"]:.0f}$ MPa'
                else:
                    ax.set_xlim(bounds[:,0].min()-.4,bounds[:,0].max()+.4);ax.set_ylim(bounds[:,1].min()-.25,bounds[:,1].max()+.25)
                    if zoom_window_m is not None:ax.add_patch(Rectangle(low,*(high-low),fill=False,edgecolor='#353F45',lw=.7,zorder=3))
                    title=_text(lang,'全齿面定位','Location on full flank')
                    if len(data['models'])>1:title=model['label']+'\n'+title
                    add_pressure_domain_legend(ax,domains[i],norm=norm,cmap=cmap,lang=lang)
                ax.set_title(title,loc='left',fontsize=7,pad=5)
                ax.set_aspect('equal',adjustable='box');ax.set_ylabel('h (mm)' if display=='continuous' and not j else
                    _text(lang,'节锥法向高度 h (mm)','Pitch-cone normal\nheight, h (mm)'))
                ax.annotate(chr(97+i*2+j),(0,1),xycoords='axes fraction',xytext=(-6,18),textcoords='offset points',
                    ha='right',va='bottom',fontsize=9,weight='bold',annotation_clip=False)
        for ax in axes[-1]:ax.set_xlabel(_text(lang,r'齿宽位置 $w_\mathrm{f}$ (mm)','Face-width position\n'+r'$w_\mathrm{f}$ (mm)'))
        fig.suptitle(_text(lang,'整轮FE接触压力','Full-gear FE contact pressure')+'\n'+
            rf'$\theta_1={data["driving_angle_rad"]*180/math.pi:.2f}^\circ$, '+rf'$T_1={data["pinion_torque_nm"]:g}$ N·m, k={tooth_offset}',fontsize=9)
        bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),pad=.025,shrink=.88,
            extend='max' if scale['over_limit_cell_count'] else 'neither')
        bar.set_label(_text(lang,'接触压力 p (MPa)','Contact pressure, p (MPa)'))
        bar.ax.set_title(_text(lang,'线性上限','Linear upper limit')+f'\n{scale["display_pressure_limit_mpa"]:.0f} MPa',fontsize=6,pad=8)
        caption=_text(lang,
            '同一物理小轮齿k、相位和载荷的真实整轮FE节点—面集总接触。每个承载节点按其原三角形中位对偶格直接着色；节点数与显示格数分别标注，vertex-lumped nodal pressure不作节点间插值。灰色为同请求完整真实裁切工作面，右列只作物理裁窗/等比例放大；黑框给出左图中裁窗的位置，不人为扩宽接触斑点。蓝色低压、红色高压，各图共用线性色标；每行标出该齿原峰值。w从外端向内端，h为相对节锥的法向高度，虚线h=0。子午投影不是等面积展开，力、面积与力作用点保留原FE定义；仍须独立验证压力网格和物理精度。',
            'Actual full-gear FE node-to-surface lumped contact at the same physical pinion tooth k, phase and load. Each loaded node is rendered over its original triangular median-dual cells; nodes and display cells are counted separately. Vertex-lumped nodal pressure is not interpolated between nodes. Grey is the complete actual trimmed flank from the same request. Right panels only crop physical coordinates at equal aspect; black boxes locate explicit crops in the full-flank panels. Patches are not widened. Blue denotes low and red high pressure with one shared linear scale; each row labels the raw peak on that tooth. w runs from outer to inner end and h is normal height relative to the pitch cone, dashed at h=0. This meridional projection is not area-preserving; force, integration area and force application points retain their FE definitions. Pressure mesh and physical accuracy require independent validation.')
        if display=='continuous':caption=_text(lang,
            '同一物理齿和相位的整轮FE压力。沿原真实三角形以原节点压力作连续P1标量云图，原求解压力、节点力矩积分点及承载面积均不改；灰色为完整真实裁切齿面。右图显示完整斑点或明确给定的等比例物理裁窗。蓝低红高，共用线性色标。显示场标量面积积分与原节点法向标量合力守恒；连续显示支撑面积与原对偶承载面积分别保存，不能混作新的接触求解或压力收敛验证。',
            'Full-gear FE pressure on the same physical tooth and phase. A continuous P1 scalar field interpolates original nodal pressures on actual FE triangles. Original solver values, nodal moment points and loaded areas remain unchanged. Grey is the actual complete trimmed flank; the local panel retains equal physical aspect. Blue denotes low and red high pressure on one linear scale. The display-field scalar integral conserves the original nodal normal resultant; display support and original dual loaded areas remain separate. This reconstruction is not a new contact solve or pressure-convergence validation.')
        caption+=_text(lang,' 零压力采用色标p=0的深蓝；已证明的求解覆盖／距离排除分片显示蓝底，未决分片单独保留灰色网纹。部分证书不等于完整面通过，未决压力不补零；蓝底严格服从实际CAD裁切边界。',
            ' Zero pressure uses the exact p=0 blue of the colour scale. Proven solved-domain or distance-excluded pieces have a blue base; unresolved pieces remain individually grey and hatched. A partial certificate does not certify the full flank or replace unknown pressures with zero. The base follows the actual CAD trim.')
        caption=caption.replace('灰色为同请求完整真实裁切工作面，','完整背景为同请求真实裁切工作面，').replace('灰色为完整真实裁切齿面。','背景为完整真实裁切齿面。').replace('Grey is the complete actual trimmed flank from the same request.','The background is the complete actual trimmed flank from the same request.').replace('Grey is the actual complete trimmed flank;','The background is the actual complete trimmed flank;')
        return _save(fig,output,stem,lang,data['rows'],caption,height_mm=max(90,58*len(data['models'])+25),
            extra=dict(**scale,contact_representation='nodal_dual_cells',interpolation='C0_P1_scalar' if display=='continuous' else 'none',geometry_aspect='equal',
                rendering=display,original_pressure_values_modified=False,
                pressure_domain_background=[pressure_domain_metadata(domain) for domain in domains],
                display_field_diagnostics=[f['diagnostics'] for f in fields],
                tooth_offset=tooth_offset,driving_angle_rad=data['driving_angle_rad'],pinion_torque_nm=data['pinion_torque_nm'],
                projection_definition=data['projection_definition'],physical_area_definition=data['physical_area_definition'],
                face_identification=data['face_identification'],zoom_windows=zoom_records,
                models=[{k:v for k,v in m.items() if k!='cells'} for m in data['models']]))


def plot_refined_fe_method_comparison(fe_result,surface_result,output,*,pinion_request,
        flank_mesh,lang='zh',tooth_offset=0,pressure_limit_mpa=None,stem='15_refined_contact_methods',display='continuous',
        fe_domain_certificate=None,surface_domain_certificate=None):
    """同输入相位/载荷下比较两种连续标量压力；cells保留离散诊断模式。"""
    import numpy as np
    from .halfspace_plotting import halfspace_contact_pattern_data,pressure_colour_scale
    from .workflow_plotting import _paper,_save,_text
    from .pressure_field import reconstruct_fe_pressure_field,reconstruct_halfspace_pressure_field,project_pressure_field,draw_pressure_field,pressure_field_rows
    from .pressure_background import pressure_domain_background,pressure_background_rows,draw_pressure_background,add_pressure_domain_legend,pressure_domain_metadata
    if display not in ('continuous','cells'):raise ValueError('display须为continuous或cells')
    fe=refined_fe_contact_data([fe_result],flank_mesh,pinion_request,labels=['FE'],tooth_offset=tooth_offset)
    hs=halfspace_contact_pattern_data([surface_result],flank_mesh,pinion_request,phase_indices=[0],tooth_offset=tooth_offset)
    phase=hs['phases'][0];fem=fe['models'][0]
    if not math.isclose(phase['driver_rad'],fe['driving_angle_rad'],rel_tol=0.,abs_tol=1e-12) or not math.isclose(
            phase['pinion_torque_nm'],fe['pinion_torque_nm'],rel_tol=1e-12,abs_tol=1e-10):
        raise ValueError('两种实际压力模型须使用同一输入相位与有符号扭矩')
    domains=[pressure_domain_background(result,fe['background_polygons_m'],method=method,tooth_offset=tooth_offset,
        certificate=certificate,pinion_request=pinion_request) for result,method,certificate in ((fe_result,'full_gear_fe',fe_domain_certificate),
        (surface_result,'local_halfspace',surface_domain_certificate))]
    models=[dict(method='full_gear_fe',polygons=[c['polygon_width_height_m'] for c in fem['cells']],
        pressures_pa=[c['pressure_pa'] for c in fem['cells']],force_n=fem['normal_force_n'],
        active_pressure_node_count=fem['active_pressure_node_count'],active_cell_count=fem['loaded_dual_cell_count'],
        area_m2=fem['integration_area_m2'],peak_pressure_mpa=fem['peak_pressure_mpa'],force_symbol='n',
        pressure_definition=fem['pressure_definition'],area_definition='原三维FE三角形中位对偶面积'),
        dict(method='local_halfspace',polygons=phase['polygons_m'],pressures_pa=phase['pressure_pa'],
        force_n=phase['plane_normal_force_n'],active_pressure_node_count=None,active_cell_count=phase['active_cell_count'],
        area_m2=phase['reference_plane_area_m2']*phase['active_cell_count'],
        peak_pressure_mpa=max(phase['pressure_pa'])*1e-6,force_symbol='n,plane',
        pressure_definition=hs['pressure_definition'],area_definition='原参考接触切平面矩形面积')]
    fields=[]
    if display=='continuous':
        fields=[project_pressure_field(reconstruct_fe_pressure_field(fe_result,tooth_offset=tooth_offset),pinion_request,fe_result['mesh_reference_frames'])]
        hs_patch=next(p for p in surface_result['geometry']['patches'] if p['tooth_offset']==tooth_offset)
        fields.append(project_pressure_field(reconstruct_halfspace_pressure_field(surface_result,tooth_offset=tooth_offset),pinion_request,hs_patch['frames']))
    norm,cmap,scale=pressure_colour_scale([p for m in models for p in m['pressures_pa']],pressure_limit_mpa)
    background=np.asarray(fe['background_polygons_m'])*1000;bounds=background.reshape(-1,2)
    all_points=np.concatenate([np.concatenate([np.asarray(p)*1000 for p in m['polygons']]) for m in models])
    if display=='continuous':all_points=np.concatenate([f['points_width_height_m']*1000 for f in fields])
    margin=np.maximum(np.ptp(all_points,axis=0)*.10,.08)
    low=all_points.min(axis=0)-margin;high=all_points.max(axis=0)+margin
    rows=[dict(row,method='full_gear_fe') for row in fe['rows']]
    rows.extend(dict(row,method='local_halfspace') for row in hs['rows']+hs['load_share_rows'])
    for model,field in zip(models,fields):rows.extend(pressure_field_rows(field,label=model['method']))
    for model,domain in zip(models,domains):rows.extend(pressure_background_rows(domain,label=model['method']))
    neighbors=[]
    for offset,field in zip(surface_result['selected_tooth_offsets'],surface_result['fields']):
        boundary=next((p for p in surface_result.get('contact_boundary_checks',[]) if p['tooth_offset']==offset),{})
        neighbors.append(dict(tooth_offset=offset,active_cell_count=int(np.count_nonzero(field['contact_mask'])),
            reference_plane_force_n=float(np.sum(field['forces_n'])),boundary_checks=boundary))
    edge_entries=[p for p in neighbors if p['active_cell_count'] and p['boundary_checks'].get('touches_trimmed_or_model_validity_edge')]
    with _paper(lang) as plt:
        from matplotlib.collections import PolyCollection
        from matplotlib.cm import ScalarMappable
        # 每行左侧小图定位，右侧大图展示完整云图，保持真实坐标比例。
        fig,axes=plt.subplots(2,2,gridspec_kw=dict(width_ratios=[.65,1.8]))
        fig.subplots_adjust(left=.085,right=.855,bottom=.22 if edge_entries else .19,top=.82,wspace=.50,hspace=1.25)
        for i,model in enumerate(models):
            polygons=[np.asarray(p)*1000 for p in model['polygons']]
            known=phase['valid_polygons_m'] if model['method']=='local_halfspace' else model['polygons']
            if model['method']=='full_gear_fe' and display=='continuous':known=np.asarray(fields[i]['points_width_height_m'])[fields[i]['triangles']]
            if model['method']=='full_gear_fe':
                name=_text(lang,'整轮FE','Full-gear FE')
                count=_text(lang,f'{model["active_pressure_node_count"]}节点 / {model["active_cell_count"]}对偶格',
                    f'{model["active_pressure_node_count"]} nodes / {model["active_cell_count"]} dual cells')
            else:
                name=_text(lang,'局部半解析','Local semi-analytical')
                count=_text(lang,f'{model["active_cell_count"]}个承载格',f'{model["active_cell_count"]} loaded cells')
            for j in range(2):
                ax=axes[i,j]
                draw_pressure_background(ax,fe['background_polygons_m'],norm=norm,cmap=cmap,info=domains[i],known_polygons_m=known)
                if display=='continuous':draw_pressure_field(ax,fields[i],norm=norm,cmap=cmap)
                else:ax.add_collection(PolyCollection(polygons,array=np.asarray(model['pressures_pa'])*1e-6,cmap=cmap,norm=norm,
                    edgecolors='none',antialiaseds=False,zorder=2))
                ax.axhline(0.,color='#8C989F',ls=(0,(3,3)),lw=.5,zorder=1)
                if j:
                    ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1])
                    title=name+'\n'+rf'$p_{{\max}}={model["peak_pressure_mpa"]:.0f}$ MPa; '+rf'$F_{{\mathrm{{{model["force_symbol"]}}}}}={model["force_n"]:.0f}$ N'
                else:
                    ax.set_xlim(bounds[:,0].min()-.4,bounds[:,0].max()+.4);ax.set_ylim(bounds[:,1].min()-.25,bounds[:,1].max()+.25)
                    title=_text(lang,'全齿面定位','Location on full flank')
                    add_pressure_domain_legend(ax,domains[i],norm=norm,cmap=cmap,lang=lang)
                ax.set_title(title,loc='left',fontsize=7,pad=5);ax.set_aspect('equal',adjustable='box')
                ax.set_ylabel(_text(lang,'节锥法向高度 h (mm)','Pitch-cone normal\nheight, h (mm)') if j else 'h (mm)')
                ax.annotate('abcd'[i*2+j],(0,1),xycoords='axes fraction',xytext=(-6,18),textcoords='offset points',
                    ha='right',va='bottom',fontsize=9,weight='bold',annotation_clip=False)
        for ax in axes[-1]:ax.set_xlabel(_text(lang,r'齿宽位置 $w_\mathrm{f}$ (mm)','Face-width position\n'+r'$w_\mathrm{f}$ (mm)'))
        fig.suptitle(_text(lang,'弧齿锥齿轮瞬时接触压力','Spiral bevel gear instantaneous contact pressure')+'\n'+
            rf'$\theta_1={fe["driving_angle_rad"]*180/math.pi:.2f}^\circ$, '+rf'$T_1={fe["pinion_torque_nm"]:g}$ N·m, k={tooth_offset}',fontsize=9)
        bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),pad=.025,shrink=.87,
            extend='max' if scale['over_limit_cell_count'] else 'neither')
        bar.set_label(_text(lang,'接触压力 p (MPa)','Contact pressure, p (MPa)'))
        edge_note=''
        for edge in edge_entries:
            edge_note+=_text(lang,
                f'半解析邻齿k={edge["tooth_offset"]}：{edge["active_cell_count"]}个承载格，承载域触及实际裁切边缘；整副压力仍需验证。',
                f'Half-space tooth k={edge["tooth_offset"]}: {edge["active_cell_count"]} loaded cells; support reaches an actual trimmed edge. Full-pair pressure remains unvalidated.')
        if edge_note:fig.text(.03,.025,edge_note,fontsize=6,ha='left',va='bottom',wrap=True)
        caption=_text(lang,
            '同一显式展示几何、输入相位和有符号扭矩下，分别显示整轮FE节点—面集总压力与局部半空间常压力。θ1为小轮输入角，T1为有符号小轮扭矩，k为物理齿对偏移。左列灰色为完整实际裁切齿面，右列使用两方法承载域的共同物理坐标窗显示整个斑点，横纵始终等比例。w沿外端到内端，h为相对节锥的法向高度，h=0为虚线。Fn为FE该齿法向合力，Fn,plane为半空间该齿参考平面合力。两行共用一个线性Pa→MPa色标，蓝低红高，标题标出各自原峰值。对偶格和原半解析格直接着色，不插值、不扩宽、不造椭圆。FE保留整轮多齿结构耦合及原节点力作用点，半空间不含有限齿体/轮体/跨齿结构柔度；两者的原积分面积、齿对载荷分配和压力定义不同，不能把图形相似称为压力精度对标。新展示几何不继承旧案例的强度或外部软件验证结果。',
            'Actual full-gear FE node-to-surface lumped pressure and local half-space constant pressure use the same explicit display geometry, input phase and signed torque. θ1 is pinion input angle, T1 signed pinion torque and k physical tooth-pair offset. Grey in the left column is the complete actual trimmed flank. The right column shows each complete patch within one shared physical coordinate window, always at equal aspect. w runs outer to inner end and h is normal height relative to the pitch cone, dashed at h=0. Fn is the FE tooth normal resultant; Fn,plane is the half-space reference-plane resultant. Both rows share one linear Pa-to-MPa map, blue low and red high, with individual raw peaks labelled. Original dual and half-space cells are directly colored without interpolation, widening or an assumed ellipse. FE retains full-gear multi-tooth structural coupling and original nodal force points; the half-space excludes finite tooth/body/cross-tooth structural compliance. Their original integration areas, load sharing and pressure definitions differ; visual similarity is not validated pressure agreement. This display geometry does not inherit strength or external-software validation from the baseline case.')
        if display=='continuous':caption=_text(lang,
            '同一几何、输入相位和有符号扭矩下的瞬时受载压力。左侧小图定位完整实际齿面中的位置，右侧大图为完整斑点的等比例连续标量云图。FE沿真实原三角形以原节点压力作C0-P1插值，半解析沿实际中心网格作P1插值；不滤波、不拟合斑点或改变原求解值。蓝低红高，两方法共用线性色标；Fn与Fn,plane分别保留原节点和原参考平面合力定义。w_f为从外端向内端的齿宽位置，h为节锥法向高度，虚线h=0。原承载面积与显示场正压支撑面积分列，不用显示重构重新定义力矩或接触求解。两方法柔度不同，图形相似不等于压力精度验证。',
            'Instantaneous loaded pressure at the same geometry, phase and signed torque. Small left panels locate each patch on the actual complete flank; large right panels show the entire patch at equal physical aspect. FE uses C0-P1 interpolation of original nodal pressures on actual FE triangles; the half-space uses P1 interpolation on actual surface-center grids. No filtering, patch fitting or solver-value changes are applied. Blue denotes low and red high pressure on a common linear scale. Fn and Fn,plane retain original nodal and reference-plane resultant definitions. w_f runs from outer to inner end; h is pitch-cone normal height, dashed at h=0. Original loaded areas and display support areas remain distinct; reconstruction does not redefine moments or contact equilibrium. The different compliance models are not mutually validated by similar plots.')
        caption+=' '+edge_note
        caption+=_text(lang,' 零压力使用同一色标p=0的深蓝；覆盖审核已证明的实际分片显示零值蓝底，部分证书中的未决分片单独灰纹，仍保留未完整认证的状态。未知压力不写入0，原压力数组不改。',
            ' Zero pressure is the same p=0 blue as the colour scale. Proven actual pieces have a zero-blue base; a partial certificate keeps unresolved pieces individually grey and hatched and retains incomplete status. Unknown pressures are not set to zero, and source pressure arrays are unchanged.')
        caption=caption.replace('左列灰色为完整实际裁切齿面，','左列为完整实际裁切齿面，').replace('Grey in the left column is the complete actual trimmed flank.','The left column shows the complete actual trimmed flank.')
        return _save(fig,output,stem,lang,rows,caption,height_mm=100,
            extra=dict(**scale,interpolation='C0_P1_scalar' if display=='continuous' else 'none',geometry_aspect='equal',tooth_offset=tooth_offset,
                rendering=display,original_pressure_values_modified=False,
                pressure_domain_background=[pressure_domain_metadata(domain) for domain in domains],
                display_field_diagnostics=[dict(method=m['method'],**f['diagnostics']) for m,f in zip(models,fields)],
                driving_angle_rad=fe['driving_angle_rad'],pinion_torque_nm=fe['pinion_torque_nm'],
                complete_patch_shared_width_height_bounds_m=[(low*.001).tolist(),(high*.001).tolist()],
                models=[{k:v for k,v in m.items() if k not in ('polygons','pressures_pa')} for m in models],
                halfspace_all_tooth_context=neighbors,fe_pair_normal_force_shares_n=fe_result['pair_force_shares'],
                comparison_scope='不同物理柔度/离散模型的同输入展示，非压力精度对标'))
