"""在完整裁切工作齿面上显示真实FE接触斑点，不重新求解或平滑压力。

齿宽/齿高坐标为本轮子午投影：s=r sinδ+z cosδ，w=R_e-s，
h=r cosδ-z sinδ。投影不是等面积展开，法向力仍按原三维面片面积积分。
"""
from __future__ import annotations

import json
import math


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name+'须为有限实数')
    return float(value)


def pinion_flank_coordinates(points_m, pitch_cone_angle_rad, outer_cone_distance_m, *,
                             rotation_rows=None, translation_m=(0., 0., 0.)):
    """将局部点或给定刚体位姿的世界点变为小轮齿宽w、节锥法向高度h，单位m。

    rotation_rows为局部→世界的正交旋转；先逆变换，禁止将世界z误作材料z。
    以外端参考截面w=0，向内端增大；h=0为节锥，正向朝齿顶。
    """
    delta=_number(pitch_cone_angle_rad,'pitch_cone_angle_rad')
    reference=_number(outer_cone_distance_m,'outer_cone_distance_m')
    if not 0<delta<math.pi/2 or reference<=0:
        raise ValueError('子午投影需要正锐角节锥和正外锥距')
    if rotation_rows is None:
        rotation_rows=((1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    if len(rotation_rows)!=3 or any(len(row)!=3 for row in rotation_rows) or len(translation_m)!=3:
        raise ValueError('刚体位姿须含3×3旋转和三维平移')
    rotation=[[_number(x,'rotation_rows') for x in row] for row in rotation_rows]
    translation=[_number(x,'translation_m') for x in translation_m]
    for i in range(3):
        for j in range(3):
            value=sum(rotation[k][i]*rotation[k][j] for k in range(3))
            if abs(value-float(i==j))>1e-9:
                raise ValueError('rotation_rows须为正交旋转，不能夹带缩放')
    determinant=(rotation[0][0]*(rotation[1][1]*rotation[2][2]-rotation[1][2]*rotation[2][1])
        -rotation[0][1]*(rotation[1][0]*rotation[2][2]-rotation[1][2]*rotation[2][0])
        +rotation[0][2]*(rotation[1][0]*rotation[2][1]-rotation[1][1]*rotation[2][0]))
    if abs(determinant-1.)>1e-9:
        raise ValueError('rotation_rows须为右手正交旋转')
    local=[];width=[];height=[]
    for point in points_m:
        if len(point)!=3:raise ValueError('points_m须为三维点列表')
        relative=[_number(point[k],'point_m')-translation[k] for k in range(3)]
        value=[sum(rotation[j][i]*relative[j] for j in range(3)) for i in range(3)]
        radius=math.hypot(value[0],value[1]);z=value[2]
        local.append(value)
        width.append(reference-radius*math.sin(delta)-z*math.cos(delta))
        height.append(radius*math.cos(delta)-z*math.sin(delta))
    return dict(local_points_m=local,width_position_m=width,normal_height_m=height,
                projection_definition='w=R_e-r*sin(delta)-z*cos(delta); h=r*cos(delta)-z*sin(delta)')


def select_pinion_left_flank(mesh, request):
    """按展成点的实际左右侧几何残差识别CAD工作面，不猜测BRep面号顺序。

    输入必须为同案例裁切侧面的真实三角网格；每个BRep面取分布样点，
    比较该点半径/机床截面处的两侧展成坐标。该门限只用于面身份识别。
    """
    from .generated_surface import _request
    from .generated_cutter import evaluate_generated_point
    request=_request(request)
    if mesh.get('units')!='m':raise ValueError('工作面网格须明确使用m')
    vertices=mesh['vertices_m'];triangles=mesh['triangles'];faces=mesh.get('face_indices',[])
    if not triangles or len(faces)!=len(triangles):
        raise ValueError('完整裁切齿面需要每个实际三角形的BRep面号')
    coordinates=pinion_flank_coordinates(vertices,request.machine.pitch_cone_angle_rad,
                                        request.machine.reference_cone_distance_m)
    tolerance=max(2e-7,request.machine.reference_cone_distance_m*1e-6)
    groups={}
    for index,(triangle,face) in enumerate(zip(triangles,faces)):
        if isinstance(face,bool) or not isinstance(face,int) or face<0 or len(triangle)!=3:
            raise ValueError('BRep面号和三角形索引不合法')
        if any(isinstance(i,bool) or not isinstance(i,int) or not 0<=i<len(vertices) for i in triangle):
            raise ValueError('齿面三角形顶点索引越界')
        groups.setdefault(face,[]).append(index)
    modified_surfaces = None
    if request.tip_relief_m:
        # 显式修缘后不能再以未修缘刀盘点识别齿侧；比较同一实际采样的NURBS。
        from .generated_surface import build_generated_surface_grid
        from .cad import build_generated_shell
        from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
        from OCP.gp import gp_Pnt
        sampling=mesh.get('metadata',{}).get('source_surface_sampling')
        if not isinstance(sampling,dict) or set(sampling)!={'width_count','profile_count'}:
            raise ValueError('修缘工作面识别需要实际CAD的source_surface_sampling，不能使用未修缘点代替')
        shell=build_generated_shell(build_generated_surface_grid(request,**sampling))
        modified_surfaces={side:shell.surfaces[side+'_flank'] for side in ('left','right')}
    selection=[];diagnostics=[]
    for face,indices in sorted(groups.items()):
        unique=sorted({i for index in indices for i in triangles[index]})
        count=min(9,len(unique))
        samples=[unique[round(i*(len(unique)-1)/(count-1))] for i in range(count)] if count>1 else unique
        errors={'left':[],'right':[]}
        for index in samples:
            point=vertices[index];radius=math.hypot(point[0],point[1])
            station=request.machine.reference_cone_distance_m-coordinates['width_position_m'][index]
            value=None if modified_surfaces is not None else evaluate_generated_point(request.machine,2*radius,station)
            for side in errors:
                if modified_surfaces is None:
                    expected=value['flanks'][side]['point_m']
                    error=math.sqrt(sum((point[k]-expected[k])**2 for k in range(3)))
                else:
                    projection=GeomAPI_ProjectPointOnSurf(gp_Pnt(*(x*1000 for x in point)),modified_surfaces[side])
                    error=projection.LowerDistance()*.001 if projection.NbPoints() else math.inf
                errors[side].append(error)
        means={side:sum(values)/len(values) for side,values in errors.items()}
        side=min(means,key=means.get)
        if max(errors[side])>tolerance or means['right' if side=='left' else 'left']<=means[side]+tolerance:
            raise ValueError('无法以真实展成点唯一识别工作面，未按面号补造分类')
        diagnostics.append(dict(face_index=face,identified_side=side,sampled_vertices=len(samples),
            maximum_matching_distance_m=max(errors[side]),mean_opposite_distance_m=means['right' if side=='left' else 'left']))
        if side=='left':selection.extend(indices)
    if not selection:raise ValueError('实际CAD选区中未识别到小轮左工作面')
    return dict(triangle_indices=selection,coordinates=coordinates,face_identification=diagnostics,
                geometric_identification_tolerance_m=tolerance)


def contact_pattern_data(bundle, mesh, *, phase_indices=None, tooth_offset=0):
    """提取同一物理齿在实际已解姿态的面片与完整齿面，保留三维面积和载荷。"""
    from .workflow_plotting import _validate,_loaded_results,_contact_values
    _validate(bundle)
    if bundle.get('kind',bundle.get('case',{}).get('kind'))!='spiral':
        raise ValueError('当前完整齿面斑点入口需要显式弧齿制造输入')
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):
        raise ValueError('tooth_offset须为整数物理齿对偏移')
    request=bundle['case']['manufacturing']['pinion']['request']
    machine=request['machine']
    delta=machine['pitch_cone_angle_rad'];reference=machine['reference_cone_distance_m']
    selected=select_pinion_left_flank(mesh,request)
    results=_loaded_results(bundle)
    if phase_indices is None:phase_indices=(0,len(results)//2,len(results)-1)
    phase_indices=list(phase_indices)
    if len(phase_indices)!=3 or len(set(phase_indices))!=3 or any(isinstance(i,bool) or not isinstance(i,int) or not 0<=i<len(results) for i in phase_indices):
        raise ValueError('接触斑点需要三个不同且实际存在的相位索引')
    rows=[];background=[];phases=[]
    coordinates=selected['coordinates']
    for index in selected['triangle_indices']:
        triangle=mesh['triangles'][index]
        polygon=[[coordinates['width_position_m'][i],coordinates['normal_height_m'][i]] for i in triangle]
        background.append(polygon)
        rows.append(dict(row_kind='actual_flank_triangle',cad_triangle_index=index,
            cad_face_index=mesh['face_indices'][index],polygon_width_height_m_json=json.dumps(polygon),
            polygon_local_m_json=json.dumps([mesh['vertices_m'][i] for i in triangle])))
    for index in phase_indices:
        result=results[index]
        contacts=_contact_values(result,result.get('contact_representation'))
        if result.get('contact_representation')!='facets':
            raise ValueError('完整齿面斑点仅接受具有真实边界的FE面片，积分点不涂成整面')
        frame=result['mesh_reference_frames']
        facets=[]
        for number,(contact,original) in enumerate(zip(contacts,result['contacts'])):
            if not contact['active'] or contact['force_n']<=0 or contact['tooth_offset']!=tooth_offset:continue
            world=json.loads(contact['polygon_world_m_json'])
            converted=pinion_flank_coordinates(world,delta,reference,
                rotation_rows=frame['pinion_rotation_rows'],translation_m=frame.get('translation_m',(0.,0.,0.)))
            polygon=list(map(list,zip(converted['width_position_m'],converted['normal_height_m'])))
            row=dict(row_kind='loaded_facet',phase_index=index,driver_rad=result['driving_angle_rad'],
                pinion_torque_nm=result['pinion_torque_nm'],tooth_offset=tooth_offset,contact_index=number,
                pressure_pa=contact['pressure_pa'],physical_area_m2=contact['physical_area_or_weight_m2'],
                force_n=contact['force_n'],polygon_width_height_m_json=json.dumps(polygon),
                polygon_world_m_json=json.dumps(world),polygon_local_m_json=json.dumps(converted['local_points_m']),
                active=True,first_facet_index=original['first_facet_index'])
            rows.append(row);facets.append(row)
        total=sum(row['force_n'] for row in facets)
        expected=result['pair_force_shares'].get(str(tooth_offset),result['pair_force_shares'].get(tooth_offset,0.))
        if not math.isclose(total,expected,rel_tol=1e-9,abs_tol=1e-8):
            raise ValueError('所选物理齿接触力与实际齿对分担不一致')
        phase=dict(phase_index=index,driver_rad=result['driving_angle_rad'],pinion_torque_nm=result['pinion_torque_nm'],
            tooth_offset=tooth_offset,loaded_facet_count=len(facets),normal_force_n=total,facets=facets)
        phases.append(phase)
        rows.append({key:value for key,value in dict(row_kind='phase_state',**phase).items() if key!='facets'})
    if not any(phase['facets'] for phase in phases):
        raise ValueError('三个姿态所选齿均无承载面片，未生成接触斑点占位图')
    return dict(rows=rows,background_polygons_m=background,phases=phases,
        tooth_offset=tooth_offset,outer_cone_distance_m=reference,pitch_cone_angle_rad=delta,
        face_width_m=request['face_width_m'],face_identification=selected['face_identification'],
        geometric_identification_tolerance_m=selected['geometric_identification_tolerance_m'],
        projection_definition=coordinates['projection_definition'],
        physical_area_definition='原三维FE面片面积；子午投影不重新定义接触面积')


def plot_workflow_contact_pattern(bundle, output, *, lang='zh', data_root=None, phase_indices=None, tooth_offset=0):
    """完整实际工作面上的同齿接触斑点三姿态图；不平滑、不补造连续斑点边界。"""
    from .workflow_plotting import _mesh,_paper,_panel,_save,_text
    import numpy as np
    mesh=_mesh(bundle['cad']['flank_mesh'],data_root or output)
    data=contact_pattern_data(bundle,mesh,phase_indices=phase_indices,tooth_offset=tooth_offset)
    all_facets=[row for phase in data['phases'] for row in phase['facets']]
    maximum=max(row['pressure_pa'] for row in all_facets)*1e-6
    background=np.asarray(data['background_polygons_m'])*1000
    bounds=background.reshape(-1,2)
    with _paper(lang) as plt:
        from matplotlib.collections import PolyCollection
        from .halfspace_plotting import pressure_colour_scale
        fig,axes=plt.subplots(3,1,sharex=True,sharey=True,layout='constrained')
        norm,pressure_cmap,scale=pressure_colour_scale([row['pressure_pa'] for row in all_facets])
        for ax,phase,letter in zip(axes,data['phases'],'abc'):
            ax.add_collection(PolyCollection(background,facecolors='#EDF0F2',edgecolors='#CBD2D6',linewidths=.13,zorder=0))
            facets=phase['facets']
            polygons=[np.asarray(json.loads(row['polygon_width_height_m_json']))*1000 for row in facets]
            if facets:
                artist=PolyCollection(polygons,array=np.asarray([row['pressure_pa']*1e-6 for row in facets]),
                    cmap=pressure_cmap,norm=norm,edgecolors='#2D3A40',linewidths=.38,zorder=2)
                ax.add_collection(artist)
            else:
                ax.text(.5,.45,_text(lang,'该齿没有承载面片','No loaded facet on this tooth'),ha='center',transform=ax.transAxes)
            ax.axhline(0.,color='#86939D',ls=(0,(3,3)),lw=.55,zorder=1)
            angle=phase['driver_rad']*180/math.pi
            facet_label=f'{len(facets)} loaded '+('facet' if len(facets)==1 else 'facets')
            title=(rf'$\theta_1={angle:+.2f}^\circ$'+f'   k={tooth_offset}   '+
                   rf'$F_{{\mathrm{{n}},{tooth_offset}}}={phase["normal_force_n"]:.1f}$ N'+
                   _text(lang,f'   {len(facets)}个承载面片','   '+facet_label))
            ax.set_title(title,loc='left',x=.02,fontsize=8,pad=5)
            ax.set_xlim(bounds[:,0].min()-.6,bounds[:,0].max()+.6)
            ax.set_ylim(bounds[:,1].min()-.4,bounds[:,1].max()+.4)
            ax.set_aspect('equal',adjustable='box')
            ax.set_ylabel(_text(lang,'节锥法向高度 h (mm)','Pitch-cone normal\nheight, h (mm)'))
            _panel(ax,letter)
        axes[-1].set_xlabel(_text(lang,r'齿宽位置 $w_\mathrm{f}$：外端 → 内端 (mm)',
            r'Face-width position, $w_\mathrm{f}$: outer → inner end (mm)'))
        from matplotlib.cm import ScalarMappable
        colorbar=fig.colorbar(ScalarMappable(norm=norm,cmap=pressure_cmap),ax=axes,shrink=.92,pad=.035,aspect=35)
        colorbar.set_label(_text(lang,'面片平均接触压力 p̄ (MPa)','Facet-mean contact pressure, p̄ (MPa)'))
        caption=_text(lang,
            '固定同一小轮物理齿对k的完整裁切左工作面接触斑点。灰色是实际CAD工作面，着色仅为该齿实际承载FE面片；邻齿不折叠到本图。w沿外端到内端，h相对节锥；h=0虚线是节锥参考。子午投影不是等面积展开，压力与力仍用原三维面积。未平滑或绘制假想椭圆；压力未完成网格/实验验证。',
            'Contact patterns of the same physical pinion tooth pair k on its complete trimmed left flank. Grey is the actual CAD flank; only loaded FE facets on that tooth are colored. Neighboring teeth are not folded onto it. w runs from the outer to the inner end; h is relative to the pitch cone, with h=0 dashed. This meridional projection is not area-preserving; pressure and force retain the original 3D areas. No smoothing or assumed ellipse is used. Pressure has not established mesh or experimental validation.')
        torque=data['phases'][0]['pinion_torque_nm']
        caption+=_text(lang,f' 所选齿对k={tooth_offset}，小轮有符号扭矩T1={torque:g}N·m。',
            f' Selected tooth pair k={tooth_offset}; signed pinion torque T1={torque:g}N·m.')
        return _save(fig,output,'09_contact_pattern',lang,data['rows'],caption,height_mm=150,
            extra=dict(**scale,case_fingerprint=bundle.get('case_fingerprint'),contact_representation='facets',
                phase_indices=[phase['phase_index'] for phase in data['phases']],
                driving_angles_rad=[phase['driver_rad'] for phase in data['phases']],
                tooth_offset=tooth_offset,loaded_facet_counts=[len(phase['facets']) for phase in data['phases']],
                projection_definition=data['projection_definition'],physical_area_definition=data['physical_area_definition'],
                face_identification=data['face_identification'],
                geometric_identification_tolerance_m=data['geometric_identification_tolerance_m']))
