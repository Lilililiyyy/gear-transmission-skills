"""细网格半解析表面接触的双语科研图；原求解值与连续显示重构分开。

原压力格及实际四角作为源数据保留；默认沿实际中心网格重构连续标量，
背景采用相同请求的完整裁切工作面。
局部半空间不含齿体／轮体／轴承及跨齿结构柔度，闭合量不能冒充全局TE。
"""
from __future__ import annotations

import json
import math


def _np():
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError('半空间接触图需要可选NumPy与Matplotlib') from error
    return np


def _phase_indices(count, indices):
    if indices is None:
        indices = tuple(range(count)) if count <= 3 else (0, count//2, count-1)
    indices = tuple(indices)
    if not 1 <= len(indices) <= 3 or len(set(indices)) != len(indices) or any(
            isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < count for i in indices):
        raise ValueError('斑点图需要一至三个不同且实际存在的相位索引')
    return indices


def _boundary_segments(mask, corners):
    """按单元边界原样提取边线，不用插值等值线或凸包补齐。"""
    np = _np()
    segments = []
    ny, nx = mask.shape
    for j, i in np.argwhere(mask):
        for dj, di, a, b in ((-1, 0, (j, i), (j, i+1)),
                              (0, 1, (j, i+1), (j+1, i+1)),
                              (1, 0, (j+1, i+1), (j+1, i)),
                              (0, -1, (j+1, i), (j, i))):
            y, x = j+dj, i+di
            if not (0 <= y < ny and 0 <= x < nx and mask[y, x]):
                segments.append([corners[a], corners[b]])
    return segments


def _project_corners(points, frames, delta, reference):
    from .contact_pattern import pinion_flank_coordinates
    np = _np()
    points = np.asarray(points, dtype=float)
    finite = np.isfinite(points).all(axis=-1)
    result = np.full(points.shape[:-1]+(2,), np.nan)
    value = pinion_flank_coordinates(points[finite], delta, reference,
        rotation_rows=frames['pinion_rotation_rows'], translation_m=frames.get('translation_m', (0., 0., 0.)))
    result[finite] = np.column_stack((value['width_position_m'], value['normal_height_m']))
    return result


def halfspace_contact_pattern_data(results, flank_mesh, pinion_request, *, phase_indices=None, tooth_offset=0):
    """校验已收敛真实压力格，返回同一物理齿一至三相位的真实多边形及源数据。

    results须保留geometry、selected_tooth_offsets和fields；四角必须来自
    geometry中的first_corner_points_world_m，不以单元中心推测物理边界。
    """
    from .contact_pattern import select_pinion_left_flank
    from .generated_surface import _request
    np = _np()
    if isinstance(tooth_offset, bool) or not isinstance(tooth_offset, int):
        raise ValueError('tooth_offset须为整数物理齿偏移')
    if not isinstance(results, (list, tuple)):
        raise ValueError('results须为实际半空间求解结果列表')
    indices = _phase_indices(len(results), phase_indices)
    request = _request(pinion_request)
    selected = select_pinion_left_flank(flank_mesh, request)
    coordinates = selected['coordinates']
    background, rows, shares, phases = [], [], [], []
    for index in selected['triangle_indices']:
        triangle = flank_mesh['triangles'][index]
        polygon = [[coordinates['width_position_m'][i], coordinates['normal_height_m'][i]] for i in triangle]
        background.append(polygon)
        rows.append(dict(row_kind='actual_flank_triangle', cad_triangle_index=index,
            polygon_width_height_m_json=json.dumps(polygon)))
    angles, torques = [], []
    for index, result in enumerate(results):
        if not isinstance(result, dict) or not result.get('converged') or not result.get('fields'):
            raise ValueError('失败或候选压力不能作为物理斑点图')
        geometry = result.get('geometry')
        if not geometry or 'selected_tooth_offsets' not in result:
            raise ValueError('需要实际有限曲面半空间结果及其geometry，不能传入旧FE面片')
        angle, torque = float(geometry['driving_angle_rad']), float(geometry['pinion_torque_nm'])
        if not math.isfinite(angle) or not math.isfinite(torque):
            raise ValueError('相位与扭矩须有限')
        angles.append(angle); torques.append(torque)
        offsets = result['selected_tooth_offsets']
        if len(offsets) != len(result['fields']) or len(set(offsets)) != len(offsets):
            raise ValueError('齿对偏移与压力场布局不一致')
        total = sum(float(np.asarray(field['forces_n']).sum()) for field in result['fields'])
        for offset, field in zip(offsets, result['fields']):
            force = float(np.asarray(field['forces_n']).sum())
            shares.append(dict(row_kind='pair_load_share', phase_index=index, driver_rad=angle,
                pinion_torque_nm=torque, tooth_offset=offset, plane_normal_force_n=force,
                plane_normal_force_share=force/total if total else 0., active_cell_count=field['active_cell_count'],
                local_halfspace_closure_rad=result.get('theta_rad'),
                local_halfspace_incremental_stiffness_nm_per_rad=result.get('incremental_stiffness_nm_per_rad')))
    if len(set(angles)) != len(angles):
        raise ValueError('主动相位重复，不能将重复结果画成相位推进')
    if not all(math.isclose(value, torques[0], rel_tol=1e-12, abs_tol=1e-12) for value in torques):
        raise ValueError('多个相位须使用同一有符号小轮扭矩')
    for index in indices:
        result = results[index]; geometry = result['geometry']
        if tooth_offset not in result['selected_tooth_offsets']:
            raise ValueError('所选同一物理齿没有实际已求解压力场')
        field = result['fields'][result['selected_tooth_offsets'].index(tooth_offset)]
        candidates = [patch for patch in geometry['patches'] if patch['tooth_offset'] == tooth_offset]
        if len(candidates) != 1:
            raise ValueError('所选物理齿几何patch缺失或重复')
        patch = candidates[0]
        if 'first_corner_points_world_m' not in patch or 'corner_valid_mask' not in patch:
            raise ValueError('缺少实际曲面单元四角，拒绝由中心半径补造压力边界')
        pressure = np.asarray(field['pressure_pa'], dtype=float)
        force = np.asarray(field['forces_n'], dtype=float)
        valid, active = np.asarray(field['valid_mask']), np.asarray(field['contact_mask'])
        shape = pressure.shape
        if pressure.ndim != 2 or any(value.shape != shape for value in (force, valid, active)):
            raise ValueError('压力／力／掩码的二维布局不一致')
        if valid.dtype != np.bool_ or active.dtype != np.bool_ or np.any(active & ~valid):
            raise ValueError('接触掩码须为有效域内的布尔单元')
        if not np.array_equal(valid, np.asarray(patch['valid_mask'])):
            raise ValueError('求解有效域与实际裁切几何不一致')
        points = np.asarray(patch['first_corner_points_world_m'], dtype=float)
        corner_valid = np.asarray(patch['corner_valid_mask'])
        expected = (shape[0]+1, shape[1]+1)
        if points.shape != expected+(3,) or corner_valid.shape != expected or corner_valid.dtype != np.bool_:
            raise ValueError('实际曲面四角布局须为(ny+1,nx+1,3)及布尔掩码')
        full_corners = corner_valid[:-1, :-1] & corner_valid[:-1, 1:] & corner_valid[1:, 1:] & corner_valid[1:, :-1]
        finite_corners = np.isfinite(points).all(axis=-1)
        full_finite = finite_corners[:-1, :-1] & finite_corners[:-1, 1:] & finite_corners[1:, 1:] & finite_corners[1:, :-1]
        if np.any(valid & ~(full_corners & full_finite)):
            raise ValueError('有效压力格含未知／裁切域外的真实四角')
        area = float(field['cell_area_m2'])
        if not math.isfinite(area) or area <= 0 or not np.isfinite(pressure[valid]).all() or np.any(pressure[valid] < 0):
            raise ValueError('有效压力与参考面积须为有限非负量')
        if not np.isfinite(force).all() or np.any(force < 0) or np.any(force[~valid] != 0):
            raise ValueError('无效域不得承载，单元力须有限非负')
        if not np.allclose(pressure[valid]*area, force[valid], rtol=2e-12, atol=1e-10):
            raise ValueError('压力乘参考切平面面积与单元力不一致')
        projected = _project_corners(points, patch['frames'], request.machine.pitch_cone_angle_rad,
                                     request.machine.reference_cone_distance_m)
        polygons, values, valid_polygons = [], [], []
        used_corners = np.zeros(expected, dtype=bool)
        for j, i in np.argwhere(valid):
            quad_index = ((j, i), (j, i+1), (j+1, i+1), (j+1, i))
            quad = [projected[q].tolist() for q in quad_index]
            valid_polygons.append(quad)
            for q in quad_index:
                used_corners[q] = True
            row = dict(row_kind='pressure_cell', phase_index=index, driver_rad=angles[index],
                pinion_torque_nm=torques[index], tooth_offset=tooth_offset, grid_y=int(j), grid_x=int(i),
                pressure_pa=float(pressure[j, i]), force_n=float(force[j, i]),
                reference_plane_area_m2=area, active=bool(active[j, i]), valid=True,
                domain_state='known_positive' if pressure[j,i]>0 else 'known_zero')
            rows.append(row)
            if active[j, i] and pressure[j, i] > 0:
                polygons.append(quad); values.append(float(pressure[j, i]))
        # 共享四角只存一次，避免细网格源CSV重复四份坐标超过仓库文件上限。
        for j, i in np.argwhere(used_corners):
            rows.append(dict(row_kind='actual_surface_corner', phase_index=index, tooth_offset=tooth_offset,
                grid_y=int(j), grid_x=int(i), world_x_m=float(points[j, i, 0]),
                world_y_m=float(points[j, i, 1]), world_z_m=float(points[j, i, 2]),
                width_position_m=float(projected[j, i, 0]), normal_height_m=float(projected[j, i, 1])))
        boundary = next((entry for entry in result.get('contact_boundary_checks', [])
                         if entry['tooth_offset'] == tooth_offset), {})
        phase = dict(phase_index=index, driver_rad=angles[index], pinion_torque_nm=torques[index],
            polygons_m=polygons, valid_polygons_m=valid_polygons, pressure_pa=values, active_cell_count=int(active.sum()),
            valid_cell_count=int(valid.sum()), plane_normal_force_n=float(force.sum()),
            window_boundary_m=_boundary_segments(valid, projected),
            contact_boundary_m=_boundary_segments(active, projected), boundary_checks=boundary,
            grid_shape=list(shape), reference_plane_area_m2=area)
        rows.append(dict(row_kind='phase_state', phase_index=index, driver_rad=angles[index],
            plane_normal_force_n=phase['plane_normal_force_n'], active_cell_count=phase['active_cell_count'],
            valid_cell_count=phase['valid_cell_count'], grid_shape_yx_json=json.dumps(list(shape)),
            boundary_checks_json=json.dumps(boundary)))
        phases.append(phase)
    if not any(phase['active_cell_count'] for phase in phases):
        raise ValueError('所选齿在给定相位均未承载，未生成占位斑点图')
    return dict(rows=rows, load_share_rows=shares, phases=phases, background_polygons_m=background,
        tooth_offset=tooth_offset, projection_definition=coordinates['projection_definition'],
        pressure_definition='参考接触切平面上的单元常压力，积分面积dx*dy；不是子午投影面积',
        model='局部半空间半解析表面接触；不含齿体／轮体／跨齿结构柔度')


def pressure_colour_scale(pressures_pa, pressure_limit_mpa=None):
    """线性色标及显式超量程颜色；仅改变显示，不裁剪原压力数组。"""
    from matplotlib.colors import Normalize
    from .publication_style import PRESSURE_BLUE, PRESSURE_RED, PRESSURE_COLOUR_POLICY, contact_pressure_colormap
    np=_np();values=np.asarray(pressures_pa,dtype=float)
    if values.size==0 or not np.isfinite(values).all() or (values<0).any() or values.max()<=0:
        raise ValueError('色标需要实际有限非负且至少一个为正的压力，单位Pa')
    maximum=float(values.max()*1e-6)
    if pressure_limit_mpa is None:limit=maximum
    else:
        if isinstance(pressure_limit_mpa,bool) or not isinstance(pressure_limit_mpa,(int,float)) or not math.isfinite(pressure_limit_mpa) or pressure_limit_mpa<=0:
            raise ValueError('显示压力上限须为正有限MPa，不代表物理压力截断')
        limit=float(pressure_limit_mpa)
    # 压力始终线性编码为蓝低红高；不通过色带或显示范围修改原数值。
    cmap=contact_pressure_colormap()
    return Normalize(0.,limit,clip=False),cmap,dict(display_pressure_limit_mpa=limit,
        raw_maximum_pressure_mpa=maximum,over_limit_cell_count=int((values*1e-6>limit).sum()),
        colour_mapping='linear',colour_map=cmap.name,low_colour=PRESSURE_BLUE,high_colour=PRESSURE_RED,
        over_colour=PRESSURE_RED,pressure_values_clipped=False,colour_policy=PRESSURE_COLOUR_POLICY)


def plot_halfspace_contact_pattern(results, output, *, pinion_request, flank_mesh, lang='zh',
                                   phase_indices=None, tooth_offset=0, pressure_limit_mpa=None, zoom=False,display='continuous',
                                   domain_certificates=None):
    """连续显示原压力；完整蓝底须有与结果逐一对应的覆盖证书。"""
    from .workflow_plotting import _paper, _panel, _save, _text
    from .pressure_field import reconstruct_halfspace_pressure_field,project_pressure_field,draw_pressure_field,pressure_field_rows
    from .pressure_background import pressure_domain_background,pressure_background_rows,draw_pressure_background,add_pressure_domain_legend,pressure_domain_metadata
    np = _np()
    if display not in ('continuous','cells'):raise ValueError('display须为continuous或cells')
    data = halfspace_contact_pattern_data(results, flank_mesh, pinion_request,
                                         phase_indices=phase_indices, tooth_offset=tooth_offset)
    certificates=[None]*len(results) if domain_certificates is None else list(domain_certificates)
    if len(certificates)!=len(results):raise ValueError('覆盖证书须与半解析结果逐一对应')
    domains=[pressure_domain_background(results[phase['phase_index']],data['background_polygons_m'],method='local_halfspace',
        tooth_offset=tooth_offset,certificate=certificates[phase['phase_index']],pinion_request=pinion_request) for phase in data['phases']]
    for phase,domain in zip(data['phases'],domains):data['rows'].extend(pressure_background_rows(domain,label=f'phase_{phase["phase_index"]:04d}'))
    fields=[]
    if display=='continuous':
        for phase in data['phases']:
            result=results[phase['phase_index']]
            patch=next(p for p in result['geometry']['patches'] if p['tooth_offset']==tooth_offset)
            field=project_pressure_field(reconstruct_halfspace_pressure_field(result,tooth_offset=tooth_offset),pinion_request,patch['frames'])
            fields.append(field);data['rows'].extend(pressure_field_rows(field,label=f'phase_{phase["phase_index"]:04d}'))
    background = np.asarray(data['background_polygons_m'])*1000
    bounds = background.reshape(-1, 2)
    if not isinstance(zoom,bool):raise ValueError('zoom须为布尔值')
    norm,cmap,scale=pressure_colour_scale([p for phase in data['phases'] for p in phase['pressure_pa']],pressure_limit_mpa)
    with _paper(lang) as plt:
        from matplotlib.collections import PolyCollection, LineCollection
        from matplotlib.cm import ScalarMappable
        count = len(data['phases'])
        fig, axes = plt.subplots(count, 2 if zoom else 1, squeeze=False,
            gridspec_kw=dict(width_ratios=[.75,1.8]) if zoom else None)
        fig.subplots_adjust(left=.085,right=.855,bottom=.19,top=.81,wspace=.50,hspace=1.15)
        fig.suptitle(_text(lang, '半解析瞬时接触压力',
            'Semi-analytical instantaneous contact pressure'), fontsize=9)
        for number,phase in enumerate(data['phases']):
            polygons = np.asarray(phase['polygons_m'])*1000
            angle = phase['driver_rad']*180/math.pi
            suffix = _text(lang, f"  {phase['active_cell_count']}个承载单元", f"  {phase['active_cell_count']} loaded cells")
            peak=max(phase['pressure_pa'],default=0.)*1e-6
            if not zoom:suffix+=rf'  $p_{{\max}}={peak:.0f}$ MPa'
            for column,ax in enumerate(axes[number]):
                draw_pressure_background(ax,data['background_polygons_m'],norm=norm,cmap=cmap,info=domains[number],
                    known_polygons_m=phase['valid_polygons_m'])
                if len(polygons):
                    if display=='continuous':draw_pressure_field(ax,fields[number],norm=norm,cmap=cmap)
                    else:
                        ax.add_collection(PolyCollection(polygons,array=np.asarray(phase['pressure_pa'])*1e-6,
                            cmap=cmap,norm=norm,edgecolors='none',antialiaseds=False,zorder=2))
                        ax.add_collection(LineCollection(np.asarray(phase['contact_boundary_m'])*1000,
                            colors='#1E373F',linewidths=.3,zorder=3))
                if phase['window_boundary_m'] and display=='cells':
                    ax.add_collection(LineCollection(np.asarray(phase['window_boundary_m'])*1000,
                        colors='#87939A',linewidths=.45,zorder=1))
                ax.axhline(0.,color='#A3ADB3',ls=(0,(3,3)),lw=.45,zorder=0)
                if column and len(polygons):
                    points=polygons.reshape(-1,2);span=np.ptp(points,axis=0);margin=np.maximum(span*.10,.10)
                    if display=='continuous':
                        points=fields[number]['points_width_height_m']*1000;span=np.ptp(points,axis=0);margin=np.maximum(span*.10,.10)
                    ax.set_xlim(points[:,0].min()-margin[0],points[:,0].max()+margin[0])
                    ax.set_ylim(points[:,1].min()-margin[1],points[:,1].max()+margin[1])
                    title=_text(lang,'完整接触斑点（等比例）','Complete contact patch (equal aspect)')+'\n'+rf'$p_{{\max}}={peak:.0f}$ MPa'
                    if display=='continuous':title+='; '+rf'$F_{{\mathrm{{n,plane}}}}={phase["plane_normal_force_n"]:.0f}$ N'
                    ax.set_title(title,fontsize=7,loc='left',pad=5)
                else:
                    ax.set_xlim(bounds[:,0].min()-.4,bounds[:,0].max()+.4)
                    ax.set_ylim(bounds[:,1].min()-.25,bounds[:,1].max()+.25)
                    title=rf'$\theta_1={angle:+.2f}^\circ$'+f'  k={tooth_offset}\n'+rf'$F_{{\mathrm{{n,plane}},{tooth_offset}}}={phase["plane_normal_force_n"]:.1f}$ N'+suffix
                    if display=='continuous' and zoom:title=_text(lang,'全齿面定位','Location on full flank')+'\n'+rf'$\theta_1={angle:+.2f}^\circ$, k={tooth_offset}'
                    ax.set_title(title,loc='left',x=.04,fontsize=7,pad=5)
                    add_pressure_domain_legend(ax,domains[number],norm=norm,cmap=cmap,lang=lang)
                ax.set_aspect('equal',adjustable='box')
                ax.set_ylabel('h (mm)' if display=='continuous' and zoom and not column else
                    _text(lang,'节锥法向高度 h (mm)','Pitch-cone normal\nheight, h (mm)'))
                ax.annotate('abcdef'[number*axes.shape[1]+column],(0,1),xycoords='axes fraction',
                    xytext=(-6,18),textcoords='offset points',ha='right',va='bottom',fontsize=9,weight='bold',annotation_clip=False)
        for ax in axes[-1]:ax.set_xlabel(_text(lang,r'齿宽位置 $w_\mathrm{f}$ (mm)','Face-width position\n'+r'$w_\mathrm{f}$ (mm)'))
        colorbar = fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),shrink=.90,pad=.03,aspect=36,
            extend='max' if scale['over_limit_cell_count'] else 'neither')
        colorbar.set_label(_text(lang,'接触压力 p (MPa)','Contact pressure, p (MPa)'))
        if scale['over_limit_cell_count']:
            colorbar.ax.set_title(_text(lang,f'红色 >{scale["display_pressure_limit_mpa"]:g}',
                f'Red >{scale["display_pressure_limit_mpa"]:g}'),fontsize=6,pad=7)
        caption = _text(lang, f'固定同一小轮物理齿k的{count}个已接受相位接触压力。',
            f'Contact pressure on the same physical pinion tooth k at {count} accepted '+('phase. ' if count == 1 else 'phases. '))
        caption += _text(lang,
            '灰色为完整实际裁切工作面；颜色为在真实双面射线求交间隙上重新求解的细网格常压力。深色阶梯线是承载单元边界，灰色边线是严格有效压力域边界；没有插值旧FE压力、Gouraud平滑或假想椭圆。齿宽坐标w沿外端到内端，高度h相对节锥；投影不是等面积展开。压力按参考接触切平面面积积分，Fn,plane,k为所选齿对参考平面上的法向合力。模型不含齿体、轮体、轴承与跨齿结构柔度，局部闭合量不能解释为整轮FE总TE。实际案例的压力网格及实验验证尚未完成。',
            'Grey is its complete actual trimmed flank. Fine constant-cell pressures are newly solved from ray intersections with both actual surfaces. Dark step lines bound loaded cells; grey lines bound the strict valid pressure domain. No interpolation of former FE pressures, Gouraud smoothing or assumed ellipse is used. w runs outer to inner end; h is relative to the pitch cone. The projection is not area-preserving. Pressure integrates over the reference contact plane; Fn,plane,k is the selected pair\'s reference-plane normal resultant. The model excludes tooth/body/bearing and cross-tooth structural compliance; local closure is not the full-gear FE total TE. Actual-case pressure mesh and experimental validation remain open.')
        caption+=_text(lang,
            f' 本图为瞬时接触压力，不是整周期累计着色印迹。显示上限{scale["display_pressure_limit_mpa"]:g}MPa，超过上限的格以红色及色标上端箭头明示；原始峰值保留在各姿态标题和源表中。局部放大保持横纵等比例，未加宽压力条带。',
            f' These are instantaneous pressures, not a cycle-accumulated marking pattern. The display limit is {scale["display_pressure_limit_mpa"]:g}MPa; over-range cells are red with an extended colorbar. Raw peaks remain in panel titles and source data. Local views retain equal aspect without widening the patch.')
        if display=='continuous':caption=_text(lang,
            f'同一物理小轮齿的{count}个已接受相位瞬时压力，非整周期累计着色印迹。灰色为实际完整裁切齿面；沿原实际曲面中心网格以原压力作连续C0-P1标量云图，蓝低红高。原求解常压力、参考平面积分及力矩定义不改；显示场支撑面积和原承载面积分别保存。真实边界承载时不外推缺失域，显示积分差单列。右图为整个斑点的等比例放大，不扩宽坐标、不作RGB滤波或斑点拟合；显示重构不等于压力物理验证。',
            f'Instantaneous pressures on the same physical pinion tooth at {count} accepted phase(s), not a cycle-accumulated marking pattern. Grey shows the actual complete trimmed flank. A continuous C0-P1 scalar field interpolates original pressures on actual surface-center grids, blue low to red high. Original constant-cell solver values, reference-plane integration and moment definitions remain unchanged. Display support and original loaded areas are recorded separately. Missing domains are not extrapolated at true loaded edges and display-integral differences remain explicit. Local views show the entire patch at equal aspect without coordinate widening, RGB filtering or fitted patches. Reconstruction is not physical pressure validation.')
        caption+=_text(lang,' 零压力采用色标p=0深蓝；有效求解域和距离界已排除的实际分片均显示蓝底。部分覆盖证书的未决分片单独灰纹，未决面积与状态保留；完整证书才认证全工作面。未计算和无效压力不改写为0。',
            ' Zero pressure is the exact p=0 blue of the colour scale. Solved-domain and distance-excluded actual pieces have a blue base. Partial coverage retains individually grey-hatched unresolved pieces, their area and incomplete status; only a complete certificate certifies the full flank. Uncomputed or invalid pressures are never rewritten as zero.')
        caption=caption.replace('灰色为完整实际裁切工作面；','背景为完整实际裁切工作面；').replace('灰色为实际完整裁切齿面；','背景为实际完整裁切齿面；').replace('Grey is its complete actual trimmed flank.','The background is its complete actual trimmed flank.').replace('Grey shows the actual complete trimmed flank.','The background shows the actual complete trimmed flank.')
        return _save(fig, output, '10_halfspace_contact_pattern', lang, data['rows']+data['load_share_rows'],
            caption, height_mm=30+40*count, extra=dict(model=data['model'], contact_representation='actual_surface_rectangular_cells',
                interpolation='C0_P1_scalar' if display=='continuous' else 'none',rendering=display,
                original_pressure_values_modified=False,display_field_diagnostics=[f['diagnostics'] for f in fields],
                pressure_domain_background=[pressure_domain_metadata(domain) for domain in domains],
                **scale,zoom_views=zoom,geometry_aspect='equal',
                over_limit_cell_counts=[sum(p*1e-6>scale['display_pressure_limit_mpa'] for p in phase['pressure_pa']) for phase in data['phases']],
                tooth_offset=tooth_offset, phase_indices=[p['phase_index'] for p in data['phases']],
                driving_angles_rad=[p['driver_rad'] for p in data['phases']],
                active_cell_counts=[p['active_cell_count'] for p in data['phases']],
                valid_cell_counts=[p['valid_cell_count'] for p in data['phases']],
                grid_shapes_yx=[p['grid_shape'] for p in data['phases']],
                boundary_checks=[p['boundary_checks'] for p in data['phases']],
                source_cell_corner_indices='(j,i), (j,i+1), (j+1,i+1), (j+1,i); actual_surface_corner行共享',
                projection_definition=data['projection_definition'], pressure_definition=data['pressure_definition']))


def hertz_validation_data(*, grid_counts=(32, 64, 128), radius_m=.02, load_n=1000.,
                           young_modulus_pa=210e9, poisson_ratio=.3):
    """实际独立求解三级轴对称Hertz基准；这是解析验证题而非齿轮或厂商对标。"""
    from .halfspace_contact import RectangularContactGrid, prepare_halfspace_operator, solve_halfspace_contact
    np = _np()
    counts = tuple(grid_counts)
    if len(counts) != 3 or len(set(counts)) != 3 or sorted(counts) != list(counts) or any(
            isinstance(n, bool) or not isinstance(n, int) or not 8 <= n <= 512 for n in counts):
        raise ValueError('Hertz图需要三个严格递增的8至512之间网格数')
    values = (radius_m, load_n, young_modulus_pa, poisson_ratio)
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in values):
        raise ValueError('Hertz材料、载荷与半径须为有限实数')
    if min(radius_m, load_n, young_modulus_pa) <= 0 or not -1 < poisson_ratio < .5:
        raise ValueError('Hertz半径／力／弹性模量须正，泊松比须在(-1,.5)')
    effective = young_modulus_pa/(2*(1-poisson_ratio**2))
    radius = (3*load_n*radius_m/(4*effective))**(1/3)
    depth = radius**2/radius_m
    peak = 3*load_n/(2*math.pi*radius**2)
    stiffness = 2*effective*radius
    reference = dict(curvature_radius_m=radius_m, load_n=load_n, young_modulus_pa=young_modulus_pa,
        poisson_ratio=poisson_ratio, effective_modulus_pa=effective, contact_radius_m=radius,
        indentation_m=depth, peak_pressure_pa=peak, contact_area_m2=math.pi*radius**2,
        incremental_stiffness_n_per_m=stiffness)
    levels, rows = [], []
    for count in counts:
        side = 4*radius; spacing = side/count
        grid = RectangularContactGrid(count, count, spacing, spacing, (-side/2+spacing/2,)*2)
        x, y = grid.centres_m()
        operator = prepare_halfspace_operator(grid, young_modulus_pa, poisson_ratio, young_modulus_pa, poisson_ratio)
        result = solve_halfspace_contact((x*x+y*y)/(2*radius_m), operator, load_n=load_n)
        if not result['converged']:
            raise RuntimeError('Hertz基准未收敛，不能绘制候选压力：'+str(result['certificate']))
        stiffness_value = result.get('incremental_stiffness_n_per_m')
        if stiffness_value is None:
            if result.get('stiffness_status') == 'requires_scipy':
                raise ImportError('Hertz切线刚度验证需要可选SciPy，请安装surface依赖')
            raise RuntimeError('Hertz切线刚度未通过，不能用未知值绘制误差：'+str(result.get('stiffness_status')))
        if not math.isfinite(stiffness_value) or stiffness_value <= 0:
            raise RuntimeError('Hertz切线刚度须为已接受的有限正值')
        field = result['fields'][0]
        exact = peak*np.sqrt(np.maximum(1-(x*x+y*y)/radius**2, 0.))
        area = result['active_cell_count']*spacing**2
        measured = dict(peak_pressure_pa=result['maximum_pressure_pa'], indentation_m=result['normal_approach_m'],
            contact_area_m2=area, incremental_stiffness_n_per_m=stiffness_value)
        errors = {key: value/reference[key]-1 for key, value in measured.items()}
        edges = -side/2+np.arange(count+1)*spacing
        levels.append(dict(count=count, spacing_m=spacing, x_centres_m=x, y_centres_m=y, edges_m=edges,
            pressure_pa=field['pressure_pa'], analytic_pressure_at_centres_pa=exact,
            measured=measured, relative_errors=errors, result=result))
        rows.append(dict(row_kind='grid_summary', grid_count=count, cell_spacing_m=spacing,
            resultant_force_n=result['resultant_load'], target_force_n=load_n,
            active_cell_count=result['active_cell_count'], **measured,
            **{key+'_relative_error': value for key, value in errors.items()},
            **{'kkt_'+key: value for key, value in result['certificate'].items()}))
        for j, i in np.ndindex(count, count):
            rows.append(dict(row_kind='hertz_pressure_cell', grid_count=count, grid_x=i, grid_y=j,
                x_m=float(x[j, i]), y_m=float(y[j, i]), cell_spacing_m=spacing,
                pressure_pa=float(field['pressure_pa'][j, i]), analytic_pressure_at_centre_pa=float(exact[j, i]),
                force_n=float(field['forces_n'][j, i]), active=bool(field['contact_mask'][j, i])))
    rows.append(dict(row_kind='analytic_reference', **reference))
    return dict(schema='halfspace_hertz_validation_v1', reference=reference, levels=levels, rows=rows)


def plot_halfspace_hertz_validation(data, output, *, lang='zh'):
    """Hertz真实压力格、归一化径向采样及三级误差图，统一183mm／600dpi。"""
    from .workflow_plotting import _paper, _panel, _save, _text, _COLORS
    from .publication_style import contact_pressure_colormap
    np = _np()
    if not isinstance(data, dict) or data.get('schema') != 'halfspace_hertz_validation_v1' or len(data['levels']) != 3:
        raise ValueError('需要实际hertz_validation_data三级求解结果')
    if any(not level['result']['converged'] for level in data['levels']):
        raise ValueError('失败Hertz解不得用于验证图')
    reference, finest = data['reference'], data['levels'][-1]
    a, p0 = reference['contact_radius_m'], reference['peak_pressure_pa']
    with _paper(lang) as plt:
        fig, axes = plt.subplots(2, 2, layout='constrained')
        ax = axes[0, 0]
        field = ax.pcolormesh(finest['edges_m']/a, finest['edges_m']/a, finest['pressure_pa']/p0,
            cmap=contact_pressure_colormap(), shading='flat', vmin=0., vmax=1.)
        from matplotlib.patches import Circle
        ax.add_patch(Circle((0, 0), 1., fill=False, color='white', linewidth=.75, linestyle=(0, (3, 2))))
        ax.set(xlim=(-1.35, 1.35), ylim=(-1.35, 1.35), aspect='equal', xlabel=r'$x/a$', ylabel=r'$y/a$')
        ax.set_title(_text(lang, f"{finest['count']}²格实际压力", f"Solved pressure, {finest['count']}² cells"), loc='left', x=.08)
        colorbar = fig.colorbar(field, ax=ax, pad=.025, fraction=.047)
        colorbar.set_label(_text(lang, r'归一化压力 $p/p_0$', r'Normalized pressure, $p/p_0$'))
        _panel(ax, 'a')
        ax = axes[0, 1]
        rho = np.linspace(0., 1.15, 400)
        ax.plot(rho, np.sqrt(np.maximum(1-rho*rho, 0)), color='#333D44', lw=1.2,
            label=_text(lang, 'Hertz解析解', 'Hertz analytical'))
        for level, color in zip(data['levels'], _COLORS):
            j = level['count']//2; x, y = level['x_centres_m'][j], level['y_centres_m'][j]
            chosen = x >= 0
            radial = np.sqrt(x[chosen]**2+y[chosen]**2)/a
            ax.plot(radial, level['pressure_pa'][j, chosen]/p0, marker='o', markersize=2.2,
                lw=.6, color=color, label=f"{level['count']} × {level['count']}")
        ax.set(xlim=(0, 1.18), ylim=(-.025, 1.04), xlabel=_text(lang, r'采样半径 $r/a$', r'Sample radius, $r/a$'),
            ylabel=_text(lang, r'归一化压力 $p/p_0$', r'Normalized pressure, $p/p_0$'))
        ax.legend(frameon=False, fontsize=7, loc='lower left')
        _panel(ax, 'b')
        ax = axes[1, 0]
        metrics = (('peak_pressure_pa', r'$p_\mathrm{max}$', _COLORS[0], 'o'),
                   ('indentation_m', r'$\delta$', _COLORS[1], 's'),
                   ('contact_area_m2', r'$A$', _COLORS[2], '^'),
                   ('incremental_stiffness_n_per_m', r'$k_\mathrm{n}$', _COLORS[3], 'D'))
        for key, label, color, marker in metrics:
            ax.plot([level['count'] for level in data['levels']],
                [abs(level['relative_errors'][key])*100 for level in data['levels']],
                label=label, marker=marker, markersize=3.5, color=color)
        ax.set_xscale('log', base=2); ax.set_yscale('log')
        ax.set_xticks([level['count'] for level in data['levels']], [str(level['count']) for level in data['levels']])
        ax.set_xlabel(_text(lang, '每方向网格数 N', 'Cells per direction, N'))
        ax.set_ylabel(_text(lang, '相对解析解的绝对误差 (%)', 'Absolute relative error (%)'))
        ax.legend(frameon=False, ncol=2, fontsize=7)
        _panel(ax, 'c')
        ax = axes[1, 1]; ax.axis('off')
        values = [
            _text(lang, 'Hertz解析基准；非厂商对标', 'Hertz benchmark; no vendor data'),
            rf'$F={reference["load_n"]:g}$ N;  $R={reference["curvature_radius_m"]*1000:g}$ mm',
            rf'$E_1=E_2={reference["young_modulus_pa"]*1e-9:g}$ GPa;  $\nu_1=\nu_2={reference["poisson_ratio"]:g}$',
            rf'$a={a*1000:.5f}$ mm;  $p_0={p0*1e-6:.2f}$ MPa',
            rf'$\delta={reference["indentation_m"]*1e6:.5f}$ μm',
            rf'$k_\mathrm{{n}}={reference["incremental_stiffness_n_per_m"]*1e-6:.3f}$ MN/m',
            _text(lang, '最细网格误差：', 'Finest-grid errors:'),
            rf'$p_\mathrm{{max}}$: {finest["relative_errors"]["peak_pressure_pa"]*100:+.5f}%;  '+
            rf'$\delta$: {finest["relative_errors"]["indentation_m"]*100:+.5f}%',
            _text(lang, '合力平衡相对残差：', 'Relative force-balance residual:'),
            f"{abs(finest['result']['balance_residual'])/reference['load_n']:.2e}",
        ]
        for index, label in enumerate(values):
            ax.text(.02, .98-index*.095, label, transform=ax.transAxes, va='top', fontsize=7.5)
        _panel(ax, 'd')
        caption = _text(lang,
            '轴对称抛物面Hertz解析验证：两体弹性半空间，常压力矩形格，32²、64²和128²各自重新求解。a为解析接触半径，p0为解析峰压，δ为法向压陷，A为承载格面积，kn=dF/dδ为切线法向刚度。a图白虚圆是解析边界；b图取实际单元中心的径向坐标，不把相邻中心插值为测量数据。c图是确定性离散误差，没有统计误差棒。压力图以flat方式显示已解单元值，色标不平滑。此图验证接触内核，不是实际锥齿面精度证明或外部LTCA数值对标。',
            'Axisymmetric paraboloid Hertz benchmark for two elastic half-spaces. Constant-pressure rectangular grids of 32², 64² and 128² cells are solved independently. a is analytical contact radius, p0 peak pressure, δ normal indentation, A loaded-cell area and kn=dF/dδ tangent normal stiffness. The white dashed circle in a is the analytical boundary. b uses radii of actual cell centres, without interpolating samples. c shows deterministic discretization errors, without statistical error bars. Flat cell colors show solved pressures without smoothing. This validates the contact kernel, not actual bevel-flank accuracy or agreement with external LTCA results.')
        return _save(fig, output, '11_halfspace_hertz_validation', lang, data['rows'], caption, height_mm=150,
            extra=dict(benchmark='Hertz轴对称解析算例', vendor_comparison=False,
                grid_counts=[level['count'] for level in data['levels']], analytic_reference=reference,
                relative_errors_by_grid={str(level['count']):level['relative_errors'] for level in data['levels']},
                interpolation='none; flat constant-pressure cells'))
