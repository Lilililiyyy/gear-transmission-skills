"""全啮合分相最大压力的双语论文图；颜色表示压力，不表示出现次数。"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path


def _path(geometry, scale=1000.):
    """真实GEOS区域转复合路径，保留所有孔洞及分离区域。"""
    import numpy as np
    from matplotlib.path import Path as PlotPath
    from shapely.geometry.polygon import orient
    paths = []
    if geometry.is_empty:
        return None
    parts = [geometry] if geometry.geom_type == 'Polygon' else list(geometry.geoms)
    for part in parts:
        if part.geom_type != 'Polygon':
            continue
        polygon = orient(part, sign=1.)
        for ring in [polygon.exterior, *polygon.interiors]:
            paths.append(PlotPath(np.asarray(ring.coords)*scale, closed=True))
    return PlotPath.make_compound_path(*paths) if paths else None


def export_mesh_cycle_pressure_data(data, output, *, stem=None):
    """导出逐点最大值、达到最大值的相位及未知标志；压缩不删除数值行。"""
    import numpy as np
    from shapely.geometry import mapping
    from .refined_fe_workflow import _write
    from .contact_cycle import _geometry_union
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    stem = stem or data['method']+'_cycle_pressure'
    path = output/(stem+'.csv.gz')
    x, y = data['x_m'], data['y_m']
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['width_position_m', 'normal_height_m', 'sampled_maximum_pressure_pa',
            'maximum_known', 'argmax_phase_index', 'argmax_tooth_offset'])
        for j, i in zip(*np.nonzero(data['inside_flank_mask'])):
            writer.writerow([format(float(x[i]), '.17g'), format(float(y[j]), '.17g'),
                format(float(data['pressure_pa'][j, i]), '.17g'), int(data['known_mask'][j, i]),
                int(data['maximum_phase_index'][j, i]), int(data['maximum_tooth_offset'][j, i])])
    arrays = output/(stem+'.npz')
    np.savez_compressed(arrays, **{key: data[key] for key in ('x_m', 'y_m', 'pressure_pa',
        'inside_flank_mask', 'known_mask', 'maximum_phase_index', 'maximum_tooth_offset')})
    metadata = {key: value for key, value in data.items() if key not in
        ('x_m', 'y_m', 'pressure_pa', 'inside_flank_mask', 'known_mask', 'maximum_phase_index',
         'maximum_tooth_offset', 'unknown_geometry', 'background_polygons_m')}
    metadata.update(outline_geometry=mapping(_geometry_union(data['background_polygons_m'])),
        unknown_geometry=mapping(data['unknown_geometry']),
        data_files={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (path, arrays)},
        query_grid='仅显示网格；原压力峰值另存，不将采样峰值冒充原峰值',
        unknown_pressure_interpretation='maximum_known=0的数值仅是已求解分相场的下界；未认证域不填作确定零压')
    _write(metadata, output/(stem+'.json'))
    return metadata


def load_mesh_cycle_pressure_data(metadata_path):
    """只凭公开小汇总及压缩SI数组重绘，不需要CAD或原FE运行环境。"""
    import numpy as np
    from shapely.geometry import shape
    from .refined_fe_postprocessing import _relative_file, _read, _sha
    path = Path(metadata_path).resolve(); meta = _read(path)
    if meta.get('schema') != 'bevel_mesh_cycle_pressure_data_v1':
        raise ValueError('需要全周期压力数据schema')
    for name, digest in meta['data_files'].items():
        if _sha(_relative_file(path.parent, name)) != digest:
            raise ValueError('压力包络SI数据SHA256已变化')
    filenames = [name for name in meta['data_files'] if name.endswith('.npz')]
    if len(filenames) != 1:
        raise ValueError('压力包络须唯一绑定NPZ数组')
    with np.load(_relative_file(path.parent, filenames[0]), allow_pickle=False) as arrays:
        data = dict(meta, **{key: arrays[key].copy() for key in arrays.files})
    data['unknown_geometry'] = shape(meta['unknown_geometry'])
    data['outline_geometry'] = shape(meta['outline_geometry'])
    return data


def _contour_supported_mask(valid, x=None, y=None):
    """返回corner_mask=True的连续支撑；排除孤立点及浮点坐标尺度的退化微格。"""
    import numpy as np
    valid = np.asarray(valid, dtype=bool)
    if valid.ndim != 2 or min(valid.shape) < 2:
        raise ValueError('连续云图至少需要二维2×2查询网格')
    count = (valid[:-1, :-1].astype(np.uint8) + valid[1:, :-1] +
        valid[:-1, 1:] + valid[1:, 1:])
    cells = count >= 3
    if x is not None or y is not None:
        x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
        if x.shape != (valid.shape[1],) or y.shape != (valid.shape[0],):
            raise ValueError('连续云图坐标轴须与二维掩码同形')
        # 原峰坐标与规则查询点可能只差舍入误差；不将这种微格视作可分辨面积。
        tolerance_x = 64*np.finfo(float).eps*max(float(np.abs(x).max()), np.finfo(float).tiny)
        tolerance_y = 64*np.finfo(float).eps*max(float(np.abs(y).max()), np.finfo(float).tiny)
        cells &= (np.diff(y) > tolerance_y)[:, None] & (np.diff(x) > tolerance_x)[None, :]
    support = np.zeros_like(valid)
    for y, x in ((slice(None, -1), slice(None, -1)), (slice(1, None), slice(None, -1)),
            (slice(None, -1), slice(1, None)), (slice(1, None), slice(1, None))):
        support[y, x] |= cells & valid[y, x]
    return support


def plot_mesh_cycle_pressure(data_sets, output, *, lang='zh', stem='16_mesh_cycle_maximum_pressure', shared_pressure_scale=False,
        zoom_window_m=None):
    """各方法给出最大压力及原峰曲线；默认各行独立线性色标并明确MPa数值。"""
    import numpy as np
    from .workflow_plotting import _paper, _save, _text
    from .halfspace_plotting import pressure_colour_scale
    from .contact_cycle import _geometry_union
    from .publication_style import COLORS
    if not data_sets or any(d.get('schema') != 'bevel_mesh_cycle_pressure_data_v1' for d in data_sets):
        raise ValueError('需要真实全周期压力数据，不能画占位图')
    torques = [d['signed_pinion_torque_nm'] for d in data_sets]
    if any(abs(t-torques[0]) > 1e-10 for t in torques):
        raise ValueError('同一最大压力比较图必须采用相同有符号扭矩')
    maximum = max(d['source_peak_pressure_pa'] for d in data_sets)
    if maximum <= 0:
        raise ValueError('全周期没有受载压力，不能生成压力包络图')
    visible_peaks = []
    for data in data_sets:
        visible = _contour_supported_mask(np.asarray(data['inside_flank_mask']) & np.asarray(data['known_mask']), data['x_m'], data['y_m'])
        values = np.asarray(data['pressure_pa'])[visible]
        visible_peaks.append(float(values.max(initial=0.)))
    # 灰域和没有连续着色支撑的孤立点不决定色标；原数值仍留在数组与峰压曲线。
    visible_maximum = max(visible_peaks)
    norm, cmap, scale = pressure_colour_scale([visible_maximum if visible_maximum > 0 else maximum])
    if not isinstance(shared_pressure_scale, bool):
        raise ValueError('共用绝对压力色标开关须为布尔值')
    zoom_windows = ({'full_gear_fe': (.00008, .000025), 'local_halfspace': (.0016, .0005)}
        if zoom_window_m is None else zoom_window_m)
    def window_for(method):
        window = zoom_windows[method] if isinstance(zoom_windows, dict) else zoom_windows
        if len(window) != 2 or any(not np.isfinite(v) or v <= 0 for v in window):
            raise ValueError('高压放大窗须显式给出正的齿宽／齿高尺寸，单位m')
        return np.asarray(window, dtype=float)
    for data in data_sets:
        window_for(data['method'])
    rows = []; method_scales = []; zoom_records = []
    with _paper(lang) as plt:
        from matplotlib.cm import ScalarMappable
        from matplotlib.patches import PathPatch, Patch, Rectangle
        fig, axes = plt.subplots(len(data_sets), 3, squeeze=False,
            gridspec_kw=dict(width_ratios=[1.9, 1.1, 1.2]))
        fig.subplots_adjust(left=.07, right=.87, bottom=.20, top=.83,
            wspace=.65, hspace=.85)
        for row_index, data in enumerate(data_sets):
            window = window_for(data['method'])
            if not shared_pressure_scale:
                bound = visible_peaks[row_index] if visible_peaks[row_index] > 0 else data['source_peak_pressure_pa']
                norm, cmap, row_scale = pressure_colour_scale([bound])
            else:
                row_scale = dict(scale)
            method_scales.append(dict(method=data['method'], **row_scale))
            ax, zoom, line = axes[row_index]
            outline = data.get('outline_geometry')
            if outline is None:
                outline = _geometry_union(data['background_polygons_m'])
            clip = _path(outline)
            xx, yy = np.meshgrid(data['x_m']*1000., data['y_m']*1000.)
            visible = _contour_supported_mask(data['inside_flank_mask'] & data['known_mask'], data['x_m'], data['y_m'])
            values = np.ma.array(data['pressure_pa']*1e-6, mask=~visible)
            unknown = _path(data['unknown_geometry'])
            j, i = np.unravel_index(np.argmax(np.where(visible, data['pressure_pa'], -1.)), visible.shape)
            centre = np.array([data['x_m'][i], data['y_m'][j]])
            low = centre-window*.5; high = centre+window*.5
            for target in (ax, zoom):
                # 无连续云图覆盖的位置不能借用零压蓝色；已绘制的真实零值由contourf编码。
                target.add_patch(PathPatch(clip, facecolor='#F3F4F5', edgecolor='#535A60', lw=.55, zorder=0))
                if target is zoom:
                    margin = np.asarray(data['actual_display_spacing_m'])*2
                    ix = np.flatnonzero((data['x_m'] >= low[0]-margin[0]) & (data['x_m'] <= high[0]+margin[0]))
                    iy = np.flatnonzero((data['y_m'] >= low[1]-margin[1]) & (data['y_m'] <= high[1]+margin[1]))
                    xx_part, yy_part = np.meshgrid(data['x_m'][ix]*1000., data['y_m'][iy]*1000.)
                    values_part = values[np.ix_(iy, ix)]
                else:
                    xx_part, yy_part, values_part = xx, yy, values
                contour = target.contourf(xx_part, yy_part, values_part, levels=np.linspace(0., norm.vmax, 257),
                    norm=norm, cmap=cmap, antialiased=False, corner_mask=True, zorder=2)
                artists = [contour] if hasattr(contour, 'set_clip_path') else contour.collections
                for collection in artists:
                    collection.set_clip_path(clip, target.transData); collection.set_clip_box(target.bbox)
                    # 密集压力面按600dpi输出；文字、坐标与曲线保持矢量。
                    collection.set_rasterized(True)
                if unknown is not None:
                    artist = PathPatch(unknown, facecolor='#D5DADD', edgecolor='#8C979E', hatch='////', linewidth=0., zorder=4)
                    artist.set_clip_path(clip, target.transData); artist.set_clip_box(target.bbox); target.add_patch(artist)
                target.set_aspect('equal', adjustable='box')
            xmin, ymin, xmax, ymax = outline.bounds
            ax.set_xlim(xmin*1000-.3, xmax*1000+.3)
            ax.set_ylim(ymin*1000-.18, ymax*1000+.18)
            ax.add_patch(Rectangle(low*1000., *(window*1000.),
                fill=False, edgecolor='#303840', lw=.65, zorder=5))
            zoom.set_xlim(low[0]*1000., high[0]*1000.); zoom.set_ylim(low[1]*1000., high[1]*1000.)
            zoom.set_title(_text(lang, '高压局部\n（等比例）', 'Peak-region detail\n(equal aspect)'), loc='left', fontsize=7, pad=6)
            zoom.set_xlabel(r'$w_\mathrm{f}$ (mm)'); zoom.set_ylabel('h (mm)')
            from matplotlib.ticker import MaxNLocator, ScalarFormatter
            for axis in (zoom.xaxis, zoom.yaxis):
                axis.set_major_locator(MaxNLocator(nbins=3))
                axis.set_major_formatter(ScalarFormatter(useOffset=False))
            zoom_records.append(dict(method=data['method'], centre_width_height_m=centre.tolist(),
                width_height_m=window.tolist(), geometry_aspect='equal', source='已审核且有连续云图支撑的查询峰值位置'))
            name = (_text(lang, '整轮FE', 'Full-gear FE') if data['method'] == 'full_gear_fe'
                else _text(lang, '局部半解析', 'Local semi-analytical'))
            cycle = data['cycle']
            title = name+'\n'+_text(lang, f'{cycle["interval_count"]}等分相位',
                f'{cycle["interval_count"]} phase intervals')
            ax.set_title(title, loc='left', fontsize=8, pad=7)
            ax.set_xlabel(_text(lang, r'齿宽位置 $w_\mathrm{f}$ (mm)', r'Face-width position, $w_\mathrm{f}$ (mm)'))
            ax.set_ylabel(_text(lang, '节锥法向高度 h (mm)', 'Pitch-cone normal\nheight, h (mm)'))
            phase = data['phase_rows']
            fraction = [r['mesh_fraction'] for r in phase]
            peaks = [r['peak_pressure_pa']*1e-6 for r in phase]
            line.plot(fraction, peaks, color=COLORS[row_index % 2], marker='o', ms=2.5, lw=1.1)
            line.set_xlim(0., 1.); line.set_xticks([0., .5, 1.]); line.set_ylim(bottom=0.)
            line.set_xlabel(_text(lang, '啮合周期位置 ξ (−)', 'Mesh-cycle position, ξ (−)'))
            line.set_ylabel(_text(lang, '瞬时峰压 (MPa)', 'Instantaneous peak\npressure (MPa)'))
            line.set_title(_text(lang, '各相位原始峰压', 'Raw phase peaks'), loc='left', fontsize=7, pad=7)
            line.grid(axis='y', color='#D9DEE2', lw=.4, alpha=.7)
            line.text(.98, .96, rf'$p_\mathrm{{peak}}={data["source_peak_pressure_pa"]*1e-6:.0f}$ MPa',
                transform=line.transAxes, va='top', ha='right', fontsize=7,
                bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=1.5))
            for col, target in enumerate((ax, zoom, line)):
                target.annotate('abcdefghijklmnopqrstuvwxyz'[3*row_index+col], (0., 1.), xycoords='axes fraction',
                    xytext=(-8, 12), textcoords='offset points', ha='right', weight='bold', fontsize=9,
                    annotation_clip=False)
            for item in phase:
                rows.append(dict(row_kind='mesh_phase_peak', method=data['method'],
                    **{key: value for key, value in item.items() if key != 'pairs'}))
                for pair in item['pairs']:
                    values = {k: v for k, v in pair.items() if k != 'display_field_diagnostics'}
                    values['display_field_diagnostics_json'] = json.dumps(pair.get('display_field_diagnostics', {}),
                        ensure_ascii=False, separators=(',', ':'))
                    rows.append(dict(row_kind='tooth_pair_phase', method=data['method'], **values))
            if not shared_pressure_scale:
                position = line.get_position()
                bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                    cax=fig.add_axes([.91, position.y0+.1*position.height, .018, .8*position.height]))
                bar.set_label(_text(lang, '已审核域压力 (MPa)', 'Audited-domain pressure (MPa)'))
        if shared_pressure_scale:
            bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes[:, 0].tolist(),
                cax=fig.add_axes([.91, .27, .018, .48]))
            bar.set_label(_text(lang, '已审核域最大压力 (MPa)', 'Audited-domain maximum pressure (MPa)'))
        handles = [Patch(facecolor=cmap(norm(0.)), label=_text(lang, '最大压力为零', 'Zero maximum pressure'))]
        if any(not d['unknown_geometry'].is_empty for d in data_sets):
            handles.append(Patch(facecolor='#D5DADD', edgecolor='#8C979E', hatch='////',
                label=_text(lang, '任一相位未认证', 'Uncertified in any phase')))
        handles.append(Patch(facecolor='#F3F4F5', edgecolor='#A0A6AB',
            label=_text(lang, '显示无连续支撑', 'No continuous display support')))
        fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.46, .025), ncol=len(handles), fontsize=6)
        full_window = all(d['cycle']['single_tooth_path_coverage_confirmed'] for d in data_sets)
        title = (_text(lang, '全啮合过程 · 同一参考齿面的最大压力包络',
            'Full mesh cycle · maximum-pressure envelope on one reference flank') if full_window else
            _text(lang, '一齿距分相 · 所选齿对窗口内的最大压力',
                'One-pitch sampling · maximum pressure within the selected pair window'))
        fig.suptitle(title+'\n'+rf'$T_1={torques[0]:g}$ N·m', fontsize=9, y=.985)
        fig.text(.47, .11, r'$\xi=(\theta_1-\theta_\mathrm{start})/(2\pi/z_1)$',
            ha='center', fontsize=6)
        caption = _text(lang,
            '每个实际分相、每个周期等价齿对的原P1标量场映回同一材料齿面后逐点取最大值。'
            '左列是离散相位最大压力，中列按原坐标等比例放大真实高压区，右列是各相位瞬时原峰压；不是压力相加或同时受载面积。'
            '全齿面以真实坐标等比例显示，零压蓝、高压红，未知域为各状态未知域的并集并以灰色网纹覆盖。'
            '图中压力峰值取原求解节点/单元，密集显示网格可能略低，不作重标定。'
            '重复端点用于周期检查；有限相位与有限网格未证明连续过程峰值收敛，也不构成真实压力验证。',
            'Original P1 scalar fields from every sampled phase and periodically equivalent tooth pair are mapped '
            'to the same material flank before taking their pointwise maximum. Left: the sampled pressure envelope; '
            'middle: an equal-aspect zoom of the actual high-pressure region; '
            'right: original instantaneous peaks. This is neither a pressure sum nor simultaneous loaded area. '
            'True coordinate aspect is retained. Zero is blue and high pressure red; the union of all uncertified '
            'domains is grey and hatched. Labelled peaks retain original nodal/cell values; the display grid may '
            'slightly under-sample them and is not rescaled. Repeated endpoints check periodicity. '
            'Neither temporal nor pressure-mesh convergence, nor physical validation, is claimed.')
        if not full_window:
            caption += _text(lang, ' 外邻齿窗口尚未认证；数值仅覆盖所选齿对，不将窗口外未知解释为零压。',
                ' Outer-pair coverage remains uncertified; the result applies only to the selected pairs and '
                'does not assign zero to unknown pairs outside the window.')
        caption += (_text(lang, ' 各行采用独立线性色标，按标出的MPa读数比较压力。',
            ' Rows use separate linear colorbars; compare their labelled MPa values.') if not shared_pressure_scale else
            _text(lang, ' 所有行共用同一绝对线性压力色标。', ' All rows share one absolute linear pressure scale.'))
        caption += _text(lang, ' 色标按具有连续着色支撑的已审核域自动定量；灰域、孤立查询点和浮点尺度退化微格不参与量程，原始峰值仍保留在右侧曲线和SI数组。显示网格没有连续支撑的位置为浅灰底，不冒充零压蓝色。',
            ' Color limits use the audited domain with continuous contour support; grey regions, isolated query points and floating-point-degenerate cells do not set the range, '
            'and their raw peaks remain in the right-hand curves. Locations without continuous display support have a neutral light-grey background rather than zero-pressure blue.')
        summaries = [{key: d[key] for key in ('method', 'cycle', 'source_peak_pressure_pa',
            'display_grid_peak_pressure_pa', 'display_peak_relative_loss', 'actual_display_spacing_m',
            'unresolved_projected_area_fraction', 'complete_zero_domain_certified',
            'force_integral_defined', 'pressure_mesh_converged', 'physical_validation')} for d in data_sets]
        for summary, data in zip(summaries, data_sets):
            summary.update({key: data[key] for key in ('phase_sampling_study', 'endpoint_periodicity', 'solver_discretization')
                if key in data})
        for summary, peak, data in zip(summaries, visible_peaks, data_sets):
            summary['contour_supported_peak_pressure_pa'] = peak
            known = data['inside_flank_mask'] & data['known_mask']
            support = _contour_supported_mask(known, data['x_m'], data['y_m'])
            summary['known_display_peak_pressure_pa'] = float(data['pressure_pa'][known].max(initial=0.))
            summary['known_query_points_without_contour_support'] = int(np.count_nonzero(known & ~support))
        figure_scale = dict(scale) if shared_pressure_scale else {}
        directory = Path(__file__).resolve().parent
        render_source = {'bevel_gear/'+name+'.py': hashlib.sha256((directory/(name+'.py')).read_bytes()).hexdigest()
            for name in ('contact_cycle_plotting', 'workflow_plotting', 'halfspace_plotting', 'publication_style')}
        return _save(fig, output, stem, lang, rows, caption, width_mm=190,
            height_mm=78*len(data_sets)+10, extra=dict(**figure_scale, models=summaries,
                shared_absolute_pressure_scale=shared_pressure_scale, method_colour_scales=method_scales,
                colour_scale_domain='inside_flank_and_known_mask_with_contour_support',
                uncovered_display_background='neutral_light_grey_not_zero_blue',
                high_pressure_zoom=zoom_records,
                rendering_source_identity=render_source,
                original_pressure_values_modified=False, geometry_aspect='equal',
                dense_scalar_rasterized=True, scalar_export_dpi=600,
                cycle_definition='pointwise maximum, not simultaneous pressure or pressure sum'))
