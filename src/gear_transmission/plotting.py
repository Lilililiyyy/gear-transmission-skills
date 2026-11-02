"""新写通用绘图适配：复用新版案例后处理，任意工况不混入参照响应。"""
from pathlib import Path

import numpy as np
from gear_contact.io import write_json
from gear_contact.plotting import _pyplot, _save, _csv
from gear_contact.contact_display import (contact_pass_data, contact_pressure_colormap,
    draw_contact_pass, plot_contact_pass_comparison)
from gear_ltca_snapshot.examples.case_contact_plotting import (
    repeat_cycle_for_display, _project_line_samples, _endpoint_display, _band_values, _draw_band)
from gear_ltca_snapshot.examples.case_geometry_plotting import plot_case_geometry


def plot_cycles(results, output):
    """原样重复三个稳态周期，端点与中心峰压分别绘制，保留物理起采相位。"""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.4))
    fig.subplots_adjust(left=.10, right=.97, bottom=.19, top=.87, hspace=.46, wspace=.30)
    raw, shown = [], []
    for label, result in results.items():
        centre = np.array([np.max(r.pressure_pa, initial=0.) for r in result.rotations])
        endpoint = np.array([np.max(r.endpoint_pressure['pressure_pa'], initial=0.)
                             for r in result.rotations]) if all(r.endpoint_pressure is not None for r in result.rotations) else None
        force = np.array([r.solution.forces_n.sum() for r in result.rotations])
        cycle = np.column_stack((np.arange(len(result.phase))/len(result.phase), result.phase,
            result.normal_approach_m, result.mesh_stiffness_n_per_m, centre,
            centre if endpoint is None else endpoint, force))
        display = repeat_cycle_for_display(cycle, 3)
        raw.extend((label, *row) for row in cycle)
        shown.extend((label, *row, i % len(cycle), i == len(display)-1) for i, row in enumerate(display))
        for ax, column, scale, title, ylabel in (
            (axes[0, 0], 2, 1e6, 'a  承载法向 TE', '接近量 / μm'),
            (axes[0, 1], 3, 1e-6, 'b  法向刚度', '刚度 / (N/μm)'),
            (axes[1, 0], 4, 1e-6, 'c  每相位峰值压力', '压力 / MPa'),
            (axes[1, 1], 6, 1., 'd  总法向载荷', '载荷 / N')):
            ax.plot(display[:, 0], display[:, column]*scale, marker='.', markersize=3, label=label)
            ax.set(title=title, xlabel='相对起采位置 / 基节周期', ylabel=ylabel, xlim=(0, 3), xticks=(0, 1, 2, 3))
            ax.grid(alpha=.6)
        if endpoint is not None:
            axes[1, 0].plot(display[:, 0], display[:, 5]*1e-6, '--', label=label+'：端点')
    axes[0, 0].legend(loc='upper left', fontsize=8)
    axes[1, 0].legend(fontsize=8)
    fig.text(.5, .045, '同一个稳态周期原样重复三次；末端增加起点作为闭合点，不是连续多周期时域求解。\n'
             '中心峰压与独立端点峰压分开保留，未计算端点时不将中心压力标为端点压力。', ha='center', fontsize=8.3)
    _save(fig, output, '03_ltca_cycle')
    plt.close(fig)
    columns = ('design', 'cycle_fraction', 'physical_phase', 'normal_approach_m', 'mesh_stiffness_n_per_m',
               'centre_peak_pressure_pa', 'display_peak_pressure_pa', 'normal_force_n')
    _csv(output/'03_ltca_cycle.csv', columns, raw)
    _csv(output/'03_ltca_cycle_display.csv', (columns[0], 'display_cycle_position', *columns[2:],
         'source_sample_index', 'periodic_closure'), shown)


def plot_result(result, output):
    """输出输入几何、三周期、瞬时压力带、材料坐标及整齿通过程图和 SI 数据。"""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    plot_case_geometry(result.pair, output, modification_pinion=result.inputs['modification_pinion'],
                       modification_wheel=result.inputs['modification_wheel'])
    plot_cycles({'当前方案': result}, output)
    phase_index = 0
    rotation = result.rotations[phase_index]
    pairs = np.unique(rotation.grid.tooth_pair_index)
    pair_index = int(max(pairs, key=lambda k: rotation.solution.forces_n[rotation.grid.tooth_pair_index == k].sum()))
    line = _project_line_samples(rotation, pair_index, np.cos(result.pair.base_helix_angle_rad))
    endpoints = line['endpoints']
    peak = line['peak_pa'] if endpoints is None else endpoints['pressure_pa']
    half_width = line['half_width_m'] if endpoints is None else endpoints['half_width_m']
    extent = max(float(np.max(half_width, initial=0.))*1.25, 1e-6)
    normal = np.linspace(-extent, extent, 122)
    intervals = line['intervals_m']
    if endpoints is not None:
        edges, peak, half_width = _endpoint_display(line['edges_m'], peak, half_width)
        intervals = np.column_stack((edges[:-1], edges[1:]))
    band = _band_values(peak, half_width, normal)
    from matplotlib.colors import Normalize
    fig, ax = plt.subplots(figsize=(10.7, 4.5))
    fig.subplots_adjust(bottom=.25, top=.85, left=.10, right=.88)
    image = _draw_band(ax, intervals, normal, band, cmap=plt.get_cmap('coolwarm'),
                       norm=Normalize(0., max(float(np.max(peak, initial=0.))*1e-6, 1.)))
    ax.set_facecolor('#E7EBEF')
    ax.set(xlim=(intervals[:, 0].min()*1e3, intervals[:, 1].max()*1e3), ylim=(-extent*1e6, extent*1e6),
           xlabel='沿接触线 s / mm', ylabel='垂直接触线 n / μm',
           title=f'同相位、单齿对瞬时压力：相位 {result.phase[phase_index]:.6g}，齿对 {pair_index:+d}')
    fig.colorbar(image, ax=ax, label='接触压力 / MPa')
    fig.text(.5, .055, ('中心条带压力带' if endpoints is None else '独立端点压力及半宽插值的 Hertz 带')+
             '；两轴显示比例不同。\n显示后处理不修改原压力、条带力或评价指标。', ha='center', fontsize=8.5)
    _save(fig, output, '04_contact_pattern')
    plt.close(fig)
    data = contact_pass_data(result)
    points = data['samples']
    np.savez_compressed(output/'contact_plot_data.npz', **{k: v for k, v in data.items() if isinstance(v, np.ndarray)},
        instant_intervals_m=intervals, instant_normal_edges_m=normal, instant_band_pressure_pa=band,
        instant_raw_cell_peak_pa=line['peak_pa'], instant_raw_cell_force_n=line['force_n'],
        instant_raw_endpoint_peak_pa=np.array([]) if endpoints is None else endpoints['pressure_pa'])
    _csv(output/'07_project_contact_pass_samples.csv', ('face_coordinate_from_left_m', 'roll_pinion_m',
        'pressure_pa', 'extended_pinion', 'extended_wheel'), points)
    cmap = contact_pressure_colormap()
    fig, ax = plt.subplots(figsize=(8., 6.2))
    maximum = max(float(points[:, 2].max())*1e-6, 1.)
    image = draw_contact_pass(ax, data, maximum_pa=maximum*1e6, title='小轮整齿通过程接触压力')
    fig.colorbar(image, ax=ax, label='接触压力 / MPa', fraction=.045, pad=.03)
    fig.subplots_adjust(left=.11, right=.84, bottom=.23, top=.91)
    fig.text(.5, .045, data['definition']+'；重复材料点取周期最大值。\n'
             '显示域为两轮实际啮合区；灰色为剩余未覆盖，延伸分支另存。\n'+data['support_definition'],
             ha='center', fontsize=8.2)
    _save(fig, output, '07_project_contact_pass')
    plt.close(fig)
    # 材料坐标图保留所有相位和延伸分支，不把通过程包络当作同时承载面积。
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.7), sharex=True, sharey=True)
    first_count = len(rotation.endpoint_pressure['pressure_pa']) if endpoints is not None else len(rotation.pressure_pa)
    for ax, own, title in ((axes[0], points[:first_count], 'a  第一采样相位'), (axes[1], points, 'b  一个周期的采样位置')):
        extended = np.any(own[:, 3:5].astype(bool), axis=1)
        for mask, marker in ((~extended, 'o'), (extended, '^')):
            if mask.any():
                image = ax.scatter(own[mask, 0]*1e3, own[mask, 1]*1e3, c=own[mask, 2]*1e-6,
                                   s=6, marker=marker, cmap=cmap, vmin=0., vmax=maximum)
        for bound in result.pair.pinion.involute_roll_bounds_m:
            ax.axhline(bound*1e3, ls='--', color='#667A88', lw=.65)
        ax.set(title=title, xlabel='小轮左端起齿宽坐标 / mm', ylabel='小轮展长 / mm')
    fig.colorbar(image, ax=list(axes), label='压力 / MPa', shrink=.7)
    _save(fig, output, '05_contact_coordinates')
    plt.close(fig)
    metrics = {'pressure_definition': data['definition'], 'instant_phase_index': phase_index,
               'instant_tooth_pair_index': pair_index, 'raw_sample_count': len(points),
               'raw_sample_peak_pa': float(points[:, 2].max()), 'interpolated': bool(data['interpolated']),
               'uncovered_display_count': int(np.count_nonzero(~np.isfinite(data['pressure_pa_x_y']))),
               'display_grid_shape_x_y': list(data['pressure_pa_x_y'].shape),
               'contact_roll_bounds_m': data['contact_roll_bounds_m'],
               'pinion_involute_roll_bounds_m': data['pinion_involute_roll_bounds_m'],
               'boundary_display_sample_count': len(data['boundary_display_samples']),
               'boundary_display_definition': data['support_definition'],
               'style_source': 'gear_contact.contact_display：与源工程案例共用',
               'reference_used': False, 'complete_2d_contact_accuracy_verified': False,
               'display_definition': '两轮实际啮合区内的材料坐标显示；原始采样和有限单元支撑以外不填零'}
    write_json(output/'plot_manifest.json', metrics)
    return metrics
