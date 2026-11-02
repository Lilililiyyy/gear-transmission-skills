"""新写接触图后处理：统一材料坐标、有限单元显示支撑与源工程色带。

本模块不求解载荷、修改压力或补造独立端点。中心模型只将首尾单元的
常值显示支撑扩到该单元真实边界；实际啮合域、原始采样与显示网格分开保存。
"""
from pathlib import Path

import numpy as np


def contact_pressure_colormap():
    """源工程整齿压力图统一使用蓝、青、绿、黄、红色带。"""
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('pressure_blue_to_red',
        ['#0000ff', '#00ffff', '#00ff00', '#ffff00', '#ff0000'], N=256)
    cmap.set_bad('#e5e9ee')
    return cmap


def project_contact_pass_samples(result):
    """保留独立端点压力与延伸标记；坐标为小轮左端起齿宽及小轮展长，单位m。"""
    samples = []
    for rotation in result.rotations:
        geometry, pressure = rotation.contact_geometry_data, rotation.endpoint_pressure
        if geometry is None or pressure is None:
            raise ValueError('端点通过程需要完整接触线几何和独立端点压力')
        ends = geometry.endpoints
        p = np.asarray(pressure['pressure_pa'])
        if p.shape != ends.face_coordinate_m.shape or np.any(~np.isfinite(p)) or np.any(p < 0):
            raise ValueError('端点压力与几何布局不一致')
        samples.append(np.column_stack((ends.face_coordinate_m+result.pair.pinion.face_width_m/2,
            ends.roll_pinion_m, p, ends.extended_pinion, ends.extended_wheel)))
    return np.concatenate(samples)


def contact_pass_data(result, *, samples=(201, 241), face_coordinate_m=None, roll_distance_m=None):
    """输出实际啮合区的显示数据；指定坐标时可复用已绑定参照的显示网格。

    两轮接触展长交集不等于小轮完整渐开线。中心模型的边界支撑为原积分
    单元内常值显示延拓，不能解释为新增求解点；凸包外继续保留NaN。
    """
    boundary = []
    endpoints = [r.endpoint_pressure is not None for r in result.rotations]
    if all(endpoints):
        points = project_contact_pass_samples(result)
        definition = '独立端点压力'
        support_definition = '原独立端点，无额外边界显示延拓'
    elif any(endpoints) or any(r.contact_geometry_data is not None for r in result.rotations):
        raise ValueError('不能混合缺少端点压力或接触几何的相位')
    else:
        points = np.concatenate([np.column_stack((
            r.grid.face_coordinate_m+result.pair.pinion.face_width_m/2, r.grid.roll_pinion_m,
            r.pressure_pa, np.zeros((len(r.pressure_pa), 2)))) for r in result.rotations])
        definition = '条带中心压力（当前模型无独立端点）'
        cosine = np.cos(result.pair.base_helix_angle_rad)
        slope = np.tan(result.pair.base_helix_angle_rad)
        for rotation in result.rotations:
            grid = rotation.grid
            for pair_index in np.unique(grid.tooth_pair_index):
                indices = np.flatnonzero(grid.tooth_pair_index == pair_index)
                indices = indices[np.argsort(grid.face_coordinate_m[indices])]
                for index, sign in ((indices[0], -1), (indices[-1], 1)):
                    # widths_m是接触线长；先投影到轴向，再沿真实接触线到边界。
                    delta = sign*grid.widths_m[index]*cosine/2
                    boundary.append((grid.face_coordinate_m[index]+delta+result.pair.pinion.face_width_m/2,
                        grid.roll_pinion_m[index]+delta*slope, rotation.pressure_pa[index], 0., 0.))
        support_definition = '首尾中心值仅延拓到各自原积分单元边界；不是独立端点解'
    boundary = np.asarray(boundary, dtype=float).reshape(-1, 5)
    ordinary = points[~np.any(points[:, 3:5].astype(bool), axis=1)]
    display_samples = np.concatenate((ordinary, boundary))
    unique, inverse = np.unique(np.round(display_samples[:, :2], 12), axis=0, return_inverse=True)
    values = np.full(len(unique), -np.inf)
    np.maximum.at(values, inverse, display_samples[:, 2])
    left, right = result.pair.face_bounds_m
    x = (np.linspace(left+result.pair.pinion.face_width_m/2, right+result.pair.pinion.face_width_m/2, samples[0])
         if face_coordinate_m is None else np.asarray(face_coordinate_m, dtype=float))
    y = (np.linspace(*result.pair.contact_roll_bounds_m, samples[1])
         if roll_distance_m is None else np.asarray(roll_distance_m, dtype=float))
    for axis in (x, y):
        if axis.ndim != 1 or len(axis) < 2 or np.any(~np.isfinite(axis)) or np.any(np.diff(axis) <= 0):
            raise ValueError('显示坐标须为至少两个有限且严格递增的SI位置')
    field = np.full((len(x), len(y)), np.nan)
    interpolated = len(unique) >= 3 and np.linalg.matrix_rank(unique-unique.mean(axis=0), tol=1e-12) == 2
    if interpolated:
        import matplotlib.tri as tri
        triangulation = tri.Triangulation(unique[:, 0], unique[:, 1])
        xx, yy = np.meshgrid(x, y, indexing='ij')
        field = np.asarray(tri.LinearTriInterpolator(triangulation, values)(
            np.round(xx, 12), np.round(yy, 12)).filled(np.nan))
    return {'face_coordinate_m': x, 'roll_distance_m': y, 'pressure_pa_x_y': field,
            'samples': points, 'ordinary_samples': ordinary, 'boundary_display_samples': boundary,
            'contact_roll_bounds_m': np.asarray(result.pair.contact_roll_bounds_m),
            'pinion_involute_roll_bounds_m': np.asarray(result.pair.pinion.involute_roll_bounds_m),
            'definition': definition, 'support_definition': support_definition, 'interpolated': bool(interpolated)}


def draw_contact_pass(ax, data, *, maximum_pa, title=None):
    """源工程单图、参照对比和候选比较共用的填色、等压线、坐标和背景。"""
    cmap = contact_pressure_colormap()
    maximum = max(float(maximum_pa)*1e-6, 1.)
    if data['interpolated']:
        field = np.ma.masked_invalid(data['pressure_pa_x_y'].T)*1e-6
        image = ax.pcolormesh(data['face_coordinate_m']*1e3, data['roll_distance_m']*1e3, field,
                             shading='gouraud', cmap=cmap, vmin=0., vmax=maximum, rasterized=True)
        if field.count() and float(field.max()) > float(field.min()):
            ax.contour(data['face_coordinate_m']*1e3, data['roll_distance_m']*1e3, field,
                       levels=np.linspace(0., maximum, 31)[1:-1], colors='black', linewidths=.35)
    else:
        own = data['ordinary_samples']
        image = ax.scatter(own[:, 0]*1e3, own[:, 1]*1e3, c=own[:, 2]*1e-6, cmap=cmap, vmin=0., vmax=maximum)
    ax.set(xlim=data['face_coordinate_m'][[0, -1]]*1e3, ylim=data['roll_distance_m'][[0, -1]]*1e3,
           xlabel='齿宽坐标 / mm', ylabel='小轮展长 / mm')
    if title is not None:
        ax.set_title(title)
    ax.set_facecolor('#e5e9ee')
    return image


def plot_contact_pass_comparison(results, output):
    """同一齿轮副的候选通过程图：统一坐标及压力色标，指标仍取原数组。"""
    from .io import write_json
    from .plotting import _pyplot, _save
    if not results:
        raise ValueError('接触比较至少需要一个方案')
    first = next(iter(results.values()))
    if any(r.pair != first.pair for r in results.values()):
        raise ValueError('统一材料坐标比较要求相同齿轮副几何')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    data = {name: contact_pass_data(result) for name, result in results.items()}
    maximum = max(float(d['samples'][:, 2].max()) for d in data.values())
    fig, axes = plt.subplots(1, len(results), figsize=(max(8., len(results)*4.4), 5.6), squeeze=False,
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=.08, right=.86, bottom=.24, top=.84, wspace=.25)
    records, arrays = [], {}
    for index, ((name, result), ax) in enumerate(zip(results.items(), axes.flat)):
        d = data[name]
        peak_name = '端点峰压' if d['definition'] == '独立端点压力' else '中心峰压'
        image = draw_contact_pass(ax, d, maximum_pa=maximum,
            title=f"{name}\n{peak_name} {d['samples'][:, 2].max()*1e-6:.1f} MPa；PPTE {result.summary['te_peak_to_peak_m']*1e6:.3f} μm")
        for key, value in d.items():
            if isinstance(value, np.ndarray):
                arrays[f'case_{index+1:03d}_{key}'] = value
        records.append({'design': name, 'pressure_definition': d['definition'],
                        'support_definition': d['support_definition'], 'summary': result.summary,
                        'uncovered_display_count': int(np.count_nonzero(~np.isfinite(d['pressure_pa_x_y'])))})
    fig.colorbar(image, cax=fig.add_axes((.90, .30, .018, .48)), label='接触压力 / MPa')
    fig.suptitle('小轮整齿通过程接触压力比较', fontsize=13)
    fig.text(.5, .065, '坐标限于两轮实际啮合区；重复材料点取周期最大值，灰色为剩余未覆盖区。\n'
             '中心模型仅在原积分单元内提供边界显示支撑；本图不是瞬时接触面积。', ha='center', fontsize=8.2)
    files = _save(fig, output, '07_contact_pass_comparison')
    plt.close(fig)
    np.savez_compressed(output/'contact_pass_comparison_data.npz', **arrays)
    record = {'implementation': '源工程统一接触后处理', 'shared_maximum_pa': maximum,
              'raw_results_modified': False, 'records': records, 'files': [Path(p).name for p in files]}
    write_json(output/'contact_pass_comparison.json', record)
    return record


def load_plot_result(directory):
    """校验已导出的文件并构造只供重绘的视图；不重新求解，不作为计算响应。

    保留每相位中心／端点及接触线offsets布局。返回的是SimpleNamespace
    绘图视图，不冒充AnalysisResult或新增数值验证。
    """
    import hashlib
    import json
    from types import SimpleNamespace
    from .geometry import Gear, GearPair
    from .elasticity import Material
    from .config import modification_from_config
    directory = Path(directory).resolve()
    read = lambda name: json.loads((directory/name).read_text(encoding='utf-8-sig'))
    manifest = read('manifest.json')
    required = {'arrays.npz', 'input.json', 'summary.json', 'diagnostics.json'}
    if manifest.get('module') != 'gear_contact' or not required <= set(manifest.get('files', {})):
        raise ValueError('缺少有效的内核导出文件清单')
    for name, digest in manifest['files'].items():
        path = (directory/name).resolve()
        if not path.is_relative_to(directory) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('绘图来源文件身份不一致：'+name)
    inputs = read('input.json')
    pair_data = dict(inputs.pop('pair'))
    pair = GearPair(pinion=Gear(**pair_data.pop('pinion')), wheel=Gear(**pair_data.pop('wheel')), **pair_data)
    for role in ('pinion', 'wheel'):
        inputs['material_'+role] = Material(**inputs['material_'+role])
        inputs['modification_'+role] = modification_from_config(inputs['modification_'+role])
    diagnostics = read('diagnostics.json')
    with np.load(directory/'arrays.npz', allow_pickle=False) as source:
        arrays = {key: source[key].copy() for key in source.files}
    offsets, phases = arrays['cell_offsets'], arrays['phase']
    if (offsets.dtype.kind not in 'iu' or offsets.shape != (len(phases)+1,) or offsets[0] != 0
            or offsets[-1] != len(arrays['pressure_pa']) or np.any(np.diff(offsets) <= 0)):
        raise ValueError('原数组的相位／单元offsets无效')
    rotations = []
    for index, (start, end) in enumerate(zip(offsets[:-1], offsets[1:])):
        section = slice(start, end)
        grid_names = ('face_coordinate_m', 'roll_pinion_m', 'roll_wheel_m', 'widths_m', 'tooth_pair_index',
                      'curvature_per_m', 'profile_factor_pinion', 'profile_factor_wheel',
                      'lead_factor_pinion', 'lead_factor_wheel')
        grid = SimpleNamespace(**{key: arrays[key][section] for key in grid_names})
        solution = SimpleNamespace(**{key: arrays[key][section] for key in ('forces_n', 'residual_gaps_m', 'active_mask')},
                                   diagnostics=diagnostics['contact'][index])
        geometry, endpoint = None, None
        if 'rotation_endpoint_offsets' in arrays:
            e0, e1 = arrays['rotation_endpoint_offsets'][index:index+2]
            l0, l1 = arrays['rotation_line_offsets'][index:index+2]
            ends = SimpleNamespace(**{key: arrays['endpoint_'+key][e0:e1] for key in (
                'face_coordinate_m', 'roll_pinion_m', 'roll_wheel_m', 'tooth_pair_index',
                'extended_pinion', 'extended_wheel')},
                line_offsets=arrays['line_endpoint_offsets'][l0:l1+1]-e0)
            geometry = SimpleNamespace(endpoints=ends, line_offsets=arrays['line_cell_offsets'][l0:l1+1]-start,
                centre_extended_pinion=arrays['extended_pinion'][section], centre_extended_wheel=arrays['extended_wheel'][section])
            endpoint = {key: arrays['endpoint_'+key][e0:e1] for key in ('pressure_pa', 'half_width_m', 'line_load_n_per_m')}
        rotations.append(SimpleNamespace(grid=grid, solution=solution, contact_geometry_data=geometry,
            endpoint_pressure=endpoint, elastic_diagnostics=diagnostics['elasticity'][index],
            **{key: arrays[key][section] for key in ('pressure_pa', 'half_width_m', 'line_load_n_per_m', 'gaps_m')}))
    return SimpleNamespace(pair=pair, inputs=inputs, summary=read('summary.json'), rotations=rotations,
        plot_only=True, **{key: arrays[key] for key in ('phase', 'normal_approach_m', 'angular_te_rad', 'mesh_stiffness_n_per_m')})
