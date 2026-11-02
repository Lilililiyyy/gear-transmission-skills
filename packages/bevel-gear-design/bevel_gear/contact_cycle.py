"""同一参考工作齿面上的分相最大压力包络；不求和、不改变原接触解。

一个主动齿距内各齿对的压力先按实际刚体位姿映回材料坐标，再逐点取
最大值。相位和空间均为离散近似；包络不能积分成某个时刻的接触力。
未认证区域取各状态未知域的并集，不因另一相位已有正压力而改判已知。
"""
from __future__ import annotations

import math
from collections.abc import Sequence


def cycle_state_sequence_fingerprint(method, angles_rad, state_sha256, tooth_offsets):
    """把周期窗口证据绑定到确切方法、升序相位、原状态SHA及请求齿对。"""
    import hashlib
    import json
    if method not in ('full_gear_fe', 'local_halfspace') or len(angles_rad) != len(state_sha256):
        raise ValueError('窗口证据方法或原状态数量不一致')
    if any(not isinstance(v, str) or len(v) != 64 or any(c not in '0123456789abcdef' for c in v) for v in state_sha256):
        raise ValueError('窗口证据须提供实际状态的SHA256')
    return hashlib.sha256(json.dumps(dict(method=method, tooth_offsets=list(tooth_offsets),
        states=[dict(driving_angle_rad=float(a), source_state_sha256=s) for a, s in zip(angles_rad, state_sha256)]),
        sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()


def validate_cycle_window_coverage(evidence, *, method, angles_rad, state_sha256, tooth_offsets):
    """核查外邻齿审核的原状态身份与逐相位记录，拒绝手填passed代替证据。"""
    if evidence is None:
        return dict(passed=False, reason='未提供绑定原周期状态的外邻齿审核')
    if not tooth_offsets or list(tooth_offsets) != list(range(min(tooth_offsets), max(tooth_offsets)+1)):
        raise ValueError('完整齿对窗口须按连续整数升序给出，不能遗漏内部齿对')
    digest = cycle_state_sequence_fingerprint(method, angles_rad, state_sha256, tooth_offsets)
    if (not isinstance(evidence, dict) or evidence.get('schema') != 'bevel_mesh_cycle_window_coverage_v1' or
            evidence.get('method') != method or evidence.get('state_sequence_sha256') != digest or
            evidence.get('requested_tooth_offsets') != list(tooth_offsets)):
        raise ValueError('周期齿对窗口证据与方法、原状态或请求窗口不一致')
    required = ('all_requested_states_accepted', 'outer_pair_search_covered',
        'remaining_gap_nonnegative', 'source_states_unchanged')
    checks = evidence.get('checks', {})
    records = evidence.get('records', [])
    if not isinstance(checks, dict) or any(checks.get(k) is not True for k in required):
        raise ValueError('周期窗口缺少原状态接受、外邻齿覆盖或剩余间隙核查')
    if not isinstance(records, list) or len(records) != len(angles_rad):
        raise ValueError('周期窗口审核记录未覆盖全部实际相位')
    outside = [min(tooth_offsets)-1, max(tooth_offsets)+1]
    for i, record in enumerate(records):
        if (record.get('phase_index') != i or record.get('source_state_sha256') != state_sha256[i] or
                not math.isclose(record.get('driving_angle_rad', float('nan')), angles_rad[i], rel_tol=0., abs_tol=1e-12)):
            raise ValueError('周期外齿记录与原相位/状态SHA不一致')
        pairs = record.get('outer_pairs', [])
        if sorted(p.get('tooth_offset') for p in pairs) != outside or any(p.get('passed') is not True for p in pairs):
            raise ValueError('周期外邻齿排除未完整通过')
        for pair in pairs:
            gap = pair.get('minimum_remaining_gap_m')
            if (isinstance(gap, bool) or not isinstance(gap, (int, float)) or not math.isfinite(gap) or gap < 0. or
                    not isinstance(pair.get('evidence_kind'), str) or not pair['evidence_kind'].strip()):
                raise ValueError('周期外齿须保留实际非负剩余间隙及证据种类')
    if evidence.get('passed') is not True:
        raise ValueError('周期窗口证据仍未通过，不能认证完整单齿历程')
    return dict(evidence, verified_for_current_cycle=True)


def compact_cycle_window_coverage(evidence):
    """只导出方法、原状态摘要和外齿余隙，不传播开发目录或几何过程缓存。"""
    if not evidence or evidence.get('verified_for_current_cycle') is not True:
        return dict(passed=False, verified_for_current_cycle=False)
    keys = ('schema', 'method', 'requested_tooth_offsets', 'outer_tooth_offsets',
        'state_sequence_sha256', 'checks', 'passed', 'verified_for_current_cycle', 'interpretation')
    output = {k: evidence[k] for k in keys if k in evidence}
    output['records'] = [dict(phase_index=r['phase_index'], driving_angle_rad=r['driving_angle_rad'],
        source_state_sha256=r['source_state_sha256'], outer_pairs=[
            {k: p[k] for k in ('tooth_offset', 'passed', 'minimum_remaining_gap_m', 'evidence_kind')}
            for p in r['outer_pairs']]) for r in evidence['records']]
    return output


def validate_mesh_cycle_angles(angles_rad, pinion_teeth):
    """校验含重复端点的一齿距均匀分相；缺相或不完整周期显式报错。"""
    if isinstance(pinion_teeth, bool) or not isinstance(pinion_teeth, int) or pinion_teeth < 3:
        raise ValueError('小轮齿数须为不小于3的整数')
    if not isinstance(angles_rad, (list, tuple)) or len(angles_rad) < 5:
        raise ValueError('完整分相周期至少需要4个区间及重复端点')
    if any(isinstance(a, bool) or not isinstance(a, (int, float)) or not math.isfinite(a) for a in angles_rad):
        raise ValueError('啮合相位须为有限弧度值')
    angles = [float(a) for a in angles_rad]
    pitch = 2 * math.pi / pinion_teeth
    step = pitch / (len(angles) - 1)
    if not math.isclose(angles[-1] - angles[0], pitch, rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError('相位起止必须覆盖恰好一个主动齿距，包含重复端点')
    if any(not math.isclose(b-a, step, rel_tol=1e-10, abs_tol=1e-12) for a, b in zip(angles, angles[1:])):
        raise ValueError('全周期相位须按均匀升序给出，不能缺相、重复或重排')
    return dict(pinion_teeth=pinion_teeth, pinion_pitch_rad=pitch,
        start_angle_rad=angles[0], end_angle_rad=angles[-1], phase_step_rad=step,
        interval_count=len(angles)-1, solved_phase_count=len(angles),
        repeated_endpoint=True, phase_resolution_converged=False)


def evaluate_p1_pressure(field, query_points_m):
    """在给定二维m坐标精确求原三角形P1值；支撑域外返回零及独立掩码。

    此处域外零仅用于最大值算术，不能作为非接触证明。调用方必须另外
    传播覆盖证书的已知域。源场极值、原节点力及积分面积均不修改。
    """
    import numpy as np
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    points = np.asarray(query_points_m, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError('显示查询点须为有限二维m坐标')
    xy = np.asarray(field['points_width_height_m'], dtype=float)
    p = np.asarray(field['pressure_pa'], dtype=float)
    triangles = np.asarray(field['triangles'])
    if (xy.ndim != 2 or xy.shape[1] != 2 or p.shape != (len(xy),) or
            not np.isfinite(xy).all() or not np.isfinite(p).all() or np.any(p < 0)):
        raise ValueError('原连续压力场须有有限坐标和非负Pa标量')
    if (triangles.ndim != 2 or triangles.shape[1] != 3 or not len(triangles) or
            triangles.dtype.kind not in 'iu' or np.any(triangles < 0) or np.any(triangles >= len(xy))):
        raise ValueError('原P1压力三角形索引无效')
    mesh = Triangulation(xy[:, 0], xy[:, 1], triangles)
    interpolator = LinearTriInterpolator(mesh, p)
    values = interpolator(points[:, 0], points[:, 1])
    valid = ~np.ma.getmaskarray(values)
    pressure = np.asarray(values.filled(0.), dtype=float)
    # 仅容许浮点舍入量级负值，不接受插值越过原场极值。
    tolerance = max(float(p.max()), 1.) * 1e-10
    if np.any(pressure < -tolerance) or np.any(pressure > p.max()+tolerance):
        raise ArithmeticError('P1求值超出原非负压力范围')
    return np.maximum(pressure, 0.), valid


def sampled_pressure_maximum(fields, query_points_m, *, known_masks=None, sampling_strides=()):
    """逐点取max并保存达到该值的相位/齿对；未知域始终按并集传播。"""
    import numpy as np
    if not isinstance(fields, (list, tuple)) or not fields:
        raise ValueError('压力包络需要非空的真实分相场')
    points = np.asarray(query_points_m, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError('最大值查询点须为有限二维m坐标')
    if known_masks is not None and len(known_masks) != len(fields):
        raise ValueError('覆盖掩码须与原相位/齿对场逐一对应')
    maximum = np.zeros(len(points))
    phase = np.full(len(points), -1, dtype=int)
    tooth = np.full(len(points), 0, dtype=int)
    known = np.ones(len(points), dtype=bool)
    support = np.zeros(len(points), dtype=bool)
    identities = set()
    if any(isinstance(s, bool) or not isinstance(s, int) or s < 2 for s in sampling_strides):
        raise ValueError('相位降采样步长须为不小于2的整数')
    coarse = {s: np.zeros(len(points)) for s in sampling_strides}
    for index, item in enumerate(fields):
        identity = (item['phase_index'], item['tooth_offset'])
        if (any(isinstance(v, bool) or not isinstance(v, int) for v in identity) or
                identity[0] < 0 or identity in identities):
            raise ValueError('相位/齿对场身份须唯一且相位索引非负')
        identities.add(identity)
        field = item.get('field')
        if field is None:
            pressure = np.zeros(len(points)); inside = np.zeros(len(points), dtype=bool)
        else:
            pressure, inside = evaluate_p1_pressure(field, points)
        select = pressure > maximum
        maximum[select] = pressure[select]
        for stride, values in coarse.items():
            if identity[0] % stride == 0:
                np.maximum(values, pressure, out=values)
        phase[select] = identity[0]; tooth[select] = identity[1]
        support |= inside
        if known_masks is None:
            # 没有外域证据时仅承认该场的实际P1支撑，不能把域外算术零当物理解。
            known &= inside
        else:
            mask = np.asarray(known_masks[index])
            if mask.dtype != np.bool_ or mask.shape != (len(points),):
                raise ValueError('覆盖掩码须为逐查询点的布尔数组')
            known &= mask
    return dict(pressure_pa=maximum, known_mask=known, positive_support_mask=support,
        maximum_phase_index=phase, maximum_tooth_offset=tooth,
        sampling_maxima=coarse,
        definition='p_env(w,h)=max_{离散相位i、周期等价齿k} p_i,k(w,h)',
        force_integral_defined=False, original_pressure_values_modified=False)


def _geometry_union(polygons_m):
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    return unary_union([Polygon(p) for p in polygons_m])


def _covered(geometry, points):
    """GEOS批量覆盖判断；保留孔洞，不缓冲、不吸附、不用凸包补齐。"""
    import numpy as np
    import shapely
    if geometry.is_empty:
        return np.zeros(len(points), dtype=bool)
    shapely.prepare(geometry)
    return shapely.covers(geometry, shapely.points(points))


def build_mesh_cycle_pressure_data(results, *, method, pinion_request, flank_mesh,
        domain_certificates=None, spacing_m=2e-5, window_coverage=None, phase_angles_rad=None, progress=None):
    """把完整一齿距各齿对映回同一材料齿面，形成有未知标志的最大压力图。

    results按主动角升序含重复端点。证书映射键为"相位索引:齿对偏移"。
    spacing_m只控制显示采样，不改变FE节点或半空间求解格距。
    """
    import numpy as np
    from .generated_surface import _request
    from .refined_fe_plotting import _background
    from .pressure_field import reconstruct_fe_pressure_field, reconstruct_halfspace_pressure_field, project_pressure_field
    from .pressure_background import pressure_domain_background, pressure_domain_metadata
    if method not in ('full_gear_fe', 'local_halfspace'):
        raise ValueError('方法须为full_gear_fe或local_halfspace')
    if progress is not None and not callable(progress):
        raise ValueError('进度回调须可调用')
    if (isinstance(spacing_m, bool) or not isinstance(spacing_m, (int, float)) or
            not math.isfinite(spacing_m) or spacing_m <= 0):
        raise ValueError('显示格距须为正有限m')
    if not isinstance(results, Sequence) or isinstance(results, (str, bytes)) or not results:
        raise ValueError('需要完整周期真实结果，不能传入空列表')
    request = _request(pinion_request)
    angles = (list(phase_angles_rad) if phase_angles_rad is not None else
        [r['driving_angle_rad'] if method == 'full_gear_fe' else r['geometry']['driving_angle_rad'] for r in results])
    if len(angles) != len(results):
        raise ValueError('相位元数据数量与实际状态序列不一致')
    cycle = validate_mesh_cycle_angles(angles, request.tooth_count)
    _, background, _ = _background(flank_mesh, pinion_request)
    outline = _geometry_union(background)
    xmin, ymin, xmax, ymax = outline.bounds
    nx = math.ceil((xmax-xmin)/spacing_m)+1; ny = math.ceil((ymax-ymin)/spacing_m)+1
    if nx*ny > 12_000_000:
        raise ValueError('压力显示网格超过1200万点；请显式增大显示格距')
    fields = []; phase_rows = []; audit_rows = []; unknown_geometries = []; peak_coordinates = []
    all_source_peaks = []; offsets_reference = None; signed_torque = None
    for index, result in enumerate(results):
        if progress is not None: progress(dict(stage='cycle_pressure_fields', method=method, phase_index=index, completed=False))
        if not isinstance(result, dict) or not result.get('converged'):
            raise ValueError('完整周期包含未接受相位，不能丢弃或以零压代替')
        actual_angle = result['driving_angle_rad'] if method == 'full_gear_fe' else result['geometry']['driving_angle_rad']
        if not math.isclose(actual_angle, angles[index], rel_tol=0., abs_tol=1e-12):
            raise ValueError('相位元数据与实际压力状态不一致')
        if method == 'full_gear_fe':
            offsets = sorted(int(k) for k in result['pair_force_shares'])
            torque = result['pinion_torque_nm']
        else:
            offsets = sorted(result['geometry'].get('build_parameters', {}).get('tooth_offsets',
                [p['tooth_offset'] for p in result['geometry']['patches']]))
            torque = result['geometry']['pinion_torque_nm']
        if offsets_reference is None:
            offsets_reference = offsets; signed_torque = torque
        if (offsets != offsets_reference or isinstance(torque, bool) or not isinstance(torque, (int, float)) or
                not math.isfinite(torque) or not math.isclose(torque, signed_torque, rel_tol=1e-12, abs_tol=1e-10)):
            raise ValueError('周期中的实际齿对窗口和有符号扭矩须一致')
        pair_rows = []
        for offset in offsets:
            if method == 'full_gear_fe':
                loaded = any(n['tooth_offset'] == offset and n.get('active') and n['force_n'] > 0 for n in result['nodal_contacts'])
                reconstruct = reconstruct_fe_pressure_field
                frames = result['mesh_reference_frames']
                force = result['pair_force_shares'].get(str(offset), result['pair_force_shares'].get(offset, 0.))
            else:
                field_index = result['selected_tooth_offsets'].index(offset) if offset in result['selected_tooth_offsets'] else None
                loaded = field_index is not None and bool(np.any(np.asarray(result['fields'][field_index]['forces_n']) > 0))
                reconstruct = reconstruct_halfspace_pressure_field
                frames = next(p['frames'] for p in result['geometry']['patches'] if p['tooth_offset'] == offset)
                force = 0. if field_index is None else float(np.asarray(result['fields'][field_index]['forces_n']).sum())
            field = project_pressure_field(reconstruct(result, tooth_offset=offset), pinion_request, frames) if loaded else None
            key = f'{index}:{offset}'
            certificate = domain_certificates.get(key) if domain_certificates is not None else None
            domain = pressure_domain_background(result, background, method=method, tooth_offset=offset,
                certificate=certificate, pinion_request=pinion_request)
            if domain['full_flank_zero_certified']:
                pass
            elif domain['domain_state'] == 'partitioned_pressure_domain':
                unknown = _geometry_union(domain['_unresolved_polygons_m'])
                unknown_geometries.append(unknown)
            else:
                if field is None:
                    unknown_geometries.append(outline)
                else:
                    support_geometry = _geometry_union(np.asarray(field['points_width_height_m'])[field['triangles']])
                    unknown_geometries.append(outline.difference(support_geometry))
            fields.append(dict(phase_index=index, tooth_offset=offset, field=field))
            if method == 'full_gear_fe':
                peak = max((float(n['pressure_pa']) for n in result['nodal_contacts']
                    if n['tooth_offset'] == offset and n.get('active')), default=0.)
            else:
                if field_index is None:
                    peak = 0.
                else:
                    original_field = result['fields'][field_index]
                    peak = float(np.asarray(original_field['pressure_pa'], dtype=float)[
                        np.asarray(original_field['valid_mask'], dtype=bool)].max(initial=0.))
            reconstructed_peak = 0. if field is None else float(np.max(field['pressure_pa']))
            if peak > reconstructed_peak * (1+1e-10) + 1e-6:
                raise ValueError('原受载峰值不在有效P1支撑内，不能将不完整重构冒充完整压力包络')
            if field is not None:
                peak_coordinates.append(np.asarray(field['points_width_height_m'])[
                    int(np.argmax(field['pressure_pa']))])
            all_source_peaks.append(peak)
            row = dict(phase_index=index, tooth_offset=offset, driver_angle_rad=angles[index],
                equivalent_reference_tooth_angle_rad=angles[index]+offset*cycle['pinion_pitch_rad'],
                peak_pressure_pa=peak, reconstructed_peak_pressure_pa=reconstructed_peak,
                normal_force_n=force, loaded=loaded,
                field_vertex_count=0 if field is None else len(field['pressure_pa']),
                display_field_diagnostics={} if field is None else field.get('diagnostics', {}))
            pair_rows.append(row)
            audit_rows.append(dict(pressure_domain_metadata(domain), phase_index=index, tooth_offset=offset))
        phase_rows.append(dict(phase_index=index, driver_angle_rad=angles[index],
            mesh_fraction=(angles[index]-angles[0])/cycle['pinion_pitch_rad'],
            peak_pressure_pa=max(r['peak_pressure_pa'] for r in pair_rows),
            loaded_tooth_pair_count=sum(r['loaded'] for r in pair_rows), pairs=pair_rows))
        if progress is not None: progress(dict(stage='cycle_pressure_fields', method=method, phase_index=index, completed=True))
    # 每个原场的峰值坐标进入查询轴，避免均匀显示网格漏掉狭窄原节点峰。
    peaks_xy = np.asarray(peak_coordinates).reshape(-1, 2)
    x = np.unique(np.concatenate((np.linspace(xmin, xmax, nx), peaks_xy[:, 0])))
    y = np.unique(np.concatenate((np.linspace(ymin, ymax, ny), peaks_xy[:, 1])))
    x = x[(x >= xmin) & (x <= xmax)]; y = y[(y >= ymin) & (y <= ymax)]
    if len(x)*len(y) > 12_000_000:
        raise ValueError('包含原峰值坐标的显示网格过大，请增大基础显示格距')
    xx, yy = np.meshgrid(x, y); points = np.column_stack((xx.ravel(), yy.ravel()))
    in_flank = _covered(outline, points); query = points[in_flank]
    from shapely.ops import unary_union
    if progress is not None: progress(dict(stage='cycle_unknown_union', method=method, geometry_count=len(unknown_geometries), completed=False))
    unknown_union = unary_union(unknown_geometries).intersection(outline)
    if progress is not None: progress(dict(stage='cycle_unknown_union', method=method, completed=True))
    known = ~_covered(unknown_union, query)
    strides = (4, 2) if cycle['interval_count'] >= 16 and cycle['interval_count'] % 4 == 0 else ()
    envelope = sampled_pressure_maximum(fields, query, known_masks=[np.ones(len(query), bool)]*len(fields),
        sampling_strides=strides)
    if progress is not None: progress(dict(stage='cycle_pressure_maximum', method=method, query_points=len(query), completed=True))
    envelope['known_mask'] = known
    shape = xx.shape
    arrays = {}
    for key, outside in (('pressure_pa', float('nan')), ('known_mask', False),
            ('maximum_phase_index', -1), ('maximum_tooth_offset', 0)):
        source = envelope[key]; array = np.full(len(points), outside, dtype=source.dtype)
        array[in_flank] = source; arrays[key] = array.reshape(shape)
    peak = max(all_source_peaks)
    display_peak = float(envelope['pressure_pa'].max())
    if peak > 0 and (peak-display_peak)/peak > .01:
        raise ValueError('显示域遗漏原压力峰值超过1%；请核查实际CAD裁切边界与原场，不生成误导的全蓝图')
    evidence = compact_cycle_window_coverage(window_coverage)
    # end(k)与start(k+1)才是同一周期等价接触位置，不能比较同名k。
    start_pairs = {r['tooth_offset']: r for r in phase_rows[0]['pairs']}
    end_pairs = {r['tooth_offset']: r for r in phase_rows[-1]['pairs']}
    endpoint_pairs = []
    for k, end in end_pairs.items():
        if k+1 not in start_pairs:
            continue
        start = start_pairs[k+1]
        reference_peak = max(start['peak_pressure_pa'], end['peak_pressure_pa'])
        endpoint_pairs.append(dict(start_tooth_offset=k+1, end_tooth_offset=k,
            normal_force_difference_n=end['normal_force_n']-start['normal_force_n'],
            peak_pressure_difference_pa=end['peak_pressure_pa']-start['peak_pressure_pa'],
            peak_relative_difference=0. if reference_peak == 0 else
                abs(end['peak_pressure_pa']-start['peak_pressure_pa'])/reference_peak))
    endpoint = dict(pair_mapping='end(k) <-> start(k+1)', pairs=endpoint_pairs,
        mapped_pair_force_l1_difference_n=sum(abs(r['normal_force_difference_n']) for r in endpoint_pairs),
        raw_global_peak_difference_pa=phase_rows[-1]['peak_pressure_pa']-phase_rows[0]['peak_pressure_pa'],
        unmatched_start_tail_force_n=start_pairs[min(start_pairs)]['normal_force_n'],
        unmatched_end_tail_force_n=end_pairs[max(end_pairs)]['normal_force_n'],
        pressure_periodicity_tolerance_certified=False)
    # 同一实际解集内抽取嵌套16/32/64等分，比较时间离散，不混用不同求解格距。
    def dual_width(axis):
        spacing = np.diff(axis)
        return np.r_[spacing[0]*.5, (spacing[:-1]+spacing[1:])*.5, spacing[-1]*.5]
    weights = np.outer(dual_width(y), dual_width(x)).ravel()[in_flank]
    weights = weights * envelope['known_mask']
    reference_integral = float(np.sum(weights*envelope['pressure_pa']))
    sampling_study = []
    for stride, coarse_pressure in envelope['sampling_maxima'].items():
        source_peak = max((row['peak_pressure_pa'] for row in phase_rows if row['phase_index'] % stride == 0), default=0.)
        sampling_study.append(dict(interval_count=cycle['interval_count']//stride,
            reference_interval_count=cycle['interval_count'], raw_sampled_peak_pressure_pa=source_peak,
            relative_peak_increase_to_reference=0. if peak == 0 else (peak-source_peak)/peak,
            known_domain_normalized_mean_absolute_difference=None if reference_integral <= 0 else
                float(np.sum(weights*np.abs(envelope['pressure_pa']-coarse_pressure))/reference_integral),
            same_solved_states_at_common_phases=True, solver_discretization_controls_changed=False,
            phase_convergence_certified=False))
    cycle.update(tooth_offsets=offsets_reference,
        single_tooth_path_coverage_confirmed=evidence.get('passed') is True and evidence.get('verified_for_current_cycle') is True,
        window_coverage=evidence, periodic_mesh_discretization_exact=False)
    return dict(schema='bevel_mesh_cycle_pressure_data_v1', method=method, cycle=cycle,
        x_m=x, y_m=y, inside_flank_mask=in_flank.reshape(shape), **arrays,
        background_polygons_m=background, unknown_geometry=unknown_union,
        phase_rows=phase_rows, domain_audits=audit_rows, signed_pinion_torque_nm=signed_torque,
        source_peak_pressure_pa=peak, display_grid_peak_pressure_pa=display_peak,
        display_peak_relative_loss=0. if peak == 0 else (peak-display_peak)/peak,
        actual_display_spacing_m=[float(np.max(np.diff(x))), float(np.max(np.diff(y)))],
        display_grid='均匀基础网格加原场峰值坐标；实际最大片距单列',
        phase_sampling_study=sampling_study,
        endpoint_periodicity=endpoint,
        unresolved_projected_area_fraction=float(unknown_union.area/outline.area),
        complete_zero_domain_certified=unknown_union.is_empty,
        numerical_definition=envelope['definition'], force_integral_defined=False,
        original_pressure_values_modified=False, pressure_mesh_converged=False,
        physical_validation=False)
