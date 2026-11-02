"""实际弧齿制造面→严格几何检查→局部半空间压力的可重放业务入口。

导入本模块只需要标准库。NumPy、OCP、SciPy和图形依赖均在实际计算
或绘图时加载。每个相位独立保存完整间隙、真实四角、压力与有效掩码；
汇总只引用相对文件和标量，不把局部半空间闭合量解释为整轮FE总TE。
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import time


SURFACE_CASE_SCHEMA = 'bevel_surface_contact_case_v1'


def _real(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name+'须为有限实数')
    if positive and value <= 0:
        raise ValueError(name+'须为正数')
    return float(value)


def surface_json_safe(value):
    """完整结果转标准JSON：数组转列表，grid转字典，未知NaN转null。

    无效几何仍由原始valid_mask标识；转换不将未知间隙改成零。
    不接受任意CAD对象或用字符串代替未序列化的物理数据。
    """
    if isinstance(value, dict):
        return {str(key): surface_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [surface_json_safe(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if callable(getattr(value, 'as_dict', None)):
        return surface_json_safe(value.as_dict())
    if callable(getattr(value, 'tolist', None)):
        return surface_json_safe(value.tolist())
    if callable(getattr(value, 'item', None)):
        return surface_json_safe(value.item())
    raise TypeError('接触结果含不能序列化的对象：'+type(value).__name__)


def _write(value, path):
    """同目录临时文件完成后替换，避免中断留下半份压力JSON。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(surface_json_safe(value), ensure_ascii=False,
        allow_nan=False, separators=(',', ':'))+'\n', encoding='utf-8', newline='\n')
    temporary.replace(path)
    return path


def _fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
        allow_nan=False, separators=(',', ':')).encode('utf-8')).hexdigest()


def _source_identity(*, nominal=False, plots=False):
    """只记录实际计算源码的公开相对名和SHA256，不导出本机路径。"""
    names = {'surface_workflow', 'halfspace_contact', 'halfspace_boundary', 'generated_halfspace',
        'generated_tca', 'generated_surface', 'generated_cutter', 'generated_solid',
        'cad', 'cad_import', 'rigid_contact', 'flank_relief', 'nurbs_batch'}
    if nominal:
        names.update(('workflow', 'gear_pair', 'geometry', 'forces', 'iso2014_nominal',
            'manufacturing_binding', 'strength_assessment'))
    if plots:
        names.update(('halfspace_plotting', 'pressure_field', 'contact_pattern', 'workflow_plotting'))
    directory = Path(__file__).resolve().parent
    return {'bevel_gear/'+name+'.py': hashlib.sha256((directory/(name+'.py')).read_bytes()).hexdigest()
            for name in sorted(names)}


def load_surface_case(path):
    """读取独立surface-case及其相对base_case，展开显式制造输入。

    原几何/强度/FE案例不作修改；半空间参数只存在本surface-case。
    """
    from .workflow import load_bevel_case
    path = Path(path).resolve()
    case = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(case, dict) or case.get('schema') != SURFACE_CASE_SCHEMA:
        raise ValueError('surface-case schema须为'+SURFACE_CASE_SCHEMA)
    if not isinstance(case.get('base_case'), str) or not case['base_case']:
        raise ValueError('surface-case须显式引用base_case文件')
    base_path = path.parent/case['base_case']
    base = load_bevel_case(base_path)
    if base.get('kind') != 'spiral':
        raise ValueError('实际制造面半空间入口当前仅支持spiral；直齿尚未接入')
    for member in ('pinion', 'wheel'):
        manufacturing = base.get('manufacturing', {}).get(member, {})
        if not isinstance(manufacturing.get('request'), dict) or not isinstance(manufacturing.get('blank'), dict):
            raise ValueError('base_case须显式提供双轮制造request和blank')
    return dict(surface_case=deepcopy(case), base_case=base,
        surface_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        base_input_sha256=hashlib.sha256(base_path.read_bytes()).hexdigest(),
        expanded_input_fingerprint=_fingerprint(dict(surface_case=case, base_case=base)))


def _configuration(loaded, phase_indices, spacing_m, allow_initial_penetration):
    """在创建实体之前检查相位、SI网格和严格非线性验收策略。"""
    case, base = loaded['surface_case'], loaded['base_case']
    analysis = deepcopy(case.get('analysis', {}))
    phases = analysis.pop('driver_angles_rad', None)
    if not isinstance(phases, list) or not phases:
        raise ValueError('analysis.driver_angles_rad须为非空实际相位列表，单位rad')
    phases = [_real(value, '主动相位rad') for value in phases]
    if len(set(phases)) != len(phases):
        raise ValueError('主动相位不能重复')
    selected = list(range(len(phases))) if phase_indices is None else list(phase_indices)
    if not selected or len(set(selected)) != len(selected) or any(
            isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(phases) for index in selected):
        raise ValueError('phase_indices须为不重复且存在的相位索引')
    allowed = {'tooth_offsets', 'grid_shape', 'spacing_m', 'position_tolerance_m',
        'maximum_surface_slope', 'maximum_normal_angle_rad', 'check_cell_corners',
        'cover_complete_working_flank', 'follower_correction_rad', 'ray_backend'}
    if set(analysis)-allowed:
        raise ValueError('不支持的表面几何参数：'+','.join(sorted(set(analysis)-allowed)))
    if spacing_m is not None:
        analysis['spacing_m'] = list(spacing_m)
    spacings = analysis.get('spacing_m', (1e-4, 2.5e-5))
    if not isinstance(spacings, (list, tuple)) or len(spacings) != 2:
        raise ValueError('spacing_m须为(dx,dy)，单位m')
    analysis['spacing_m'] = [_real(value, '接触格距m', positive=True) for value in spacings]
    if analysis.get('check_cell_corners', True) is not True:
        raise ValueError('业务入口要求真实中心及四角均通过裁切域检查')
    analysis['check_cell_corners'] = True
    analysis.setdefault('cover_complete_working_flank', True)
    if not isinstance(analysis['cover_complete_working_flank'], bool):
        raise ValueError('cover_complete_working_flank须为布尔值')
    if 'grid_shape' in analysis:
        shape = analysis['grid_shape']
        if not isinstance(shape, (list, tuple)) or len(shape) != 2 or any(
                isinstance(n, bool) or not isinstance(n, int) or n < 3 for n in shape):
            raise ValueError('grid_shape须为(nx,ny)两个不小于3的整数')
    if 'tooth_offsets' in analysis:
        offsets = analysis['tooth_offsets']
        if not isinstance(offsets, (list, tuple)) or not offsets or len(set(offsets)) != len(offsets) or any(
                isinstance(k, bool) or not isinstance(k, int) for k in offsets):
            raise ValueError('tooth_offsets须为不重复整数')
    for name in ('position_tolerance_m', 'maximum_surface_slope', 'maximum_normal_angle_rad'):
        if name in analysis:
            _real(analysis[name], name, positive=True)
    if 'follower_correction_rad' in analysis:
        # 显式共同刚体相位改变后，正式几何入口重新求双面射线；不夹零间隙。
        _real(analysis['follower_correction_rad'], '共同从动轮相位修正rad')
    solver = deepcopy(case.get('solver', {}))
    supported = {'maximum_iterations', 'gap_tolerance_m', 'relative_tolerance',
        'allow_initial_penetration', 'maximum_initial_penetration_m',
        'maximum_outer_iterations', 'maximum_linearization_error_m'}
    if set(solver)-supported:
        raise ValueError('不支持的半空间求解参数：'+','.join(sorted(set(solver)-supported)))
    solver.setdefault('gap_tolerance_m', 1e-10)
    solver.setdefault('maximum_linearization_error_m', 1e-10)
    for name in ('gap_tolerance_m', 'maximum_linearization_error_m'):
        if _real(solver[name], name, positive=True) > 1e-10:
            raise ValueError(name+'不得宽于1e-10m严格验收门限')
    solver.setdefault('allow_initial_penetration', False)
    if allow_initial_penetration is not None:
        solver['allow_initial_penetration'] = allow_initial_penetration
    if not isinstance(solver['allow_initial_penetration'], bool):
        raise ValueError('allow_initial_penetration须为显式布尔值')
    for name in ('maximum_iterations', 'maximum_outer_iterations'):
        if name in solver and (isinstance(solver[name], bool) or not isinstance(solver[name], int) or solver[name] < 1):
            raise ValueError(name+'须为正整数')
    for name in ('relative_tolerance', 'maximum_initial_penetration_m'):
        if name in solver:
            _real(solver[name], name, positive=True)
    torque = _real(base['load']['pinion_torque_n_m'], '主动扭矩幅值Nm', positive=True)
    sign = base['load']['ltca_torque_sign']
    if isinstance(sign, bool) or sign not in (-1, 1):
        raise ValueError('ltca_torque_sign须为-1或+1')
    materials = {}
    for member in ('pinion', 'wheel'):
        values = base['materials'][member]
        materials[member+'_young_modulus_pa'] = _real(values['young_modulus_pa'], member+'弹性模量Pa', positive=True)
        poisson = _real(values['poisson_ratio'], member+'泊松比')
        if not -1 < poisson < .5:
            raise ValueError('泊松比须在(-1,.5)')
        materials[member+'_poisson_ratio'] = poisson
    return phases, selected, analysis, solver, torque*sign, materials


def _accepted(result, source_unchanged, maximum_pinion_torque_relative_residual):
    history = result.get('geometry_outer_history') or []
    residual = result.get('pinion_torque_relative_residual')
    return bool(result.get('converged') and result.get('fields') and source_unchanged
        and result.get('geometry_linearization_checked') is True and history
        and history[-1].get('geometry_check', {}).get('passed') is True
        and result.get('active_contact_domain_accepted') is True
        and isinstance(residual, (int, float)) and math.isfinite(residual)
        and abs(residual) <= maximum_pinion_torque_relative_residual)


def _flank_mesh(prepared, *, linear_deflection_m=2e-5, angular_deflection_rad=.08):
    """只三角化同一真实裁切小轮工作面compound，不重造背景齿形。"""
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    from .cad_import import ImportedCadModel, triangulate_imported_model
    compound = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(compound)
    for face in prepared.trimmed['pinion'].faces:
        builder.Add(compound, face)
    mesh = triangulate_imported_model(ImportedCadModel(compound, dict(units='m',
        interpretation='同一实际闭实体中slot0小轮左裁切工作面')),
        linear_deflection_m=linear_deflection_m, angular_deflection_rad=angular_deflection_rad)
    model=prepared.tca.cad_models['pinion']
    source=model.grid['patches']['left_flank']
    mesh['metadata']['source_surface_sampling']=dict(width_count=source['width_count'],
        profile_count=source['profile_count'])
    return mesh


def run_surface_workflow(surface_case_path, output_dir, *, phase_indices=None, spacing_m=None,
        plots=False, nominal=False, reference=None, cad_models=None,
        allow_initial_penetration=None, progress=None):
    """运行显式弧齿实际制造面的独立半空间业务；逐相位立即落盘。

    cad_models仅允许调用方显式传入同请求的真实实体对象；公开CLI始终
    从request+blank新建实体，不读取私有缓存。reference只用于可选的
    through=strength名义几何/强度流程；不参与表面压力求解。
    """
    for value, name in ((plots, 'plots'), (nominal, 'nominal')):
        if not isinstance(value, bool):
            raise ValueError(name+'须为布尔值')
    if reference is not None and not nominal:
        raise ValueError('reference只用于显式开启的nominal几何/强度流程')
    loaded = load_surface_case(surface_case_path)
    phases, selected, analysis, solver, torque, materials = _configuration(
        loaded, phase_indices, spacing_m, allow_initial_penetration)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    maximum_residual = _real(loaded['base_case'].get('loaded_contact', {}).get(
        'maximum_pinion_torque_relative_residual', .02), '主动扭矩相对残差门限', positive=True)
    identity = _source_identity(nominal=nominal, plots=plots)
    effective_fingerprint = _fingerprint(dict(base_case=loaded['base_case'], analysis=analysis,
        driver_angles_rad=phases, solver=solver, pinion_signed_torque_nm=torque, materials=materials))
    summary = dict(schema='bevel_surface_workflow_v1', analysis_completed=False,
        name=loaded['surface_case'].get('name', '实际弧齿半空间表面接触'), kind='spiral',
        input_file='surface_inputs.json', input_fingerprint=effective_fingerprint,
        source_identity_before=identity, source_identity_after=None, source_identity_unchanged=None,
        cad_mode='explicit_same_request_objects' if cad_models is not None else 'built_from_explicit_request_and_blank',
        closure_definition='局部弹性半空间的从动参数轴闭合量rad；不含齿体/轮体/跨齿结构柔度，不是整轮FE总TE',
        selected_phase_indices=selected, states=[], plots_status='not_requested', nominal_file=None,
        physical_validation=False, pressure_mesh_converged=None)
    inputs = dict(**loaded, effective_analysis=analysis, effective_solver=solver,
        effective_input_fingerprint=effective_fingerprint,
        pinion_signed_torque_nm=torque, materials=materials,
        maximum_pinion_torque_relative_residual=maximum_residual)
    _write(inputs, output/'surface_inputs.json')
    _write(summary, output/'surface_summary.json')

    def emit(stage, **values):
        if progress is not None:
            progress(dict(stage=stage, **values))

    if nominal:
        from .workflow import analyze_bevel_case
        if isinstance(reference, (str, Path)):
            reference = json.loads(Path(reference).read_text(encoding='utf-8-sig'))
        analyze_bevel_case(loaded['base_case'], output/'nominal', through='strength',
            reference=reference, progress=progress)
        summary['nominal_file'] = 'nominal/workflow.json'
    from .generated_halfspace import (prepare_generated_halfspace,
        build_generated_halfspace_patches, solve_generated_halfspace_contact)
    manufacturing = loaded['base_case']['manufacturing']
    emit('prepare_actual_solids', completed=False)
    prepared = prepare_generated_halfspace(manufacturing['pinion']['request'],
        manufacturing['wheel']['request'], cad_models=cad_models,
        pinion_blank_settings=manufacturing['pinion']['blank'],
        wheel_blank_settings=manufacturing['wheel']['blank'],
        **loaded['base_case'].get('cad_sampling', {}))
    _write(prepared.diagnostics, output/'prepared_geometry.json')
    summary['prepared_geometry_file'] = 'prepared_geometry.json'
    emit('prepare_actual_solids', completed=True)
    accepted_results = []
    for index in selected:
        tick = time.perf_counter()
        emit('surface_state', phase_index=index, driving_angle_rad=phases[index], completed=False)
        state = dict(schema='bevel_surface_state_v1', phase_index=index,
            driving_angle_rad=phases[index], input_fingerprint=effective_fingerprint,
            source_identity_before=identity, accepted=False)
        geometry = None
        try:
            geometry = build_generated_halfspace_patches(prepared, phases[index], torque, **analysis)
            result = solve_generated_halfspace_contact(geometry, prepared=prepared, **materials, **solver)
            state['result'] = result
        except (ArithmeticError, ValueError, RuntimeError) as error:
            result = dict(converged=False, status='surface_computation_failed', fields=[], geometry=geometry)
            state.update(result=result, failure=dict(type=type(error).__name__, message=str(error)))
        after = _source_identity(nominal=nominal, plots=plots)
        unchanged = after == identity
        state.update(source_identity_after=after, source_identity_unchanged=unchanged,
            seconds=time.perf_counter()-tick,
            accepted=_accepted(result, unchanged, maximum_residual))
        relative = f'states/phase_{index:04d}.json'
        _write(state, output/relative)
        row = dict(phase_index=index, driving_angle_rad=phases[index], state_file=relative,
            accepted=state['accepted'], converged=bool(result.get('converged')),
            status=result.get('status'), seconds=state['seconds'],
            local_halfspace_closure_rad=result.get('theta_rad') if state['accepted'] else None,
            local_halfspace_stiffness_nm_per_rad=result.get('incremental_stiffness_nm_per_rad') if state['accepted'] else None,
            maximum_pressure_pa=result.get('maximum_pressure_pa') if state['accepted'] else None,
            active_cell_count=result.get('active_cell_count') if state['accepted'] else None,
            pinion_torque_relative_residual=result.get('pinion_torque_relative_residual'),
            wheel_balance_residual_nm=result.get('balance_residual'))
        summary['states'].append(row)
        summary.update(source_identity_after=after, source_identity_unchanged=unchanged)
        _write(summary, output/'surface_summary.json')
        emit('surface_state', **row, completed=True)
        if state['accepted'] and plots:
            accepted_results.append(result)
        if not unchanged:
            raise RuntimeError('计算过程中公开计算源码发生变化；已保留当前state，不能验收本次业务')
    summary['analysis_completed'] = all(row['accepted'] for row in summary['states'])
    if plots:
        if not accepted_results or not summary['analysis_completed']:
            summary['plots_status'] = 'requires_accepted_states'
        else:
            from .halfspace_plotting import plot_halfspace_contact_pattern
            mesh = _flank_mesh(prepared, **loaded['surface_case'].get('flank_mesh', {}))
            _write(mesh, output/'actual_pinion_working_flank.json')
            summary['actual_flank_mesh_file'] = 'actual_pinion_working_flank.json'
            offset = loaded['surface_case'].get('plot_tooth_offset', 0)
            for language in ('zh', 'en'):
                plot_halfspace_contact_pattern(accepted_results, output/'figures',
                    pinion_request=manufacturing['pinion']['request'], flank_mesh=mesh,
                    tooth_offset=offset, lang=language, **loaded['surface_case'].get('plot', {}))
            summary.update(plots_status='completed',
                figure_metadata_files=['figures/zh/10_halfspace_contact_pattern.json',
                    'figures/en/10_halfspace_contact_pattern.json'])
    after = _source_identity(nominal=nominal, plots=plots)
    summary.update(source_identity_after=after, source_identity_unchanged=after == identity)
    summary['analysis_completed'] = summary['analysis_completed'] and summary['source_identity_unchanged']
    _write(summary, output/'surface_summary.json')
    return summary
