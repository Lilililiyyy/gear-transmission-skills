"""直齿/弧齿同一案例的几何→实体→强度→TCA→LTCA业务链。

配置使用SI；机床、材料、设计限值和孔壁支承显式输入。标准评级使用
名义几何，制造绑定另列实际两侧差异；迭代成功、设计限值满足、网格
收敛和外部对标分别报告，不以单一passed覆盖这些不同结论。
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sys

from .gear_pair import StraightBevelPair, SpiralBevelPair
from .forces import calculate_mesh_forces
from .iso2014_nominal import NominalMaterial
from .strength_assessment import assess_nominal_strength
from .io import write_json


STAGES = ('geometry', 'strength', 'cad', 'tca', 'ltca')


def _calculation_identity():
    """只记录参与业务计算的公开源码相对名与哈希，不输出本机路径。"""
    names = ('geometry','forces','gear_pair','manufacturing_binding','iso2014_nominal','strength_assessment',
        'generated_cutter','generated_surface','generated_solid','cad','cad_import','rigid_contact',
        'straight_cad','straight_tca','generated_tca','elasticity','contact_mesh','contact_equilibrium',
        'straight_ltca','generated_ltca','workflow')
    folder = Path(__file__).resolve().parent
    return {name:hashlib.sha256((folder/(name+'.py')).read_bytes()).hexdigest() for name in names}


def load_bevel_case(path):
    """读取案例，将相对于配置文件的制造请求展开，返回可独立重放的字典。"""
    path = Path(path).resolve()
    case = json.loads(path.read_text(encoding='utf-8-sig'))
    case = deepcopy(case)
    for member in case.get('manufacturing', {}).values():
        if 'request_file' in member:
            if 'request' in member:
                raise ValueError('制造输入不能同时包含request和request_file')
            member['request'] = json.loads((path.parent / member.pop('request_file')).read_text(encoding='utf-8-sig'))
    return case


def case_fingerprint(case):
    """绑定全部数学输入；名称不参与数值身份，未忽略材料/载荷/制造条件。"""
    values = {key: value for key, value in case.items() if key not in ('name', 'description')}
    return hashlib.sha256(json.dumps(values, sort_keys=True, separators=(',', ':'), ensure_ascii=False,
                                    allow_nan=False).encode('utf-8')).hexdigest()


def reference_input_fingerprint(case):
    """几何/强度对标的实际数学输入身份；FE网格变化不改变标准评级工况。"""
    keys = ('kind', 'geometry', 'materials', 'load', 'strength', 'force_convention')
    return case_fingerprint({key: case.get(key) for key in keys})


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(name + '须为正有限实数')
    return float(value)


def _selected_tca(results, sign):
    """选择有符号加载方向最先闭合的有效齿对，始终使用同一整轮轴相位。"""
    candidates = [value for value in results if value['converged']]
    if not candidates:
        return None
    key = 'follower_increment_angle_rad' if 'follower_increment_angle_rad' in candidates[0] else 'follower_angle_rad'
    return min(candidates, key=lambda value: sign * value[key])


def _straight_tca_record(value):
    """为统一后处理补等义字段；保留直齿原始法向误差的定义。"""
    return {**value, 'tooth_offset':value['frames']['tooth_offset'],
        'follower_increment_angle_rad':value['follower_angle_rad'],
        'unloaded_te_rad':value['unloaded_transmission_error_rad'],
        'pinion_contact_point_world_m':value['mean_contact_point_world_m'],
        'normal_alignment_residual':value['opposing_normal_residual'],
        'normal_alignment_definition':'norm(n1+n2)，两轮外法向反向残差'}


def _assembly(models, frames):
    """将实际材料实体按相切姿态装配；CAD内部mm，公开位移m。"""
    from OCP.TopoDS import TopoDS_Compound
    from OCP.BRep import BRep_Builder
    from .rigid_contact import RigidTransform, _transform
    from .cad_import import ImportedCadModel
    compound = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(compound)
    for member in ('pinion', 'wheel'):
        pose = RigidTransform(frames[member + '_rotation_rows'], frames.get('translation_m', (0., 0., 0.)))
        builder.Add(compound, _transform(models[member].shape, pose, .001))
    return ImportedCadModel(compound, {'units': 'm', 'cad_internal_length_unit': 'mm', 'valid': True,
        'pose': frames, 'interpretation': '真实闭实体在无载相切姿态的装配；未作全周期干涉认证'})


def _surface_subset(mesh, face_ids):
    """按真实BRep面号选三角网格，不拉伸、不补面或改变原始顶点位置。"""
    selected = [index for index, face in enumerate(mesh['face_indices']) if face in face_ids]
    if not selected:
        raise ValueError('指定实际工作面/根区没有三角化结果')
    vertices = sorted({vertex for index in selected for vertex in mesh['triangles'][index]})
    mapping = {value:index for index,value in enumerate(vertices)}
    return dict(units='m', vertices_m=[mesh['vertices_m'][index] for index in vertices],
        triangles=[[mapping[vertex] for vertex in mesh['triangles'][index]] for index in selected],
        face_indices=[mesh['face_indices'][index] for index in selected],
        metadata={'closed_solid':False,'interpretation':'真实裁切材料表面的选区；不是另建完整材料体'})


def compare_workflow_reference(result, reference):
    """对外部已采集参考逐量比较；工况身份不匹配时不判定对标通过。"""
    def lookup(value, path):
        for key in path.split('.'):
            value = value[int(key)] if isinstance(value, (list, tuple)) else value[key]
        return value
    keys = ('kind', 'geometry', 'materials', 'load', 'strength', 'force_convention')
    actual_inputs = {key: result['case'].get(key) for key in keys}
    saved_inputs = reference.get('reference_inputs')
    identity = (saved_inputs == actual_inputs and
                reference.get('input_fingerprint') == case_fingerprint(saved_inputs) == result['reference_input_fingerprint'])
    groups = {}
    for group in ('geometry', 'strength'):
        rows = []
        for item in reference.get(group + '_rows', []):
            for key in ('absolute_tolerance','relative_tolerance'):
                tolerance_value = item.get(key, 0.)
                if isinstance(tolerance_value,bool) or not isinstance(tolerance_value,(int,float)) or not math.isfinite(tolerance_value) or tolerance_value < 0:
                    raise ValueError('逐量比较容差须为非负有限实数，不接受无限容差或布尔值')
            if isinstance(item['reference'],bool):
                raise ValueError('数值参考不能用布尔值代替')
            actual = float(lookup(result, item['result_path']))
            expected = float(item['reference'])
            if not math.isfinite(actual) or not math.isfinite(expected):
                raise ValueError('对标量须为有限数；缺失参照不能写成0')
            tolerance = max(float(item.get('absolute_tolerance', 0.)),
                            abs(expected) * float(item.get('relative_tolerance', 0.)))
            difference = actual - expected
            rows.append({**item, 'calculated': actual, 'signed_error': difference, 'absolute_error': abs(difference),
                'relative_error': None if expected == 0 else difference / abs(expected),
                'passed': abs(difference) <= tolerance})
        groups[group + '_rows'] = rows
    rows = groups['geometry_rows'] + groups['strength_rows']
    return dict(**groups, same_case_inputs=identity, compared_fields=len(rows),
        matched_fields=sum(row['passed'] for row in rows),
        comparison_passed=bool(identity and rows and all(row['passed'] for row in rows)),
        ltca_reference_used=False, scope='同工况名义几何与标准应力；LTCA为独立计算结果')


def analyze_bevel_case(case, directory, *, through='ltca', include_fields=False, reference=None,
                       render3d=False, progress=None):
    """运行同一案例到指定阶段并导出输入与结果；支持直齿和显式制造弧齿。

    through可取geometry/strength/cad/tca/ltca。负设计安全裕度是合法结果；
    求解失败单独保留status并令analysis_completed=False，不以空数组当成功。
    绘图由独立后处理入口调用，计算函数不依赖Matplotlib。
    """
    if through not in STAGES:
        raise ValueError('业务阶段须为geometry、strength、cad、tca或ltca')
    if not isinstance(include_fields, bool) or not isinstance(render3d, bool):
        raise ValueError('输出/渲染开关须为布尔值')
    case = deepcopy(case)
    if case.get('schema') != 'bevel_analysis_case_v1':
        raise ValueError('案例schema须为bevel_analysis_case_v1')
    kind = case.get('kind')
    if kind not in ('straight', 'spiral'):
        raise ValueError('对象kind仅支持straight或spiral')
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    write_json(case, directory / 'inputs.json')
    load = case['load']
    torque = _positive(load['pinion_torque_n_m'], '主动扭矩幅值Nm')
    speed = _positive(load['pinion_speed_rpm'], '主动转速rpm')
    sign = load['ltca_torque_sign']
    if isinstance(sign, bool) or sign not in (-1, 1):
        raise ValueError('承载参数轴的扭矩符号须为-1或+1')
    if kind == 'straight' and sign != 1:
        raise ValueError('本版理想直齿承载入口仅开放right工作侧的正向参数轴扭矩')
    config = {'geometry': case['geometry'], 'load': {'pinion_torque_n_m': torque, 'pinion_speed_rpm': speed}}
    pair = (StraightBevelPair(config) if kind == 'straight' else
            SpiralBevelPair(config, binding_tolerances=case.get('binding_tolerances')))
    geometry = pair.geometry()
    result = dict(schema='bevel_workflow_v1', case=case, case_fingerprint=case_fingerprint(case),
        kind=kind, completed_stages=[], analysis_completed=False, geometry=asdict(geometry),
        reference_input_fingerprint=reference_input_fingerprint(case),
        forces=asdict(calculate_mesh_forces(geometry, torque, **case.get('force_convention', {}))),
        standard_rating=None, strength_assessment=None, manufacturing_binding=None,
        tca={'rows': [], 'closure_direction': sign}, ltca={'results': [], 'load_sweep': []}, cad={},
        quality={'discretization_converged': None, 'physical_validation': False})
    from . import __version__
    from importlib.metadata import version, PackageNotFoundError
    versions = {}
    for distribution in ('numpy','scipy','gmsh','cadquery-ocp','vtk','matplotlib'):
        try:versions[distribution]=version(distribution)
        except PackageNotFoundError:versions[distribution]=None
    result['runtime'] = {'package_version':__version__, 'python_version':sys.version.split()[0], 'dependencies':versions}
    result['source_identity_before'] = _calculation_identity()
    target = STAGES.index(through)

    def finished(stage):
        result['completed_stages'].append(stage)
        result['analysis_completed'] = STAGES.index(stage) >= target
        result['source_identity_after'] = _calculation_identity()
        result['source_identity_unchanged'] = result['source_identity_before'] == result['source_identity_after']
        if not result['source_identity_unchanged']:
            result['analysis_completed'] = False
            write_json(result, directory / 'workflow.json')
            raise RuntimeError('计算期间参与计算的源码发生变化，本次结果保留但不判定可复现业务完成')
        write_json(result, directory / 'workflow.json')
        if progress:
            progress({'stage': stage, 'completed': True})
        return result['analysis_completed']

    if kind == 'spiral':
        manufacturing = case['manufacturing']
        requests = {m: manufacturing[m]['request'] for m in ('pinion', 'wheel')}
        result['manufacturing_binding'] = pair.manufacturing_binding(requests['pinion'], requests['wheel'])
    if finished('geometry'):
        return result
    materials = {m: NominalMaterial(**case['materials'][m]) for m in ('pinion', 'wheel')}
    result['standard_rating'] = pair.nominal_strength_2014(**materials, **case['strength'])
    result['strength_assessment'] = assess_nominal_strength(result['standard_rating'], case['stress_limits'])
    result['quality']['design_meets_stress_targets'] = result['strength_assessment']['meets_all_targets']
    if reference is not None:
        result['masta_benchmark'] = compare_workflow_reference(result, reference)
    if finished('strength'):
        return result

    if kind == 'spiral':
        models = pair.build_generated_models(requests['pinion'], requests['wheel'],
            manufacturing['pinion']['blank'], manufacturing['wheel']['blank'], **case.get('cad_sampling', {}))
        prepared_tca = pair.prepare_tca(requests['pinion'], requests['wheel'], cad_models=models)
        solver_tca = lambda angle, offset: pair.tca(angle, prepared=prepared_tca, tooth_offset=offset)
    else:
        models = pair.build_models(case['cad']['pinion'], case['cad']['wheel'], full_gears=True)
        solver_tca = lambda angle, offset: _straight_tca_record(pair.tca(angle, tooth_offset=offset))
    from .cad_import import ImportedCadModel, export_model_step, export_model_brep, triangulate_imported_model
    cad_dir = directory / 'cad'
    for member, model in models.items():
        imported = ImportedCadModel(model.shape, model.diagnostics)
        export_model_step(imported, cad_dir / (member + '.step'))
        export_model_brep(imported, cad_dir / (member + '.brep'))
        result['cad'][member] = dict(step='cad/' + member + '.step', brep='cad/' + member + '.brep',
                                   diagnostics=model.diagnostics)
        mesh = triangulate_imported_model(imported, linear_deflection_m=1e-4)
        write_json(mesh, cad_dir / (member + '_mesh.json'))
        result['cad'][member + '_mesh'] = 'cad/' + member + '_mesh.json'
        if kind == 'spiral' and member == 'pinion':
            groups = model.slot_face_groups[0]
            flank = _surface_subset(mesh, {face for name in ('left_flank','right_flank') for face in groups[name]})
            root = _surface_subset(mesh, {face for name in ('left_fillet','root_bottom','right_fillet') for face in groups[name]})
            for name, detail in (('flank',flank),('root',root)):
                write_json(detail, cad_dir / (name + '_mesh.json'))
                result['cad'][name + '_mesh'] = 'cad/' + name + '_mesh.json'
        elif kind == 'straight' and member == 'pinion':
            groups = model.face_groups
            flank_ids = {face for name, ids in groups.items() if name.startswith('involute_') for face in ids}
            root_ids = {face for name, ids in groups.items() if 'root' in name for face in ids}
            for name, ids in (('flank',flank_ids),('root',root_ids)):
                write_json(_surface_subset(mesh, ids), cad_dir / (name + '_mesh.json'))
                result['cad'][name + '_mesh'] = 'cad/' + name + '_mesh.json'
    zero_candidates = [solver_tca(0., k) for k in (-1, 0, 1)]
    zero = _selected_tca(zero_candidates, sign)
    if zero is None:
        raise ArithmeticError('零主动角没有有效相切姿态，不能按理想齿数比补画装配')
    assembly = _assembly(models, zero['frames'])
    export_model_step(assembly, cad_dir / 'pair.step')
    result['cad']['pair'] = dict(step='cad/pair.step', frames=zero['frames'], driving_angle_rad=0.)
    pair_mesh = triangulate_imported_model(assembly, linear_deflection_m=1e-4)
    write_json(pair_mesh, cad_dir / 'pair_mesh.json')
    result['cad']['pair_mesh'] = 'cad/pair_mesh.json'
    if render3d:
        from .cad_preview import render_cad_mesh, export_vtp
        export_vtp(pair_mesh, cad_dir / 'pair.vtp')
        for lang in ('zh', 'en'):
            render_cad_mesh(pair_mesh, cad_dir / ('pair_' + lang + '.png'), lang=lang)
    if finished('cad'):
        return result

    angles = list(case['analysis']['driver_angles_rad'])
    offsets = tuple(case['analysis'].get('tooth_offsets', (-1, 0, 1)))
    if not angles or any(isinstance(a, bool) or not isinstance(a, (int, float)) or not math.isfinite(a) for a in angles):
        raise ValueError('主动角采样须为非空有限弧度数列')
    if len(set(angles)) != len(angles) or any(b <= a for a, b in zip(angles, angles[1:])):
        raise ValueError('主动角采样须严格递增且不得重复')
    for angle in angles:
        candidates = zero_candidates if angle == 0. and offsets == (-1, 0, 1) else [solver_tca(angle, k) for k in offsets]
        selected = _selected_tca(candidates, sign)
        result['tca']['rows'].append(dict(driving_angle_rad=angle, results=candidates, selected=selected,
                                         converged=selected is not None))
    result['tca']['all_positions_tangent'] = all(row['converged'] for row in result['tca']['rows'])
    write_json(result['tca'], directory / 'tca.json')
    if finished('tca'):
        result['analysis_completed'] = result['tca']['all_positions_tangent']
        write_json(result, directory / 'workflow.json')
        return result

    elastic = {m: {key: case['materials'][m][key] for key in ('young_modulus_pa', 'poisson_ratio')} for m in ('pinion', 'wheel')}
    if kind == 'spiral':
        prepared = pair.prepare_ltca(requests['pinion'], requests['wheel'],
            pinion_blank_settings=manufacturing['pinion']['blank'], wheel_blank_settings=manufacturing['wheel']['blank'],
            pinion_material=elastic['pinion'], wheel_material=elastic['wheel'], cad_models=models,
            directory=directory / 'finite_element', **case['finite_element'])
    else:
        prepared = pair.prepare_ltca(pinion_cad_settings=case['cad']['pinion'], wheel_cad_settings=case['cad']['wheel'],
            pinion_material=elastic['pinion'], wheel_material=elastic['wheel'],
            directory=directory / 'finite_element', **case['finite_element'])
    result['ltca']['preparation'] = prepared.diagnostics
    contact_options = dict(case['loaded_contact'])
    contact_options.setdefault('tooth_offsets', offsets)
    def solve(angle, magnitude, label):
        try:
            value = pair.solve_ltca(prepared, angle, sign*magnitude, include_fields=include_fields, **contact_options)
        except (ValueError, ArithmeticError) as error:
            value = dict(converged=False, status='failed', error=str(error), driving_angle_rad=angle,
                         pinion_torque_nm=sign*magnitude, contacts=[])
        write_json(value, directory / 'ltca' / (label + '.json'))
        if progress:
            progress({'stage': 'ltca', 'label': label, 'angle_rad': angle, 'torque_nm': sign*magnitude,
                      'converged': value['converged'], 'status': value['status']})
        return value
    zero_loaded = solve(0., torque, 'reference')
    magnitudes = case['analysis'].get('load_sweep_torques_nm', [torque])
    for index, magnitude in enumerate(magnitudes):
        _positive(magnitude, '载荷扫描扭矩幅值Nm')
        value = zero_loaded if magnitude == torque else solve(0., magnitude, f'load_{index:02}')
        result['ltca']['load_sweep'].append(value)
    prepared._contact_cache.clear()
    for index, angle in enumerate(angles):
        value = zero_loaded if angle == 0. else solve(angle, torque, f'phase_{index:02}')
        result['ltca']['results'].append(value)
        prepared._contact_cache.clear()
    all_loaded = result['ltca']['results'] + result['ltca']['load_sweep']
    result['ltca']['all_positions_converged'] = bool(all_loaded and all(row['converged'] for row in all_loaded))
    finished('ltca')
    result['analysis_completed'] = result['ltca']['all_positions_converged'] and result['tca']['all_positions_tangent']
    write_json(result, directory / 'workflow.json')
    return result
