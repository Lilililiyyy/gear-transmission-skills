"""新写锥齿轮Skill适配：调用原接口、保护输出并实际重算；不重写计算公式。"""
from copy import deepcopy
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path
import platform


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def identity():
    """逐文件核对实际导入内核；代码变化不能沿用旧快照身份。"""
    import bevel_gear
    manifest = read_json(Path(__file__).with_name('bevel_snapshot.json'))
    package = Path(bevel_gear.__file__).resolve().parent
    for relative, expected in manifest['files'].items():
        if relative.startswith('bevel_gear/'):
            path = package / relative.removeprefix('bevel_gear/')
            if not path.is_file() or digest(path) != expected:
                raise ValueError('锥齿轮内核与已登记快照不同：'+relative)
    versions = {'python': platform.python_version(), 'kernel': bevel_gear.__version__}
    for name in ('gear-transmission-skills', 'numpy', 'scipy', 'matplotlib', 'cadquery-ocp', 'gmsh', 'vtk', 'shapely'):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {'source_commit': manifest['source_commit'], 'versions': versions,
            'snapshot_sha256': digest(Path(__file__).with_name('bevel_snapshot.json')),
            'adapter_sha256': digest(__file__)}


def prepare_case(case, through):
    """检查阶段相关的必需输入与原生几何类型；不执行接触数值求解。"""
    from bevel_gear.workflow import STAGES
    from bevel_gear.gear_pair import StraightBevelPair, SpiralBevelPair
    from bevel_gear.forces import calculate_mesh_forces
    if through not in STAGES:
        raise ValueError('阶段须为geometry、strength、cad、tca或ltca')
    if case.get('schema') != 'bevel_analysis_case_v1' or case.get('kind') not in ('straight', 'spiral'):
        raise ValueError('只接受bevel_analysis_case_v1的straight或spiral案例')
    if case.get('geometry', {}).get('offset_m', 0.) != 0.:
        raise ValueError('当前仅支持零偏置相交轴锥齿轮，不能以此分析准双曲面齿轮')
    load = case['load']
    for name in ('pinion_torque_n_m', 'pinion_speed_rpm'):
        value = load[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(name+'须为正有限幅值，单位由字段规定')
    sign = load['ltca_torque_sign']
    if isinstance(sign, bool) or sign not in (-1, 1):
        raise ValueError('ltca_torque_sign须显式为-1或+1')
    if case['kind'] == 'straight' and sign != 1:
        raise ValueError('本版直齿统一入口仅支持right侧正向参数轴扭矩')
    config = {'geometry': case['geometry'], 'load': load}
    if case['kind'] == 'straight':
        pair = StraightBevelPair(config)
    else:
        pair = SpiralBevelPair(config, binding_tolerances=case.get('binding_tolerances'))
        requests = case['manufacturing']
        binding = pair.manufacturing_binding(requests['pinion']['request'], requests['wheel']['request'])
        if not binding['accepted']:
            fields = sorted({row['field'] for row in binding['failures']})
            raise ValueError('名义几何与实际制造绑定未接受：'+', '.join(fields))
    calculate_mesh_forces(pair.geometry(), load['pinion_torque_n_m'], **case.get('force_convention', {}))
    stage = STAGES.index(through)
    if stage >= 1:
        from bevel_gear.iso2014_nominal import NominalMaterial
        from bevel_gear.strength_assessment import assess_nominal_strength
        materials = {member: NominalMaterial(**case['materials'][member]) for member in ('pinion', 'wheel')}
        # 标准评级为代数校核；调用原函数验证系数与限值，不将其视作接触求解。
        rating = pair.nominal_strength_2014(**materials, **case['strength'])
        assess_nominal_strength(rating, case['stress_limits'])
    if stage >= 2:
        key = 'cad' if case['kind'] == 'straight' else 'manufacturing'
        for member in ('pinion', 'wheel'):
            if member not in case[key]:
                raise ValueError('缺少'+key+'.'+member)
    if stage >= 3:
        angles = case['analysis']['driver_angles_rad']
        if not angles or any(isinstance(a, bool) or not isinstance(a, (int, float)) or not math.isfinite(a) for a in angles):
            raise ValueError('主动角须为非空有限弧度数组')
        if any(b <= a for a, b in zip(angles, angles[1:])):
            raise ValueError('主动角须严格递增')
    if stage >= 4:
        for key in ('finite_element', 'loaded_contact'):
            if not isinstance(case.get(key), dict):
                raise ValueError('缺少'+key+'设置')
    return deepcopy(case)


def _new_directory(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    return path


def _files(directory):
    return {p.relative_to(directory).as_posix(): digest(p) for p in sorted(directory.rglob('*')) if p.is_file() and p.name != 'run.json'}


def _report(result, through):
    assessment = result.get('strength_assessment')
    geometry = result['geometry']
    lines = [
        '# 锥齿轮设计分析结果', '',
        '实现性质：新写Skill适配；调用原字节保留的公开重建计算内核。', '',
        f"对象：{result['kind']}；请求阶段：{through}；实际阶段：{', '.join(result['completed_stages'])}。",
        f"业务完成：{result['analysis_completed']}；设计应力目标：{None if assessment is None else assessment['meets_all_targets']}。", '',
        '| 主要量 | 数值 |', '|---|---:|',
        f"| 齿数比 z₂/z₁ | {geometry['ratio']:.8g} |",
        f"| 外端端面模数 / mm | {geometry['outer_transverse_module_m']*1000:.8g} |",
        f"| 中点法向模数 / mm | {geometry['mean_normal_module_m']*1000:.8g} |",
        f"| 小轮切向力 / N | {result['forces']['pinion_tangential_n']:.8g} |",
        f"| 小轮径向力 / N | {result['forces']['pinion_radial_n']:.8g} |",
        f"| 小轮轴向力 / N | {result['forces']['pinion_axial_n']:.8g} |", '']
    if assessment is not None:
        labels = {'H': '接触', 'F1': '小轮齿根', 'F2': '大轮齿根'}
        lines += ['| 指标 | 应力 / MPa | 安全系数 | 满足给定目标 |', '|---|---:|---:|---|']
        for row in assessment['rows']:
            lines.append(f"| {labels[row['quantity']]} | {row['calculated_stress_pa']/1e6:.8g} | {row['safety_factor']:.8g} | {row['meets_target']} |")
        lines += ['', '设计限值依据：'+str(assessment['limits']['basis']), '']
    lines += ['输入说明：'+str(result['case'].get('description', '未另附说明')), '',
        '用户声明的假设：'+json.dumps(result['case'].get('assumptions', []), ensure_ascii=False), '',
        '实际采用的完整几何参数（含内核默认值；SI）：', '', '```json',
        json.dumps(geometry['inputs'], ensure_ascii=False, indent=2), '```', '',
        '原输入见 inputs.json；实际计算与具名SI数据见 calculation/workflow.json。',
        '名义应力与LTCA空间压力分别解释；阶段完成不代表设计合格或离散收敛。',
        '统一LTCA使用原节点FE路线；未自动切换为原生分布牵引FE或半解析模型。',
        '数值复现须另行执行replay；当前未据此宣称网格／相位收敛、实验验证或整机寿命。', '']
    return '\n'.join(lines)


def run_case(case, output, *, through='strength', plots=False, reference=None):
    """运行所选阶段；所有结果和失败信息保留在新目录，LTCA保存原场数组。"""
    from bevel_gear.workflow import analyze_bevel_case
    case = prepare_case(case, through)
    before = identity()
    output = _new_directory(output)
    run = dict(schema=1, workflow='bevel-gear-design', status='running', through=through,
               include_fields=through == 'ltca', identity=before, plots=plots)
    write_json(output/'inputs.json', case)
    if reference is not None:
        write_json(output/'reference.json', reference)
    write_json(output/'run.json', run)
    try:
        result = analyze_bevel_case(case, output/'calculation', through=through,
            include_fields=run['include_fields'], reference=reference)
        if not result['analysis_completed']:
            raise ArithmeticError('所选阶段未全部完成；请读取workflow.json中的逐相位状态与残差')
        if identity() != before:
            raise ValueError('运行期间源码或数值环境发生变化')
        if plots:
            from bevel_gear.workflow_plotting import plot_spiral_workflow
            figures = plot_spiral_workflow(result, output/'figures', lang='zh', data_root=output/'calculation')
            if not figures['complete']:
                raise ArithmeticError('部分图表未完成，详见figures/zh/manifest.json')
        (output/'report.md').write_text(_report(result, through), encoding='utf-8')
        run.update(status='completed', analysis_completed=True)
    except Exception as error:
        run.update(status='failed', error_type=type(error).__name__, message=str(error))
        raise
    finally:
        run['files'] = _files(output)
        write_json(output/'run.json', run)
    return run


def compare_values(old, new, *, rtol=1e-8, atol=1e-12, path='', differences=None):
    """递归核对结构、状态、完整数值列表与非数值量；未知值不能补零。"""
    if differences is None:
        differences = []
    if isinstance(old, bool) or isinstance(new, bool):
        same = type(old) is type(new) and old == new
    elif isinstance(old, (int, float)) and isinstance(new, (int, float)):
        same = (old == new) or (math.isfinite(old) and math.isfinite(new) and math.isclose(old, new, rel_tol=rtol, abs_tol=atol))
    elif isinstance(old, dict) and isinstance(new, dict):
        same = old.keys() == new.keys()
        if same:
            for key in old:
                compare_values(old[key], new[key], rtol=rtol, atol=atol, path=path+'.'+key, differences=differences)
            return differences
    elif isinstance(old, list) and isinstance(new, list):
        same = len(old) == len(new)
        if same:
            for index, (a, b) in enumerate(zip(old, new)):
                compare_values(a, b, rtol=rtol, atol=atol, path=f'{path}[{index}]', differences=differences)
            return differences
    else:
        same = type(old) is type(new) and old == new
    if not same:
        differences.append(path)
    return differences


def replay_run(original, output, *, plots=False):
    """先验原文件、源码和环境，再独立求解并比较数学结果；不重放冻结数组。"""
    original = Path(original).resolve()
    saved = read_json(original/'run.json')
    if saved.get('workflow') != 'bevel-gear-design' or saved.get('status') != 'completed':
        raise ValueError('仅可复算已完成的锥齿轮Skill结果')
    for relative, expected in saved['files'].items():
        path = (original/relative).resolve()
        if not path.is_relative_to(original) or not path.is_file() or digest(path) != expected:
            raise ValueError('原结果缺失、越界或SHA256不符：'+relative)
    if identity() != saved['identity']:
        raise ValueError('当前源码或数值环境与原结果不同，不能登记为同环境复算')
    reference = read_json(original/'reference.json') if (original/'reference.json').exists() else None
    run_case(read_json(original/'inputs.json'), output, through=saved['through'], plots=plots, reference=reference)
    old = read_json(original/'calculation/workflow.json')
    new = read_json(Path(output)/'calculation/workflow.json')
    keys = ('case', 'case_fingerprint', 'kind', 'completed_stages', 'analysis_completed', 'geometry',
            'forces', 'standard_rating', 'strength_assessment', 'manufacturing_binding', 'tca', 'ltca', 'quality')
    keys = keys + (('masta_benchmark',) if 'masta_benchmark' in old else ())
    differences = compare_values({k: old[k] for k in keys}, {k: new[k] for k in keys})
    result = dict(passed=not differences, compared_groups=list(keys), differences=differences,
                  rtol=1e-8, atol=1e-12, scope='重新求解后的完整数学结果；不比较STEP等CAD文件的序列化字节或图片')
    write_json(Path(output)/'replay.json', result)
    # 复算报告也登记到文件清单，保证后续复算仍可完整校验。
    run = read_json(Path(output)/'run.json')
    run['files'] = _files(Path(output))
    write_json(Path(output)/'run.json', run)
    return result
