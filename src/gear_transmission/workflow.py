"""新写工作流：调用保留内核、保存证据、生成中文报告并独立复算全部原数组。"""
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import platform

import numpy as np
import scipy
import gear_contact
from gear_contact import analyze_pair, compare_designs, export_result
from gear_contact.io import write_json

from . import __version__
from .protocol import prepare_request, read_json


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity():
    """验证实际加载的内核文件，防止同名安装包或后续改动被误认作本快照。"""
    here = Path(__file__).resolve().parent
    manifest = read_json(here / 'kernel_snapshot.json')
    kernel = Path(gear_contact.__file__).resolve().parent
    expected = {k.removeprefix('gear_contact/'): v for k, v in manifest['files'].items()
                if k.startswith('gear_contact/') and k.endswith('.py')}
    actual = {p.relative_to(kernel).as_posix(): sha256(p) for p in sorted(kernel.rglob('*.py'))}
    if actual != expected:
        raise RuntimeError('实际加载的gear_contact与登记快照不一致，请重新安装本项目并核对内核来源')
    return {'kernel_files': actual, 'snapshot_sha256': sha256(here / 'kernel_snapshot.json'),
            'workflow_files': {p.name: sha256(p) for p in sorted(here.glob('*.py'))},
            'versions': {'workflow': __version__, 'kernel': gear_contact.__version__,
                         'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}}


def artifact_hashes(directory):
    directory = Path(directory)
    return {p.relative_to(directory).as_posix(): sha256(p)
            for p in sorted(directory.rglob('*')) if p.is_file() and p != directory / 'run.json'}


def _label(value):
    return str(value).replace('|', '／').replace('\n', ' ').replace('\r', ' ')


def _report(directory, manifest, request, study):
    rows = ['# 圆柱齿轮承载接触分析', '',
            '本报告由新写Skill适配流程调用保留的重建计算内核生成。内部采用m、N、Pa、rad和N·m。', '',
            '|方案|工况|扭矩 N·m|TE峰峰值 μm|峰压 MPa|平均法向刚度 N/μm|收敛|',
            '|---|---|---:|---:|---:|---:|---|']
    for record in manifest['cases']:
        s = record['summary']
        rows.append(f"|{_label(record['design'])}|{_label(record['operating_point'])}|{s['torque_n_m']:.8g}|"
                    f"{s['te_peak_to_peak_m']*1e6:.8g}|{s['max_pressure_pa']*1e-6:.8g}|"
                    f"{s['mean_mesh_stiffness_n_per_m']*1e-6:.8g}|{'是' if s['all_converged'] else '否'}|")
    if study:
        rows += ['', '给定候选与工况集合中的非支配方案：' + '、'.join(map(_label, study.pareto_names)) + '。',
                 '目标为最坏工况TE峰峰值与最大Hertz压力；没有隐含综合权重，不宣称连续全局最优。']
    rows += ['', '## 输入与默认值', '',
             '原请求见request.json，换算后的SI请求见normalized_request.json；每个cases子目录的input.json记录实际求解输入及展开的默认值。', '',
             '本次单位换算：' + str(len(manifest['unit_conversions'])) + '项，具体路径和倍率保存在run.json。', '',
             '内核默认字段：' + ('、'.join(manifest['defaulted_fields']) or '无') + '。']
    rows += ['', '用户提供或明确接受的假设：'] + ['- ' + _label(a) for a in request.get('assumptions', [])]
    rows += ['', '## 结果含义与验证范围', '',
             '- 本次求解完成且各相位达到内核退出条件；这不等于离散收敛或工程精度通过。',
             '- 当前采用' + _label(manifest['cases'][0]['summary']['model_scope']) + '。',
             '- TE为法向接近量，角度TE另存rad；法向刚度为N/m，不能当作扭转刚度。',
             '- 压力为局部Hertz线接触结果；瞬时接触图固定相位与齿对，周期包络不能积分为同时承载面积。',
             '- 本次未执行空间／相位收敛研究、原组件对标或物理试验验证；不输出标准强度或疲劳寿命结论。',
             '- 标准评级需要另外提供对应标准因子及材料限值，不能以LTCA峰压替代标准接触应力。',
             '- 原数组NPZ不含pickle，cell_offsets保留每个相位的变长单元索引；未覆盖区域不能补作确定零压。', '',
             '## 复算', '',
             '运行`gear-design replay 本结果目录 --output 新目录`。程序先验证保存文件和代码身份，再重新求解，',
             '逐项比较全部NPZ数组、活动掩码和汇总指标；replay.json只说明数值复现，不证明物理精度。', '']
    (directory / 'report.md').write_text('\n'.join(rows), encoding='utf-8')


def run_request(request, output, *, plots=False):
    """执行任意已支持参数的请求；输出目录须为新目录，异常保留失败状态。"""
    prepared = prepare_request(request)
    before = identity()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / 'request.json', request)
    write_json(output / 'normalized_request.json', prepared['normalized'])
    manifest = {'schema': 'gear_ltca_run_v1', 'status': 'running', 'analysis_completed': False,
                'created_utc': datetime.now(timezone.utc).isoformat(), 'identity': before, 'plots': plots,
                'unit_conversions': prepared['conversions'], 'defaulted_fields': prepared['defaults'],
                'validation': {'numerical_replay': '尚未复算', 'discretization_converged': None,
                               'vendor_comparison': '本次未执行', 'physical_validation': False}, 'cases': []}
    write_json(output / 'run.json', manifest)
    try:
        study = None
        if prepared['designs'] is None:
            result = analyze_pair(prepared['pair'], **prepared['kwargs'])
            results = [('当前方案', '给定工况', result)]
        else:
            study = compare_designs(prepared['pair'], prepared['designs'],
                    operating_points=prepared['points'], pressure_limit_pa=prepared['pressure_limit_pa'],
                    **prepared['kwargs'])
            write_json(output / 'study.json', study.to_dict())
            results = [(design, point, result) for design, points in study.results.items()
                       for point, result in points.items()]
        for index, (design, point, result) in enumerate(results, 1):
            relative = f'cases/case_{index:03d}'
            export_result(result, output / relative)
            manifest['cases'].append({'directory': relative, 'design': design, 'operating_point': point,
                                      'summary': result.summary})
            if not result.summary['all_converged']:
                raise ArithmeticError(f'第{index}个求解含未收敛相位，结果不能标记完成')
        if plots:
            from gear_contact.plotting import plot_analysis, plot_comparison
            if study is None:
                plot_analysis(results[0][2], output / 'plots')
            else:
                for index, point in enumerate(prepared['points'], 1):
                    plot_comparison({name: values[point.name] for name, values in study.results.items()},
                                    output / f'plots/point_{index:03d}')
        if identity() != before:
            raise RuntimeError('运行期间源码或依赖版本发生变化，本次结果不能标记可复现完成')
        manifest.update(status='completed', analysis_completed=True)
        _report(output, manifest, request, study)
    except Exception as error:
        manifest.update(status='failed', analysis_completed=False,
                        error={'type': type(error).__name__, 'message': str(error)})
        raise
    finally:
        manifest['artifacts'] = artifact_hashes(output)
        write_json(output / 'run.json', manifest)
    return manifest


def _verified_run(directory):
    directory = Path(directory).resolve()
    manifest = read_json(directory / 'run.json')
    if manifest.get('schema') != 'gear_ltca_run_v1' or manifest.get('status') != 'completed':
        raise ValueError('仅允许复算已完成的gear_ltca_run_v1结果')
    if not manifest.get('artifacts') or artifact_hashes(directory) != manifest['artifacts']:
        raise ValueError('保存的结果文件已变化或缺失，不能将其作为有效复算参照')
    if manifest['identity'] != identity():
        raise ValueError('当前源码／数值环境与原运行身份不一致，不能声明同环境复现')
    if not manifest.get('cases'):
        raise ValueError('原运行没有求解记录')
    for case in manifest['cases']:
        relative = Path(case['directory'])
        if relative.is_absolute() or not (directory / relative).resolve().is_relative_to(directory):
            raise ValueError('原求解记录路径超出结果目录')
    return directory, manifest


def _compare_arrays(reference, current):
    """比较解压后的全部数组，忽略NPZ容器时间戳；离散索引与掩码要求精确一致。"""
    rows = []
    with np.load(reference, allow_pickle=False) as a, np.load(current, allow_pickle=False) as b:
        if set(a.files) != set(b.files):
            raise ValueError('复算数组字段集合变化')
        for name in sorted(a.files):
            left, right = a[name], b[name]
            atol = 1e-14
            if name.endswith('_pa'):
                atol = 1e-3
            elif name.endswith('_n_per_m'):
                atol = 1e-3
            elif name.endswith('_n'):
                atol = 1e-9
            exact = left.dtype.kind in 'biu'
            passed = left.shape == right.shape and left.dtype == right.dtype
            passed = bool(passed and (np.array_equal(left, right) if exact else
                                      np.allclose(left, right, rtol=1e-10, atol=atol, equal_nan=False)))
            delta = float(np.max(np.abs(left.astype(float) - right.astype(float)))) if left.shape == right.shape else None
            rows.append({'name': name, 'shape': list(left.shape), 'passed': passed,
                         'maximum_absolute_difference': delta, 'exact_required': exact,
                         'absolute_tolerance': 0 if exact else atol, 'relative_tolerance': 0 if exact else 1e-10})
    return rows


def replay_run(directory, output, *, plots=False):
    """先校验已保存输入与证据，再实际重新调用内核；不读取冻结数组作为求解响应。"""
    directory, reference = _verified_run(directory)
    output = Path(output).resolve()
    if output.is_relative_to(directory) or directory.is_relative_to(output):
        raise ValueError('复算目录与原结果目录不能互相包含')
    current = run_request(read_json(directory / 'request.json'), output, plots=plots)
    records = []
    if len(current['cases']) != len(reference['cases']):
        raise ValueError('复算求解数量与原记录不同')
    for previous, repeated in zip(reference['cases'], current['cases']):
        relative = previous['directory']
        rows = _compare_arrays(directory / relative / 'arrays.npz', output / relative / 'arrays.npz')
        summary_passed = previous['summary'].keys() == repeated['summary'].keys()
        for key, value in previous['summary'].items():
            other = repeated['summary'].get(key)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                summary_passed &= bool(np.isclose(value, other, rtol=1e-10, atol=1e-14))
            else:
                summary_passed &= value == other
        input_passed = read_json(directory / relative / 'input.json') == read_json(output / relative / 'input.json')
        records.append({'directory': relative, 'input_equal': input_passed, 'summary_passed': bool(summary_passed),
                        'arrays': rows, 'passed': bool(input_passed and summary_passed and all(r['passed'] for r in rows))})
    study_passed = (read_json(directory / 'study.json') == read_json(output / 'study.json')
                    if (directory / 'study.json').exists() else not (output / 'study.json').exists())
    result = {'schema': 'gear_ltca_replay_v1', 'reference_run_sha256': sha256(directory / 'run.json'),
              'passed': all(r['passed'] for r in records) and study_passed, 'study_equal': study_passed,
              'cases': records, 'scope': '同源码与数值环境的独立重新求解和全数组复现，不构成工程精度验证'}
    write_json(output / 'replay.json', result)
    current['validation']['numerical_replay'] = '通过' if result['passed'] else '失败'
    current['artifacts'] = artifact_hashes(output)
    write_json(output / 'run.json', current)
    return result
