"""新写重绘入口：验证原结果，只读原数组，以源工程统一风格输出新图。"""
from collections import defaultdict
from pathlib import Path

from gear_contact.contact_display import load_plot_result
from gear_contact.io import write_json

from .plotting import plot_result, plot_cycles, plot_contact_pass_comparison
from .protocol import read_json
from .workflow import artifact_hashes, identity, sha256


def replot_run(directory, output):
    """允许旧版本结果按新版显示规则重绘；代码版本变化不被冒充数值复算。"""
    directory, output = Path(directory).resolve(), Path(output).resolve()
    run = read_json(directory/'run.json')
    if (run.get('schema') != 'gear_ltca_run_v1' or run.get('status') != 'completed'
            or not run.get('cases') or artifact_hashes(directory) != run.get('artifacts')):
        raise ValueError('原运行未完成，或保存的输入／结果文件已变化，不能重绘')
    if output.is_relative_to(directory) or directory.is_relative_to(output):
        raise ValueError('重绘目录与原结果目录不能互相包含')
    views = []
    for case in run['cases']:
        relative = Path(case['directory'])
        if relative.is_absolute() or not (directory/relative).resolve().is_relative_to(directory):
            raise ValueError('原求解目录超出结果范围')
        views.append(load_plot_result(directory/relative))
    before = identity()
    output.mkdir(parents=True, exist_ok=False)
    record = {'schema': 'gear_ltca_plot_v1', 'status': 'running', 'solver_executed': False,
              'numerical_replay_performed': False, 'source_run_sha256': sha256(directory/'run.json'),
              'source_identity': run['identity'], 'plotting_identity': before, 'cases': []}
    groups = defaultdict(dict)
    try:
        for index, (case, view) in enumerate(zip(run['cases'], views), 1):
            relative = f'cases/case_{index:03d}/plots'
            metrics = plot_result(view, output/relative)
            record['cases'].append({'design': case['design'], 'directory': relative, 'summary': view.summary,
                'source_arrays_sha256': sha256(directory/case['directory']/'arrays.npz'), 'plot_definitions': metrics})
            groups[case['operating_point']][case['design']] = view
        for index, candidates in enumerate(groups.values(), 1):
            plot_cycles(candidates, output/f'plots/point_{index:03d}')
            plot_contact_pass_comparison(candidates, output/f'plots/point_{index:03d}')
        if identity() != before or artifact_hashes(directory) != run['artifacts']:
            raise RuntimeError('重绘期间源码或原始结果发生变化')
        record['status'] = 'completed'
        (output/'report.md').write_text('# 接触图重绘\n\n原输入、数组、指标及原目录保持不变。'
            '本次只重绘，未执行求解或数值复算。\n\n'
            '比较图位于 plots/point_NNN/07_contact_pass_comparison.png；单方案图位于 cases/case_NNN/plots/。\n', encoding='utf-8')
    except Exception as error:
        record.update(status='failed', error={'type': type(error).__name__, 'message': str(error)})
        raise
    finally:
        record['artifacts'] = artifact_hashes(output)
        write_json(output/'plot.json', record)
    return record
