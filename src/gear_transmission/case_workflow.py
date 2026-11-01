"""新写适配：直接运行保留的新版案例入口，保存代码身份并支持全数组复算。"""
from datetime import datetime, timezone
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from gear_contact.io import write_json

from .protocol import read_json
from .upstream import snapshot_root


def prepare_case_request(request):
    """复用原案例载荷记账和因子数据类，校验实际 LTCA 输入但不求解。"""
    from importlib import import_module
    from .protocol import prepare_case
    previous_path = sys.path[:]
    try:
        module = import_module('gear_ltca_snapshot.examples.run_cylindrical_workflow')
        config, strength = module.prepare_case(request)
    finally:
        sys.path[:] = previous_path
    pair, kwargs, defaults = prepare_case(config)
    for name in ('label', 'case_id'):
        if not isinstance(request.get(name), str) or not request[name].strip():
            raise ValueError('完整案例需要非空字段：'+name)
    if not isinstance(request.get('instant'), dict):
        raise ValueError('完整案例需要 instant，指定同相位同齿对的选择')
    return {'pair': pair, 'kwargs': kwargs, 'defaults': defaults, 'strength': strength,
            'conversions': [], 'solve_count': 1}


def run_case_request(request, output, *, plots=False):
    """案例协议保持上游 SI 字段；强度因子和材料应力限均为显式输入。"""
    from .workflow import artifact_hashes, identity
    before = identity()
    prepared = prepare_case_request(request)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'request.json', request)
    manifest = {'schema': 'gear_ltca_run_v1', 'workflow_kind': 'cylindrical_case_v1',
                'status': 'running', 'analysis_completed': False, 'plots': plots,
                'created_utc': datetime.now(timezone.utc).isoformat(), 'identity': before,
                'validation': {'numerical_replay': '尚未复算', 'physical_validation': False}, 'cases': []}
    write_json(output/'run.json', manifest)
    try:
        # 输入放在独立临时目录，遵守原入口对输入与输出相互包含的保护。
        with tempfile.TemporaryDirectory(prefix='gear-case-', dir=output.parent) as temporary:
            case_path = Path(temporary)/'case.json'
            write_json(case_path, request)
            command = [sys.executable, '-B', str(snapshot_root()/'examples/run_cylindrical_workflow.py'),
                       str(case_path), '--output', str(output)]
            if plots:
                command.append('--plots')
            process = subprocess.run(command, capture_output=True, text=True, encoding='utf-8',
                                     env={**os.environ, 'PYTHONUTF8': '1', 'PYTHONDONTWRITEBYTECODE': '1'})
        (output/'workflow.log').write_text(process.stdout+process.stderr, encoding='utf-8')
        if process.returncode:
            raise RuntimeError('独立案例入口执行失败，详情见 workflow.log：'+process.stderr[-1200:])
        summary = read_json(output/'workflow_summary.json')
        if plots:
            from gear_ltca_snapshot.examples.case_geometry_plotting import plot_case_geometry
            summary['geometry_figures'] = plot_case_geometry(prepared['pair'], output,
                modification_pinion=prepared['kwargs']['modification_pinion'],
                modification_wheel=prepared['kwargs']['modification_wheel'])
            write_json(output/'workflow_summary.json', summary)
        for array in sorted(output.rglob('arrays.npz')):
            relative = array.parent.relative_to(output).as_posix()
            case_summary = read_json(array.parent/'summary.json')
            if not case_summary['all_converged']:
                raise ArithmeticError('案例含未收敛求解：'+relative)
            manifest['cases'].append({'directory': relative, 'design': relative,
                                      'operating_point': request['label'], 'summary': case_summary})
        if not manifest['cases']:
            raise RuntimeError('独立案例入口未导出原始求解数组')
        manifest['strength'] = summary['strength']
        manifest['validation'].update({key: summary[key] for key in (
            'reference_binding_verified', 'same_inputs_as_reference', 'geometry_comparison_passed',
            'strength_comparison_passed', 'phase_geometry_registered', 'ltca_response_accuracy_verified')})
        if identity() != before:
            raise RuntimeError('运行期间源码或资源发生变化')
        manifest.update(status='completed', analysis_completed=True)
        rows = ['# 圆柱齿轮新版案例工作流', '',
                '直接调用随库保留的独立项目案例入口；所有本项目响应均重新计算。', '',
                f"案例：{request['label']}。主工况 TE 峰峰值：{summary['ltca_summary']['te_peak_to_peak_m']*1e6:.8g} μm。", '',
                '完整图表、几何／显式因子强度对比及限制见 workflow_summary.json；原始 LTCA 数据见 ltca/。',
                '启用图表时，06 文件为原入口固定的附加失准／鼓形敏感性案例，不表示自动优化。', '',
                '标准强度使用输入的评级因子及材料应力限；没有独立预测完整标准因子或材料寿命链。',
                'MASTA 资源仅用于已绑定案例对标，不作为求解响应；接触图为显示后处理，不是完整二维接触求解。',
                '几何与显式因子强度对比、LTCA 响应精度、数值复算及物理验证分别记录。', '',
                '运行 gear-design replay 本目录 --output 新目录，重新求解并比较全部原数组和强度结果。', '']
        (output/'report.md').write_text('\n'.join(rows), encoding='utf-8')
    except Exception as error:
        manifest.update(status='failed', analysis_completed=False,
                        error={'type': type(error).__name__, 'message': str(error)})
        raise
    finally:
        manifest['artifacts'] = artifact_hashes(output)
        write_json(output/'run.json', manifest)
    return manifest
