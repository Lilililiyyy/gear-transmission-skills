"""新写命令适配：校验请求、运行内核或复算已保存结果，错误返回非零退出码。"""
import argparse
import json
import sys

from .protocol import prepare_request, read_json
from .workflow import identity, replay_run, run_request


def main(argv=None):
    parser = argparse.ArgumentParser(description='齿轮传动设计Skill库：圆柱齿轮LTCA试点')
    commands = parser.add_subparsers(dest='command', required=True)
    check = commands.add_parser('check', help='校验请求、显式单位与内核身份，不执行数值求解')
    check.add_argument('request', help='请求JSON文件')
    run = commands.add_parser('run', help='实际求解并保存可复算结果')
    run.add_argument('request', help='请求JSON文件')
    case = commands.add_parser('case', help='直接运行新版完整案例：显式因子强度、LTCA及接触图对标')
    case.add_argument('request', nargs='?', help='独立项目案例JSON；省略时使用随附18/77案例')
    replay = commands.add_parser('replay', help='验证保存文件后，重新求解并比较全部原数组')
    replay.add_argument('result', help='已完成的原结果目录')
    for command in (run, replay, case):
        command.add_argument('--output', required=True, help='尚不存在的新结果目录')
        command.add_argument('--plots', action='store_true', help='同时生成中文PNG／SVG图表')
    args = parser.parse_args(argv)
    try:
        if args.command == 'check':
            request = read_json(args.request)
            if request.get('schema_version') == 1 and 'analysis' in request and 'loading' in request:
                from .case_workflow import prepare_case_request
                prepared = prepare_case_request(request)
            else:
                prepared = prepare_request(request)
            result = {'request_valid': True, 'solve_count': prepared['solve_count'],
                      'defaulted_fields': prepared['defaults'], 'unit_conversions': prepared['conversions'],
                      'identity': identity(), 'scope': '结构与基本参数校验；不证明网格可用或数值求解成功'}
        elif args.command in ('run', 'case'):
            from .upstream import snapshot_root
            path = args.request or snapshot_root()/'examples/helical_18_77_case.json'
            run = run_request(read_json(path), args.output, plots=args.plots)
            result = {'status': run['status'], 'solve_count': len(run['cases']), 'output': args.output}
        else:
            result = replay_run(args.result, args.output, plots=args.plots)
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 1 if result.get('passed') is False else 0
    except Exception as error:
        print(json.dumps({'status': 'failed', 'error': type(error).__name__, 'message': str(error)},
                         ensure_ascii=False), file=sys.stderr)
        return 1
