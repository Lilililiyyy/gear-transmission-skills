"""新写锥齿轮命令入口，使用原生案例协议及SI单位。"""
import argparse
import json
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description='锥齿轮几何、强度与分阶段接触分析Skill')
    commands = parser.add_subparsers(dest='command', required=True)
    check = commands.add_parser('check', help='校验输入并执行代数几何／评级检查，不求接触')
    run = commands.add_parser('run', help='实际计算所选阶段，保存输入、身份和原始结果')
    replay = commands.add_parser('replay', help='校验原文件后重新计算并比较数学结果')
    for command in (check, run):
        command.add_argument('request', help='原生bevel_analysis_case_v1 JSON')
        command.add_argument('--through', choices=('geometry', 'strength', 'cad', 'tca', 'ltca'), default='strength', help='终止阶段，默认strength')
    for command in (run, replay):
        command.add_argument('--output', required=True, help='尚不存在的新结果目录')
        command.add_argument('--plots', action='store_true', help='另存原内核的中文图表和SI源表')
    run.add_argument('--reference', help='与同工况输入绑定的外部几何／强度参照JSON')
    replay.add_argument('result', help='已完成结果目录')
    args = parser.parse_args(argv)
    try:
        from .bevel_workflow import identity, prepare_case, read_json, run_case, replay_run
        if args.command in ('check', 'run'):
            from bevel_gear.workflow import load_bevel_case
            case = load_bevel_case(args.request)
            if args.command == 'check':
                prepare_case(case, args.through)
                result = dict(request_valid=True, identity=identity(), through=args.through,
                              scope='输入、代数几何及所选评级分支检查；未验证CAD或接触数值收敛')
            else:
                result = run_case(case, args.output, through=args.through, plots=args.plots,
                                  reference=read_json(args.reference) if args.reference else None)
                result = {key: result[key] for key in ('status', 'through', 'analysis_completed')}
        else:
            result = replay_run(args.result, args.output, plots=args.plots)
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 1 if result.get('passed') is False else 0
    except Exception as error:
        print(json.dumps(dict(status='failed', error=type(error).__name__, message=str(error)), ensure_ascii=False), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
