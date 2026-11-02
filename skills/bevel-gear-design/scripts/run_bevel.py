"""新写Skill入口：使用已安装的完整项目，不改变锥齿轮计算公式。"""
import sys

try:
    from gear_transmission.bevel_cli import main
except ModuleNotFoundError as error:
    print('缺少运行依赖，请按references/runtime.md安装完整项目并使用其Python。缺失模块：'+str(error), file=sys.stderr)
    raise SystemExit(1)


if __name__ == '__main__':
    raise SystemExit(main())
