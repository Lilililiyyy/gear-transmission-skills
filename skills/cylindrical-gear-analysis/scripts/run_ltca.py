"""新写Skill入口：依赖已安装的gear-transmission-skills，不改变计算公式。"""
import sys

try:
    from gear_transmission.cli import main
except ModuleNotFoundError as error:
    print('缺少运行依赖：请按references/runtime.md在当前Python环境安装完整项目。缺失模块：' + str(error),
          file=sys.stderr)
    raise SystemExit(1)

if __name__ == '__main__':
    raise SystemExit(main())
