"""由命名JSON执行可复现计算的简洁命令入口。"""
import argparse
import json
from pathlib import Path

from .config import analyze_config
from .io import export_result


def main():
    parser = argparse.ArgumentParser(description="渐开线圆柱齿轮准静态承载齿面接触分析（LTCA）")
    parser.add_argument("input", type=Path, help="命名JSON输入配置，内部采用SI单位")
    parser.add_argument("--output", required=True, type=Path, help="结果目录")
    parser.add_argument("--plots", action="store_true", help="同时生成中文PNG/SVG图表")
    args = parser.parse_args()
    result = analyze_config(json.loads(args.input.read_text(encoding="utf-8-sig")))
    export_result(result, args.output)
    if args.plots:
        from .plotting import plot_analysis
        plot_analysis(result, args.output)
    print(json.dumps(result.summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
