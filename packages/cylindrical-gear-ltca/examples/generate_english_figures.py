"""新写英文图表入口：仅翻译图中文字，复用同一案例、计算函数和数值输出。"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
from functools import partial
from pathlib import Path
import re
import sys
import textwrap

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p/"pyproject.toml").is_file() and (p/"gear_contact").is_dir())
for folder in (PROJECT, PROJECT/"examples"):
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))


TRANSLATIONS = {
    "a  关键几何参数": "a  Key geometry",
    "b  ISO 6336:2006 显式因子应力": "b  ISO 6336:2006 stresses",
    "c  同工况的两种载荷层级": "c  LTCA and nominal rating loads",
    "d  外部材料／寿命应力限下的安全系数": "d  Safety factors",
    "中心距 a": "Centre distance a", "重叠齿宽 b": "Overlapping face width b",
    "小轮分度圆直径": "Pinion reference diameter", "大轮分度圆直径": "Wheel reference diameter",
    "小轮基圆直径": "Pinion base diameter", "大轮基圆直径": "Wheel base diameter",
    "小轮名义齿根": "Pinion nominal root", "大轮名义齿根": "Wheel nominal root",
    "小轮评级接触": "Pinion rated contact", "大轮评级接触": "Wheel rated contact",
    "小轮评级齿根": "Pinion rated root", "大轮评级齿根": "Wheel rated root",
    "名义接触": "Nominal contact", "小轮接触": "Pinion contact", "大轮接触": "Wheel contact",
    "小轮弯曲": "Pinion bending", "大轮弯曲": "Wheel bending",
    "LTCA：评级后总法向 Fn": "LTCA: normal Fn",
    "强度：名义参考圆 Ft": "Rating: nominal Ft",
    "18/77斜齿轮副：关键几何尺寸、标准强度与安全系数对比。":
        "18/77 helical gear pair: key geometry, rated stresses and safety factors.",
    "强度计算采用相同评级因子与材料应力限；标准强度与LTCA局部接触压力分别评价。":
        "The same rating factors and material stress limits are used; rated stresses and local LTCA pressure are assessed separately.",
    "a  承载法向TE": "a  Loaded normal TE", "b  法向刚度": "b  Normal stiffness",
    "c  每相位峰值接触压力": "c  Peak pressure at each mesh position", "d  总法向载荷": "d  Total normal load",
    "18/77斜齿轮副：相同载荷与工作面下的LTCA周期响应。":
        "18/77 helical gear pair: LTCA response at the same load and working flank.",
    "展示3个稳态啮合周期，每周期32个采样位置。":
        "Three steady-state mesh periods are displayed, with 32 samples per period.",
    "同相位、同齿对的瞬时接触压力": "Instantaneous contact pressure: same phase and tooth pair",
    "沿接触线长度单位为mm，垂直接触线宽度单位为μm，两轴显示比例不同。":
        "Distance along the contact line is in mm; transverse width is in μm. The two axes use different scales.",
    "两图采用相同压力色标；灰色为未覆盖区域，蓝色为零压力。":
        "Both plots use the same pressure scale. Grey denotes uncovered regions; blue denotes zero pressure.",
    "圆点为普通接触，三角为延伸接触，虚线为有效齿廓边界。":
        "Dots: ordinary contact; triangles: extended contact; dashed lines: active profile limits.",
    "上方为单相位接触位置，下方为一个啮合周期的接触覆盖。":
        "Top: contact positions at one phase. Bottom: contact coverage over one mesh period.",
    "18/77斜齿轮副：保持几何、材料与载荷，比较附加失准和齿向鼓形的影响。":
        "18/77 helical gear pair: geometry, material and load are fixed while additional misalignment and lead crowning are compared.",
    "本项目计算结果，展示三个稳态啮合周期。":
        "Results from this project, displayed over three steady-state mesh periods.",
    "18/77斜齿轮副：相同载荷与工作面下，一个齿通过程的接触压力分布。":
        "18/77 helical gear pair: tooth-passage contact pressure at the same load and working flank.",
    "横轴为齿宽坐标，纵轴为小轮展长；灰色为未覆盖区域。":
        "Horizontal axis: face-width coordinate. Vertical axis: pinion roll distance. Grey denotes uncovered regions.",
    "本项目：小轮整齿通过程接触压力": "This project: pinion tooth-passage contact pressure",
    "附加1 μm失准": "Additional 1 μm misalignment",
    "失准 + 0.5 μm鼓形": "Misalignment + 0.5 μm crowning", "原案例": "Baseline",
    "a  法向TE": "a  Normal TE", "b  端点峰值接触压力": "b  Endpoint peak contact pressure",
    "相对起采位置 / 基节周期": "Mesh periods from the first sampled position",
    "条带中心线载荷 / (N/mm)": "Strip-centre line load / (N/mm)",
    "垂直接触线 n / μm": "Across contact line n / μm", "沿接触线 s / mm": "Along contact line s / mm",
    "齿宽坐标 z / mm": "Face-width coordinate z / mm", "齿宽坐标 / mm": "Face-width coordinate / mm",
    "小轮展长 q / mm": "Pinion roll distance q / mm", "小轮展长 / mm": "Pinion roll distance / mm",
    "端点峰压 / MPa": "Endpoint peak pressure / MPa", "法向接近量 / μm": "Normal approach / μm",
    "局部接触压力 / MPa": "Local contact pressure / MPa", "接触压力 / MPa": "Contact pressure / MPa",
    "接近量 / μm": "Normal approach / μm", "刚度 / (N/μm)": "Stiffness / (N/μm)",
    "压力 / MPa": "Pressure / MPa", "载荷 / N": "Load / N", "数值 / ": "Value / ",
    "周期离散采样位置": "sampled cycle coverage", "同一瞬时": "instantaneous",
    "本项目": "This project", "齿对": "tooth pair ", "条带力": "line force ", "峰压": "peak pressure ",
    "：": ": ", "；": "; ",
}


def translate_text(value):
    """翻译显示标签并检查遗漏；不处理数据、单位、数字和MASTA原始位图。"""
    translated = value
    for original in sorted(TRANSLATIONS, key=len, reverse=True):
        translated = translated.replace(original, TRANSLATIONS[original])
    if re.search(r"[\u3400-\u9fff]", translated):
        raise ValueError("英文图表存在未翻译标签："+value)
    return translated


@contextmanager
def english_figure_labels():
    """在保存图表时替换文字，退出后恢复Matplotlib原保存方法。"""
    from matplotlib.figure import Figure
    from matplotlib.text import Text
    from matplotlib.ticker import FixedFormatter
    original_save = Figure.savefig

    def save(figure, *args, **kwargs):
        # 类别刻度在draw时会由Formatter重新生成，须同步翻译其标签来源。
        for axes in figure.axes:
            for axis in (axes.xaxis, axes.yaxis):
                formatter = axis.get_major_formatter()
                if isinstance(formatter, FixedFormatter):
                    formatter.seq = [translate_text(str(label)) for label in formatter.seq]
                elif isinstance(getattr(formatter, "func", None), partial):
                    # set_ticks(labels=...)使用带字典的FuncFormatter，绘制前更新字典。
                    for argument in formatter.func.args:
                        if isinstance(argument, dict):
                            argument.update({key: translate_text(str(label)) for key, label in argument.items()})
        for artist in figure.findobj(match=Text):
            translated = translate_text(artist.get_text())
            if artist in figure.texts and artist.get_position()[1] < .2:
                translated = "\n".join(textwrap.fill(line, width=112) for line in translated.splitlines())
                artist.set_fontsize(7.3)
            artist.set_text(translated)
            artist.set_fontfamily("DejaVu Sans")
        for legend in figure.legends:
            for artist in legend.get_texts():
                artist.set_fontsize(8.)
        return original_save(figure, *args, **kwargs)

    Figure.savefig = save
    try:
        yield
    finally:
        Figure.savefig = original_save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    from run_cylindrical_workflow import run
    with english_figure_labels():
        summary = run(PROJECT/"examples/helical_18_77_case.json", args.output,
            reference_path=PROJECT/"examples/helical_18_77_reference.json", plots=True)
    print("英文案例图已生成："+str(args.output.resolve()))


if __name__ == "__main__":
    main()
