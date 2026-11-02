"""线图统一红蓝分类色，接触压力用常规色阶；不改变几何或原数值。"""
from __future__ import annotations


BLUE = "#2166AC"
RED = "#B2182B"
NEUTRAL = "#616970"
COLORS = (BLUE, RED, "#67A9CF", "#D6604D", NEUTRAL)
COLOUR_POLICY = "red_blue_publication_v1"
PRESSURE_COLOUR_POLICY = "blue_cyan_green_yellow_red_pressure_v1"
PRESSURE_BLUE = "#000080"
PRESSURE_RED = "#800000"


def red_blue_colormap():
    """返回低值蓝、高值红的连续线性色图，未知值保留灰色。"""
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list(
        "bevel_red_blue", ((0.0, BLUE), (0.5, "#F7F7F7"), (1.0, RED)), N=1024)
    cmap.set_under(BLUE)
    cmap.set_over("#67001F")
    cmap.set_bad("#B5BBC0")
    return cmap


def contact_pressure_colormap():
    """接触压力采用常规蓝—青—绿—黄—红色阶，零压蓝、高压红。

    只控制非负压力的显示；线图仍采用红蓝分类色。数值归一化和各图
    共用的绝对压力上限由调用者给定，不对原压力作平滑或截断。
    """
    from matplotlib import colormaps

    cmap = colormaps["jet"].resampled(1024)
    cmap.set_under(PRESSURE_BLUE)
    cmap.set_over(PRESSURE_RED)
    cmap.set_bad("#B5BBC0")
    return cmap


def figure_style():
    """论文图共用样式：白底、细轴、可编辑文字与红蓝分类色循环。"""
    from cycler import cycler

    return {"axes.prop_cycle": cycler(color=COLORS),
        "image.cmap": "RdBu_r", "axes.unicode_minus": False,
        "legend.frameon": False, "xtick.direction": "out", "ytick.direction": "out",
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"}
