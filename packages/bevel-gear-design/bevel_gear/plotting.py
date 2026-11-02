"""将本次计算结果绘制为有源数据的中文论文图。

统一单／双栏物理尺寸与字体。默认保留PNG/SVG/CSV接口；
publication=True额外导出PDF。图表不补造未知计算量或统计误差线。
"""
from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass
import math
from pathlib import Path

import numpy as np


from .publication_style import COLORS as _PUBLICATION_COLORS, NEUTRAL as _NEUTRAL, figure_style
_COLORS = _PUBLICATION_COLORS[:4]
_SINGLE_MM = 89.0
_DOUBLE_MM = 183.0


def _mapping(value):
    if callable(getattr(value, "to_dict", None)):
        return value.to_dict()
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return value
    raise TypeError("图表输入需要结果字典、数据类或to_dict接口")


def _pyplot():
    """仅调用时加载后端；在最终物理尺寸下采用7–9pt中文无衬线字体。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    fonts = [name for name in ("Microsoft YaHei", "Noto Sans CJK SC", "SimHei",
                              "Source Han Sans SC", "DejaVu Sans") if name in available]
    plt.rcParams.update({**figure_style(),
        "font.family": "sans-serif", "font.sans-serif": fonts, "font.size": 7.5,
        "axes.labelsize": 7.5, "axes.titlesize": 8, "axes.titleweight": "normal",
        "axes.unicode_minus": False, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .6, "axes.edgecolor": "#333333",
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "xtick.major.width": .6, "ytick.major.width": .6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "lines.linewidth": 1, "lines.markersize": 3,
        "legend.fontsize": 7, "legend.frameon": False, "legend.borderpad": .2,
        "legend.handlelength": 1.6, "legend.labelspacing": .3,
        "legend.columnspacing": 1, "legend.handletextpad": .4,
        "grid.color": "#C6CDD2", "grid.linewidth": .4, "grid.alpha": .55,
        "svg.fonttype": "none", "pdf.fonttype": 42, "pdf.compression": 6,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "figure.constrained_layout.h_pad": .045,
        "figure.constrained_layout.w_pad": .045,
        "figure.constrained_layout.hspace": .08,
        "figure.constrained_layout.wspace": .08,
    })
    return plt


def _size(width_mm, height_mm):
    return width_mm / 25.4, height_mm / 25.4


def _panel(ax, letter, heading):
    """小写面板标记与短标题共用左对齐基线，避免大型总标题。"""
    ax.set_title(r"$\bf{" + letter + "}$  " + heading, loc="left", pad=6)


def _grid(ax, axis="y"):
    ax.set_axisbelow(True)
    ax.grid(axis=axis)


def _save(fig, directory: Path, stem: str, *, publication=False, title=""):
    """保持实际画布尺寸；不以tight裁剪改变单／双栏输出宽度。"""
    if not isinstance(publication, bool):
        raise TypeError("publication须为布尔值")
    files = []
    description = "确定性工程计算；合成案例需结合输入和模型范围解释；无统计误差线"
    for extension in ("png", "svg", *(("pdf",) if publication else ())):
        path = directory / f"{stem}.{extension}"
        metadata = ({"Title": title, "Description": description, "Date": None}
                    if extension == "svg" else
                    {"Title": title, "Subject": description, "CreationDate": None}
                    if extension == "pdf" else
                    {"Title": title, "Description": description})
        fig.savefig(path, dpi=600, facecolor="white", metadata=metadata)
        files.append(str(path.resolve()))
    return files


def _csv(path: Path, header, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    return str(path.resolve())


def plot_geometry(geometry, output_dir, *, title="锥齿轮宏观节锥剖面", publication=False):
    """由已计算锥角与锥距画宏观剖面；非实际制造齿形。"""
    g = _mapping(geometry)
    pinion, wheel = _mapping(g["pinion"]), _mapping(g["wheel"])
    delta1, delta2 = (float(gear["pitch_cone_angle_rad"]) for gear in (pinion, wheel))
    shaft_angle = delta1 + delta2
    outer, inner, mean = (float(g[key]) for key in
                          ("outer_cone_distance_m", "inner_cone_distance_m", "mean_cone_distance_m"))
    if not 0 < inner <= mean <= outer:
        raise ValueError("节锥剖面要求正的内、平均、外锥距并保持顺序")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    from matplotlib.patches import Arc
    fig, ax = plt.subplots(figsize=_size(_SINGLE_MM, 84), constrained_layout=True)
    rows = []
    for name, axis, delta, color in (("小轮", 0., delta1, _COLORS[0]),
                                     ("大轮", shaft_angle, delta2, _COLORS[1])):
        directions = [np.array([math.cos(axis + side * delta), math.sin(axis + side * delta)])
                      for side in (-1, 1)]
        points = np.array([inner * directions[0], outer * directions[0],
                           outer * directions[1], inner * directions[1]]) * 1e3
        ax.fill(points[:, 0], points[:, 1], color=color, alpha=.08)
        ax.plot(np.r_[points[:, 0], points[0, 0]], np.r_[points[:, 1], points[0, 1]],
                color=color, linewidth=1.1)
        midpoints = np.array([mean * direction for direction in directions]) * 1e3
        ax.plot(midpoints[:, 0], midpoints[:, 1], "--", color=color, linewidth=.8)
        for side, direction in zip((-1, 1), directions):
            for position, radius in (("内端", inner), ("平均截面", mean), ("外端", outer)):
                point = radius * direction
                rows.append((name, side, position, point[0], point[1]))
        direction = np.array([math.cos(axis), math.sin(axis)])
        ax.plot([0, 1.08 * outer * direction[0] * 1e3],
                [0, 1.08 * outer * direction[1] * 1e3], color=color, linestyle="-.", linewidth=.7)
        ax.annotate(f"{name}轴线", xy=direction * outer * 1.02e3,
                    xytext=(3, 3), textcoords="offset points", color=color, fontsize=7)
    shared = np.array([math.cos(delta1), math.sin(delta1)])
    ax.plot([0, outer * shared[0] * 1e3], [0, outer * shared[1] * 1e3],
            color=_NEUTRAL, linestyle=":", linewidth=.7)
    radius_mm = outer * 1e3 * .24
    ax.add_patch(Arc((0, 0), 2 * radius_mm, 2 * radius_mm,
                     theta1=0, theta2=math.degrees(shaft_angle), color=_NEUTRAL, linewidth=.6))
    middle = shaft_angle / 2
    ax.text(radius_mm * 1.23 * math.cos(middle), radius_mm * 1.23 * math.sin(middle),
            f"Σ={math.degrees(shaft_angle):.1f}°", ha="center", va="center", fontsize=7)
    ax.scatter([0], [0], color=_NEUTRAL, s=9)
    ax.annotate("公共锥顶", (0, 0), xytext=(-3, -12), textcoords="offset points", ha="center", fontsize=7)
    ax.text(.98, .98, "宏观节锥\n非加工齿面", transform=ax.transAxes,
            ha="right", va="top", color=_NEUTRAL, fontsize=7)
    ax.set(xlabel="剖面坐标 x / mm", ylabel="剖面坐标 y / mm")
    _panel(ax, "a", "节锥与齿宽")
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(.16)
    _grid(ax, "both")
    files = _save(fig, output, "cone_section", publication=publication, title=title)
    plt.close(fig)
    files.append(_csv(output / "cone_section.csv",
                      ("gear", "side", "section", "x_m", "y_m"), rows))
    return files


def plot_forces(forces, output_dir, *, title="两轮局部坐标下的啮合力", publication=False):
    """保留有符号分量；两个局部坐标系不能直接相加为全局力。"""
    f = _mapping(forces)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=_size(_SINGLE_MM, 70), constrained_layout=True)
    rows, locations = [], np.arange(2)
    for i, (key, label, hatch) in enumerate((("tangential", "切向 Ft", ""),
                                            ("radial", "径向 Fr", "//"),
                                            ("axial", "轴向 Fa", ".."))):
        values = [float(f[f"{gear}_{key}_n"]) for gear in ("pinion", "wheel")]
        ax.bar(locations + (i - 1) * .23, values, width=.21, label=label,
               color=_COLORS[i], edgecolor="white", linewidth=.4, hatch=hatch)
        rows.extend((gear, key, value) for gear, value in zip(("小轮", "大轮"), values))
    ax.axhline(0, color=_NEUTRAL, linewidth=.6)
    ax.set(xticks=locations, xticklabels=("小轮", "大轮"), ylabel="局部力分量 / N")
    _panel(ax, "a", "局部啮合力")
    ax.legend(ncol=3, loc="upper center")
    ax.margins(y=.32)
    _grid(ax)
    files = _save(fig, output, "mesh_forces", publication=publication, title=title)
    plt.close(fig)
    files.append(_csv(output / "mesh_forces.csv", ("gear", "component", "force_n"), rows))
    return files


_RATING_KEYS = (
    "contact_stress_pa", "pinion_bending_stress_pa", "wheel_bending_stress_pa",
    "pinion_allowable_contact_stress_pa", "wheel_allowable_contact_stress_pa",
    "pinion_allowable_bending_stress_pa", "wheel_allowable_bending_stress_pa",
    "pinion_contact_capacity_pa", "wheel_contact_capacity_pa",
    "pinion_bending_capacity_pa", "wheel_bending_capacity_pa",
    "pinion_contact_safety", "wheel_contact_safety", "pinion_bending_safety", "wheel_bending_safety",
    "required_contact_safety", "required_bending_safety",
)


def _rating(rating):
    values = _mapping(rating)
    for name in _RATING_KEYS:
        if name not in values or not isinstance(values[name], (int, float)) or not math.isfinite(values[name]):
            raise ValueError(f"强度图缺少有效计算字段：{name}")
    return values


def plot_strength(rating, output_dir, *, title="给定评级条件下的强度与安全系数", publication=False):
    """分开呈现接触／弯曲应力与各自目标；不推断未知许用值。"""
    r = _rating(rating)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    from matplotlib.patches import Patch
    fig, axes = plt.subplots(2, 2, figsize=_size(_DOUBLE_MM, 120), constrained_layout=True)
    x = np.arange(2)
    for ax, mode, letter, heading in ((axes[0, 0], "contact", "a", "接触应力"),
                                      (axes[0, 1], "bending", "b", "齿根弯曲应力")):
        actual = ([r["contact_stress_pa"] / 1e6] * 2 if mode == "contact" else
                  [r[f"{gear}_bending_stress_pa"] / 1e6 for gear in ("pinion", "wheel")])
        allowable = [r[f"{gear}_allowable_{mode}_stress_pa"] / 1e6 for gear in ("pinion", "wheel")]
        ax.bar(x - .17, actual, width=.29, color=_COLORS[:2], linewidth=.6)
        ax.bar(x + .17, allowable, width=.29, facecolor="white", edgecolor=_COLORS[:2],
               hatch="///", linewidth=.7)
        ax.set(xticks=x, xticklabels=("小轮", "大轮"), ylabel="应力 / MPa",
               ylim=(0, max(actual + allowable) * 1.35))
        _panel(ax, letter, heading)
        _grid(ax)
    axes[0, 0].legend(handles=(Patch(facecolor=_NEUTRAL, label="计算应力"),
                              Patch(facecolor="white", edgecolor=_NEUTRAL, hatch="///", label="目标许用")),
                       loc="upper center", ncol=2)
    for ax, mode, letter, symbol, heading in ((axes[1, 0], "contact", "c", "SH", "接触安全系数"),
                                              (axes[1, 1], "bending", "d", "SF", "弯曲安全系数")):
        values = [r[f"{gear}_{mode}_safety"] for gear in ("pinion", "wheel")]
        target = r[f"required_{mode}_safety"]
        ax.bar(x, values, width=.52, color=_COLORS[:2], edgecolor="white", linewidth=.4)
        ax.axhline(target, color=_NEUTRAL, linestyle="--", linewidth=.8,
                   label=f"目标 {symbol}={target:g}")
        ax.set(xticks=x, xticklabels=("小轮", "大轮"), ylabel=f"{symbol} / 1",
               ylim=(0, max(values + [target]) * 1.32))
        _panel(ax, letter, heading)
        ax.legend(loc="upper right")
        _grid(ax)
    files = _save(fig, output, "strength_and_safety", publication=publication, title=title)
    plt.close(fig)
    files.append(_csv(output / "strength_and_safety.csv",
                      ("gear", "contact_stress_pa", "allowable_contact_stress_pa", "contact_safety",
                       "bending_stress_pa", "allowable_bending_stress_pa", "bending_safety",
                       "contact_capacity_pa", "bending_capacity_pa",
                       "required_contact_safety", "required_bending_safety"),
                      ((gear, r["contact_stress_pa"], r[f"{key}_allowable_contact_stress_pa"],
                        r[f"{key}_contact_safety"], r[f"{key}_bending_stress_pa"],
                        r[f"{key}_allowable_bending_stress_pa"], r[f"{key}_bending_safety"],
                        r[f"{key}_contact_capacity_pa"], r[f"{key}_bending_capacity_pa"],
                        r["required_contact_safety"], r["required_bending_safety"])
                       for gear, key in (("小轮", "pinion"), ("大轮", "wheel")))))
    return files


def _coefficient_label(name):
    """显示紧凑数学下标，CSV仍保留完整协议名称。"""
    if name == "Ifactor":
        return "$I$"
    if name.startswith("Jfactor_"):
        return "$J_" + ("1" if name.endswith("P") else "2") + "$"
    if name == "ZLZvZR":
        return "$Z_LZ_vZ_R$"
    if "_" in name:
        family, suffix = name.split("_", 1)
        suffix = suffix.replace("deltarel", r"\delta rel").replace("epsilon", r"\epsilon").replace("beta", r"\beta").replace("alpha", r"\alpha")
        return "$" + family + "_{" + suffix + "}$"
    return name


def plot_coefficients(coefficients, output_dir, *, publication=False):
    """用分组点线图展示无量纲系数；有量纲弹性系数单独保留。"""
    coefficients = _mapping(coefficients)
    groups = (
        ("载荷系数", ("K_A", "K_O", "K_v", "K_m", "K_Hbeta", "K_Halpha", "K_Fbeta", "K_Falpha")),
        ("接触系数", ("Z_H", "Z_MB", "Z_LS", "Z_beta", "Z_K", "ZLZvZR", "Z_W", "Ifactor", "C_s", "C_xc", "C_H", "K_T", "C_R")),
        ("齿根系数", ("Y_Fa1", "Y_Fa2", "Y_Sa1", "Y_Sa2", "Y_A1", "Y_A2", "Y_J1", "Y_J2", "Y_epsilon", "Y_K", "Y_LS",
                        "Y_deltarelT1", "Y_deltarelT2", "Y_RrelT1", "Y_RrelT2", "Y_X1", "Y_X2", "Jfactor_P", "Jfactor_G", "K_s", "K_x", "K_R")),
    )
    selected = [(label, [(name, float(coefficients[name])) for name in names
                          if name in coefficients and isinstance(coefficients[name], (int, float))
                          and math.isfinite(coefficients[name])]) for label, names in groups]
    selected = [(label, values) for label, values in selected if values]
    if not selected:
        raise ValueError("没有可绘制的已计算无量纲系数")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    height = max(82, 30 + 3.3 * max(len(values) for _, values in selected))
    fig, axes = plt.subplots(1, len(selected), figsize=_size(_DOUBLE_MM, height),
                             constrained_layout=True, squeeze=False)
    for index, (ax, (label, values)) in enumerate(zip(axes.ravel(), selected)):
        names, numbers = zip(*values)
        positions = np.arange(len(values))
        ax.hlines(positions, 0, numbers, color="#BCC6CD", linewidth=.7)
        ax.plot(numbers, positions, "o", color=_COLORS[0], markersize=3)
        ax.axvline(1, color=_NEUTRAL, linestyle=":", linewidth=.6)
        ax.set(yticks=positions, yticklabels=[_coefficient_label(name) for name in names],
               xlabel="系数 / 1", ylim=(len(values) - .4, -.6))
        ax.tick_params(axis="y", length=0)
        ax.margins(x=.12)
        _panel(ax, chr(97 + index), label)
        _grid(ax, "x")
    files = _save(fig, output, "rating_coefficients", publication=publication,
                  title="当前评级的无量纲系数；虚线为数值1参考，不是安全阈值")
    plt.close(fig)
    files.append(_csv(output / "rating_coefficients.csv", ("group", "coefficient", "value", "unit"),
                      ((label, name, value, "1") for label, values in selected for name, value in values)))
    return files


def plot_design(result, output_dir, *, title="合成锥齿轮设计案例", publication=False):
    """真实设计结果生成几何、力、强度和系数图；不重复运行计算内核。"""
    result = _mapping(result)
    files = plot_geometry(result["geometry"], output_dir, title=title + "：节锥剖面", publication=publication)
    files.extend(plot_forces(result["forces"], output_dir, title=title + "：啮合力", publication=publication))
    if result.get("rating") is not None:
        files.extend(plot_strength(result["rating"], output_dir, title=title + "：强度与安全系数", publication=publication))
        files.extend(plot_coefficients(result["rating"]["coefficients"], output_dir, publication=publication))
    return files


def plot_sweep(results, x_values, output_dir, *, xlabel, x_field="x", stem="parameter_sweep",
               title="参数扫描", publication=False):
    """读取逐点重算结果，连接线只表示点序；CSV保留实际采样量。"""
    if len(results) != len(x_values) or len(results) < 2:
        raise ValueError("参数扫描至少需要两个结果，并与横坐标逐点对应")
    x = np.asarray(x_values, dtype=float)
    if not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0):
        raise ValueError("扫描横坐标必须为严格递增有限数")
    ratings = [_rating(_mapping(result)["rating"]) for result in results]
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=_size(_DOUBLE_MM, 126), constrained_layout=True)
    axes[0, 0].plot(x, [r["contact_stress_pa"] / 1e6 for r in ratings], "o-",
                    color=_COLORS[0], label="计算")
    axes[0, 0].plot(x, [r["pinion_allowable_contact_stress_pa"] / 1e6 for r in ratings], "--",
                    color=_COLORS[0], label="小轮目标许用")
    axes[0, 0].plot(x, [r["wheel_allowable_contact_stress_pa"] / 1e6 for r in ratings], ":",
                    color=_COLORS[1], label="大轮目标许用")
    axes[0, 0].set(ylabel="接触应力 / MPa")
    for gear, name, color, marker in (("pinion", "小轮", _COLORS[0], "o"),
                                     ("wheel", "大轮", _COLORS[1], "s")):
        axes[0, 1].plot(x, [r[f"{gear}_bending_stress_pa"] / 1e6 for r in ratings],
                        marker=marker, color=color, label=name)
        axes[1, 0].plot(x, [r[f"{gear}_contact_safety"] for r in ratings],
                        marker=marker, color=color, label=name)
        axes[1, 1].plot(x, [r[f"{gear}_bending_safety"] for r in ratings],
                        marker=marker, color=color, label=name)
    axes[0, 1].set(ylabel="弯曲应力 / MPa")
    axes[1, 0].set(ylabel="SH / 1")
    axes[1, 1].set(ylabel="SF / 1")
    for ax, key in ((axes[1, 0], "required_contact_safety"), (axes[1, 1], "required_bending_safety")):
        ax.plot(x, [r[key] for r in ratings], "--", color=_NEUTRAL, linewidth=.8, label="输入目标")
    for ax, letter, heading in zip(axes.ravel(), "abcd",
                                    ("接触应力", "齿根弯曲应力", "接触安全系数", "弯曲安全系数")):
        ax.set_xlabel(xlabel)
        ax.margins(y=.20)
        _panel(ax, letter, heading)
        _grid(ax)
        ax.legend(loc="best")
    files = _save(fig, output, stem, publication=publication, title=title)
    plt.close(fig)
    files.append(_csv(output / f"{stem}.csv", (x_field, *_RATING_KEYS),
                      ((value, *(r[key] for key in _RATING_KEYS)) for value, r in zip(x, ratings))))
    return files


def plot_comparison(results, labels, output_dir, *, title="合成设计案例的安全系数比较", publication=False):
    """保留各案例两轮与各自目标，不将多参数差异归因为单一参数。"""
    if len(results) != len(labels) or not results:
        raise ValueError("案例结果与名称必须非空且逐一对应")
    ratings = [_rating(_mapping(result)["rating"]) for result in results]
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=_size(_DOUBLE_MM, 78), constrained_layout=True)
    x = np.arange(len(results))
    for ax, mode, letter, heading in ((axes[0], "contact", "a", "接触安全系数 SH"),
                                      (axes[1], "bending", "b", "弯曲安全系数 SF")):
        for offset, gear, name, color, hatch in ((-.17, "pinion", "小轮", _COLORS[0], ""),
                                                 (.17, "wheel", "大轮", _COLORS[1], "//")):
            ax.bar(x + offset, [r[f"{gear}_{mode}_safety"] for r in ratings],
                   width=.31, color=color, edgecolor="white", linewidth=.4, hatch=hatch, label=name)
        ax.plot(x, [r[f"required_{mode}_safety"] for r in ratings], "D--", color=_NEUTRAL,
                linewidth=.8, markersize=2.8, label="各例目标")
        ax.set(xticks=x, xticklabels=labels, ylabel="安全系数 / 1")
        _panel(ax, letter, heading)
        _grid(ax)
        ax.legend(ncol=3, loc="upper center")
        ax.margins(y=.32)
    files = _save(fig, output, "case_comparison", publication=publication, title=title)
    plt.close(fig)
    files.append(_csv(output / "case_comparison.csv", ("case", *_RATING_KEYS),
                      ((label, *(r[key] for key in _RATING_KEYS)) for label, r in zip(labels, ratings))))
    return files






def plot_generated_geometry(grid,output_dir,*,publication=False):
    """实际展成截面与刀尖切换界；只读取网格与边界，不以图形推断接触性能。"""
    if grid.get('units')!='m' or not grid.get('sections') or not grid.get('patches'):
        raise ValueError('展成图需要单位m的五面网格和真实切换边界')
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    plt=_pyplot()
    fig,axes=plt.subplots(1,2,figsize=_size(_DOUBLE_MM,81),constrained_layout=True)
    width_index=len(grid['sections'])//2
    station=grid['sections'][width_index]['station_m']
    profiles=[]
    for name,label,color,linestyle in (
        ('left_flank','左侧有效齿面',_COLORS[0],'-'),('left_fillet','左刀尖根过渡',_COLORS[0],'--'),
        ('root_bottom','槽底连接圆弧',_NEUTRAL,':'),('right_fillet','右刀尖根过渡',_COLORS[1],'--'),
        ('right_flank','右侧有效齿面',_COLORS[1],'-')):
        patch=grid['patches'][name]
        points=[patch['points_m'][v*patch['width_count']+width_index] for v in range(patch['profile_count'])]
        azimuth=np.unwrap([math.atan2(p[1],p[0]) for p in points])
        diameters=[2*math.hypot(p[0],p[1]) for p in points]
        axes[0].plot(np.degrees(azimuth),np.array(diameters)*1000,color=color,linestyle=linestyle,label=label)
        for angle,diameter in zip(azimuth,diameters):profiles.append((name,station,angle,diameter))
    axes[0].set(xlabel='局部方位角 / °',ylabel='截面直径 / mm')
    axes[0].legend(fontsize=7,loc='best')
    _panel(axes[0],'a',f'截面位置 {station*1000:.2f} mm')
    boundary=[]
    for flank,label,color,line in (('left','左侧',_COLORS[0],'-'),('right','右侧',_COLORS[1],'--')):
        distances=[section['station_m'] for section in grid['sections']]
        diameters=[section['transitions'][flank]['last_straight_diameter_m'] for section in grid['sections']]
        axes[1].plot(np.array(distances)*1000,np.array(diameters)*1000,color=color,linestyle=line,label=label)
        boundary.extend((flank,r,d) for r,d in zip(distances,diameters))
    axes[1].set(xlabel='宽向截面位置 / mm',ylabel='直刃与刀尖切换直径 / mm')
    axes[1].legend(loc='best');_panel(axes[1],'b','实际包络区域切换')
    for ax in axes:_grid(ax)
    files=_save(fig,output,'generated_sections',publication=publication,title='显式机床设置的展成侧面与根过渡')
    plt.close(fig)
    files.append(_csv(output/'generated_profiles.csv',('patch','station_m','azimuth_rad','diameter_m'),profiles))
    files.append(_csv(output/'generated_transition.csv',('flank','station_m','transition_diameter_m'),boundary))
    return files
