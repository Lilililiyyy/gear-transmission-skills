"""同一齿轮案例的双语科研图；只读取计算结果，不生成几何或接触解。

输入为bevel_workflow_v1。积分点、面片平均压力和实际CAD分别保持原
含义，缺少数据时拒绝绘图。输出固定物理尺寸、600dpi和SI源数据。
"""
from __future__ import annotations

import csv
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path


from .publication_style import COLORS as _COLORS, COLOUR_POLICY, figure_style, red_blue_colormap
_ARCSEC = 180/math.pi*3600
_NAMES = {
    "outer_pitch_diameter_m": ("外端节圆直径", "Outer pitch diameter", "dₑ", "mm", 1000),
    "mean_pitch_diameter_m": ("中点节圆直径", "Mean pitch diameter", "dₘ", "mm", 1000),
    "inner_pitch_diameter_m": ("内端节圆直径", "Inner pitch diameter", "dᵢ", "mm", 1000),
    "pitch_cone_angle_rad": ("节锥角", "Pitch cone angle", "δ", "deg", 180/math.pi),
    "outer_cone_distance_m": ("外端锥距", "Outer cone distance", "Rₑ", "mm", 1000),
    "mean_cone_distance_m": ("中点锥距", "Mean cone distance", "Rₘ", "mm", 1000),
    "inner_cone_distance_m": ("内端锥距", "Inner cone distance", "Rᵢ", "mm", 1000),
    "face_width_m": ("齿宽", "Face width", "b", "mm", 1000),
    "outer_transverse_module_m": ("外端端面模数", "Outer transverse module", "mₜₑ", "mm", 1000),
    "mean_normal_module_m": ("中点法向模数", "Mean normal module", "mₙₘ", "mm", 1000),
    "spiral_angle_rad": ("中点名义螺旋角", "Nominal mean spiral angle", "βₘ", "deg", 180/math.pi),
    "pressure_angle_rad": ("名义法向压力角", "Nominal normal pressure angle", "αₙ", "deg", 180/math.pi),
    "tip_cone_angle_rad": ("齿顶锥角", "Tip cone angle", "δa", "deg", 180/math.pi),
    "root_cone_angle_rad": ("齿根锥角", "Root cone angle", "δf", "deg", 180/math.pi),
    "mean_spiral_angle_rad": ("中点螺旋角", "Mean spiral angle", "βₘ", "deg", 180/math.pi),
    "mean_normal_tooth_thickness_m": ("中点法向齿厚", "Mean normal tooth thickness", "s_nm", "mm", 1000),
    "shaft_angle_rad": ("轴交角", "Shaft angle", "Σ", "deg", 180/math.pi),
    "tangential": ("名义切向力", "Nominal tangential force", "Fₜ", "N", 1),
    "normal": ("名义法向力", "Nominal normal force", "Fₙ", "N", 1),
    "H0": ("名义接触应力", "Nominal contact stress", "σH0", "MPa", 1e-6),
    "H": ("修正接触应力", "Corrected contact stress", "σH", "MPa", 1e-6),
    "F01": ("小轮名义齿根应力", "Pinion nominal root stress", "σF01", "MPa", 1e-6),
    "F1": ("小轮修正齿根应力", "Pinion corrected root stress", "σF1", "MPa", 1e-6),
    "F02": ("大轮名义齿根应力", "Wheel nominal root stress", "σF02", "MPa", 1e-6),
    "F2": ("大轮修正齿根应力", "Wheel corrected root stress", "σF2", "MPa", 1e-6),
}


def _text(lang, zh, en):
    if lang not in ("zh", "en"):
        raise ValueError("lang须为zh或en")
    return zh if lang == "zh" else en


def _finite(value, field):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{field}须为有限实数；未知或失败值不得补零")
    return float(value)


def _label(field, lang, *, member=None, unit=True):
    key = field.rsplit(".", 1)[-1]
    if key not in _NAMES:
        raise ValueError(f"{field}缺少明确的双语变量定义")
    zh, en, symbol, units, _ = _NAMES[key]
    # 在刻度格式器接管文字之前排好下标，避免其draw过程恢复缺字Unicode。
    symbol = {"dₑ": r"$d_{\mathrm{e}}$", "dₘ": r"$d_{\mathrm{m}}$", "dᵢ": r"$d_{\mathrm{i}}$",
        "Rₑ": r"$R_{\mathrm{e}}$", "Rₘ": r"$R_{\mathrm{m}}$", "Rᵢ": r"$R_{\mathrm{i}}$",
        "mₜₑ": r"$m_{\mathrm{te}}$", "mₙₘ": r"$m_{\mathrm{nm}}$", "βₘ": r"$\beta_{\mathrm{m}}$",
        "αₙ": r"$\alpha_{\mathrm{n}}$", "s_nm": r"$s_{\mathrm{nm}}$",
        "δa": r"$\delta_{\mathrm{a}}$", "δf": r"$\delta_{\mathrm{f}}$"}.get(symbol, symbol)
    prefix = ""
    if member:
        prefix = _text(lang, "小轮" if member == "pinion" else "大轮", "Pinion " if member == "pinion" else "Wheel ")
    return f"{prefix}{_text(lang, zh, en)}, {symbol}"+(f" ({units})" if unit else "")


def _value(geometry, field, member=None):
    source = geometry[member] if member else geometry
    return _finite(source[field], field)


def _validate(bundle):
    if not isinstance(bundle, dict) or bundle.get("schema") != "bevel_workflow_v1":
        raise ValueError("需要schema=bevel_workflow_v1的实际工作流结果")
    return bundle


@contextmanager
def _paper(lang):
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib import font_manager
    except ImportError as error:
        raise ImportError("科研图需要可选Matplotlib依赖") from error
    _text(lang, "", "")
    families = ["DejaVu Sans"]
    if lang == "zh":
        available = {font.name for font in font_manager.fontManager.ttflist}
        families = [name for name in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Noto Sans CJK JP") if name in available]+families
        if len(families) == 1:
            raise RuntimeError("未找到可用中文字体，停止生成缺字图")
    with mpl.rc_context({**figure_style(), "font.family": "sans-serif", "font.sans-serif": families,
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "axes.titlesize": 8, "axes.linewidth": .6,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"}):
        yield plt


def _panel(ax, letter):
    writer = getattr(ax, "text2D", ax.text)
    writer(-.02, 1.04, letter, transform=ax.transAxes, fontsize=9, weight="bold", ha="left", va="bottom")


def _normalize_svg_text(value):
    """规范Matplotlib导出的换行与行尾空白；保留行内文字和路径数据。"""
    value = value.replace("\r\n", "\n").replace("\r", "\n")
    return "\n".join(line.rstrip(" \t") for line in value.split("\n")).rstrip("\n")+"\n"


def _save(fig, output, stem, lang, rows, caption, *, width_mm=183, height_mm=120, extra=None):
    output = Path(output)
    directory, data = output/lang, output/"data"
    directory.mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("图形没有实际源数据，未保存占位图")
    columns = list(dict.fromkeys(key for row in rows for key in row))
    csv_path = data/f"{stem}.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    fig.set_size_inches(width_mm/25.4, height_mm/25.4)
    # 中文字体未必包含Unicode下标；统一用Matplotlib数学字体排物理符号。
    from matplotlib.text import Text
    symbols = {"Δθ₂,₀": r"$\Delta\theta_{2,0}$", "Kθ₂": r"$K_{\theta_2}$", "p̄": r"$\bar{p}$",
        "Fₙ,k/ΣFₙ": r"$F_{\mathrm{n},k}/\sum F_{\mathrm{n}}$", "θ₁": r"$\theta_1$", "θ₂": r"$\theta_2$",
        "mₜₑ": r"$m_{\mathrm{te}}$", "mₙₘ": r"$m_{\mathrm{nm}}$", "dₑ": r"$d_{\mathrm{e}}$",
        "dₘ": r"$d_{\mathrm{m}}$", "dᵢ": r"$d_{\mathrm{i}}$", "Rₑ": r"$R_{\mathrm{e}}$",
        "Rₘ": r"$R_{\mathrm{m}}$", "Rᵢ": r"$R_{\mathrm{i}}$", "βₘ": r"$\beta_{\mathrm{m}}$",
        "αₙ": r"$\alpha_{\mathrm{n}}$", "Fₜ": r"$F_{\mathrm{t}}$", "Fₙ": r"$F_{\mathrm{n}}$",
        "s_nm": r"$s_{\mathrm{nm}}$", "δa": r"$\delta_{\mathrm{a}}$", "δf": r"$\delta_{\mathrm{f}}$",
        "σH0": r"$\sigma_{\mathrm{H0}}$", "σH": r"$\sigma_{\mathrm{H}}$",
        "σF01": r"$\sigma_{\mathrm{F01}}$", "σF02": r"$\sigma_{\mathrm{F02}}$",
        "σF1": r"$\sigma_{\mathrm{F1}}$", "σF2": r"$\sigma_{\mathrm{F2}}$"}
    def converted(label):
        for old, new in symbols.items():
            label = label.replace(old, new)
        return label
    for artist in fig.findobj(match=Text):
        artist.set_text(converted(artist.get_text()))
    # FixedFormatter会在draw时重设刻度文字，须同时转换其原始字符串。
    from matplotlib.ticker import FixedFormatter
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            formatter = axis.get_major_formatter()
            if isinstance(formatter, FixedFormatter):
                axis.set_major_formatter(FixedFormatter([converted(str(value)) for value in formatter.seq]))
    files = []
    for suffix in ("png", "svg", "pdf"):
        target = directory/f"{stem}.{suffix}"
        fig.savefig(target, dpi=600, metadata={"Title": stem, "Description": caption} if suffix == "svg" else
                    {"Title": stem, "Subject": caption} if suffix == "pdf" else {"Description": caption})
        if suffix == "svg":
            target.write_text(_normalize_svg_text(target.read_text(encoding="utf-8")), encoding="utf-8", newline="\n")
        files.append(target.relative_to(output).as_posix())
    record = dict(figure=stem, language=lang, caption=caption, width_mm=width_mm, height_mm=height_mm,
        png_dpi=600, files=files, file_sha256={name:hashlib.sha256((output/name).read_bytes()).hexdigest() for name in files},
        source_csv=csv_path.relative_to(output).as_posix(), source_csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        statistics="none", **{"colour_policy":COLOUR_POLICY, **(extra or {})})
    (directory/f"{stem}.json").write_text(json.dumps(record, ensure_ascii=False, indent=2)+"\n", encoding="utf-8", newline="\n")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return record


def plot_workflow_geometry(bundle, output, *, lang="zh"):
    """真实宏观节锥、明确模数与已有名义力，不在绘图层补算分力。"""
    bundle = _validate(bundle)
    g = bundle["geometry"]
    inputs = g["inputs"]
    forces = bundle.get("nominal_forces") or bundle.get("forces") or (bundle.get("standard_rating") or {}).get("forces_n")
    if not forces:
        raise ValueError("宏观图需要实际名义力输出")
    rows = []
    with _paper(lang) as plt:
        fig = plt.figure(layout="constrained")
        grid = fig.add_gridspec(2, 2, width_ratios=(1.15, 1.))
        ax = fig.add_subplot(grid[:, 0])
        sigma = _finite(inputs["shaft_angle_rad"], "shaft_angle_rad")
        for member, color in zip(("pinion", "wheel"), _COLORS):
            delta = _value(g, "pitch_cone_angle_rad", member)
            axis = (0., 1.) if member == "pinion" else (math.sin(sigma), math.cos(sigma))
            radial = (1., 0.) if member == "pinion" else (-math.cos(sigma), math.sin(sigma))
            outer, inner = _value(g, "outer_cone_distance_m"), _value(g, "inner_cone_distance_m")
            for side in (-1, 1):
                points = [(r*(math.cos(delta)*axis[0]+side*math.sin(delta)*radial[0])*1000,
                           r*(math.cos(delta)*axis[1]+side*math.sin(delta)*radial[1])*1000) for r in (inner, outer)]
                ax.plot(*zip(*points), color=color)
                for index, point in enumerate(points):
                    rows.append(dict(kind="pitch_cone", member=member, side=side, station=index,
                        x_m=point[0]/1000, z_m=point[1]/1000))
            ax.plot([0, 1.1*outer*axis[0]*1000], [0, 1.1*outer*axis[1]*1000], color=color, ls="--", lw=.7,
                    label=_text(lang, "小轮" if member == "pinion" else "大轮", member.capitalize()))
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel(_text(lang, "径向坐标 x (mm)", "Radial coordinate, x (mm)"))
        ax.set_ylabel(_text(lang, "轴向坐标 z (mm)", "Axial coordinate, z (mm)"))
        ax.legend(loc="best")
        _panel(ax, "a")
        table = fig.add_subplot(grid[0, 1])
        table.axis("off")
        fields = [("outer_pitch_diameter_m", member) for member in ("pinion", "wheel")]+[("face_width_m", None), ("outer_transverse_module_m", None), ("mean_normal_module_m", None)]
        for index, (field, member) in enumerate(fields):
            source = inputs if field == "face_width_m" else g
            value = _value(source, field, member)
            scale = _NAMES[field][4]
            table.text(0., .9-index*.18, _label(field, lang, member=member), va="center", fontsize=7)
            table.text(1., .9-index*.18, f"{value*scale:.5g}", ha="right", va="center", fontsize=8)
            rows.append(dict(kind="parameter", member=member or "pair", field=field, value_si=value))
        _panel(table, "b")
        force_ax = fig.add_subplot(grid[1, 1])
        if "pinion_tangential_n" in forces:
            components = (("tangential", "切向力", "Tangential force", r"$F_\mathrm{t}$"),
                          ("radial", "径向力", "Radial force", r"$F_\mathrm{r}$"),
                          ("axial", "轴向力", "Axial force", r"$F_\mathrm{a}$"))
            for member, shift, color in (("pinion", -.17, _COLORS[0]), ("wheel", .17, _COLORS[1])):
                values = [_finite(forces[f"{member}_{name}_n"], f"{member}_{name}_n") for name, _, _, _ in components]
                force_ax.bar([i+shift for i in range(3)], values, color=color, width=.32,
                    label=_text(lang, "小轮" if member == "pinion" else "大轮", member.capitalize()))
                rows.extend(dict(kind="nominal_force", member=member, field=name, value_si=value)
                            for (name, _, _, _), value in zip(components, values))
            force_ax.set_xticks(range(3), [_text(lang, zh, en)+"\n"+symbol for _, zh, en, symbol in components])
            force_ax.legend(loc="best", ncol=2)
        else:
            force_keys = [key for key in ("tangential", "normal") if key in forces]
            if not force_keys:
                raise ValueError("名义力数据需要明确具名分量，不能猜测")
            values = [_finite(forces[key], key) for key in force_keys]
            force_ax.bar(range(len(values)), values, color=_COLORS[:len(values)], width=.55)
            force_ax.set_xticks(range(len(values)), [_label(key, lang, unit=False).replace(", ", "\n") for key in force_keys])
            rows.extend(dict(kind="nominal_force", field=key, value_si=value) for key, value in zip(force_keys, values))
        force_ax.set_ylabel(_text(lang, "名义局部力 F (N)", "Nominal local force, F (N)"))
        _panel(force_ax, "c")
        return _save(fig, output, "01_geometry_forces", lang, rows,
            _text(lang, "节锥与名义力来自同一案例；外端端面模数与中点法向模数分别定义，力不是节点载荷分配。",
                "Pitch cones and nominal forces use the same case. Outer transverse and mean normal modules are distinct; forces are not a nodal load distribution."))


def _benchmark_rows(bundle, kind):
    if bundle["masta_benchmark"].get("same_case_inputs") is False:
        raise ValueError("对标案例身份不一致，不能生成同工况比较图")
    rows = bundle["masta_benchmark"][kind+"_rows"]
    if not rows:
        raise ValueError("缺少同案例参照字段，未用自身值补参照")
    result = []
    for row in rows:
        record = dict(row)
        record["field"] = row.get("field", row.get("result_path"))
        if not record["field"]:
            raise ValueError("对标量须有具名原始字段")
        record["reference"] = _finite(row["reference"], "reference")
        record["project"] = _finite(row.get("project", row.get("calculated")), "project/calculated")
        record["signed_difference"] = record["project"]-record["reference"]
        record["relative_difference"] = None if record["reference"] == 0 else record["signed_difference"]/abs(record["reference"])
        result.append(record)
    return result


def _benchmark_label(row, lang):
    if row.get("label_"+lang):
        return row["label_"+lang]
    field = row["field"]
    member = "pinion" if "pinion." in field else "wheel" if "wheel." in field else None
    label = _label(field, lang, member=member, unit=False)
    if lang == "en":
        import textwrap
        label = "\n".join(textwrap.wrap(label, width=31, break_long_words=False))
    return label


def plot_workflow_benchmark(bundle, output, *, kind="geometry", lang="zh"):
    """同定义数据配准对照，保留实际非零差异；不涉及承载参照。"""
    _validate(bundle)
    source_rows = _benchmark_rows(bundle, kind)
    # 正文图显示代表性宏观量与六个应力；全部已比较字段仍写入源CSV。
    if kind == "geometry":
        primary = {"geometry.outer_cone_distance_m", "geometry.mean_cone_distance_m",
                   "geometry.outer_transverse_module_m", "geometry.mean_normal_module_m"}
        primary.update(f"geometry.{member}.{field}" for member in ("pinion", "wheel")
                       for field in ("outer_pitch_diameter_m", "mean_normal_tooth_thickness_m",
                                     "pitch_cone_angle_rad", "tip_cone_angle_rad", "root_cone_angle_rad",
                                     "mean_spiral_angle_rad"))
        rows = [row for row in source_rows if row["field"] in primary]
    elif kind == "strength":
        rows = [row for row in source_rows if row["field"].startswith("standard_rating.stress_pa.")]
    else:
        raise ValueError("kind须为geometry或strength")
    if not rows:
        raise ValueError("缺少正文图所需的明确宏观量或应力")
    for row in source_rows:
        row["plotted"] = row in rows
    units = {"m": ("mm", 1000), "rad": ("deg", 180/math.pi), "Pa": ("MPa", 1e-6), "1": ("1", 1)}
    if any(row.get("unit") not in units for row in rows):
        raise ValueError("参照字段缺少明确支持的原单位")
    with _paper(lang) as plt:
        if kind == "geometry":
            groups = list(dict.fromkeys(row["unit"] for row in rows))
            fig, axes = plt.subplots(1, len(groups), squeeze=False, layout="constrained")
            for column, unit in enumerate(groups):
                ax = axes[0, column]
                selected = [row for row in rows if row["unit"] == unit]
                shown, scale = units[unit]
                for index, row in enumerate(selected):
                    ax.plot([row["reference"]*scale, row["project"]*scale], [index]*2, color="#777777", lw=.7)
                ax.scatter([row["reference"]*scale for row in selected], range(len(selected)),
                    facecolors="none", edgecolors=_COLORS[0], s=25, label="MASTA")
                ax.scatter([row["project"]*scale for row in selected], range(len(selected)), marker="x", s=18,
                    color=_COLORS[1], label=_text(lang, "本模块", "This module"))
                ax.set_yticks(range(len(selected)), [_benchmark_label(row, lang) for row in selected])
                ax.invert_yaxis()
                ax.set_xlabel(_text(lang, f"同定义数值 ({shown})", f"Matched quantity ({shown})"))
                if column == 0:
                    ax.legend(loc="lower right")
                _panel(ax, chr(97+column))
            stem, height = "03_geometry_benchmark", max(90, min(170, 13*max(sum(row["unit"] == unit for row in rows) for unit in groups)))
        elif kind == "strength":
            if any(row["unit"] != "Pa" for row in rows):
                raise ValueError("强度对照仅接受同定义Pa应力，不混容量/寿命")
            assessment = bundle.get("strength_assessment")
            safety = assessment.get("rows", []) if assessment else []
            fig, axes = plt.subplots(3 if safety else 2, 1, layout="constrained",
                gridspec_kw={"height_ratios": [2., 1., 1.3] if safety else [2., 1.]})
            x = list(range(len(rows)))
            axes[0].bar([v-.17 for v in x], [r["reference"]*1e-6 for r in rows], width=.32, color=_COLORS[0], label="MASTA")
            axes[0].bar([v+.17 for v in x], [r["project"]*1e-6 for r in rows], width=.32, color=_COLORS[1], label=_text(lang, "本模块", "This module"))
            axes[0].set_ylabel(_text(lang, "名义／修正应力 σ (MPa)", "Nominal/corrected stress, σ (MPa)"))
            axes[0].legend(ncol=2)
            axes[1].plot(x, [math.nan if r["relative_difference"] is None else r["relative_difference"]*1e6 for r in rows], "o", color=_COLORS[0], ms=3)
            axes[1].axhline(0., color="#777777", lw=.6)
            axes[1].set_ylabel(_text(lang, "应力相对差\nΔσ/|σref| (ppm)", "Relative stress difference\nΔσ/|σref| (ppm)"))
            labels = []
            for row in rows:
                key = row["field"].rsplit(".", 1)[-1]
                zh, en, symbol, _, _ = _NAMES[key]
                if lang == "en":
                    import textwrap
                    name = "\n".join(textwrap.wrap(en, width=18, break_long_words=False))
                else:
                    name = zh
                labels.append(name+"\n"+r"$\sigma_{\mathrm{"+key+r"}}$")
            axes[0].set_xticks(x, labels)
            axes[1].set_xticks(x, [r"$\sigma_{\mathrm{"+r["field"].rsplit(".", 1)[-1]+r"}}$" for r in rows])
            if safety:
                labels, factors, targets = [], [], []
                for entry in safety:
                    quantity = entry["quantity"]
                    if quantity not in ("H", "F1", "F2"):
                        raise ValueError("设计安全系数须对应明确H/F1/F2应力")
                    stress = _finite(entry["calculated_stress_pa"], "calculated_stress_pa")
                    limit = _finite(entry["design_stress_limit_pa"], "design_stress_limit_pa")
                    factor = _finite(entry["safety_factor"], "safety_factor")
                    target = _finite(entry["minimum_safety"], "minimum_safety")
                    if min(stress, limit, factor, target) <= 0 or not math.isclose(factor, limit/stress, rel_tol=1e-8):
                        raise ValueError("给定限值与设计安全系数不一致")
                    labels.append(_text(lang, {"H": "接触", "F1": "小轮齿根", "F2": "大轮齿根"}[quantity],
                        {"H": "Contact", "F1": "Pinion root", "F2": "Wheel root"}[quantity]))
                    factors.append(factor)
                    targets.append(target)
                    source_rows.append(dict(kind="specified_design_limit", plotted=True, **entry))
                axes[2].bar(range(len(safety)), factors, color=_COLORS[0], width=.55,
                            label=_text(lang, "给定限值下安全系数", "Safety with specified limits"))
                axes[2].scatter(range(len(safety)), targets, marker="_", color="#333333", s=130,
                                label=_text(lang, "最低目标", "Minimum target"))
                axes[2].set_xticks(range(len(safety)), labels)
                axes[2].set_ylabel(_text(lang, "设计安全系数 S (1)", "Design safety factor, S (1)"))
                axes[2].set_ylim(0, max(factors+targets)*1.32)
                axes[2].legend(ncol=2, loc="upper left")
            for ax, letter in zip(axes, "abc"):
                _panel(ax, letter)
            stem, height = "04_strength_benchmark", 175 if safety else 130
        else:
            raise ValueError("kind须为geometry或strength")
        return _save(fig, output, stem, lang, source_rows,
            _text(lang, "同一案例、共同定义字段对照；原始差异保留。安全系数若显示则依据案例给定设计限值，非自动2014容量或寿命。",
                "Matched fields for the same case retain nonzero differences. Any safety factors use specified design limits, not an automatically derived 2014 capacity or life."), height_mm=height)


def _mesh(entry, data_root):
    if isinstance(entry, dict) and "vertices_m" in entry:
        result = entry
    else:
        specification = entry if isinstance(entry, dict) else {"path": entry}
        path = Path(specification["path"])
        if not path.is_absolute():
            if data_root is None:
                raise ValueError("相对CAD数据路径须显式给data_root")
            path = Path(data_root)/path
        if path.suffix.lower() == ".json":
            result = json.loads(path.read_text(encoding="utf-8-sig"))
        elif path.suffix.lower() == ".vtp":
            from .cad_preview import _vtk
            vtk = _vtk()
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(str(path))
            reader.Update()
            poly = reader.GetOutput()
            field = poly.GetFieldData().GetAbstractArray("geometry_metadata_json")
            metadata = json.loads(field.GetValue(0)) if field and field.GetNumberOfValues() else {}
            if metadata.get("length_unit", specification.get("units")) != "m":
                raise ValueError("VTP须有明确m单位，不按坐标量级猜单位")
            faces = poly.GetCellData().GetArray("face_indices")
            triangles = []
            for index in range(poly.GetNumberOfCells()):
                cell = poly.GetCell(index)
                if cell.GetNumberOfPoints() != 3:
                    raise ValueError("科研三维图只接收实际三角网格")
                triangles.append([cell.GetPointId(j) for j in range(3)])
            result = dict(units="m", vertices_m=[poly.GetPoint(i) for i in range(poly.GetNumberOfPoints())],
                triangles=triangles, face_indices=[int(faces.GetValue(i)) for i in range(poly.GetNumberOfCells())] if faces else None,
                metadata={"closed_solid": bool(metadata.get("closed_solid")), **metadata.get("geometry_metadata", {})})
            if result["face_indices"] is None:
                result.pop("face_indices")
        else:
            raise ValueError("绘图层只读SI网格JSON或VTP；不在此生成CAD")
    if result.get("units") != "m" or not result.get("vertices_m") or not result.get("triangles"):
        raise ValueError("CAD图需要明确m单位的非空实际网格")
    return result


def plot_workflow_cad(bundle, output, *, lang="zh", data_root=None, detail_layout="standard"):
    """实际CAD图版；contact布局突出单侧工作面并将根区图例放图外。

    standard保留既有布局。SVG/PDF注释为矢量，CAD为真实网格渲染；
    白边裁剪只用于contact图像排版，不缩放或改变实际几何。
    """
    _validate(bundle)
    if detail_layout not in ("standard", "contact"):
        raise ValueError("detail_layout须为standard或contact")
    from .cad_preview import render_cad_mesh
    cad = bundle["cad"]
    straight = bundle.get("kind", bundle.get("case", {}).get("kind")) == "straight"
    entries = [("pair_mesh", "双轮材料实体", "Gear-pair material solids"),
               ("flank_mesh", "球面渐开线工作面" if straight else "实际展成工作侧面",
                "Spherical-involute working flanks" if straight else "Generated working flanks"),
               ("root_mesh", "工程根连接" if straight else "刀尖根过渡与连接槽底",
                "Engineering root connection" if straight else "Cutter-tip transition and connected root bottom")]
    images, rows = [], []
    palette = {"left_flank": ((33, 102, 172), "左工作面", "Left working flank"),
        "right_flank": ((103, 169, 207), "右工作面", "Right working flank"),
        "left_fillet": ((178, 24, 43), "左刀尖根过渡", "Left cutter-tip transition"),
        "right_fillet": ((214, 96, 77), "右刀尖根过渡", "Right cutter-tip transition"),
        "root_bottom": ((97, 105, 112), "根圆连接槽底", "Root-circle connecting bottom")}
    if detail_layout == "contact":
        if straight:
            raise ValueError("contact精选布局需要实际弧齿制造选区")
        palette["adjacent_working_flanks"] = ((190, 198, 204), "邻接工作齿面（根区参照）", "Adjacent working flanks (root context)")
    temporary = Path(output)/"cad"
    temporary.mkdir(parents=True, exist_ok=True)
    for role, zh, en in entries:
        if role not in cad:
            raise ValueError(f"几何图版缺少实际{role}，未用其他工况或占位图代替")
        mesh = _mesh(cad[role], data_root)
        regions = mesh.get("metadata", {}).get("face_regions", {})
        colors, legend = None, []
        if regions:
            by_face = {}
            for name, faces in regions.items():
                if name not in palette:
                    continue
                rgb, zh_name, en_name = palette[name]
                legend.append((rgb, _text(lang, zh_name, en_name)))
                for face in faces:
                    if face in by_face:
                        raise ValueError("实际几何区域的面号不可重叠")
                    by_face[face] = rgb
            if legend:
                if "face_indices" not in mesh:
                    raise ValueError("区域着色须保留真实BRep面号")
                colors = [by_face.get(face, (170, 180, 188)) for face in mesh["face_indices"]]
        elif role in ("flank_mesh", "root_mesh"):
            # 没有子区面号时只按实际选区整体配色，不猜测STEP面序或刀具子区。
            rgb = (33, 102, 172) if role == "flank_mesh" else (178, 24, 43)
            colors = [rgb]*len(mesh["triangles"])
        target = temporary/f"{role}_{lang}.png"
        # 长时间计算进程可能已加载旧版预览函数；此时仅回退为原几何色。
        import inspect
        color_options = {"cell_colors_rgb": colors} if "cell_colors_rgb" in inspect.signature(render_cad_mesh).parameters else {}
        if not color_options:
            legend = []
        render_size = (3200, 1500) if role == "pair_mesh" else (2800, 2100)
        render_cad_mesh(mesh, target, size=render_size, lang=lang, annotate=False, **color_options)
        images.append((target, _text(lang, zh, en), legend))
        rows.append(dict(role=role, vertices=len(mesh["vertices_m"]), triangles=len(mesh["triangles"]),
            units="m", closed_solid=bool(mesh.get("metadata", {}).get("closed_solid")),
            face_regions_json=json.dumps(regions, sort_keys=True)))
    with _paper(lang) as plt:
        if detail_layout == "contact":
            # 复用真实渲染，仅裁掉白底边距；不拉伸或补造物理曲面。
            import numpy as np
            from PIL import Image
            cropped = []
            for target, title, legend in images:
                pixels = np.asarray(Image.open(target).convert("RGB"))
                ys, xs = np.where(np.any(pixels != 255, axis=2))
                if not len(xs):
                    raise ValueError("实际CAD渲染为空白，拒绝输出占位图")
                pad = max(6, int(.035*max(xs.max()-xs.min(), ys.max()-ys.min())))
                cropped.append(pixels[max(0, ys.min()-pad):min(len(pixels), ys.max()+pad+1),
                    max(0, xs.min()-pad):min(pixels.shape[1], xs.max()+pad+1)])
            fig = plt.figure()
            axes = [fig.add_axes([.045, .44, .91, .46]), fig.add_axes([.055, .11, .39, .29]),
                    fig.add_axes([.555, .11, .39, .29])]
            requests = bundle.get("case", {}).get("manufacturing", {})
            counts = [requests.get(member, {}).get("request", {}).get("tooth_count") for member in ("pinion", "wheel")]
            pair_title = _text(lang, "双轮真实材料实体", "Actual gear-pair solids")
            if all(isinstance(count, int) and count > 0 for count in counts):
                pair_title += f": {counts[0]}/{counts[1]}"
            titles = [pair_title, _text(lang, "小轮左工作齿面（单侧实际面）", "Pinion left working flank\n(actual selected surface)"),
                      _text(lang, "小轮根部及邻接工作齿面", "Pinion root region and\nadjacent working flanks")]
            for ax, pixels in zip(axes, cropped):
                ax.imshow(pixels); ax.axis("off")
            for (x, y), title, letter in zip(((.045, .928), (.055, .432), (.555, .432)), titles, "abc"):
                fig.text(x, y, letter, fontsize=9, weight="bold", va="top")
                fig.text(x+.025, y, title, fontsize=8, va="top")
            details = []
            for member in ("pinion", "wheel"):
                request = requests.get(member, {}).get("request", {})
                relief = request.get("tip_relief_m", 0.)*1e6
                start = request.get("tip_relief_start", .65)
                label = _text(lang, "小轮" if member == "pinion" else "大轮", "Pinion" if member == "pinion" else "Wheel")
                details.append(f"{label}: {relief:g} µm, v = {start:g}")
            fig.suptitle(_text(lang, "实际法向齿顶修缘；"+"；".join(details)+"\nv：归一化有效齿廓位置（齿根0 → 齿顶1）",
                "Actual normal tip relief; "+"; ".join(details)+"\nv: normalized active-profile position (root 0 → tip 1)"), y=.982, fontsize=8)
            from matplotlib.patches import Patch
            legends = list(dict.fromkeys(item for _, _, legend in images for item in legend))
            fig.legend(handles=[Patch(facecolor=tuple(channel/255 for channel in rgb), label=name) for rgb, name in legends],
                loc="lower center", bbox_to_anchor=(.5, .015), ncol=2, fontsize=6.5, columnspacing=1.7, handlelength=1.7)
            caption = _text(lang, "同一显式制造案例的实际闭实体、单侧工作面及邻面根区。蓝色为左工作面，红色为刀尖根过渡，灰色为连接槽底或邻面；均不是压力图，槽底连接不代表完整刀具包络。",
                "Actual solids, one working flank and the root region from the same explicit manufacturing case. Blue marks the left working flank, red the cutter-tip transitions, and gray the root connector or adjacent flanks. These are not pressure maps; the root connector is not a complete cutter envelope.")
            return _save(fig, output, "02_geometry_details", lang, rows, caption, height_mm=168,
                extra={"cad_rendering": "actual mesh; embedded raster with editable vector annotations", "detail_layout": "contact", "legend_outside_models": True})
        fig = plt.figure(layout="constrained")
        grid = fig.add_gridspec(2, 2, height_ratios=(1.6, 1.))
        axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
        for ax, (target, title, legend), letter in zip(axes, images, "abc"):
            ax.imshow(plt.imread(target))
            ax.axis("off")
            ax.set_title(title, loc="left", x=.045)
            _panel(ax, letter)
            if legend:
                from matplotlib.patches import Patch
                ax.legend(handles=[Patch(facecolor=tuple(channel/255 for channel in rgb), label=name)
                                   for rgb, name in legend], loc="lower left", frameon=False, fontsize=7)
        caption = (_text(lang, "同案例理想直齿CAD：球面渐开线工作面及工程根连接，根区不是实际切削圆角。局部面图不是压力图。",
            "Same-case ideal straight-gear CAD: spherical-involute flanks and an engineering root connection, not an actual cutting fillet. Detail surfaces are not pressure maps.") if straight else
            _text(lang, "图像来自同案例实际CAD。刀尖根过渡与根圆连接槽底分别保留；局部面图不是压力或完整刀具包络。",
                "Images use the actual CAD for the same case. Cutter-tip transitions and the root-circle connector remain distinct; detail surfaces are not pressure maps or a complete cutter envelope."))
        return _save(fig, output, "02_geometry_details", lang, rows, caption, height_mm=155,
            extra={"cad_rendering": "actual mesh; embedded raster with editable vector annotations"})


def _accepted(result):
    return isinstance(result, dict) and result.get("converged") is True


def _phase(result):
    return _finite(result["driving_angle_rad"], "driving_angle_rad")


def plot_workflow_tca(bundle, output, *, lang="zh"):
    """只绘实际接受相切点；失败姿态保留缺值，不作为零TE。"""
    _validate(bundle)
    dataset = bundle["tca"]
    groups = dataset.get("rows")
    results = dataset.get("results", []) if groups is None else [result for group in groups for result in group["results"]]
    if not results or not any(_accepted(result) for result in results):
        raise ValueError("无实际接受相切点，未生成成功TCA图")
    direction = dataset.get("closure_direction")
    if isinstance(direction, bool) or direction not in (-1, 1):
        raise ValueError("TCA包络须显式给closure_direction=-1或1，不能强制参考齿对0")
    rows = []
    for result in results:
        row = dict(driver_rad=_phase(result), tooth_offset=int(result.get("tooth_offset", 0)),
                   converged=_accepted(result), status=result.get("status"), selected_first_touch=False)
        if _accepted(result):
            row["unloaded_te_rad"] = _finite(result["unloaded_te_rad"], "unloaded_te_rad")
            row["follower_increment_rad"] = _finite(result["follower_increment_angle_rad"], "follower_increment_angle_rad")
            for axis, value in zip("xyz", result["pinion_contact_point_world_m"]):
                row["contact_"+axis+"_m"] = _finite(value, "contact_point")
            row["position_residual_m"] = _finite(result["position_residual_m"], "position_residual_m")
            if "normal_cross_residual" in result:
                row["normal_cross_residual"] = _finite(result["normal_cross_residual"], "normal_cross_residual")
                row["normal_residual_definition"] = "norm(cross(n1,n2)); opposite direction checked separately"
            elif "normal_alignment_residual" in result:
                row["normal_alignment_residual"] = _finite(result["normal_alignment_residual"], "normal_alignment_residual")
                row["normal_residual_definition"] = "norm(n1+n2)"
            else:
                raise ValueError("实际相切结果缺少明确法向残差及定义")
        rows.append(row)
    envelope = []
    selected_by_phase = {group["driving_angle_rad"]: group["selected"] for group in groups
                         if "selected" in group} if groups is not None else {}
    for phase in sorted({row["driver_rad"] for row in rows}):
        candidates = [row for row in rows if row["driver_rad"] == phase and row["converged"]]
        if phase in selected_by_phase:
            selected = selected_by_phase[phase]
            matches = [row for row in candidates if selected is not None and row["tooth_offset"] == selected["tooth_offset"]]
            if selected is None:
                chosen = None
            elif len(matches) != 1 or selected.get("converged") is not True:
                raise ValueError("业务层首接触选择与有效物理候选不一致")
            else:
                chosen = matches[0]
        elif candidates:
            chosen = min(candidates, key=lambda row: direction*row["follower_increment_rad"])
        else:
            chosen = None
        if chosen is not None:
            chosen["selected_first_touch"] = True
            envelope.append(chosen)
        else:
            envelope.append(dict(driver_rad=phase, converged=False))
    with _paper(lang) as plt:
        fig = plt.figure(layout="constrained")
        grid = fig.add_gridspec(2, 2, width_ratios=(1., 1.1))
        first = fig.add_subplot(grid[:, 0])
        successful = [row for row in envelope if row["converged"]]
        points = first.scatter([r["contact_x_m"]*1000 for r in successful], [r["contact_z_m"]*1000 for r in successful],
                              c=[r["driver_rad"]*180/math.pi for r in successful], cmap=red_blue_colormap(), s=14)
        first.set_xlabel(_text(lang, "世界坐标 x (mm)", "World coordinate, x (mm)"))
        first.set_ylabel(_text(lang, "世界坐标 z (mm)", "World coordinate, z (mm)"))
        first.set_aspect("equal", adjustable="box")
        fig.colorbar(points, ax=first, shrink=.6, pad=.025, label=_text(lang, "小轮增量角 θ₁ (deg)", "Pinion increment angle, θ₁ (deg)"))
        second = fig.add_subplot(grid[0, 1])
        third = fig.add_subplot(grid[1, 1], sharex=second)
        offsets = sorted({row["tooth_offset"] for row in rows})
        for number, offset in enumerate(offsets):
            ordered = sorted((row for row in rows if row["tooth_offset"] == offset), key=lambda row: row["driver_rad"])
            third.plot([r["driver_rad"]*180/math.pi for r in ordered],
                        [r.get("unloaded_te_rad", math.nan)*_ARCSEC for r in ordered], color="#A6ABB1", lw=.65,
                        ls=(":", "--", "-.")[number % 3], label=f"k={offset}")
        for ax in (second, third):
            ax.plot([r["driver_rad"]*180/math.pi for r in envelope],
                        [r.get("unloaded_te_rad", math.nan)*_ARCSEC for r in envelope], "o-", color=_COLORS[0], ms=2.6,
                        label=_text(lang, "首接触包络", "First-touch envelope") if ax is second else None)
        second.legend(loc="lower left")
        third.legend(loc="best", ncol=3)
        third.set_xlabel(_text(lang, "小轮增量角 θ₁ (deg)", "Pinion increment angle, θ₁ (deg)"))
        second.set_ylabel(_text(lang, "未加载传动误差\nΔθ₂,₀ (arcsec)", "Unloaded transmission error\nΔθ₂,₀ (arcsec)"))
        third.set_ylabel(_text(lang, "齿对候选传动误差\nΔθ₂,₀ (arcsec)", "Tooth-pair transmission error\nΔθ₂,₀ (arcsec)"))
        for ax, letter in ((first, "a"), (second, "b"), (third, "c")):
            _panel(ax, letter)
        return _save(fig, output, "05_unloaded_tca", lang, rows,
            _text(lang, "灰线为各齿对有效相切解，突出曲线为给定闭合方向的首接触包络；点为实际相切位置的世界xz投影。连线只连接计算姿态，保留原始TE零点。",
                "Grey lines show valid tangencies of individual tooth pairs; the highlighted curve is the first-touch envelope for the stated closure direction. Points are computed tangencies projected onto world x-z. Lines join solved poses and retain the original transmission-error zero."), height_mm=125,
            extra={"closure_direction": direction, "first_touch_selection": "min(closure_direction*follower_increment_rad)"})


def _loaded_results(bundle):
    results = bundle["ltca"]["results"]
    if not results:
        raise ValueError("缺少实际承载结果")
    for result in results:
        _phase(result)
        if _accepted(result) and result.get("equilibrium", {}).get("converged") is not True:
            raise ValueError("承载成功状态与互补平衡状态不一致")
    phases = [_phase(result) for result in results]
    if len(set(phases)) != len(phases):
        raise ValueError("周期结果的主动角不可重复；载荷扫描须另存load_sweep")
    torques = [_finite(result["pinion_torque_nm"], "pinion_torque_nm") for result in results]
    if any(not math.isclose(value, torques[0], rel_tol=1e-12, abs_tol=1e-12) for value in torques):
        raise ValueError("周期图只能采用同一有符号输入扭矩，不混载荷扫描")
    return sorted(results, key=_phase)


def _contact_values(result, representation):
    if not _accepted(result) or result.get("equilibrium", {}).get("converged") is not True:
        raise ValueError("拒绝将失败或试探接触数据画成物理压力图")
    if representation not in ("facets", "integration_points"):
        raise ValueError("须明确contact_representation=facets或integration_points")
    contacts = result.get("contacts", [])
    if not contacts:
        raise ValueError("没有实际接触数据，未构造椭圆或插值边界")
    rows = []
    for index, contact in enumerate(contacts):
        row = dict(contact_index=index, active=contact["active"], tooth_offset=contact["tooth_offset"])
        if not isinstance(row["active"], bool):
            raise ValueError("接触活跃状态须为明确布尔值")
        pressure = _finite(contact["mean_pressure_pa" if representation == "facets" else "pressure_pa"], "pressure_pa")
        area = _finite(contact["area_m2" if representation == "facets" else "weight_m2"], "physical_area_or_quadrature_weight")
        force = _finite(contact["force_n"], "force_n")
        if area <= 0 or pressure < 0 or force < 0 or not math.isclose(pressure*area, force, rel_tol=1e-7, abs_tol=1e-8):
            raise ValueError("压力、真实面积/积分权重和法向力不一致")
        point = contact["first_point_world_m"] if representation == "facets" else contact["point_world_m"]
        normal = contact["normal_world"]
        if len(point) != 3 or len(normal) != 3:
            raise ValueError("接触点和法向须为三维")
        row.update(pressure_pa=pressure, physical_area_or_weight_m2=area, force_n=force)
        for axis, value, direction in zip("xyz", point, normal):
            row["point_"+axis+"_m"] = _finite(value, "point_world_m")
            row["normal_"+axis] = _finite(direction, "normal_world")
        if representation == "facets":
            vertices = contact.get("polygon_world_m", contact.get("first_vertices_world_m"))
            if not vertices or len(vertices) < 3 or any(len(point) != 3 for point in vertices):
                raise ValueError("面片压力需要真实边界顶点，不能用半径/椭圆代替")
            row["polygon_world_m_json"] = json.dumps([[_finite(x, "facet_vertex_m") for x in point] for point in vertices])
        rows.append(row)
    if not any(row["active"] and row["force_n"] > 0 for row in rows):
        raise ValueError("没有实际承载接触点或面片")
    return rows


def plot_workflow_contact(bundle, output, *, lang="zh", phase_index=None):
    """真实面片或积分点的固定切平面投影，不补造连续接触斑。"""
    _validate(bundle)
    import numpy as np
    results = _loaded_results(bundle)
    if phase_index is None:
        valid = [result for result in results if _accepted(result)]
        if not valid:
            raise ValueError("没有接受承载姿态")
        result = min(valid, key=lambda result: abs(_phase(result)))
    else:
        result = results[phase_index]
    representation = result.get("contact_representation", bundle["ltca"].get("contact_representation"))
    rows = _contact_values(result, representation)
    active = [row for row in rows if row["active"]]
    forces = np.array([row["force_n"] for row in active])
    positions = np.array([[row["point_"+axis+"_m"] for axis in "xyz"] for row in active])
    normals = np.array([[row["normal_"+axis] for axis in "xyz"] for row in active])
    origin = np.average(positions, axis=0, weights=forces)
    normal = np.average(normals, axis=0, weights=forces)
    normal /= np.linalg.norm(normal)
    reference = np.array((0., 0., 1.)) if abs(normal[2]) < .9 else np.array((1., 0., 0.))
    tangent_u = np.cross(reference, normal)
    tangent_u /= np.linalg.norm(tangent_u)
    tangent_v = np.cross(normal, tangent_u)
    def project(points):
        array = np.asarray(points)-origin
        return np.column_stack((array@tangent_u, array@tangent_v))*1000
    for row in rows:
        position = project([[row["point_"+axis+"_m"] for axis in "xyz"]])[0]
        row["projection_x_mm"], row["projection_y_mm"] = map(float, position)
    with _paper(lang) as plt:
        from matplotlib.collections import PolyCollection
        from matplotlib.colors import Normalize
        fig, axes = plt.subplots(1, 2, layout="constrained", gridspec_kw={"width_ratios": [1.5, 1.]})
        pressures = [row["pressure_pa"]*1e-6 for row in active]
        from .halfspace_plotting import pressure_colour_scale
        norm,pressure_cmap,_ = pressure_colour_scale([value*1e6 for value in pressures])
        if representation == "facets":
            polygons = [project(json.loads(row["polygon_world_m_json"])) for row in active]
            background = PolyCollection([project(json.loads(row["polygon_world_m_json"])) for row in rows],
                facecolors="#E8EAED", edgecolors="#BAC1C6", linewidths=.18, zorder=0)
            axes[0].add_collection(background)
            artist = PolyCollection(polygons, array=np.asarray(pressures), cmap=pressure_cmap, norm=norm,
                                    edgecolors="#555555", linewidths=.3, zorder=1)
            axes[0].add_collection(artist)
            bounds = np.concatenate(polygons)
            margin = max(float(np.ptp(bounds[:, 0])), float(np.ptp(bounds[:, 1])))*.12
            axes[0].set_xlim(float(bounds[:, 0].min())-margin, float(bounds[:, 0].max())+margin)
            axes[0].set_ylim(float(bounds[:, 1].min())-margin, float(bounds[:, 1].max())+margin)
            pressure_label = _text(lang, "面片平均接触压力 p̄ (MPa)", "Facet-mean contact pressure, p̄ (MPa)")
        else:
            artist = axes[0].scatter([r["projection_x_mm"] for r in active], [r["projection_y_mm"] for r in active],
                c=pressures, cmap=pressure_cmap, norm=norm, s=9, linewidths=0)
            pressure_label = _text(lang, "接触积分点压力 p (MPa)", "Contact quadrature-point pressure, p (MPa)")
        axes[0].set_aspect("equal", adjustable="datalim")
        axes[0].set_xlabel(_text(lang, "固定切平面坐标 ξ (mm)", "Fixed tangent-plane coordinate, ξ (mm)"))
        axes[0].set_ylabel(_text(lang, "固定切平面坐标 η (mm)", "Fixed tangent-plane coordinate, η (mm)"))
        fig.colorbar(artist, ax=axes[0], shrink=.75, label=pressure_label)
        pairs = sorted({row["tooth_offset"] for row in rows})
        totals = [sum(row["force_n"] for row in rows if row["tooth_offset"] == pair) for pair in pairs]
        axes[1].bar(range(len(pairs)), totals, color=[_COLORS[i % len(_COLORS)] for i in range(len(pairs))], width=.6)
        axes[1].set_xticks(range(len(pairs)), [f"k={pair}" for pair in pairs])
        axes[1].set_xlabel(_text(lang, "实际齿对偏移 k (1)", "Actual tooth-pair offset, k (1)"))
        axes[1].set_ylabel(_text(lang, "齿对法向合力 Fₙ (N)", "Tooth-pair normal force, Fₙ (N)"))
        for ax, letter in zip(axes, "ab"):
            _panel(ax, letter)
        return _save(fig, output, "06_loaded_contact", lang, rows,
            _text(lang, "真实接触数据的固定切平面投影；灰色为实际候选FE面片，着色为承载面片。投影不重新定义物理面积，未插值或补造椭圆边界。",
                "Actual contact data are projected onto a fixed tangent plane. Grey shows actual candidate FE facets; color shows loaded facets. Projection does not redefine physical area; no interpolation or elliptical boundary is added."), height_mm=110,
            extra={"contact_representation": representation, "driving_angle_rad": _phase(result),
                   "projection_origin_m": origin.tolist(), "projection_u": tangent_u.tolist(), "projection_v": tangent_v.tolist()})


def _result_scalar(result, *fields):
    for field in fields:
        if field in result:
            return _finite(result[field], field)
    raise ValueError(f"实际结果缺少{fields}，未推导或补零")


def plot_workflow_cycle(bundle, output, *, lang="zh"):
    """实际计算姿态上的TE分解与从动轴增量扭转刚度；失败点断开。"""
    _validate(bundle)
    rows = []
    for result in _loaded_results(bundle):
        row = dict(driver_rad=_phase(result), converged=_accepted(result), status=result.get("status"))
        if row["converged"]:
            row.update(unloaded_te_rad=_result_scalar(result, "actual_unloaded_te_rad", "unloaded_tca_te_rad"),
                mesh_te_rad=_result_scalar(result, "mesh_correction_rad", "mesh_unloaded_correction_rad"),
                elastic_te_rad=_result_scalar(result, "elastic_te_rad", "elastic_follower_lag_rad"),
                total_te_rad=_result_scalar(result, "total_relative_te_rad", "total_follower_te_rad"),
                wheel_axis_stiffness_nm_per_rad=_result_scalar(result, "incremental_mesh_stiffness_nm_per_rad"))
            if not math.isclose(row["total_te_rad"], row["unloaded_te_rad"]+row["mesh_te_rad"]+row["elastic_te_rad"], rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError("总TE与未加载/网格/弹性三项不一致")
        rows.append(row)
    with _paper(lang) as plt:
        fig = plt.figure(layout="constrained")
        grid = fig.add_gridspec(2, 2)
        axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])]
        x = [row["driver_rad"]*180/math.pi for row in rows]
        axes[0].plot(x, [r.get("unloaded_te_rad", math.nan)*_ARCSEC for r in rows], "o-", color=_COLORS[0], ms=2.8)
        axes[0].set_ylabel(_text(lang, "未加载传动误差\nΔθ₂,₀ (arcsec)", "Unloaded transmission error\nΔθ₂,₀ (arcsec)"))
        for field, color, zh, en in (("mesh_te_rad", _COLORS[4], "网格初相位", "Mesh phase"),
            ("elastic_te_rad", _COLORS[1], "弹性传动误差", "Elastic transmission error"),
            ("total_te_rad", _COLORS[0], "总传动误差", "Total transmission error")):
            axes[1].plot(x, [r.get(field, math.nan)*_ARCSEC for r in rows], "o-", color=color, ms=2.5, label=_text(lang, zh, en))
        axes[1].set_ylabel(_text(lang, "大轮传动误差\nΔθ₂ (arcsec)", "Wheel transmission error\nΔθ₂ (arcsec)"))
        axes[1].legend(loc="best")
        axes[2].plot(x, [r.get("wheel_axis_stiffness_nm_per_rad", math.nan) for r in rows], "o-", color=_COLORS[2], ms=2.8)
        axes[2].set_ylabel(_text(lang, "从动轴增量扭转刚度\nKθ₂ (N·m/rad)", "Wheel-axis incremental\ntorsional stiffness\nKθ₂ (N·m/rad)"))
        for ax, letter in zip(axes, "abc"):
            ax.set_xlabel(_text(lang, "小轮增量角 θ₁ (deg)", "Pinion increment angle, θ₁ (deg)"))
            _panel(ax, letter)
        return _save(fig, output, "07_te_stiffness_cycle", lang, rows,
            _text(lang, "TE保留符号与实际未加载零点；从动轴刚度为N·m/rad。线只连接有效计算点，不表示插值或网格精度保证。",
                "Transmission errors retain their signs and actual unloaded zero. Stiffness is about the wheel axis in N·m/rad. Lines join valid computed points without an interpolation or mesh-accuracy guarantee."), height_mm=140)


def plot_workflow_equilibrium(bundle, output, *, lang="zh"):
    """实际齿对分担、两轴残差及互补数值门槛，不混为物理认证。"""
    _validate(bundle)
    rows = []
    for result in _loaded_results(bundle):
        phase = _phase(result)
        if not _accepted(result):
            rows.append(dict(driver_rad=phase, converged=False, kind="failed", status=result.get("status")))
            continue
        equilibrium = result["equilibrium"]
        diagnostics = equilibrium["diagnostics"]
        total = _result_scalar(result, "total_normal_force_n")
        if total <= 0:
            raise ValueError("接受承载结果的总法向力须为正")
        for offset, force in result["pair_force_shares"].items():
            rows.append(dict(driver_rad=phase, converged=True, kind="pair", tooth_offset=int(offset),
                force_n=_finite(force, "pair_force"), fraction=_finite(force, "pair_force")/total))
        target = _result_scalar(result, "wheel_target_torque_nm")
        nominal = _result_scalar(result, "pinion_torque_nm")
        rows.append(dict(driver_rad=phase, converged=True, kind="torque",
            pinion_relative_error=(_result_scalar(result, "pinion_torque_reconstructed_nm")-nominal)/abs(nominal),
            wheel_relative_error=(_result_scalar(result, "wheel_torque_reconstructed_nm")-target)/abs(target)))
        limits = (("balance", abs(equilibrium["balance_residual"]), "balance_acceptance"),
                  ("complementarity", equilibrium["complementarity_residual_nm"], "complementarity_acceptance_nm"),
                  ("negative_gap", max(0., -equilibrium["minimum_remaining_gap_m"]), "gap_acceptance_m"))
        for field, value, limit_key in limits:
            limit = _finite(diagnostics[limit_key], limit_key)
            if limit <= 0:
                raise ValueError("数值接受门槛必须为正")
            rows.append(dict(driver_rad=phase, converged=True, kind="certificate", field=field,
                value=value, acceptance=limit, normalized_residual=value/limit))
    with _paper(lang) as plt:
        fig, axes = plt.subplots(3, 1, sharex=True, layout="constrained")
        offsets = sorted({row["tooth_offset"] for row in rows if row["kind"] == "pair"})
        phases = sorted({_phase(result) for result in bundle["ltca"]["results"]})
        for index, offset in enumerate(offsets):
            data = {row["driver_rad"]: row["fraction"]*100 for row in rows if row["kind"] == "pair" and row["tooth_offset"] == offset}
            axes[0].plot([x*180/math.pi for x in phases], [data.get(x, math.nan) for x in phases], "o-", ms=2.4,
                         color=_COLORS[index % len(_COLORS)], label=f"k={offset}")
        axes[0].set_ylabel(_text(lang, "齿对法向力占比\nFₙ,k/ΣFₙ (%)", "Tooth-pair normal-force share\nFₙ,k/ΣFₙ (%)"))
        axes[0].legend(ncol=max(1, len(offsets)))
        for field, zh, en, color in (("pinion_relative_error", "小轮轴重构", "Pinion reconstruction", _COLORS[1]),
            ("wheel_relative_error", "大轮轴约束", "Wheel constraint", _COLORS[0])):
            data = {r["driver_rad"]: r[field]*100 for r in rows if r["kind"] == "torque"}
            axes[1].plot([x*180/math.pi for x in phases], [data.get(x, math.nan) for x in phases], "o-", ms=2.4, label=_text(lang, zh, en), color=color)
        axes[1].set_ylabel(_text(lang, "力矩相对残差\nΔT/|T| (%)", "Relative torque residual\nΔT/|T| (%)"))
        axes[1].legend(ncol=2)
        for index, (field, zh, en) in enumerate((("balance", "力矩平衡", "Torque balance"),
            ("complementarity", "互补条件", "Complementarity"), ("negative_gap", "负间隙", "Negative gap"))):
            data = {r["driver_rad"]: r["normalized_residual"] for r in rows if r["kind"] == "certificate" and r["field"] == field}
            axes[2].plot([x*180/math.pi for x in phases], [data.get(x, math.nan) for x in phases], "o-", ms=2.4,
                         color=_COLORS[index], label=_text(lang, zh, en))
        axes[2].axhline(1., color="#777777", ls="--", lw=.7)
        axes[2].set_ylabel(_text(lang, "残差/接受门槛\nr/r_lim (1)", "Residual/acceptance threshold\nr/r_lim (1)"))
        axes[2].set_yscale("symlog", linthresh=1e-8, linscale=.3)
        axes[2].set_yticks([0., 1e-8, 1e-6, 1e-4, 1e-2, 1.])
        axes[2].set_ylim(0., 1.5)
        axes[2].legend(ncol=3)
        axes[2].set_xlabel(_text(lang, "小轮增量角 θ₁ (deg)", "Pinion increment angle, θ₁ (deg)"))
        for ax, letter in zip(axes, "abc"):
            _panel(ax, letter)
        return _save(fig, output, "08_load_sharing_equilibrium", lang, rows,
            _text(lang, "分担来自实际接触力。两轴力矩定义与数值接受门槛分别显示，数值平衡不代替网格或物理验证。",
                "Load shares use actual contact forces. Torque definitions and numerical acceptance thresholds are shown separately; numerical equilibrium does not establish mesh or physical validation."), height_mm=155)


def plot_spiral_workflow(bundle, output, *, lang="zh", data_root=None, require_all=True,
        include_contact_diagnostics=False):
    """输出实际完成的业务图；粗接触面片仅在显式诊断模式绘制。"""
    if not isinstance(include_contact_diagnostics,bool):raise ValueError('接触诊断开关须为布尔值')
    _validate(bundle)
    from .contact_pattern import plot_workflow_contact_pattern
    functions = [("geometry", lambda: plot_workflow_geometry(bundle, output, lang=lang)),
        ("cad", lambda: plot_workflow_cad(bundle, output, lang=lang, data_root=data_root)),
        ("geometry_benchmark", lambda: plot_workflow_benchmark(bundle, output, kind="geometry", lang=lang)),
        ("strength_benchmark", lambda: plot_workflow_benchmark(bundle, output, kind="strength", lang=lang)),
        ("tca", lambda: plot_workflow_tca(bundle, output, lang=lang)),
        ("cycle", lambda: plot_workflow_cycle(bundle, output, lang=lang)),
        ("equilibrium", lambda: plot_workflow_equilibrium(bundle, output, lang=lang))]
    if include_contact_diagnostics:
        functions.extend([("contact",lambda: plot_workflow_contact(bundle,output,lang=lang)),
            ("contact_pattern",lambda: plot_workflow_contact_pattern(bundle,output,lang=lang,data_root=data_root))])
    stages = set(bundle.get("completed_stages", []))
    benchmark = bundle.get("masta_benchmark", {})
    applicable = {"geometry": bool(bundle.get("geometry")),
        "cad": "cad" in stages or bool(bundle.get("cad", {}).get("pair_mesh")),
        "geometry_benchmark": bool(benchmark.get("geometry_rows")),
        "strength_benchmark": bool(benchmark.get("strength_rows")),
        "tca": "tca" in stages or bool(bundle.get("tca", {}).get("rows", bundle.get("tca", {}).get("results"))),
        "contact_pattern": bundle.get("kind",bundle.get("case",{}).get("kind")) == "spiral" and
            bool(bundle.get("cad",{}).get("flank_mesh")) and len(bundle.get("ltca",{}).get("results",[]))>=3,
        **{name: "ltca" in stages or bool(bundle.get("ltca", {}).get("results"))
           for name in ("contact", "cycle", "equilibrium")}}
    records, failures, omitted = [], [], []
    if not include_contact_diagnostics:
        omitted.extend(dict(figure=name,status='not_requested',reason='粗面片仅用于显式接触诊断；压力结果使用连续齿面云图')
            for name in ('contact','contact_pattern'))
    for name, function in functions:
        if not applicable[name]:
            omitted.append(dict(figure=name, status="not_applicable",
                reason="未提供外部参照" if name.endswith("benchmark") else
                    "完整制造工作面或三个实际姿态尚未提供" if name=="contact_pattern" else "尚无该阶段实际结果"))
            continue
        try:
            record = function()
            record["case_fingerprint"] = bundle.get("case_fingerprint")
            (Path(output)/lang/(record["figure"]+".json")).write_text(
                json.dumps(record, ensure_ascii=False, indent=2)+"\n", encoding="utf-8", newline="\n")
            records.append(record)
        except (KeyError, ValueError, FileNotFoundError) as error:
            failures.append(dict(figure=name, status="not_rendered", reason=str(error)))
            if require_all:
                raise
    directory = Path(output)/lang
    directory.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema="bevel_workflow_figures_v1", language=lang, figures=records, failures=failures,
                    omitted=omitted, complete=not failures, rendered_count=len(records),
                    case_schema=bundle["schema"], case_fingerprint=bundle.get("case_fingerprint"),
                    figure_metadata_sha256={record["figure"]:hashlib.sha256((directory/(record["figure"]+".json")).read_bytes()).hexdigest()
                                            for record in records})
    (directory/"manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8", newline="\n")
    return manifest
