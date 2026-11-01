"""新写工程对标绘图工具：只展示双方原标量，不配准或插值周期数组。

MASTA仅用于对标，另一侧由本项目实际求解。每个字段独占一个
量纲明确的面板；绘图不改变逐字段容差判断，也不构成精度认证。
"""
from __future__ import annotations

import math
from pathlib import Path
import textwrap


_LABELS = {
    "te_peak_to_peak_m": "承载 TE 峰峰值",
    "max_pressure_pa": "最大 Hertz 压力",
    "maximum_pressure_pa": "最大 Hertz 压力",
    "mean_normal_approach_m": "平均法向接近量",
    "mean_mesh_stiffness_n_per_m": "平均法向刚度（定义见说明）",
    "normal_force_n": "法向载荷",
    "tangential_force_n": "名义切向力",
    "normal_module_m": "法向模数",
    "transverse_module_m": "端面模数",
    "centre_distance_m": "中心距",
    "center_distance_m": "中心距",
    "pressure_angle_rad": "法向压力角",
    "transverse_pressure_angle_rad": "端面压力角",
    "helix_angle_rad": "螺旋角",
    "contact_ratio": "端面重合度",
    "transverse_contact_ratio": "端面重合度",
    "overlap_ratio": "轴向重合度",
    "total_contact_ratio": "总重合度",
    "pinion_pitch_radius_m": "小轮节圆半径",
    "wheel_pitch_radius_m": "大轮节圆半径",
    "pinion_base_radius_m": "小轮基圆半径",
    "wheel_base_radius_m": "大轮基圆半径",
}


def _display(unit, field):
    """由 SI 转为显示单位；刚度与线载荷虽同量纲仍保留各自名称。"""
    if unit == "m":
        return (1e6, "μm") if "te_" in field or "approach" in field else (1e3, "mm")
    if unit == "rad":
        return (1e6, "μrad") if "te_" in field else (180 / math.pi, "°")
    if unit in ("N_per_m", "N/m"):
        return (1e-6, "N/μm") if "stiffness" in field else (1e-3, "N/mm")
    conversions = {"Pa": (1e-6, "MPa"), "N": (1., "N"), "1": (1., "1"),
                   "N_m": (1., "N·m"), "N*m": (1., "N·m"),
                   "m_per_N": (1e6, "μm/N"), "Hz": (1., "Hz")}
    if unit not in conversions:
        raise ValueError(f"对标图尚未定义显示换算的单位：{unit}")
    return conversions[unit]


def _validated_rows(rows):
    """缺失双方数据的字段不画；拒绝非有限值、重复项和不一致单位。"""
    selected, identities, cases, units = [], set(), set(), {}
    for source in rows:
        row = dict(source)
        for key in ("case_id", "label", "field", "unit", "group"):
            if not isinstance(row.get(key), str) or not row[key]:
                raise ValueError(f"对标图字段 {key} 须为非空文字")
        if row["group"] not in ("common", "ltca"):
            raise ValueError("对标图 group 只接受 common 或 ltca")
        if row.get("project") is None or row.get("reference") is None:
            continue
        for key in ("project", "reference"):
            value = row[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"对标图 {key} 须为有限实数；缺失值不得填零")
        identity = (row["group"], row["case_id"], row["field"])
        if identity in identities:
            raise ValueError(f"对标图存在重复项：{identity}")
        identities.add(identity)
        quantity = (row["group"], row["field"])
        if quantity in units and units[quantity] != row["unit"]:
            raise ValueError(f"同一对标字段的单位不一致：{row['field']}")
        units[quantity] = row["unit"]
        _display(row["unit"], row["field"])
        # 从双方原值求差；不依赖展示字段被四舍五入后的相对差。
        row["relative_error"] = (None if row["reference"] == 0 else
                                 (row["project"] - row["reference"]) / abs(row["reference"]))
        cases.add(row["case_id"])
        selected.append(row)
    if len(cases) > 9:
        raise ValueError("单次工程对标图最多接受 9 个案例，请分批调用")
    return selected


def plot_comparison(rows, output):
    """按字段分面保存 PNG、SVG、PDF，返回绝对路径列表。

    rows 包含 case_id、label、field、unit、project、reference、
    relative_error、matched、group。两侧标量使用 SI 原值，group 为
    common 或 ltca；缺少刚度字段时不会补造刚度面板。此函数不导出
    CSV，源表由调用者统一保存。每页最多六个字段，不连接独立案例。
    """
    selected = _validated_rows(rows)
    if not selected:
        return []
    from gear_contact.plotting import _pyplot

    plt = _pyplot()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for group in ("common", "ltca"):
        grouped = {}
        for row in selected:
            if row["group"] == group:
                grouped.setdefault(row["field"], []).append(row)
        fields = list(grouped)
        for page, start in enumerate(range(0, len(fields), 6), 1):
            page_fields = fields[start:start + 6]
            columns = min(2, len(page_fields))
            panel_rows = math.ceil(len(page_fields) / columns)
            largest = max(len(grouped[field]) for field in page_fields)
            height = .65 + panel_rows * max(2.25, .24 * largest + .7) + .85
            with plt.rc_context({"font.size": 8, "axes.titlesize": 9,
                                 "axes.labelsize": 8, "xtick.labelsize": 7.5,
                                 "ytick.labelsize": 7.5}):
                fig, axes = plt.subplots(panel_rows, columns, squeeze=False,
                                         figsize=(10.6 if columns == 2 else 6.0, height))
                fig.subplots_adjust(left=.26 if columns == 1 else .17,
                                    right=.91, top=1 - .65 / height,
                                    bottom=1.05 / height, hspace=.6, wspace=.8)
                for index, field in enumerate(page_fields):
                    ax = axes.flat[index]
                    data = grouped[field]
                    scale, display_unit = _display(data[0]["unit"], field)
                    reference = [row["reference"] * scale for row in data]
                    project = [row["project"] * scale for row in data]
                    for y, (a, b) in enumerate(zip(reference, project)):
                        ax.plot((a, b), (y, y), color="#BAC5CB", linewidth=.8, zorder=1)
                    ax.scatter(reference, range(len(data)), s=27, marker="s", facecolors="white",
                               edgecolors="#D07537", linewidths=1., label="MASTA", zorder=3)
                    ax.scatter(project, range(len(data)), s=18, marker="o", color="#176B91",
                               label="本项目", zorder=4)
                    labels = [textwrap.fill(row["label"], 15) for row in data]
                    ax.set_yticks(range(len(data)), labels)
                    ax.set_ylim(len(data) - .5, -.5)
                    title = _LABELS.get(field, field.replace("_", " "))
                    ax.set(title=f"{chr(97 + index)}  {title}", xlabel=f"数值 / {display_unit}")
                    ax.grid(axis="x")
                    ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 4), useMathText=True)
                    ax.margins(x=.15)
                    ax.text(1.03, 1.015, "相对差 / %", transform=ax.transAxes, fontsize=7.5,
                            ha="left", va="bottom", color="#526470")
                    for y, row in enumerate(data):
                        error = row["relative_error"]
                        text = "未定义*" if error is None else f"{100 * error:+.3g}"
                        ax.text(1.03, y, text, transform=ax.get_yaxis_transform(), va="center",
                                color="#934F35" if error and abs(error) > .01 else "#526470",
                                fontsize=7.5)
                for ax in list(axes.flat)[len(page_fields):]:
                    ax.set_visible(False)
                handles, labels = axes.flat[0].get_legend_handles_labels()
                fig.legend(handles, labels, loc="upper center", ncol=2,
                           bbox_to_anchor=(.5, 1), fontsize=8.5)
                if group == "common":
                    note = "共同几何与名义载荷；逐字段容差和同输入状态见源表。差异 = (当前−参照)/|参照|。"
                else:
                    note = "LTCA 响应比较；齿顶延伸、结构模型和边界不完全等效，响应差异不构成精度认证。"
                fig.text(.5, .40 / height, note + "\n独立案例不插值、不连为周期曲线；* 零参照的相对差未定义，原值保留。",
                         ha="center", va="center", fontsize=7.5, color="#526470")
                stem = f"masta_{group}_comparison" + (f"_{page}" if len(fields) > 6 else "")
                for extension in ("png", "svg", "pdf"):
                    path = output / f"{stem}.{extension}"
                    fig.savefig(path, dpi=300, facecolor="white")
                    paths.append(str(path.resolve()))
                plt.close(fig)
    return paths
