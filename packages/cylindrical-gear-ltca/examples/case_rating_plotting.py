"""新写单案例几何、名义载荷与显式因子强度对比图；数值来自实际复算源表。"""
from pathlib import Path


LABELS = {
    "center_distance_m": "中心距 a", "base_pitch_m": "基节 pb", "face_overlap_m": "重叠齿宽 b",
    "pinion.reference_radius_m": "小轮参考半径", "wheel.reference_radius_m": "大轮参考半径",
    "pinion.base_radius_m": "小轮基圆半径", "wheel.base_radius_m": "大轮基圆半径",
    "pinion.root_form_radius_m": "小轮根形半径", "wheel.root_form_radius_m": "大轮根形半径",
    "pinion.tip_form_radius_m": "小轮顶形半径", "wheel.tip_form_radius_m": "大轮顶形半径",
    "nominal_contact_stress_pa": "名义接触 σH0",
    "pinion.nominal_root_stress_pa": "小轮名义齿根 σF01", "wheel.nominal_root_stress_pa": "大轮名义齿根 σF02",
    "pinion.contact_stress_pa": "小轮评级接触 σH1", "wheel.contact_stress_pa": "大轮评级接触 σH2",
    "pinion.root_stress_pa": "小轮评级齿根 σF1", "wheel.root_stress_pa": "大轮评级齿根 σF2",
    "pinion.contact_safety_factor": "小轮接触 SH1", "wheel.contact_safety_factor": "大轮接触 SH2",
    "pinion.bending_safety_factor": "小轮弯曲 SF1", "wheel.bending_safety_factor": "大轮弯曲 SF2",
    "normal_force_n": "LTCA：评级后总法向 Fn", "nominal_tangential_force_n": "强度：名义参考圆 Ft",
}


def plot_case_ratings(rows, output):
    """同单位分面，保持标准评级应力与LTCA局部峰压的物理含义分开。"""
    import numpy as np
    from gear_contact.plotting import _pyplot
    from gear_contact.io import write_json
    plt = _pyplot()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    selected_geometry = ("center_distance_m", "face_overlap_m", "pinion.reference_radius_m",
                         "wheel.reference_radius_m", "pinion.base_radius_m", "wheel.base_radius_m")
    lengths = [dict(r) for r in rows if r["group"] == "geometry" and r["field"] in selected_geometry]
    for row in lengths:
        if row["field"].endswith("radius_m"):
            row["project"] *= 2
            row["reference"] *= 2
    stresses = [r for r in rows if r["group"] == "strength" and r["unit"] == "Pa"]
    safeties = [r for r in rows if r["group"] == "strength" and r["field"].endswith("_safety_factor")]
    forces = [r for r in rows if r["unit"] == "N"]
    groups = [(lengths, 1e3, "mm", "a  关键几何参数"),
              (stresses, 1e-6, "MPa", "b  ISO 6336:2006 显式因子应力"),
              (forces, 1., "N", "c  同工况的两种载荷层级"),
              (safeties, 1., "1", "d  外部材料／寿命应力限下的安全系数")]
    with plt.rc_context({"font.size": 8.5, "axes.titlesize": 10, "ytick.labelsize": 8.5}):
        fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.7))
        fig.subplots_adjust(left=.18, right=.95, top=.90, bottom=.15, wspace=.80, hspace=.66)
        for ax, (items, scale, unit, title) in zip(axes.flat, groups):
            if not items:
                ax.set_visible(False)
                continue
            actual, reference = [r["project"]*scale for r in items], [r["reference"]*scale for r in items]
            y = np.arange(len(items))
            ax.barh(y-.16, reference, height=.30, facecolor="white", edgecolor="#CC793C", linewidth=1.0,
                    label="MASTA")
            ax.barh(y+.16, actual, height=.28, color="#1D6D91", alpha=.85, label="本项目")
            ax.set_yticks(y, [(LABELS.get(r["field"], r["field"]).replace("参考半径", "分度圆直径")
                              .replace("基圆半径", "基圆直径") if title.startswith("a")
                             else LABELS.get(r["field"], r["field"])) for r in items])
            ax.invert_yaxis()
            ax.set(title=title, xlabel=f"数值 / {unit}")
            ax.grid(axis="x", alpha=.35)
            xmax = max([*actual, *reference, 1e-30])
            ax.set_xlim(0., xmax*1.28)
            for k, value in enumerate(actual):
                ax.text(max(actual[k], reference[k])+xmax*.025, y[k], f"{value:.4g}", va="center", fontsize=8)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .98), ncol=2, frameon=False)
        fig.text(.5, .065, "18/77斜齿轮副：关键几何尺寸、标准强度与安全系数对比。\n"
                 "强度计算采用相同评级因子与材料应力限；标准强度与LTCA局部接触压力分别评价。",
                 ha="center", va="center", fontsize=8, color="#536777")
        paths = []
        for extension in ("png", "svg"):
            path = output/f"02_case_benchmarks.{extension}"
            fig.savefig(path, dpi=240, facecolor="white")
            paths.append(path.name)
        plt.close(fig)
    write_json(output/"rating_plot_data.json", rows)
    return paths
