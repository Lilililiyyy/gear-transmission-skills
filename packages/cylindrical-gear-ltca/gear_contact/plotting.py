"""计算结果的中文工程图及其原始绘图数据。

瞬时接触图按单相位、单齿对展开 Hertz 半宽；周期包络单独输出。
展示重建不参与指标计算；原始接触数据完整导出。Matplotlib 为可选依赖。
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Mapping

import numpy as np


_COLORS = ("#176B91", "#D47732", "#7259A0", "#328577", "#B55B72")


def _pyplot():
    """选择无界面后端和可用中文字体，保持 SVG 中的文字可编辑。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    fonts = [name for name in ("Microsoft YaHei", "Noto Sans CJK SC", "SimHei",
                              "Source Han Sans SC", "DejaVu Sans") if name in available]
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": fonts,
        "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
        "axes.unicode_minus": False, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .8,
        "axes.edgecolor": "#A0ACB6", "axes.labelcolor": "#334552",
        "text.color": "#223440", "xtick.color": "#526470", "ytick.color": "#526470",
        "grid.color": "#DCE3E8", "grid.linewidth": .6, "axes.axisbelow": True,
        "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    return plt


def _save(fig, output: Path, name: str) -> list[str]:
    """同时保存便于预览的 PNG 和可编辑的 SVG。"""
    files = []
    for extension in ("png", "svg"):
        path = output / f"{name}.{extension}"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        files.append(str(path.resolve()))
    return files


def _csv(path: Path, columns, rows) -> str:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)
    return str(path.resolve())


def analysis_metrics(result) -> dict[str, float]:
    """从原始求解数组计算统一指标，不从图像或插值网格反算。"""
    force_totals = [float(np.sum(r.solution.forces_n)) for r in result.rotations]
    peak_pressure = max((float(np.max(r.pressure_pa)) for r in result.rotations), default=0.)
    return {
        "te_peak_to_peak_um": float(np.ptp(result.normal_approach_m) * 1e6),
        "mean_approach_um": float(np.mean(result.normal_approach_m) * 1e6),
        "maximum_pressure_mpa": peak_pressure / 1e6,
        "mean_stiffness_n_per_um": float(np.mean(result.mesh_stiffness_n_per_m) * 1e-6),
        "normal_force_n": float(np.mean(force_totals)),
    }


def _contact_field(result, *, samples=(360, 240)):
    """重建周期扫掠包络：线性插值不超出采样凸包，非活动单元保持遮罩。

    两个通道分别为 MPa 和 N/mm；不是某一相位的瞬时 Hertz 接触带。
    最近邻活动遮罩保留未接触区和孔洞，边界精度由原始离散决定。
    """
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import Delaunay, cKDTree

    points = np.concatenate([np.column_stack((r.grid.face_coordinate_m,
                                              r.grid.roll_pinion_m)) for r in result.rotations]) * 1e3
    active = np.concatenate([r.solution.active_mask for r in result.rotations]).astype(bool)
    values = np.concatenate([np.column_stack((r.pressure_pa * 1e-6,
                                              r.line_load_n_per_m * 1e-3)) for r in result.rotations])
    x = np.linspace(*(np.asarray(result.pair.face_bounds_m) * 1e3), num=samples[0])
    y = np.linspace(*(np.asarray(result.pair.contact_roll_bounds_m) * 1e3), num=samples[1])
    xx, yy = np.meshgrid(x, y)
    if not np.any(active) or np.linalg.matrix_rank(points - points[0]) < 2:
        return xx, yy, np.ma.masked_all((*xx.shape, 2))
    triangulation = Delaunay(points)
    field = LinearNDInterpolator(triangulation, np.where(active[:, None], values, 0.))(xx, yy)
    nearest = cKDTree(points).query(np.column_stack((xx.ravel(), yy.ravel())))[1].reshape(xx.shape)
    mask = ~active[nearest] | ~np.all(np.isfinite(field), axis=2)
    return xx, yy, np.ma.array(field, mask=np.broadcast_to(mask[:, :, None], field.shape))


def _draw_contact_field(ax, field, *, channel=0, maximum, xlabel, ylabel, title):
    """连续色阶和共用量程；浅灰色表示未接触或未采样。"""
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    xx, yy, values = field
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.set_facecolor("#EDF1F4")
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_aspect("equal", adjustable="box")
    norm = Normalize(0., maximum if maximum > 0 else 1.)
    if values[:, :, channel].count():
        ax.contourf(xx, yy, values[:, :, channel], levels=np.linspace(0, norm.vmax, 160),
                    cmap="inferno", norm=norm, antialiased=False, corner_mask=False)
    else:
        ax.text(.5, .5, "本工况无承载单元\nNo loaded contact", ha="center", va="center",
                transform=ax.transAxes, color="#526470")
    return ScalarMappable(norm=norm, cmap="inferno")


def plot_contact_envelope_comparison(results: Mapping[str, object], output_dir, *,
                                     title: str = "Loaded contact envelope | One mesh cycle") -> list[str]:
    """用相同压力色标比较工况的周期接触包络，原始数据另存 CSV。"""
    if not results:
        raise ValueError("接触对比至少需要一组计算结果")
    plt = _pyplot()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    maximum = max(analysis_metrics(result)["maximum_pressure_mpa"] for result in results.values())
    fig, axes = plt.subplots(1, len(results), figsize=(4.4 * len(results), 4.1),
                             squeeze=False, constrained_layout=True)
    rows = []
    for index, (name, result) in enumerate(results.items()):
        ax = axes[0, index]
        image = _draw_contact_field(ax, _contact_field(result), maximum=maximum,
                                    xlabel="Face width coordinate / mm", ylabel="Pinion roll / mm",
                                    title=f"{chr(97 + index)}   {name}")
        ax.text(0., 1.18, f"Peak  {analysis_metrics(result)['maximum_pressure_mpa']:.0f} MPa",
                transform=ax.transAxes, fontsize=10, color=_COLORS[index % len(_COLORS)])
        for phase, rotation in zip(result.phase, result.rotations):
            rows.extend((name, phase, x * 1e3, y * 1e3, int(active), pressure * 1e-6)
                        for x, y, active, pressure in zip(rotation.grid.face_coordinate_m,
                            rotation.grid.roll_pinion_m, rotation.solution.active_mask, rotation.pressure_pa))
    bar = fig.colorbar(image, ax=list(axes.ravel()), label="Hertz pressure / MPa", shrink=.65, pad=.025)
    bar.outline.set_visible(False)
    fig.suptitle(title, fontsize=15, fontweight="bold")
    files = _save(fig, output, "contact_envelope")
    plt.close(fig)
    files.append(_csv(output / "contact_envelope_data.csv",
                      ("design", "phase", "face_coordinate_mm", "pinion_roll_mm", "active_contact", "pressure_MPa"), rows))
    return files


def _instantaneous_contact_field(result, *, phase_index=0, tooth_pair_index=None,
                                 samples=(501, 251)):
    """单齿对的条带式 Hertz 压力重建，不将其他相位或齿对叠到局部坐标。

    沿线线性插值 q 和曲率，再按 Hertz 式计算 b(s)、p0(s)，沿横向展开
    p(s,n)=p0(s) sqrt(1-(n/b(s))²)。这不是二维半空间接触边界求解。
    接触线端部延伸半个积分单元，使用端点值；不人为添加圆滑端部效应。
    """
    from .hertz import hertz_line_contact

    if isinstance(phase_index, bool) or not isinstance(phase_index, (int, np.integer)) or not 0 <= phase_index < len(result.rotations):
        raise ValueError("phase_index 必须是有效相位下标")
    rotation = result.rotations[phase_index]
    grid = rotation.grid
    pairs = np.unique(grid.tooth_pair_index)
    if tooth_pair_index is None:
        tooth_pair_index = int(max(pairs, key=lambda k: np.sum(rotation.solution.forces_n[grid.tooth_pair_index == k])))
    if tooth_pair_index not in pairs:
        raise ValueError("选定相位中不存在该齿对")
    selected = np.flatnonzero(grid.tooth_pair_index == tooth_pair_index)
    selected = selected[np.argsort(grid.face_coordinate_m[selected])]
    cosine = np.cos(result.pair.base_helix_angle_rad)
    s = grid.face_coordinate_m[selected] / cosine
    widths = grid.widths_m[selected]
    left, right = s[0] - widths[0] / 2, s[-1] + widths[-1] / 2
    along = np.unique(np.r_[np.linspace(left, right, samples[0]), s])
    q = np.interp(along, s, rotation.line_load_n_per_m[selected])
    curvature = np.interp(along, s, grid.curvature_per_m[selected])
    hertz = hertz_line_contact(q, curvature, result.inputs["material_pinion"], result.inputs["material_wheel"])
    extent = max(float(np.max(hertz.half_width_m)) * 1.25, 1e-6)
    across = np.linspace(-extent, extent, samples[1])
    distance = np.divide(across[:, None], hertz.half_width_m[None, :],
                         out=np.full((len(across), len(along)), np.inf),
                         where=hertz.half_width_m[None, :] > 0)
    pressure = hertz.maximum_pressure_pa[None, :] * np.sqrt(np.maximum(0., 1 - distance ** 2))
    return {"along_m": along, "across_m": across, "pressure_pa": pressure,
            "line_load_n_per_m": q, "half_width_m": hertz.half_width_m,
            "maximum_pressure_pa": hertz.maximum_pressure_pa,
            "phase": float(result.phase[phase_index]), "tooth_pair_index": tooth_pair_index,
            "selected": selected, "tooth_pair_force_n": float(np.sum(rotation.solution.forces_n[selected]))}


def plot_contact_comparison(results: Mapping[str, object], output_dir, *, phase_index: int = 0,
                            tooth_pair_index: int | None = None,
                            title: str = "Instantaneous contact footprint | Strip-wise Hertz reconstruction") -> list[str]:
    """固定相位和单齿对的瞬时压力图；默认选择该相位载荷最大的齿对。

    纵轴为接触线横向的齿面内距离，单位 μm；横轴沿接触线，单位 mm。
    两个方向采用不同显示尺度，保持微米级接触宽度可见。
    """
    if not results:
        raise ValueError("接触对比至少需要一组计算结果")
    plt = _pyplot()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fields = {name: _instantaneous_contact_field(result, phase_index=phase_index,
                                                 tooth_pair_index=tooth_pair_index) for name, result in results.items()}
    maximum = max(float(np.max(field["maximum_pressure_pa"])) * 1e-6 for field in fields.values())
    extent = max(float(np.max(field["across_m"])) * 1e6 for field in fields.values())
    fig, axes = plt.subplots(1, len(results), figsize=(max(8.6, 5.1 * len(results)), 4.8),
                             squeeze=False, constrained_layout=True)
    rows, source_rows = [], []
    for index, (name, field) in enumerate(fields.items()):
        ax = axes[0, index]
        image = ax.contourf(field["along_m"] * 1e3, field["across_m"] * 1e6,
                            field["pressure_pa"] * 1e-6, levels=np.linspace(0, maximum if maximum > 0 else 1., 128),
                            cmap="turbo", antialiased=False)
        ax.set_facecolor(plt.get_cmap("turbo")(0.))
        ax.set(xlabel="Along contact line / mm", ylabel="Across contact line / μm", ylim=(-extent, extent),
               title=f"{name} | phase {field['phase']:.3f} | tooth pair {field['tooth_pair_index']:+d}")
        ax.set_box_aspect(.4)
        if not np.any(field["pressure_pa"]):
            ax.text(.5, .5, "本工况无承载单元\nNo loaded contact", color="white", ha="center", va="center", transform=ax.transAxes)
        result = results[name]
        rotation = result.rotations[phase_index]
        for i in field["selected"]:
            source_rows.append((name, field["phase"], field["tooth_pair_index"],
                                rotation.grid.face_coordinate_m[i] * 1e3, rotation.grid.widths_m[i] * 1e3,
                                rotation.line_load_n_per_m[i] * 1e-3, rotation.half_width_m[i] * 1e6,
                                rotation.pressure_pa[i] * 1e-6, rotation.solution.forces_n[i]))
        rows.extend((name, field["phase"], field["tooth_pair_index"], s * 1e3, n * 1e6, p * 1e-6)
                    for n, line in zip(field["across_m"], field["pressure_pa"])
                    for s, p in zip(field["along_m"], line))
    bar = fig.colorbar(image, ax=list(axes.ravel()), label="Contact pressure / MPa", shrink=.7, pad=.025)
    bar.outline.set_visible(False)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.supxlabel("Transverse width enlarged for display; read the mm / μm axis units", fontsize=9, color="#526470")
    files = _save(fig, output, "contact_pattern")
    plt.close(fig)
    files.append(_csv(output / "contact_pattern_data.csv",
                      ("design", "phase", "tooth_pair_index", "along_line_mm", "across_line_um", "pressure_MPa"), rows))
    files.append(_csv(output / "contact_pattern_source.csv",
                      ("design", "phase", "tooth_pair_index", "face_coordinate_mm", "strip_length_mm",
                       "line_load_N_per_mm", "half_width_um", "maximum_pressure_MPa", "force_N"), source_rows))
    return files


def plot_analysis(result, output_dir, *, title: str = "圆柱齿轮承载接触分析") -> list[str]:
    """生成周期响应、接触云图和原始 CSV；展示插值不改变计算数组。"""
    plt = _pyplot()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    phase = np.asarray(result.phase, dtype=float)
    peak_pressure, loaded_length, force_total = [], [], []
    cell_rows = []
    for p, rotation in zip(phase, result.rotations):
        grid, solution = rotation.grid, rotation.solution
        active = np.asarray(solution.active_mask, dtype=bool)
        peak_pressure.append(float(np.max(rotation.pressure_pa)) / 1e6)
        loaded_length.append(float(np.sum(grid.widths_m[active])) * 1e3)
        force_total.append(float(np.sum(solution.forces_n)))
        for i in range(len(grid.widths_m)):
            cell_rows.append((p, int(grid.tooth_pair_index[i]),
                              grid.face_coordinate_m[i] * 1e3,
                              grid.roll_pinion_m[i] * 1e3, int(active[i]),
                              solution.forces_n[i], rotation.line_load_n_per_m[i] * 1e-3,
                              rotation.pressure_pa[i] * 1e-6, rotation.gaps_m[i] * 1e6,
                              solution.residual_gaps_m[i] * 1e9))
    response_rows = zip(phase, np.asarray(result.normal_approach_m) * 1e6,
                        result.angular_te_rad, np.asarray(result.mesh_stiffness_n_per_m) * 1e-6,
                        peak_pressure, loaded_length, force_total)
    files = [_csv(output / "plot_phase_data.csv",
                  ("phase", "normal_approach_um", "angular_te_rad", "tangent_stiffness_N_per_um",
                   "maximum_pressure_MPa", "loaded_contact_length_mm", "total_normal_force_N"),
                  response_rows),
             _csv(output / "plot_contact_data.csv",
                  ("phase", "tooth_pair_index", "face_coordinate_mm", "pinion_roll_mm",
                   "active_contact", "force_N", "line_load_N_per_mm", "pressure_MPa",
                   "initial_gap_um", "residual_gap_nm"), cell_rows)]

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.4), constrained_layout=True)
    panels = (
        (axes[0, 0], np.asarray(result.normal_approach_m) * 1e6, "法向接近量 / μm", "a  准静态传递误差"),
        (axes[0, 1], np.asarray(result.mesh_stiffness_n_per_m) * 1e-6, "切线啮合刚度 / (N/μm)", "b  切线刚度"),
        (axes[1, 0], peak_pressure, "峰值接触压力 / MPa", "c  各相位最大压力"),
        (axes[1, 1], loaded_length, "实际承载接触线长度 / mm", "d  有载接触范围"),
    )
    for ax, values, ylabel, panel_title in panels:
        ax.plot(phase, values, color=_COLORS[0], linewidth=1.9)
        ax.set(xlabel="啮合周期相位", ylabel=ylabel, title=panel_title)
        ax.grid(alpha=.7)
    fig.suptitle(title, fontsize=14)
    files.extend(_save(fig, output, "analysis_dashboard"))
    plt.close(fig)

    field = _contact_field(result)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5), constrained_layout=True)
    for ax, channel, label, panel_title, maximum in (
        (axes[0], 0, "接触压力 / MPa", "a  接触压力云图", max(peak_pressure)),
        (axes[1], 1, "线载荷 / (N/mm)", "b  线载荷云图",
         max(float(np.max(r.line_load_n_per_m)) * 1e-3 for r in result.rotations)),
    ):
        image = _draw_contact_field(ax, field, channel=channel, maximum=maximum,
                                    xlabel="共同轴向坐标 / mm", ylabel="小轮展长 / mm", title=panel_title)
        bar = fig.colorbar(image, ax=ax, label=label, shrink=.65, pad=.025)
        bar.outline.set_visible(False)
    fig.suptitle("一个啮合周期的承载接触包络 · 连续填色，非接触区遮罩", fontsize=13)
    files.extend(_save(fig, output, "contact_envelope"))
    plt.close(fig)
    files.extend(plot_contact_comparison({"Highest-load tooth pair": result}, output, phase_index=len(phase) * 3 // 8))

    from .analysis import transmission_error_spectrum
    spectrum = transmission_error_spectrum(result.normal_approach_m)
    orders = np.asarray(spectrum["order"])
    amplitudes = np.asarray(spectrum["amplitude"]) * 1e6
    arguments = np.asarray(spectrum["phase_rad"])
    valid_phase = np.isfinite(arguments)
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4), constrained_layout=True)
    axes[0].bar(orders, amplitudes, color=_COLORS[0], width=.65)
    axes[0].set(xlabel="啮合阶次", ylabel="法向 TE 谐波幅值 / μm", title="a  去均值后的周期激励")
    axes[0].set_xticks(orders)
    if np.any(valid_phase):
        axes[1].scatter(orders[valid_phase], arguments[valid_phase], s=35, color=_COLORS[0])
    else:
        axes[1].text(.5, .5, "幅值近零，无可报告相位", ha="center", va="center", transform=axes[1].transAxes)
    axes[1].set(xlabel="啮合阶次", ylabel="余弦相位 / rad", title="b  有效幅值的相位")
    axes[1].set_xticks(orders)
    axes[1].set_ylim(-np.pi * 1.1, np.pi * 1.1)
    for ax in axes:
        ax.grid(axis="y", alpha=.18)
    fig.suptitle("TE 谐波：首样本位于相位 0；近零幅值不报告相位", fontsize=12)
    files.extend(_save(fig, output, "te_harmonics"))
    plt.close(fig)
    files.append(_csv(output / "te_spectrum.csv", ("mesh_order", "amplitude_um", "cosine_phase_rad"),
                      ((order, amplitude, argument if np.isfinite(argument) else "")
                       for order, amplitude, argument in zip(orders, amplitudes, arguments))))
    return files


def plot_comparison(results: Mapping[str, object], output_dir, *, title: str = "修形方案性能对比") -> list[str]:
    """比较相同工况下各方案的周期曲线和主要指标。"""
    if not results:
        raise ValueError("方案对比至少需要一组计算结果")
    plt = _pyplot()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    metrics = {name: analysis_metrics(result) for name, result in results.items()}
    names = list(results)
    fig = plt.figure(figsize=(12, 7.3), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(1.2, 1))
    hero = fig.add_subplot(grid[0, :])
    for i, (name, result) in enumerate(results.items()):
        hero.plot(result.phase, np.asarray(result.normal_approach_m) * 1e6,
                  label=name, color=_COLORS[i % len(_COLORS)], linewidth=1.7)
    hero.set(xlabel="啮合周期相位", ylabel="法向接近量 / μm", title="a  各方案的有载周期响应")
    hero.legend(ncol=min(4, len(names)))
    hero.grid(alpha=.18)
    for i, (key, label, panel_title) in enumerate((
        ("te_peak_to_peak_um", "TE 峰峰值 / μm", "b  传递误差波动"),
        ("maximum_pressure_mpa", "峰值接触压力 / MPa", "c  局部压力"),
        ("mean_stiffness_n_per_um", "平均刚度 / (N/μm)", "d  平均切线刚度"),
    )):
        ax = fig.add_subplot(grid[1, i])
        ax.bar(np.arange(len(names)), [metrics[n][key] for n in names],
               color=[_COLORS[j % len(_COLORS)] for j in range(len(names))], width=.62)
        ax.set(xticks=np.arange(len(names)), xticklabels=names, ylabel=label, title=panel_title)
        ax.tick_params(axis="x", labelrotation=20)
        ax.grid(axis="y", alpha=.18)
    fig.suptitle(title, fontsize=14)
    files = _save(fig, output, "design_comparison")
    plt.close(fig)
    files.extend(_plot_response_comparisons(results, output))
    keys = list(next(iter(metrics.values())))
    files.append(_csv(output / "comparison_metrics.csv", ("design", *keys),
                      ((name, *(metrics[name][key] for key in keys)) for name in names)))
    files.append(_csv(output / "comparison_curves.csv",
                      ("design", "phase", "normal_approach_um", "tangent_stiffness_N_per_um"),
                      ((name, phase, approach * 1e6, stiffness * 1e-6)
                       for name, result in results.items() for phase, approach, stiffness in
                       zip(result.phase, result.normal_approach_m, result.mesh_stiffness_n_per_m))))
    return files


def _plot_response_comparisons(results, output):
    """独立的 TE 与刚度图，保持真实相位采样的分段线性曲线。"""
    plt = _pyplot()
    files = []
    for filename, attribute, factor, title, ylabel, metric in (
        ("transmission_error", "normal_approach_m", 1e6, "Loaded static transmission error",
         "Normal approach / μm", "te_peak_to_peak_um"),
        ("mesh_stiffness", "mesh_stiffness_n_per_m", 1e-6, "Tangent normal mesh stiffness",
         "Mesh stiffness / (N/μm)", "mean_stiffness_n_per_um"),
    ):
        fig = plt.figure(figsize=(11.4, 4.7), constrained_layout=True)
        grid = fig.add_gridspec(1, 2, width_ratios=(3.3, 1))
        ax, key = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
        key.set_axis_off()
        key.text(0, .96, "PEAK-TO-PEAK / μm" if filename == "transmission_error" else "MEAN / (N/μm)",
                 transform=key.transAxes, fontsize=9, color="#526470")
        styles = ("-", "--", "-.", ":")
        for index, (name, result) in enumerate(results.items()):
            color = _COLORS[index % len(_COLORS)]
            ax.plot(result.phase, np.asarray(getattr(result, attribute)) * factor, color=color,
                    linewidth=2.1, linestyle=styles[index % len(styles)])
            top = .82 - index * (.75 / max(3, len(results)))
            key.plot([0, .15], [top, top], color=color, linewidth=2.1,
                     linestyle=styles[index % len(styles)], transform=key.transAxes)
            key.text(.21, top, name, transform=key.transAxes, va="center", fontsize=10)
            key.text(.21, top - .1, f"{analysis_metrics(result)[metric]:.3f}",
                     transform=key.transAxes, color=color, fontsize=18, fontweight="bold")
        ax.set(xlabel="Mesh phase / base-pitch cycle", ylabel=ylabel, xlim=(0, 1))
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.grid(alpha=.8)
        ax.margins(y=.15)
        fig.suptitle(title, fontsize=16, fontweight="bold", x=.055, ha="left")
        files.extend(_save(fig, output, filename))
        plt.close(fig)
    return files


def plot_surface(modification, output_dir, *, title: str = "法向修形曲面", samples: int = 81) -> list[str]:
    """直接评价参数化修形曲面；正材料去除量表示间隙增大。"""
    if samples < 3:
        raise ValueError("修形曲面至少需要三个采样点")
    plt = _pyplot()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    u = np.linspace(0, 1, samples)
    v = np.linspace(0, 1, samples)
    uu, vv = np.meshgrid(u, v)
    removal = np.asarray(modification.evaluate(uu, vv)) * 1e6
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    image = axes[0].contourf(u, v, removal, levels=100, cmap="viridis", antialiased=False)
    axes[0].set(xlabel="齿廓因子：根部 0 → 顶部 1", ylabel="齿向因子：左端 0 → 右端 1", title="a  材料去除量")
    fig.colorbar(image, ax=axes[0], label="法向去除量 / μm")
    for lead in (0., .5, 1.):
        axes[1].plot(u, np.asarray(modification.evaluate(u, lead)) * 1e6,
                     label=f"齿向因子 {lead:g}", linewidth=1.7)
    axes[1].set(xlabel="齿廓因子", ylabel="法向去除量 / μm", title="b  三个齿向位置的齿廓截线")
    axes[1].legend()
    axes[1].grid(alpha=.18)
    fig.suptitle(title, fontsize=13)
    files = _save(fig, output, "modification_surface")
    plt.close(fig)
    files.append(_csv(output / "modification_surface.csv",
                      ("profile_factor", "lead_factor", "normal_removal_um"),
                      zip(uu.ravel(), vv.ravel(), removal.ravel())))
    return files


def plot_design_study(study, output_dir) -> list[str]:
    """以方案和工况原始结果构建最坏性能散点，不混用不同单位的加权和。"""
    plt = _pyplot()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    names = list(study.results)
    worst_te, worst_pressure = [], []
    for name in names:
        metrics = [analysis_metrics(r) for r in study.results[name].values()]
        worst_te.append(max(m["te_peak_to_peak_um"] for m in metrics))
        worst_pressure.append(max(m["maximum_pressure_mpa"] for m in metrics))
    fig, ax = plt.subplots(figsize=(7.4, 5.2), constrained_layout=True)
    for i, name in enumerate(names):
        selected = name in study.pareto_names
        ax.scatter(worst_te[i], worst_pressure[i], s=80 if selected else 45,
                   marker="D" if selected else "o", color=_COLORS[i % len(_COLORS)])
        ax.annotate(name + ("（Pareto）" if selected else ""),
                    (worst_te[i], worst_pressure[i]), xytext=(7, 7), textcoords="offset points")
    ax.set(xlabel="工况集合中的最坏 TE 峰峰值 / μm",
           ylabel="工况集合中的最高接触压力 / MPa",
           title="多工况方案取舍：两项指标越小越好")
    ax.margins(.25)
    ax.grid(alpha=.18)
    files = _save(fig, output, "robust_design_tradeoff")
    plt.close(fig)
    files.append(_csv(output / "robust_design_plot_data.csv",
                      ("design", "worst_te_peak_to_peak_um", "worst_pressure_MPa", "pareto"),
                      ((name, te, pressure, name in study.pareto_names)
                       for name, te, pressure in zip(names, worst_te, worst_pressure))))
    return files
