"""新写单案例LTCA对标后处理：原周期、同相位同齿Hertz带与材料坐标。

MASTA的端点压力和半宽保持对应采样定义；本项目保留各求解点的值。
Hertz带只是局部线接触重建，不是完整二维接触解，也不校准响应。
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from gear_contact.io import write_json
from gear_contact.plotting import _pyplot, _save


def _rows(path, columns, values):
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(values)
    return str(path.resolve())


def _band_values(peaks, half_widths, normal_edges):
    """在局部垂线坐标上评价Hertz形状；不改峰压、半宽或原条带力。"""
    peaks, half_widths = np.asarray(peaks), np.asarray(half_widths)
    centres = (normal_edges[:-1]+normal_edges[1:])/2
    ratio = np.divide(centres[None, :], half_widths[:, None],
                      out=np.full((len(peaks), len(centres)), np.inf),
                      where=half_widths[:, None] > 0)
    return peaks[:, None]*np.sqrt(np.maximum(0., 1-ratio**2))


def _endpoint_display(edges, peaks, half_widths):
    """端点间只细分显示网格，保持输入端点峰压和半宽，不作载荷归一化。"""
    fractions = np.arange(5)/4
    display_edges = np.r_[np.concatenate([a+(b-a)*fractions[:-1]
        for a, b in zip(edges[:-1], edges[1:])]), edges[-1]]
    centres = (display_edges[:-1]+display_edges[1:])/2
    return display_edges, np.interp(centres, edges, peaks), np.interp(centres, edges, half_widths)


def repeat_cycle_for_display(cycle, display_cycles=3):
    """将完整均匀稳态周期原样展开，末端追加起点作为周期闭合点。

    cycle第0列为从0起采且不重复端点的基节周期分数，其余列包括双方
    物理相位与SI响应。只平移显示坐标，不插值、拟合或重算响应；返回
    数组不替换原单周期数据，也不表示连续多周期时域求解。
    """
    if (isinstance(display_cycles, bool) or not isinstance(display_cycles, (int, np.integer))
            or display_cycles < 1):
        raise ValueError("显示周期数必须为正整数")
    cycle = np.asarray(cycle, dtype=float)
    if (cycle.ndim != 2 or cycle.shape[0] < 2 or cycle.shape[1] < 2
            or np.any(~np.isfinite(cycle))):
        raise ValueError("周期展开需要至少两个采样位置及完整有限的二维原周期数据")
    count = len(cycle)
    if not np.allclose(cycle[:, 0], np.arange(count)/count, rtol=0., atol=1e-12):
        raise ValueError("周期展开只接受从0起采且不重复端点的完整均匀周期")
    repeated = np.tile(cycle, (display_cycles, 1))
    repeated[:, 0] += np.repeat(np.arange(display_cycles), count)
    closing = cycle[0].copy()
    closing[0] = display_cycles
    return np.vstack((repeated, closing))


def _project_line_samples(rotation, pair_index, cosine):
    """整理本项目同齿中心与可选端点；延伸曲线按真实累计线长定位。

    普通域仍使用中心face/cos坐标。延伸域只用普通端点锚定累计width，
    不把整条非直接触曲线误当作face/cos，也不从MASTA坐标反标长度。
    """
    grid = rotation.grid
    mask = np.asarray(grid.tooth_pair_index) == pair_index
    indices = np.flatnonzero(mask)
    if not len(indices):
        raise ValueError("本项目在相同物理位置没有所选齿对，不能替换为其他齿对")
    geometry = getattr(rotation, "contact_geometry_data", None)
    endpoint_values = getattr(rotation, "endpoint_pressure", None)
    if geometry is None:
        order = np.argsort(np.asarray(grid.face_coordinate_m)[indices])
        width = np.asarray(grid.widths_m)[indices][order]
        centres = np.asarray(grid.face_coordinate_m)[indices][order]/cosine
        intervals = np.column_stack((centres-width/2, centres+width/2))
        edges = np.r_[intervals[:, 0], intervals[-1, 1]]
        anchor_spread = 0.
        endpoint_slice = None
        if endpoint_values is not None:
            raise ValueError("本项目端点压力缺少对应接触线几何，不能猜测端点布局")
    else:
        offsets = np.asarray(geometry.line_offsets, dtype=int)
        matches = [i for i, (start, end) in enumerate(zip(offsets[:-1], offsets[1:]))
                   if end > start and np.all(np.asarray(grid.tooth_pair_index)[start:end] == pair_index)]
        if len(matches) != 1:
            raise ValueError("本项目齿对没有唯一的中心接触线布局")
        line = matches[0]
        start, end = offsets[line:line+2]
        if not np.array_equal(indices, np.arange(start, end)):
            raise ValueError("本项目中心齿对索引与接触线offsets不一致")
        order = np.arange(len(indices))
        width = np.asarray(grid.widths_m)[indices]
        if np.any(width <= 0) or not np.all(np.isfinite(width)):
            raise ValueError("本项目真实接触线条带长度必须正且有限")
        edges = np.r_[0., np.cumsum(width)]
        endpoint_offsets = np.asarray(geometry.endpoints.line_offsets, dtype=int)
        if len(endpoint_offsets) != len(offsets):
            raise ValueError("本项目中心与端点接触线数量不一致")
        first, last = endpoint_offsets[line:line+2]
        if last-first != len(width)+1:
            raise ValueError("本项目端点布局必须保持每线N+1个位置")
        endpoint_slice = slice(first, last)
        ordinary = ~(np.asarray(geometry.endpoints.extended_pinion)[endpoint_slice]
                     | np.asarray(geometry.endpoints.extended_wheel)[endpoint_slice])
        if not np.any(ordinary):
            raise ValueError("本项目所选延伸接触线没有普通端点，不能锚定局部沿线坐标")
        # 粗网格的单元可能横跨普通/延伸边界；中心坐标未必位于真实弧长中点。
        # 使用独立端点锚定累计实际线长，不能用中心flag替代端点几何。
        anchors = edges[ordinary]-np.asarray(geometry.endpoints.face_coordinate_m)[endpoint_slice][ordinary]/cosine
        anchor_spread = float(np.ptp(anchors))
        if anchor_spread > 1e-11:
            raise ValueError("本项目普通端点累计线长与几何坐标不相容")
        edges -= np.median(anchors)
        intervals = np.column_stack((edges[:-1], edges[1:]))
    sample = {"mask": mask, "order": order, "widths_m": width, "intervals_m": intervals,
        "edges_m": edges, "peak_pa": np.asarray(rotation.pressure_pa)[indices][order],
        "half_width_m": np.asarray(rotation.half_width_m)[indices][order],
        "force_n": np.asarray(rotation.solution.forces_n)[indices][order],
        "line_load_n_per_m": np.asarray(rotation.line_load_n_per_m)[indices][order],
        "ordinary_anchor_spread_m": anchor_spread, "endpoints": None}
    if endpoint_values is not None:
        count = int(geometry.endpoints.line_offsets[-1])
        data = {}
        for name in ("pressure_pa", "half_width_m", "line_load_n_per_m"):
            values = np.asarray(endpoint_values[name], dtype=float)
            if values.shape != (count,) or np.any(~np.isfinite(values)) or np.any(values < 0):
                raise ValueError("本项目端点压力、半宽与线载荷必须保持有效全局端点布局")
            data[name] = values[endpoint_slice].copy()
        sample["endpoints"] = data
    return sample


def prepare_case_contact(result, metadata, reference_arrays):
    """核验纯几何相位后整理原值；返回数据供绘图与无图形测试。

    有N+1端点时，两侧均沿端点分段重建Hertz带；中心力和压力仍分别保留。
    仅有中心值的路线保持逐原单元常数，不跨单元或未采样区插值。
    """
    source = {key: np.asarray(value) for key, value in reference_arrays.items()}
    count = int(metadata["phase_registration"]["sample_count"])
    if len(result.rotations) != count or len(source["rotation_phase"]) != count:
        raise ValueError("同相位对标要求双方使用相同周期采样数")
    phase_error = ((np.asarray(result.phase)-source["rotation_phase"]+.5) % 1.)-.5
    if np.max(np.abs(phase_error)) > 1e-10:
        raise ValueError("本项目相位未与MASTA接触线几何配准，不能绘制逐相位对比")
    if metadata["phase_registration"]["phase_forward_direction"] != 1:
        raise NotImplementedError("本显示入口尚未实现反向周期的索引映射")
    beta = result.pair.base_helix_angle_rad
    if not math.isclose(beta, metadata["phase_registration"]["base_helix_angle_rad"], abs_tol=1e-12):
        raise ValueError("双方所选工作面的基圆螺旋角不一致")
    cosine = math.cos(beta)
    instant = metadata["suggested_instant"]
    ri, source_line = int(instant["rotation_index"]), int(instant["source_line_index"])
    rstart, rend = source["rotation_line_offsets"][ri:ri+2]
    lines = np.arange(rstart, rend, dtype=int)
    selected = lines[source["line_source_index"][lines] == source_line]
    if len(selected) != 1:
        raise ValueError("MASTA主承载齿对没有唯一接触线")
    line = int(selected[0])
    pair_index = int(source["line_project_tooth_pair_index"][line])
    if pair_index != instant["project_tooth_pair_index"]:
        raise ValueError("MASTA齿对索引和主齿选择记录不一致")
    cstart, cend = source["centre_offsets"][line:line+2]
    estart, eend = source["endpoint_offsets"][line:line+2]
    cs, es = slice(cstart, cend), slice(estart, eend)
    ds = source["centre_width_m"][cs]
    if len(ds)+1 != eend-estart or np.any(ds <= 0):
        raise ValueError("原主齿中心/端点布局或条带长度无效")
    source_s_edges = np.r_[0., np.cumsum(ds)]
    normal = ~(source["endpoint_extended_pinion"][es] | source["endpoint_extended_wheel"][es])
    if not np.any(normal):
        raise ValueError("原主齿没有普通渐开线端点，不能锚定局部沿线长度坐标")
    anchor_offsets = source_s_edges[normal]-source["endpoint_face_coordinate_m"][es][normal]/cosine
    if np.ptp(anchor_offsets) > 1e-11:
        raise ValueError("原普通接触段的累计线长与共同几何坐标不相容")
    source_s_edges -= np.median(anchor_offsets)
    rotation = result.rotations[ri]
    own = _project_line_samples(rotation, pair_index, cosine)
    own_mask, order = own["mask"], own["order"]
    own_width, own_intervals = own["widths_m"], own["intervals_m"]
    own_peaks, own_half = own["peak_pa"], own["half_width_m"]
    own_force, own_q = own["force_n"], own["line_load_n_per_m"]
    own_endpoints = own["endpoints"]
    ref_peak, ref_half = source["endpoint_pressure_pa"][es], source["endpoint_half_width_m"][es]
    ref_q = source["endpoint_line_load_n_per_m"][es]
    display_edges, source_peak_display, source_half_display = _endpoint_display(source_s_edges, ref_peak, ref_half)
    if own_endpoints is None:
        own_display_intervals, own_peak_display, own_half_display = own_intervals, own_peaks, own_half
    else:
        own_display_edges, own_peak_display, own_half_display = _endpoint_display(
            own["edges_m"], own_endpoints["pressure_pa"], own_endpoints["half_width_m"])
        own_display_intervals = np.column_stack((own_display_edges[:-1], own_display_edges[1:]))
    max_half = max(float(np.max(ref_half)), float(np.max(own_half_display)), float(np.max(own_half)))
    # 零载时只为显示选择一个有限坐标窗，不赋予有限接触面积。
    normal_edges = np.linspace(-max(max_half*1.15, 1e-8), max(max_half*1.15, 1e-8), 122)
    source_pressure = _band_values(source_peak_display, source_half_display, normal_edges)
    project_pressure = _band_values(own_peak_display, own_half_display, normal_edges)
    project_cell_pressure = _band_values(own_peaks, own_half, normal_edges)
    cycle_peak = []
    for r in range(count):
        ls, le = source["rotation_line_offsets"][r:r+2]
        endpoint_slice = slice(source["endpoint_offsets"][ls], source["endpoint_offsets"][le])
        cycle_peak.append(float(np.max(source["endpoint_pressure_pa"][endpoint_slice])))
    current_force = np.array([np.sum(r.solution.forces_n) for r in result.rotations])
    current_pressure = np.array([np.max(r.endpoint_pressure["pressure_pa"])
        if getattr(r, "endpoint_pressure", None) is not None else np.max(r.pressure_pa) for r in result.rotations])
    source_force = float(np.sum(source["centre_force_n"][cs]))
    source_band_force = float(np.trapz(ref_q, source_s_edges))
    own_total = float(np.sum(own_force))
    own_cell_band_force = float(np.sum(math.pi/2*own_peaks*own_half*own_width))
    own_band_force = (float(np.trapz(own_endpoints["line_load_n_per_m"], own["edges_m"]))
                      if own_endpoints is not None else own_cell_band_force)
    own_band_area = (float(np.trapz(2*own_endpoints["half_width_m"], own["edges_m"]))
                     if own_endpoints is not None else float(np.sum(2*own_half*own_width)))
    metrics = {
        "reference_phase_index": ri, "phase": float(result.phase[ri]), "project_tooth_pair_index": pair_index,
        "source_pinion_tooth_id": int(source["line_pinion_tooth_id"][line]),
        "source_wheel_tooth_id": int(source["line_wheel_tooth_id"][line]),
        "reference_tooth_force_n": source_force, "project_tooth_force_n": own_total,
        "reference_endpoint_peak_pressure_pa": float(np.max(ref_peak)),
        "project_cell_peak_pressure_pa": float(np.max(own_peaks)),
        "project_endpoint_peak_pressure_pa": float(np.max(own_endpoints["pressure_pa"])) if own_endpoints is not None else None,
        "project_band_evaluation": "endpoints" if own_endpoints is not None else "cell_centres",
        "reference_hertz_band_integral_n": source_band_force,
        "project_hertz_band_integral_n": own_band_force,
        "project_cell_hertz_band_integral_n": own_cell_band_force,
        "project_cell_band_minus_strip_force_n": own_cell_band_force-own_total,
        "reference_band_minus_strip_force_n": source_band_force-source_force,
        "project_band_minus_strip_force_n": own_band_force-own_total,
        "reference_hertz_band_area_m2": float(np.trapz(2*ref_half, source_s_edges)),
        "project_hertz_band_area_m2": own_band_area,
        "project_cell_hertz_band_area_m2": float(np.sum(2*own_half*own_width)),
        "project_ordinary_length_anchor_spread_m": own["ordinary_anchor_spread_m"],
        "maximum_phase_registration_error_cycles": float(np.max(np.abs(phase_error))),
        "source_ordinary_length_anchor_spread_m": float(np.ptp(anchor_offsets)),
        "reference_te_peak_to_peak_m": float(np.ptp(source["rotation_te_normal_m"])),
        "project_te_peak_to_peak_m": float(np.ptp(result.normal_approach_m)),
        "scope": "同几何相位的周期与同齿Hertz带；端点及中心指标分别保留，未作完整二维接触精度认证",
    }
    return {"source": source, "line_index": line, "centre_slice": cs, "endpoint_slice": es,
        "rotation_index": ri, "own_mask": own_mask, "own_order": order, "metrics": metrics,
        "source_s_edges_m": source_s_edges, "source_display_s_edges_m": display_edges,
        "normal_edges_m": normal_edges, "source_pressure_pa": source_pressure,
        "project_pressure_pa": project_pressure, "project_s_intervals_m": own_intervals,
        "project_cell_pressure_pa": project_cell_pressure, "project_display_intervals_m": own_display_intervals,
        "project_s_edges_m": own["edges_m"], "project_endpoints": own_endpoints,
        "source_peak_pa": ref_peak, "source_half_width_m": ref_half, "source_q_n_per_m": ref_q,
        "project_peak_pa": own_peaks, "project_half_width_m": own_half,
        "project_force_n": own_force, "project_q_n_per_m": own_q,
        "cycle": np.column_stack((np.arange(count)/count, source["rotation_phase"], result.phase,
            source["rotation_te_normal_m"], result.normal_approach_m,
            source["rotation_stiffness_n_per_m"], result.mesh_stiffness_n_per_m,
            cycle_peak, current_pressure, source["rotation_force_n"], current_force))}


def _draw_band(ax, intervals, normal_edges, values, *, cmap, norm):
    """每个纵向单元独立着色，未采样的纵向间隙保持背景，不跨空洞插值。"""
    from matplotlib.collections import PolyCollection
    intervals = np.asarray(intervals)
    left = np.repeat(intervals[:, 0]*1e3, values.shape[1])
    right = np.repeat(intervals[:, 1]*1e3, values.shape[1])
    bottom = np.tile(normal_edges[:-1]*1e6, len(intervals))
    top = np.tile(normal_edges[1:]*1e6, len(intervals))
    vertices = np.stack((np.column_stack((left, bottom)), np.column_stack((right, bottom)),
                         np.column_stack((right, top)), np.column_stack((left, top))), axis=1)
    collection = PolyCollection(vertices, array=np.asarray(values).ravel()/1e6,
                                cmap=cmap, norm=norm, edgecolors="none", rasterized=True)
    ax.add_collection(collection)
    return collection


def plot_case_contact(result, metadata, reference_arrays, output, *, display_cycles=3):
    """生成三图及SI复绘数据；周期图默认将单个稳态周期展开为三个。"""
    data = prepare_case_contact(result, metadata, reference_arrays)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    source, cycle, metrics = data["source"], data["cycle"], data["metrics"]
    displayed = repeat_cycle_for_display(cycle, display_cycles)
    source_samples = np.arange(len(displayed)) % len(cycle)
    period_indices = np.arange(len(displayed)) // len(cycle)
    closing = np.zeros(len(displayed), dtype=bool)
    closing[-1] = True
    files = []
    columns = ("cycle_fraction", "reference_phase", "project_phase", "reference_te_normal_m", "project_te_normal_m",
        "reference_stiffness_n_per_m", "project_stiffness_n_per_m", "reference_endpoint_peak_pa",
        "project_peak_pressure_pa", "reference_normal_force_n", "project_normal_force_n")
    files.append(_rows(output/"03_ltca_cycle.csv", columns, cycle))
    files.append(_rows(output/"03_ltca_cycle_display.csv",
        ("display_cycle_position", *columns[1:], "source_sample_index", "display_period_index", "periodic_closure"),
        zip(*displayed.T, source_samples, period_indices, closing)))
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.4))
    fig.subplots_adjust(left=.09, right=.97, bottom=.21, top=.89, hspace=.48, wspace=.30)
    panels = ((3, 4, 1e6, "a  承载法向TE", "接近量 / μm"),
              (5, 6, 1e-6, "b  法向刚度", "刚度 / (N/μm)"),
              (7, 8, 1e-6, "c  每相位峰值接触压力", "压力 / MPa"),
              (9, 10, 1., "d  总法向载荷", "载荷 / N"))
    for ax, (ref, own, scale, title, ylabel) in zip(axes.flat, panels):
        ax.plot(displayed[:, 0], displayed[:, ref]*scale, color="#C56F33", linewidth=1.4,
                marker="o", markersize=3, markerfacecolor="white", label="MASTA")
        ax.plot(displayed[:, 0], displayed[:, own]*scale, color="#176B91", linewidth=1.6,
                marker=".", markersize=4, label="本项目")
        ax.set(title=title, xlabel="相对起采位置 / 基节周期", ylabel=ylabel, xlim=(0, display_cycles))
        ax.set_xticks(np.arange(display_cycles+1))
        ax.grid(alpha=.7)
    force_max = max(float(np.max(cycle[:, 9])), float(np.max(cycle[:, 10])))
    axes[1, 1].set_ylim(0., max(1., force_max*1.12))
    axes[1, 1].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[0, 0].legend(loc="upper center", bbox_to_anchor=(1.15, 1.5), ncol=2, fontsize=9)
    fig.text(.5, .050, f"18/77斜齿轮副：相同载荷与工作面下的LTCA周期响应。\n"
        f"展示{display_cycles}个稳态啮合周期，每周期{len(cycle)}个采样位置。",
        ha="center", fontsize=8.5, color="#526470")
    files.extend(_save(fig, output, "03_ltca_cycle"))
    plt.close(fig)
    from matplotlib.colors import Normalize
    project_peak = (metrics["project_endpoint_peak_pressure_pa"] if metrics["project_endpoint_peak_pressure_pa"] is not None
                    else metrics["project_cell_peak_pressure_pa"])
    peak = max(metrics["reference_endpoint_peak_pressure_pa"], project_peak)/1e6
    norm, cmap = Normalize(0., max(peak, 1e-9)), plt.get_cmap("coolwarm")
    fig, axes = plt.subplots(2, 1, figsize=(10.7, 6.9), sharex=True, sharey=True)
    fig.subplots_adjust(left=.10, right=.84, bottom=.19, top=.90, hspace=.58)
    source_intervals = np.column_stack((data["source_display_s_edges_m"][:-1], data["source_display_s_edges_m"][1:]))
    _draw_band(axes[0], source_intervals, data["normal_edges_m"], data["source_pressure_pa"], cmap=cmap, norm=norm)
    image = _draw_band(axes[1], data["project_display_intervals_m"], data["normal_edges_m"], data["project_pressure_pa"], cmap=cmap, norm=norm)
    bounds = (min(source_intervals[:, 0].min(), data["project_display_intervals_m"][:, 0].min()),
              max(source_intervals[:, 1].max(), data["project_display_intervals_m"][:, 1].max()))
    for ax in axes:
        ax.set_facecolor("#E7EBEF")
        ax.set(xlim=np.asarray(bounds)*1e3, ylim=data["normal_edges_m"][[0, -1]]*1e6,
               ylabel="垂直接触线 n / μm", xlabel="沿接触线 s / mm")
        ax.axhline(0, color="white", linewidth=.45, alpha=.6)
    axes[0].set_title(f"a  MASTA：齿对{metrics['source_pinion_tooth_id']}/{metrics['source_wheel_tooth_id']}；"
                      f"条带力 {metrics['reference_tooth_force_n']:.3f} N；峰压 {metrics['reference_endpoint_peak_pressure_pa']/1e6:.1f} MPa", fontsize=10)
    axes[1].set_title(f"b  本项目：齿对{metrics['source_pinion_tooth_id']}/{metrics['source_wheel_tooth_id']}；"
                      f"条带力 {metrics['project_tooth_force_n']:.3f} N；峰压 {project_peak/1e6:.1f} MPa", fontsize=10)
    bar = fig.add_axes((.88, .27, .022, .52))
    fig.colorbar(image, cax=bar, label="局部接触压力 / MPa")
    fig.suptitle("同相位、同齿对的瞬时接触压力", fontsize=13)
    fig.text(.5, .065, "沿接触线长度单位为mm，垂直接触线宽度单位为μm，两轴显示比例不同。\n"
        "两图采用相同压力色标；灰色为未覆盖区域，蓝色为零压力。",
        ha="center", fontsize=8.2, color="#526470")
    files.extend(_save(fig, output, "04_contact_pattern"))
    plt.close(fig)
    files.append(_rows(output/"04_reference_endpoints.csv",
        ("s_m", "peak_pressure_pa", "half_width_m", "line_load_n_per_m"),
        zip(data["source_s_edges_m"], data["source_peak_pa"], data["source_half_width_m"], data["source_q_n_per_m"])))
    files.append(_rows(output/"04_project_cells.csv",
        ("s_left_m", "s_right_m", "peak_pressure_pa", "half_width_m", "line_load_n_per_m", "force_n"),
        zip(data["project_s_intervals_m"][:, 0], data["project_s_intervals_m"][:, 1], data["project_peak_pa"],
            data["project_half_width_m"], data["project_q_n_per_m"], data["project_force_n"])))
    if data["project_endpoints"] is not None:
        endpoints = data["project_endpoints"]
        files.append(_rows(output/"04_project_endpoints.csv",
            ("s_m", "peak_pressure_pa", "half_width_m", "line_load_n_per_m"),
            zip(data["project_s_edges_m"], endpoints["pressure_pa"], endpoints["half_width_m"],
                endpoints["line_load_n_per_m"])))
    own_offsets = np.r_[0, np.cumsum([len(r.grid.widths_m) for r in result.rotations])]
    project_face = np.concatenate([r.grid.face_coordinate_m for r in result.rotations])
    project_roll = np.concatenate([r.grid.roll_pinion_m for r in result.rotations])
    project_q = np.concatenate([r.line_load_n_per_m for r in result.rotations])
    project_force = np.concatenate([r.solution.forces_n for r in result.rotations])
    project_tooth = np.concatenate([r.grid.tooth_pair_index for r in result.rotations])
    project_extended = np.concatenate([
        (r.contact_geometry_data.centre_extended_pinion | r.contact_geometry_data.centre_extended_wheel)
        if getattr(r, "contact_geometry_data", None) is not None else np.zeros(len(r.grid.widths_m), dtype=bool)
        for r in result.rotations])
    ri = data["rotation_index"]
    first, last = source["rotation_line_offsets"][ri:ri+2]
    source_instant = slice(source["centre_offsets"][first], source["centre_offsets"][last])
    own_instant = slice(own_offsets[ri], own_offsets[ri+1])
    source_extended = source["centre_extended_pinion"] | source["centre_extended_wheel"]
    limit = max(float(np.max(source["centre_line_load_n_per_m"])), float(np.max(project_q)))/1e3
    qnorm = Normalize(0., max(limit, 1e-12))
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.09, right=.86, bottom=.17, top=.91, hspace=.40, wspace=.18)
    coordinate_data = ((0, 0, source["centre_face_coordinate_m"][source_instant], source["centre_roll_pinion_m"][source_instant],
                       source["centre_line_load_n_per_m"][source_instant], source_extended[source_instant], "a  MASTA：同一瞬时"),
                       (0, 1, project_face[own_instant], project_roll[own_instant], project_q[own_instant], project_extended[own_instant], "b  本项目：同一瞬时"),
                       (1, 0, source["centre_face_coordinate_m"], source["centre_roll_pinion_m"],
                        source["centre_line_load_n_per_m"], source_extended, "c  MASTA：周期离散采样位置"),
                       (1, 1, project_face, project_roll, project_q, project_extended, "d  本项目：周期离散采样位置"))
    for row, col, x, y, q, extended, title in coordinate_data:
        ax = axes[row, col]
        extended = np.zeros(len(x), dtype=bool) if extended is None else extended
        for mask, marker in ((~extended, "o"), (extended, "^")):
            if np.any(mask):
                image = ax.scatter(x[mask]*1e3, y[mask]*1e3, c=q[mask]/1e3, s=12 if row == 0 else 3,
                    marker=marker, cmap=cmap, norm=qnorm, linewidths=0, alpha=.9 if row == 0 else .7)
        for bound in result.pair.pinion.involute_roll_bounds_m:
            ax.axhline(bound*1e3, color="#667A88", linestyle="--", linewidth=.65)
        ax.set(title=title, xlabel="齿宽坐标 z / mm", ylabel="小轮展长 q / mm")
        ax.grid(alpha=.3)
    bar = fig.add_axes((.90, .25, .018, .55))
    fig.colorbar(image, cax=bar, label="条带中心线载荷 / (N/mm)")
    fig.text(.5, .065, "圆点为普通接触，三角为延伸接触，虚线为有效齿廓边界。\n"
        "上方为单相位接触位置，下方为一个啮合周期的接触覆盖。",
        ha="center", fontsize=8.3, color="#526470")
    files.extend(_save(fig, output, "05_contact_coordinates"))
    plt.close(fig)
    plot_data = dict(cycle=cycle, cycle_display=displayed,
        cycle_display_source_sample_index=source_samples, cycle_display_period_index=period_indices,
        cycle_display_periodic_closure=closing,
        source_s_edges_m=data["source_s_edges_m"], source_display_s_edges_m=data["source_display_s_edges_m"],
        normal_edges_m=data["normal_edges_m"], source_band_pressure_pa=data["source_pressure_pa"],
        project_band_pressure_pa=data["project_pressure_pa"], project_s_intervals_m=data["project_s_intervals_m"],
        project_display_intervals_m=data["project_display_intervals_m"],
        project_cell_band_pressure_pa=data["project_cell_pressure_pa"],
        source_face_coordinate_m=source["centre_face_coordinate_m"], source_roll_pinion_m=source["centre_roll_pinion_m"],
        source_line_load_n_per_m=source["centre_line_load_n_per_m"], source_extended=source_extended,
        source_centre_offsets=source["centre_offsets"], source_rotation_line_offsets=source["rotation_line_offsets"],
        project_face_coordinate_m=project_face, project_roll_pinion_m=project_roll,
        project_line_load_n_per_m=project_q, project_force_n=project_force, project_extended=project_extended,
        project_tooth_pair_index=project_tooth, project_rotation_cell_offsets=own_offsets)
    if data["project_endpoints"] is not None:
        plot_data.update(project_endpoint_s_m=data["project_s_edges_m"],
            project_endpoint_pressure_pa=data["project_endpoints"]["pressure_pa"],
            project_endpoint_half_width_m=data["project_endpoints"]["half_width_m"],
            project_endpoint_line_load_n_per_m=data["project_endpoints"]["line_load_n_per_m"])
    np.savez_compressed(output/"contact_plot_data.npz", **plot_data)
    files.append(str((output/"contact_plot_data.npz").resolve()))
    write_json(output/"contact_comparison_metrics.json", metrics)
    files.append(str((output/"contact_comparison_metrics.json").resolve()))
    return {"files": files, "metrics": metrics,
            "cycle_display": {"periods": int(display_cycles), "original_sample_count": len(cycle),
                "display_sample_count": len(displayed), "periodic_closing_point_added": True,
                "definition": "同一稳态周期原样重复，物理相位及响应不变；末端闭合点取原起点，不是多周期时域求解"},
            "reference_values_recomputed": False, "project_values_recomputed": True,
            "complete_2d_contact_accuracy_verified": False}
