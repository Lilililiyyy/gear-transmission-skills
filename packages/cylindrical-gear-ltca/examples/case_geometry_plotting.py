"""新写案例几何展示工具：从当前 Gear/GearPair 输入生成端面齿廓与输入示意。

齿廓由 tooth_thickness_m 的理想渐开线定义生成，基圆以下沿用该接口的
径向延拓；不生成刀具根部圆角、不判断干涉，也不以图形冒充接触结果。
"""
from __future__ import annotations

import csv
from dataclasses import asdict
import math
from pathlib import Path

import numpy as np

from gear_contact.io import write_json


def tooth_outline(gear, *, phase_rad=0., radial_samples=40):
    """返回单轮闭合端面齿廓的米制坐标，根圆至顶圆逐齿保留原几何。

    外齿为外边界，内齿为内圈齿形边界；内齿轮缘外径没有输入，不能由
    此多段线补造。根部间隙使用根圆圆弧连接，不表示真实过渡曲线。
    """
    if isinstance(radial_samples, bool) or not isinstance(radial_samples, int) or radial_samples < 4:
        raise ValueError("端面齿廓径向采样数须为不小于4的整数")
    if not math.isfinite(phase_rad):
        raise ValueError("几何展示相位必须为有限弧度")
    radii = np.linspace(gear.root_radius_m, gear.tip_radius_m, radial_samples)
    halves = np.asarray(gear.tooth_thickness_m(radii)) / (2 * radii)
    pitch = 2 * math.pi / gear.teeth
    points, tooth_ids, parts = [], [], []

    def append(rr, aa, tooth, part):
        rr, aa = np.broadcast_arrays(rr, aa)
        points.extend(np.column_stack((rr * np.cos(aa), rr * np.sin(aa))))
        tooth_ids.extend([tooth] * rr.size)
        parts.extend([part] * rr.size)

    for tooth in range(gear.teeth):
        centre = phase_rad + tooth * pitch
        append(radii, centre - halves, tooth, "负角侧齿廓")
        append(gear.tip_radius_m, np.linspace(centre-halves[-1], centre+halves[-1], 12), tooth, "齿顶圆弧")
        append(radii[::-1], centre + halves[::-1], tooth, "正角侧齿廓")
        append(gear.root_radius_m, np.linspace(centre+halves[0], centre+pitch-halves[0], 12), tooth, "根圆连接弧")
    points = np.asarray(points)
    points = np.vstack((points, points[0]))
    return {"xy_m": points, "tooth_index": np.r_[tooth_ids, tooth_ids[0]],
            "part": parts + [parts[0]], "phase_rad": float(phase_rad)}


def display_phases(pair):
    """在工作节圆位置对齐两条渐开线；这是几何展示位置，不是LTCA相位。

    若实际根顶范围不包含工作节圆，不强制推断几何相切，改用明确的
    任意展示相位。单个齿宽截面无法代替完整斜齿三维装配。
    """
    cosine = math.cos(pair.working_pressure_angle_rad)
    rp = pair.pinion.base_radius_m / cosine
    rw = pair.wheel.base_radius_m / cosine
    for gear, radius in ((pair.pinion, rp), (pair.wheel, rw)):
        if not min(gear.root_radius_m, gear.tip_radius_m) <= radius <= max(gear.root_radius_m, gear.tip_radius_m):
            return {"pinion_phase_rad": 0., "wheel_phase_rad": math.pi,
                    "display_method": "工作节圆超出实体根顶范围，采用任意展示相位，未配准接触",
                    "pitch_point_m": [-rp if pair.wheel.internal else rp, 0.]}
    hp = float(pair.pinion.tooth_thickness_m(rp)) / (2 * rp)
    hw = float(pair.wheel.tooth_thickness_m(rw)) / (2 * rw)
    if pair.wheel.internal:
        pp, pw, x = math.pi-hp, math.pi+hw, -rp
    else:
        pp, pw, x = -hp, math.pi-hw, rp
    return {"pinion_phase_rad": pp, "wheel_phase_rad": pw,
            "display_method": "在工作节圆点对齐理论渐开线，仅为端面几何展示",
            "pitch_point_m": [x, 0.]}


def plot_case_geometry(pair, output, *, modification_pinion=None, modification_wheel=None):
    """输出端面、局部齿形、实际齿宽和可选修形输入的PNG/SVG及SI数据。

    pair 为公开 GearPair。修形对象须有 evaluate(profile_factor,lead_factor)
    接口；其曲面若包含等效边界间隙，本图也按输入原值保留。不传入时
    明确记为未提供，不将其解释成无修形。输出不含压力、
    接触斑点或制造根部过渡。两轮孔型和轮缘外边界也不由本工具推断。
    """
    from gear_contact.plotting import _pyplot, _save

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    phases = display_phases(pair)
    curves = {}
    roles = (("pinion", "小轮", pair.pinion, 0., "#176B91"),
             ("wheel", "大轮", pair.wheel, pair.center_distance_m, "#D47732"))
    with (output / "geometry_outlines.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("gear", "vertex_index", "tooth_index", "part", "x_m", "y_m"))
        for role, label, gear, centre, color in roles:
            curve = tooth_outline(gear, phase_rad=phases[f"{role}_phase_rad"])
            curve["xy_m"][:, 0] += centre
            curves[role] = curve
            writer.writerows((role, i, int(t), part, float(x), float(y))
                             for i, (t, part, (x, y)) in enumerate(zip(curve["tooth_index"], curve["part"], curve["xy_m"])))

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.7))
    fig.subplots_adjust(left=.07, right=.97, bottom=.14, top=.91, hspace=.45, wspace=.32)
    for ax in axes[0]:
        for role, label, gear, centre, color in roles:
            xy = curves[role]["xy_m"] * 1e3
            ax.plot(xy[:, 0], xy[:, 1], color=color, linewidth=.9, label=f"{label} z={gear.teeth}")
            angle = np.linspace(0., 2*math.pi, 500)
            for radius, linestyle in ((gear.reference_radius_m, "--"), (gear.base_radius_m, ":")):
                ax.plot((centre+radius*np.cos(angle))*1e3, radius*np.sin(angle)*1e3,
                        color=color, alpha=.45, linestyle=linestyle, linewidth=.65)
            ax.plot(centre*1e3, 0., marker="+", color=color, markersize=7)
        ax.set_aspect("equal", adjustable="box")
        ax.set(xlabel="端面 x / mm", ylabel="端面 y / mm")
        ax.grid(alpha=.55)
    axes[0, 0].set_title("a  输入几何的端面齿形与参考圆／基圆")
    axes[0, 0].legend(loc="upper right", fontsize=9)
    pitchpoint = np.asarray(phases["pitch_point_m"]) * 1e3
    extent = 4.3 * pair.pinion.normal_module_m * 1e3
    axes[0, 1].set_xlim(pitchpoint[0]-extent, pitchpoint[0]+extent)
    axes[0, 1].set_ylim(-extent, extent)
    axes[0, 1].plot(*pitchpoint, "o", color="#263E4A", markersize=3)
    axes[0, 1].set_title("b  工作节圆点附近的理论齿廓（同尺度）")

    width_ax = axes[1, 0]
    for index, (role, label, gear, centre, color) in enumerate(roles):
        shift = 0. if role == "pinion" else pair.axial_offset_m
        left, right = shift-gear.face_width_m/2, shift+gear.face_width_m/2
        width_ax.plot(np.array((left, right))*1e3, (index, index), color=color, linewidth=14, alpha=.2)
        width_ax.plot(np.array((left, left, right, right))*1e3, (index-.15, index+.15, index+.15, index-.15),
                      color=color, linewidth=.9)
        width_ax.annotate(f"b={gear.face_width_m*1e3:.3f} mm；β={math.degrees(gear.helix_angle_rad):+.3f}°",
                          ((left+right)*500, index+.22), ha="center", va="bottom", fontsize=9)
    overlap_left, overlap_right = pair.face_bounds_m
    width_ax.axvspan(overlap_left*1e3, overlap_right*1e3, color="#8BB3A7", alpha=.12)
    width_ax.set(xlabel="轴向位置 / mm", title="c  实际齿宽与轴偏；竖向只区分两轮",
                 yticks=(0, 1), yticklabels=("小轮", "大轮"), ylim=(-.4, 1.6))
    width_ax.grid(axis="x")

    factors = np.linspace(0., 1., 201)
    modification_ax = axes[1, 1]
    modifications = (("pinion", "小轮", modification_pinion, "#176B91"),
                     ("wheel", "大轮", modification_wheel, "#D47732"))
    modifications_data, supplied = [], False
    for role, label, modification, color in modifications:
        if modification is None:
            continue
        supplied = True
        for axis, linestyle in (("profile", "-"), ("lead", "--")):
            p, l = (factors, np.full_like(factors, .5)) if axis == "profile" else (np.full_like(factors, .5), factors)
            values = np.broadcast_to(np.asarray(modification.evaluate(p, l), float), factors.shape)
            if not np.all(np.isfinite(values)):
                raise ValueError("修形输入示意需要有限法向去除量，单位m")
            modification_ax.plot(factors, values*1e6, linestyle=linestyle, color=color,
                                  label=f"{label}{'齿廓' if axis == 'profile' else '齿向'}中线")
            modifications_data.extend((role, axis, float(f), float(v)) for f, v in zip(factors, values))
    modification_ax.set(xlabel="对应完整齿廓／齿宽的归一化因子", ylabel="法向去除／间隙 / μm",
                        title="d  输入曲面的两条中线（含原零基准）")
    modification_ax.grid()
    if supplied:
        modification_ax.legend(fontsize=8, ncol=2)
    else:
        modification_ax.text(.5, .5, "未传入修形对象；不推断为无修形", transform=modification_ax.transAxes,
                              ha="center", va="center", fontsize=10)
    with (output / "geometry_modification_sections.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("gear", "axis", "factor", "normal_removal_m"))
        writer.writerows(modifications_data)

    fig.suptitle(f"圆柱齿轮单案例输入几何：{pair.pinion.teeth}/{pair.wheel.teeth} 齿，"
                 f"法向模数={pair.pinion.normal_module_m*1e3:.3f} mm，a={pair.center_distance_m*1e3:.3f} mm", fontsize=14)
    fig.text(.5, .045, "端面为理想渐开线；基圆以下径向延拓、根圆弧连接，不表示刀具过渡或干涉判定。\n"
             "虚线为参考圆、点线为基圆；展示相位不等于LTCA相位。曲面若含等效边界间隙，也按输入原值保留。",
             ha="center", va="center", fontsize=8.5, color="#526470")
    files = _save(fig, output, "01_case_geometry")
    plt.close(fig)
    record = {"implementation": "新写输入几何展示工具", "pair": asdict(pair), "display_phases": phases,
              "coordinates_unit": "m", "plot_coordinates_unit": "mm", "modification_unit": "m",
              "modification_sections_supplied": {"pinion": modification_pinion is not None,
                                                  "wheel": modification_wheel is not None},
              "geometry_limits": ["理想端面渐开线和给定根顶圆", "基圆以下径向延拓、根圆连接，不是真实刀具过渡",
                                  "显示相位不对应LTCA滚动相位", "未给轮体孔型及内齿圈轮缘外径", "不作根切、干涉或接触载荷判断"],
              "figures": [Path(p).name for p in files],
              "data_files": ["geometry_outlines.csv", "geometry_modification_sections.csv"]}
    write_json(output / "geometry_display.json", record)
    return record
