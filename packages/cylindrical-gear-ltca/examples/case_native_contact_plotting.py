"""新写整齿通过程压力展示：MASTA原生位图与本项目独立端点计算分别保留。"""
from __future__ import annotations
import csv
import hashlib
from pathlib import Path
import shutil
import numpy as np
from gear_contact.io import write_json
from gear_contact.plotting import _pyplot, _save
from gear_contact.contact_display import (contact_pass_data, draw_contact_pass,
    project_contact_pass_samples as _project_contact_pass_samples)


def project_contact_pass_samples(result):
    """收集单个稳态周期全部接触线端点，坐标为齿宽左端起算/小轮展长，单位m。"""
    return _project_contact_pass_samples(result)


def plot_native_contact_comparison(result, metadata, image_path, surface_path, output):
    """原图直接复制；本项目仅从自身端点压力插值，不把MASTA压力用于求解或拟合。"""
    output, image_path, surface_path = Path(output), Path(image_path), Path(surface_path)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    if (digest(image_path) != metadata["image_sha256"] or digest(surface_path) != metadata["arrays_sha256"]
            or metadata["reference_name"] != "MASTA" or metadata["gear"] != "pinion"):
        raise ValueError("MASTA原生图或显示面身份与元数据不一致")
    with np.load(surface_path, allow_pickle=False) as source:
        x, y = source["face_coordinate_m"].copy(), source["roll_distance_m"].copy()
        reference = source["pressure_pa_x_y"].copy()
    points = project_contact_pass_samples(result)
    # 同一材料点可在多个接触线出现；逐点最大值定义通过程包络，原端点值另外导出。
    # 齿顶延伸点属于独立接触分支，其展长可能与普通齿廓点重合。
    # 齿廓图仅插值双方均未延伸的点，不能把齿顶分支叠到渐开线面上。
    profile_points = points[~np.any(points[:, 3:5].astype(bool), axis=1)]
    display = contact_pass_data(result, face_coordinate_m=x, roll_distance_m=y)
    project = display['pressure_pa_x_y']
    maximum = max(float(np.nanmax(reference)), 1.)
    plt = _pyplot()

    def draw(ax):
        return draw_contact_pass(ax, display, maximum_pa=maximum)

    output.mkdir(parents=True, exist_ok=True)
    native_name = "07_masta_native_contact.png"
    shutil.copyfile(image_path, output/native_name)
    files = [native_name]
    fig, ax = plt.subplots(figsize=(8., 6.2))
    fig.subplots_adjust(left=.11, right=.84, bottom=.12, top=.91)
    image = draw(ax)
    ax.set_title("本项目：小轮整齿通过程接触压力")
    fig.colorbar(image, ax=ax, label="接触压力 / MPa", fraction=.045, pad=.03)
    files.extend(Path(p).name for p in _save(fig, output, "07_project_contact_pass"))
    plt.close(fig)
    fig = plt.figure(figsize=(14., 5.9))
    native_ax = fig.add_axes((.005, .10, .49, .83))
    native_ax.imshow(plt.imread(image_path), interpolation="none")
    native_ax.axis("off")
    ax = fig.add_axes((.55, .141, .39, .758))
    image = draw(ax)
    fig.colorbar(image, cax=fig.add_axes((.96, .27, .012, .49)), label="接触压力 / MPa")
    fig.text(.25, .96, "a  MASTA", ha="center")
    fig.text(.74, .96, "b  本项目", ha="center")
    fig.text(.5, .025, "18/77斜齿轮副：相同载荷与工作面下，一个齿通过程的接触压力分布。\n"
        "横轴为齿宽坐标，纵轴为小轮展长；灰色为未覆盖区域。",
        ha="center", fontsize=8.3, color="#526470")
    files.extend(Path(p).name for p in _save(fig, output, "07_native_contact_comparison"))
    plt.close(fig)
    with (output/"07_project_contact_pass_samples.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(("face_coordinate_from_left_m", "roll_pinion_m", "endpoint_pressure_pa", "extended_pinion", "extended_wheel"))
        writer.writerows(points)
    np.savez_compressed(output/"07_contact_pass_data.npz", face_coordinate_m=x, roll_distance_m=y,
        reference_pressure_pa_x_y=reference, project_pressure_pa_x_y=project, project_endpoint_samples=points)
    summary = {"native_image_sha256": digest(output/native_name),
        "native_image_unchanged": digest(output/native_name) == metadata["image_sha256"],
        "definition": metadata["definition"], "grid_shape_x_y": list(project.shape),
        "reference_display_peak_pa": float(np.nanmax(reference)), "project_display_peak_pa": float(np.nanmax(project)),
        "project_raw_endpoint_peak_pa": float(points[:, 2].max()), "project_raw_endpoint_count": len(points),
        "project_profile_endpoint_count": len(profile_points),
        "project_uncovered_display_count": int(np.count_nonzero(~np.isfinite(project))),
        "project_interpolation": "双方未延伸的普通齿廓端点材料坐标上Delaunay线性三角插值；重复点取通过程最大值，延伸分支另存，凸包外不填零",
        "scope": "显示面用于图义、坐标和分布比较；显示插值误差不替代原端点TE、刚度或峰压误差。",
        "files": files+["07_project_contact_pass_samples.csv", "07_contact_pass_data.npz", "07_native_contact_metrics.json"]}
    write_json(output/"07_native_contact_metrics.json", summary)
    return summary
