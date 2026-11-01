"""运行合成案例、解析基准、多工况设计和图表导出。

本文件为新写案例编排。全部案例调用同一计算内核，正式计算参数由本文件
中的构造函数给出；没有安装目录、外部工程或固定工作区路径依赖。
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, is_dataclass
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

# 仅定位随本案例分发的包，不查找宿主工作区。
_PACKAGE_DIRECTORY = Path(__file__).resolve().parent.parent
if str(_PACKAGE_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(_PACKAGE_DIRECTORY))

from gear_contact.analysis import analyze_pair
from gear_contact.design import OperatingPoint, robust_design_search
from gear_contact.elasticity import Material
from gear_contact.geometry import Gear, GearPair
from gear_contact.hertz import hertz_line_contact
from gear_contact.io import export_result
from gear_contact.plotting import (analysis_metrics, plot_analysis, plot_comparison,
                                   plot_design_study, plot_surface, _pyplot, _save)
from gear_contact.solver import solve_contact
from gear_contact.surfaces import SurfaceModification, MeasuredSurface, CompositeModification


def json_value(value):
    """将案例的输入和结果转成可移植 JSON，拒绝隐藏非有限数。"""
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def standard_pair(*, helix_degrees=0., internal=False):
    """标准合成几何；内啮合采用 28/84 齿，外啮合采用 28/56 齿。"""
    beta = math.radians(helix_degrees)
    pinion = Gear(28, .0025, .025, pressure_angle_rad=math.radians(20), helix_angle_rad=beta)
    wheel = Gear(84 if internal else 56, .0025, .025,
                 pressure_angle_rad=math.radians(20),
                 helix_angle_rad=beta if internal else -beta, internal=internal)
    return GearPair(pinion, wheel)


def run_analytic_baselines(output: Path):
    """用分段线性解析解和 Hertz 压力积分检查独立物理量。"""
    output.mkdir(parents=True, exist_ok=True)
    gaps = np.array([0., 2e-6, 5e-6])
    compliance = np.eye(3) * 1e-8
    result = solve_contact(gaps, compliance, 600.)
    expected_forces = np.array([400., 200., 0.])
    np.testing.assert_allclose(result.forces_n, expected_forces, rtol=1e-8, atol=1e-6)
    np.testing.assert_allclose(result.approach_m, 4e-6, rtol=1e-8, atol=1e-13)
    np.testing.assert_allclose(result.tangent_stiffness_n_per_m, 2e8, rtol=1e-7)
    material = Material(210e9, .3)
    q, r1, r2 = 5e4, .01, .02
    curvature = 1 / r1 + 1 / r2
    contact = hertz_line_contact(q, curvature, material, material)
    half_width = float(np.asarray(contact.half_width_m))
    pressure = float(np.asarray(contact.maximum_pressure_pa))
    effective_modulus = 1 / (2 * (1 - .3 ** 2) / 210e9)
    expected_half = math.sqrt(4 * q / (math.pi * effective_modulus * curvature))
    np.testing.assert_allclose(half_width, expected_half, rtol=1e-12)
    x = np.linspace(-half_width, half_width, 20001)
    p = pressure * np.sqrt(np.maximum(0, 1 - (x / half_width) ** 2))
    integrated_q = float(np.trapz(p, x))
    assert abs(integrated_q / q - 1) < 1e-6, "Hertz 压力积分未恢复给定线载荷"
    record = {
        "status": "passed",
        "linear_contact": {"gaps_m": gaps, "compliance_m_per_n": compliance,
                           "normal_force_n": 600., "expected_forces_n": expected_forces,
                           "computed_forces_n": result.forces_n,
                           "approach_m": result.approach_m,
                           "tangent_stiffness_n_per_m": result.tangent_stiffness_n_per_m},
        "hertz_contact": {"line_load_n_per_m": q, "radius_a_m": r1, "radius_b_m": r2,
                          "effective_modulus_pa": effective_modulus,
                          "half_width_m": half_width, "maximum_pressure_pa": pressure,
                          "integrated_line_load_n_per_m": integrated_q,
                          "integration_relative_error": abs(integrated_q / q - 1)},
    }
    write_json(output / "analytic_baselines.json", record)
    with (output / "linear_contact.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("cell_index", "initial_gap_m", "expected_force_n", "computed_force_n", "residual_gap_m"))
        writer.writerows(zip(range(3), gaps, expected_forces, result.forces_n, result.residual_gaps_m))
    with (output / "hertz_pressure_profile.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("normal_coordinate_m", "pressure_pa"))
        writer.writerows(zip(x, p))
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    positions = np.arange(3)
    axes[0].bar(positions - .16, expected_forces, width=.3, label="解析解", color="#B87936")
    axes[0].bar(positions + .16, result.forces_n, width=.3, label="计算解", color="#286B83")
    axes[0].set(xticks=positions, xticklabels=("单元1", "单元2", "单元3"),
                ylabel="法向接触力 / N", title="a  非负接触与力平衡解析基准")
    axes[0].legend()
    axes[1].plot(x * 1e6, p / 1e6, color="#286B83", linewidth=1.7)
    axes[1].fill_between(x * 1e6, p / 1e6, color="#286B83", alpha=.13)
    axes[1].set(xlabel="接触带横向坐标 / μm", ylabel="接触压力 / MPa",
                title="b  Hertz 半椭圆压力积分恢复线载荷")
    _save(fig, output, "analytic_baselines")
    plt.close(fig)
    return record


def run_examples(output_dir, *, quick=False, phases=None, strips=None):
    """运行全部案例；快速模式减小离散，仅用于安装与调用检查。"""
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    phase_count = phases if phases is not None else (8 if quick else 32)
    strip_count = strips if strips is not None else (9 if quick else 21)
    if phase_count < 4 or strip_count < 3:
        raise ValueError("案例至少需要 4 个相位和 3 个条带")
    started = time.perf_counter()
    manifest = {"module": "圆柱齿轮准静态承载齿面接触分析（LTCA）", "implementation": "新写物理模型实现",
                "input_origin": "参数化合成案例", "physical_experiments": "未进行",
                "quick_mode": quick, "phases": phase_count, "strips": strip_count,
                "cases": [], "figures": []}
    analytic = run_analytic_baselines(output / "01_analytic_baselines")
    manifest["analytic_status"] = analytic["status"]
    pair = standard_pair()
    steel = Material(210e9, .3)
    kwargs = dict(material_pinion=steel, material_wheel=steel,
                  phases=phase_count, strips=strip_count, elastic_mode="coupled", local_contact=True)

    def calculate(name, label, *, gear_pair=pair, torque=80., modification=None,
                  misalignment=0., make_plots=False, custom_strips=None):
        params = dict(kwargs)
        if custom_strips is not None:
            params["strips"] = custom_strips
        surface = SurfaceModification() if modification is None else modification
        case_input = {"pair": gear_pair, "torque_n_m": torque,
                      "material_pinion": steel, "material_wheel": steel,
                      "modification_pinion": surface, "misalignment_m": misalignment,
                      "phases": params["phases"], "strips": params["strips"],
                      "elastic_mode": params["elastic_mode"], "local_contact": True}
        result = analyze_pair(gear_pair, torque_n_m=torque, modification_pinion=surface,
                              misalignment_m=misalignment, **params)
        destination = output / name
        export_result(result, destination)
        write_json(destination / "case_input.json", case_input)
        row = {"case": name, "label": label, **analysis_metrics(result)}
        manifest["cases"].append(row)
        if make_plots:
            files = plot_analysis(result, destination, title=label)
            manifest["figures"].extend(str(Path(p).relative_to(output)) for p in files
                                       if Path(p).suffix in (".png", ".svg"))
        return result

    baseline = calculate("02_spur_baseline", "无修形直齿轮副", make_plots=True)
    crown = SurfaceModification(lead_crown_m=10e-6)
    crowned = calculate("03_lead_crown", "齿向鼓形：端部 10 μm", modification=crown)
    combined_surface = SurfaceModification(lead_crown_m=8e-6, tip_relief_m=6e-6,
                                           tip_start=.72, root_relief_m=3e-6, root_end=.18)
    combined = calculate("04_combined_modification", "齿廓与齿向组合修形",
                         modification=combined_surface, make_plots=True)
    plot_surface(combined_surface, output / "04_combined_modification")
    plot_comparison({"无修形": baseline, "齿向鼓形": crowned, "组合修形": combined}, output / "comparisons")

    misaligned = calculate("05_misalignment", "两端法向失准 12 μm", misalignment=12e-6, make_plots=True)
    compensated = calculate("06_misalignment_with_crown", "失准下的鼓形方案",
                            modification=crown, misalignment=12e-6, make_plots=True)
    plot_comparison({"对中无修形": baseline, "失准无修形": misaligned, "失准鼓形": compensated},
                    output / "misalignment_comparison", title="相同载荷下的失准与修形响应")

    loads = {}
    for torque in (20., 40., 80., 120.):
        name = f"07_load_{torque:g}Nm"
        loads[f"{torque:g} N·m"] = calculate(name, f"组合修形：扭矩 {torque:g} N·m",
                                             torque=torque, modification=combined_surface)
    plot_comparison(loads, output / "load_comparison", title="同一修形方案的载荷适应性")
    calculate("08_helical", "20° 斜齿轮副", gear_pair=standard_pair(helix_degrees=20), make_plots=True)
    calculate("09_internal", "28/84 齿内啮合", gear_pair=standard_pair(internal=True), make_plots=True)

    # 这些点是人为定义的测量残差样式，用于展示接口，不是仪器实测数据。
    measured = MeasuredSurface([0., .5, 1.], [0., .5, 1.],
                               np.array([[0., .6, 0.], [.8, -.3, .8], [0., .6, 0.]]) * 1e-6)
    measured_combination = CompositeModification((combined_surface, measured))
    calculate("12_measured_residual", "参数修形叠加合成测量残差",
              modification=measured_combination, make_plots=True)
    plot_surface(measured_combination, output / "12_measured_residual", title="参数修形与合成测量残差组合")

    convergence_rows = []
    levels = (5, 9, 17) if quick else (11, 21, 41)
    for count in levels:
        result = calculate(f"10_grid_{count}", f"条带数 {count} 的离散检查",
                           modification=combined_surface, custom_strips=count)
        convergence_rows.append({"strips": count, **analysis_metrics(result)})
    grid_metrics = ("te_peak_to_peak_um", "maximum_pressure_mpa", "mean_stiffness_n_per_um")
    grid_changes = {key: abs(convergence_rows[-1][key] - convergence_rows[-2][key])
                        / max(abs(convergence_rows[-1][key]), 1e-15) for key in grid_metrics}
    grid_passed = all(value <= .02 for value in grid_changes.values())
    write_json(output / "grid_convergence.json", {
        "status": "末两级指标变化达标" if grid_passed else "需要进一步细化条带",
        "passed": grid_passed, "relative_tolerance": .02,
        "fixed_phases": phase_count, "levels": convergence_rows,
        "last_two_relative_changes": grid_changes,
        "interpretation": "只评价所列工况与指标的空间离散敏感性，不代表物理模型误差"})
    (output / "grid_convergence").mkdir(exist_ok=True)
    with (output / "grid_convergence" / "grid_convergence.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(convergence_rows[0]))
        writer.writeheader()
        writer.writerows(convergence_rows)
    _plot_grid_convergence(convergence_rows, output / "grid_convergence")

    candidates = {"无修形": SurfaceModification(), "鼓形5": SurfaceModification(lead_crown_m=5e-6),
                  "组合修形": combined_surface}
    operating_points = [OperatingPoint("轻载对中", 40., 0.), OperatingPoint("额定失准", 80., 12e-6),
                        OperatingPoint("过载反向失准", 120., -8e-6)]
    study = robust_design_search(pair, candidates, operating_points, **kwargs)
    robust_output = output / "11_robust_design"
    robust_output.mkdir(exist_ok=True)
    write_json(robust_output / "design_study.json", study.to_dict())
    write_json(robust_output / "study_input.json", {"pair": pair, "candidates": candidates,
                                                  "operating_points": operating_points, **kwargs})
    result_index = []
    for design_index, (design_name, design_results) in enumerate(study.results.items(), start=1):
        for point_index, (point_name, result) in enumerate(design_results.items(), start=1):
            relative = f"design_{design_index:02d}/point_{point_index:02d}"
            export_result(result, robust_output / relative)
            result_index.append({"design": design_name, "operating_point": point_name,
                                 "result_directory": relative})
    write_json(robust_output / "result_index.json", result_index)
    plot_design_study(study, robust_output)
    manifest["pareto_designs"] = list(study.pareto_names)
    manifest["elapsed_seconds"] = time.perf_counter() - started
    manifest["figures"] = sorted(str(path.relative_to(output))
                                 for extension in ("*.png", "*.svg")
                                 for path in output.rglob(extension))
    manifest["status"] = "全部案例运行完成"
    manifest["model_scope"] = "固定齿根变截面梁与齿宽剪切耦合降阶近似；局部弹性线接触"
    write_json(output / "case_manifest.json", manifest)
    with (output / "case_metrics.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest["cases"][0]))
        writer.writeheader()
        writer.writerows(manifest["cases"])
    print(json.dumps({"status": manifest["status"], "gear_cases": len(manifest["cases"]),
                      "analytic_status": manifest["analytic_status"],
                      "pareto_designs": manifest["pareto_designs"],
                      "elapsed_seconds": round(manifest["elapsed_seconds"], 3)}, ensure_ascii=False))
    return manifest


def _plot_grid_convergence(rows, output):
    output.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), constrained_layout=True)
    for ax, key, label in zip(axes,
                              ("te_peak_to_peak_um", "maximum_pressure_mpa", "mean_stiffness_n_per_um"),
                              ("TE 峰峰值 / μm", "峰值压力 / MPa", "平均切线刚度 / (N/μm)")):
        ax.plot([r["strips"] for r in rows], [r[key] for r in rows],
                marker="o", color="#286B83", linewidth=1.7)
        ax.set(xlabel="每接触线条带数", ylabel=label)
        ax.grid(alpha=.18)
    fig.suptitle("空间离散敏感性：保持几何、载荷与相位数不变", fontsize=12)
    _save(fig, output, "grid_convergence")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="运行圆柱齿轮LTCA的合成案例与图表")
    parser.add_argument("--output", type=Path, required=True, help="案例结果目录")
    parser.add_argument("--quick", action="store_true", help="使用小网格检查安装与调用，不作为收敛结论")
    parser.add_argument("--phases", type=int, help="覆盖默认相位数")
    parser.add_argument("--strips", type=int, help="覆盖默认每条接触线的条带数")
    args = parser.parse_args()
    run_examples(args.output, quick=args.quick, phases=args.phases, strips=args.strips)


if __name__ == "__main__":
    main()
