"""保持空间离散不变，独立检查啮合相位采样的敏感性。"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from run_examples import standard_pair, write_json
from gear_contact import Material, SurfaceModification, analyze_pair
from gear_contact.io import export_result
from gear_contact.plotting import analysis_metrics, _pyplot, _save


def check_phase_resolution(output_dir, *, phase_counts=(16, 32, 64), strips=21,
                           relative_tolerance=.02):
    """比较末两级的指标相对变化，不假定离散误差单调或已知收敛阶。"""
    if len(phase_counts) < 2 or any(n < 4 for n in phase_counts) or any(
            b <= a for a, b in zip(phase_counts, phase_counts[1:])):
        raise ValueError("相位数应为至少两级、严格递增且不小于4的整数")
    if not 0 < relative_tolerance < 1:
        raise ValueError("相位指标相对变化阈值应在(0,1)内")
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    pair = standard_pair()
    surface = SurfaceModification(lead_crown_m=8e-6, tip_relief_m=6e-6,
                                   tip_start=.72, root_relief_m=3e-6, root_end=.18)
    material = Material(210e9, .3)
    rows = []
    responses = []
    for phases in phase_counts:
        result = analyze_pair(pair, torque_n_m=80, material_pinion=material, material_wheel=material,
                              modification_pinion=surface, phases=phases, strips=strips)
        export_result(result, output / f"phases_{phases}")
        rows.append({"phases": phases, **analysis_metrics(result),
                     "force_residual_n": result.summary["max_force_residual_n"],
                     "gap_residual_m": result.summary["max_gap_residual_m"]})
        responses.append(result)
    metrics = ("te_peak_to_peak_um", "maximum_pressure_mpa", "mean_stiffness_n_per_um")
    changes = {key: abs(rows[-1][key] - rows[-2][key]) / max(abs(rows[-1][key]), 1e-15)
               for key in metrics}
    passed = all(change <= relative_tolerance for change in changes.values())
    report = {"status": "末两级指标变化达标" if passed else "需要进一步细化相位",
              "passed": passed, "relative_tolerance": relative_tolerance,
              "fixed_strips": strips, "phase_counts": phase_counts,
              "input": {"pair": pair, "torque_n_m": 80, "material": material,
                        "modification_pinion": surface},
              "results": rows, "last_two_relative_changes": changes,
              "interpretation": "只评价所列工况与指标的相位敏感性，不代表全参数域或物理模型误差"}
    write_json(output / "phase_resolution.json", report)
    with (output / "phase_resolution.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), constrained_layout=True)
    for result in responses:
        axes[0, 0].plot(result.phase, result.normal_approach_m * 1e6,
                        label=f"{len(result.phase)} 相位", linewidth=1.5)
    axes[0, 0].set(xlabel="啮合周期相位", ylabel="法向接近量 / μm", title="a  周期响应采样")
    axes[0, 0].legend()
    for ax, key, label, title in zip(axes.ravel()[1:], metrics,
                                     ("TE 峰峰值 / μm", "峰值压力 / MPa", "平均切线刚度 / (N/μm)"),
                                     ("b  TE 波动", "c  峰值压力", "d  平均切线刚度")):
        ax.plot(phase_counts, [row[key] for row in rows], "o-", color="#286B83")
        ax.set(xlabel="每周期相位数", ylabel=label, title=title)
    for ax in axes.ravel():
        ax.grid(alpha=.18)
    fig.suptitle(f"相位离散敏感性：固定每接触线 {strips} 条带", fontsize=13)
    _save(fig, output, "phase_resolution")
    plt.close(fig)
    return report


def main():
    parser = argparse.ArgumentParser(description="检查组合修形案例的相位采样敏感性")
    parser.add_argument("--output", type=Path, required=True, help="检查结果目录")
    parser.add_argument("--phases", type=int, nargs="+", default=[16, 32, 64], help="递增的相位数")
    parser.add_argument("--strips", type=int, default=21, help="固定的每接触线条带数")
    parser.add_argument("--relative-tolerance", type=float, default=.02, help="末两级指标相对变化阈值")
    args = parser.parse_args()
    report = check_phase_resolution(args.output, phase_counts=tuple(args.phases), strips=args.strips,
                                    relative_tolerance=args.relative_tolerance)
    print(report["status"])
    print(report["last_two_relative_changes"])
    if not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
