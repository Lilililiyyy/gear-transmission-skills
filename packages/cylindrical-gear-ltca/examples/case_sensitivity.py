"""新写同一齿轮案例的附加失准/鼓形敏感性工具，计算使用现有LTCA内核。

保留案例原几何、材料、载荷、修形及给定边界，仅叠加明确的微米量级
扰动。本工具不计算新方案的标准强度，也不继承MASTA对标通过状态。
"""
from __future__ import annotations

from copy import deepcopy
import csv
from pathlib import Path
import sys
import time

import numpy as np

PROJECT = next(p for p in Path(__file__).resolve().parents
               if (p/"pyproject.toml").is_file() and (p/"gear_contact").is_dir())
for folder in (PROJECT, PROJECT/"examples"):
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))

from gear_contact import analyze_config, export_result
from gear_contact.io import write_json
from gear_contact.plotting import _pyplot, _save
from case_contact_plotting import repeat_cycle_for_display
from run_cylindrical_workflow import prepare_case


def _csv(path, columns, rows):
    """写出SI原值；显示单位转换只在绘图步骤完成。"""
    with Path(path).open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def _scenario_config(case, additional_misalignment_m, additional_pinion_lead_crown_m):
    """复用案例载荷记账，叠加扰动而不替换已有修形或给定边界。"""
    config, _ = prepare_case(deepcopy(case))
    if config.get("extended_contact_settings") is None:
        raise ValueError("此敏感性比较要求启用齿顶延伸，以分别保存并比较端点压力")
    config["misalignment_m"] = config.get("misalignment_m", 0.)+additional_misalignment_m
    if additional_pinion_lead_crown_m:
        config["modification_pinion"] = {
            "kind": "composite", "components": [deepcopy(config.get("modification_pinion", {})),
                {"kind": "parameters", "lead_crown_m": additional_pinion_lead_crown_m}]}
    return config


def run_case_sensitivity(case: dict, output: Path, *, plots=False):
    """计算原案例、附加1μm失准及再叠加0.5μm小轮齿向鼓形三个方案。

    output/sensitivity/<方案>/analysis_config.json可直接传给analyze_config；
    同目录保存完整LTCA数组。根目录06文件保存指标、单周期和三周期
    展示数据；三周期仅重复稳态周期，不表示连续瞬态求解。
    """
    specifications = (
        ("baseline", "原案例", 0., 0.),
        ("additional_misalignment", "附加1 μm失准", 1e-6, 0.),
        ("misalignment_and_crown", "失准 + 0.5 μm鼓形", 1e-6, .5e-6),
    )
    # 先构造全部配置，缺失端点模型时不创建看似完成的结果目录。
    configurations = [_scenario_config(case, misalignment, crown)
                      for _, _, misalignment, crown in specifications]
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    records, cycles, cycle_rows, display_rows = [], [], [], []
    for (scenario_id, label, misalignment, crown), config in zip(specifications, configurations):
        result = analyze_config(config)
        if (not result.summary["all_converged"] or result.summary["max_force_residual_n"] > 1e-7
                or result.summary["max_gap_residual_m"] > 3e-12):
            raise RuntimeError(f"{label}未满足接触平衡与协调要求，不作为有效敏感性结果")
        if any(rotation.endpoint_pressure is None for rotation in result.rotations):
            raise ValueError("敏感性方案缺少端点压力，不能改用中心峰压冒充相同定义")
        endpoint_peaks = np.asarray([np.max(rotation.endpoint_pressure["pressure_pa"])
                                     for rotation in result.rotations])
        record = {"scenario_id": scenario_id, "label": label,
            "additional_misalignment_m": misalignment,
            "additional_pinion_lead_crown_m": crown,
            "mean_normal_approach_m": float(np.mean(result.normal_approach_m)),
            "te_peak_to_peak_m": float(np.ptp(result.normal_approach_m)),
            "mean_mesh_stiffness_n_per_m": float(np.mean(result.mesh_stiffness_n_per_m)),
            "max_endpoint_pressure_pa": float(np.max(endpoint_peaks)),
            "normal_force_n": result.summary["normal_force_n"],
            "max_force_residual_n": result.summary["max_force_residual_n"],
            "max_gap_residual_m": result.summary["max_gap_residual_m"],
            "all_converged": True,
            "config_file": f"sensitivity/{scenario_id}/analysis_config.json",
            "ltca_directory": f"sensitivity/{scenario_id}/ltca"}
        records.append(record)
        directory = output/"sensitivity"/scenario_id
        directory.mkdir(parents=True, exist_ok=True)
        write_json(directory/"analysis_config.json", config)
        export_result(result, directory/"ltca")
        cycle = np.column_stack((np.arange(len(result.phase))/len(result.phase), result.phase,
            result.normal_approach_m, result.mesh_stiffness_n_per_m, endpoint_peaks))
        displayed = repeat_cycle_for_display(cycle, 3)
        cycles.append(displayed)
        cycle_rows.extend((scenario_id, *row) for row in cycle)
        display_rows.extend((scenario_id, *row) for row in displayed)
    metrics_file = "06_modification_sensitivity.csv"
    _csv(output/metrics_file, tuple(records[0]), ([row[key] for key in records[0]] for row in records))
    columns = ("scenario_id", "cycle_fraction", "physical_phase", "normal_approach_m",
               "mesh_stiffness_n_per_m", "endpoint_peak_pressure_pa")
    _csv(output/"06_sensitivity_cycles.csv", columns, cycle_rows)
    _csv(output/"06_sensitivity_cycle_display.csv", (columns[0], "display_cycle_position", *columns[2:]), display_rows)
    files = [metrics_file, "06_sensitivity_cycles.csv", "06_sensitivity_cycle_display.csv"]
    if plots:
        plt = _pyplot()
        fig, axes = plt.subplots(2, 1, figsize=(11.4, 7.2), sharex=True)
        colors = ("#176B91", "#D47732", "#7259A0")
        for record, displayed, color in zip(records, cycles, colors):
            axes[0].plot(displayed[:, 0], displayed[:, 2]*1e6, color=color, lw=1.6, label=record["label"])
            axes[1].plot(displayed[:, 0], displayed[:, 4]*1e-6, color=color, lw=1.6, label=record["label"])
        for axis in axes:
            axis.set_xlim(0, 3)
            axis.set_xticks(np.arange(0., 3.01, .5))
            axis.grid(True)
            for position in (1., 2.):
                axis.axvline(position, color="#A0ACB6", ls="--", lw=.8, zorder=0)
        axes[0].set(title="a  法向TE", ylabel="法向接近量 / μm")
        axes[1].set(title="b  端点峰值接触压力", ylabel="端点峰压 / MPa",
                    xlabel="相对起采位置 / 基节周期")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(.5, .99))
        fig.subplots_adjust(top=.90, bottom=.19, left=.10, right=.98, hspace=.27)
        fig.text(.5, .055, "18/77斜齿轮副：保持几何、材料与载荷，比较附加失准和齿向鼓形的影响。\n"
            "本项目计算结果，展示三个稳态啮合周期。",
            ha="center", fontsize=8.6, color="#526470")
        files.extend(Path(path).name for path in _save(fig, output, "06_modification_sensitivity"))
        plt.close(fig)
    summary = {"case_id": case.get("case_id"), "analysis": "同一案例的附加失准/齿向鼓形敏感性",
        "scenarios": records, "scenario_count": len(records), "display_cycles": 3,
        "statistics_from_single_period": True,
        "original_geometry_material_loading_modification_and_boundary_preserved": True,
        "strength_evaluated": False, "masta_comparison_performed": False,
        "physical_validation": False, "elapsed_seconds": time.perf_counter()-started,
        "files": files,
        "scope": "仅本项目LTCA敏感性分析；TE、刚度和端点峰压分别比较，不继承对标状态，不作最优性结论。"}
    write_json(output/"06_modification_sensitivity.json", summary)
    return summary
