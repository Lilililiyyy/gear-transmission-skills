"""重新计算 README 的接触、TE 和刚度示例，同时保留输入、数组和绘图数据。"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gear_contact.analysis import analyze_pair
from gear_contact.io import export_result
from gear_contact.plotting import (plot_comparison, plot_contact_comparison,
                                   plot_contact_envelope_comparison, plot_surface)
from gear_contact.surfaces import SurfaceModification
from run_examples import standard_pair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=_ROOT / "output/readme_examples")
    parser.add_argument("--assets", type=Path, default=_ROOT / "examples/results")
    parser.add_argument("--phases", type=int, default=64)
    parser.add_argument("--strips", type=int, default=41)
    args = parser.parse_args()
    pair = standard_pair()
    results = {}
    cases = (
        ("Aligned", "aligned", 0., SurfaceModification()),
        ("Misaligned", "misaligned", 12e-6, SurfaceModification()),
        ("Lead crowned", "lead_crowned", 12e-6, SurfaceModification(lead_crown_m=10e-6)),
    )
    for label, name, misalignment, modification in cases:
        print(f"计算 {label}：{args.phases} 相位／{args.strips} 条带", flush=True)
        result = analyze_pair(pair, torque_n_m=80., modification_pinion=modification,
                              misalignment_m=misalignment, phases=args.phases, strips=args.strips)
        if not result.summary["all_converged"]:
            raise RuntimeError(f"{label} 未收敛，停止生成示例图")
        export_result(result, args.output / name)
        results[label] = result
    plot_contact_envelope_comparison(results, args.output,
                                     title="Loaded contact envelope | 28/56 spur gears · 80 N·m")
    # 瞬时斑点单独选用对中鼓形案例，与周期包络及失准补偿案例明确区分。
    centered = analyze_pair(pair, torque_n_m=80., modification_pinion=SurfaceModification(lead_crown_m=10e-6),
                            phases=args.phases, strips=args.strips)
    if not centered.summary["all_converged"]:
        raise RuntimeError("对中鼓形案例未收敛")
    export_result(centered, args.output / "centered_crown")
    plot_contact_comparison({"Centered lead crown · 10 μm": centered}, args.output,
                            phase_index=args.phases * 3 // 8)
    plot_comparison(results, args.output, title="相同载荷下的失准与修形响应")
    plot_surface(SurfaceModification(lead_crown_m=8e-6, tip_relief_m=6e-6,
                                     tip_start=.72, root_relief_m=3e-6, root_end=.18), args.output)
    args.assets.mkdir(parents=True, exist_ok=True)
    for name in ("contact_pattern", "contact_envelope", "transmission_error", "mesh_stiffness", "design_comparison", "modification_surface"):
        shutil.copy2(args.output / f"{name}.png", args.assets / f"{name}.png")
    data = args.assets
    data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.output / "comparison_curves.csv", data / "example_responses.csv")
    print(f"示例图已更新：{args.assets.resolve()}", flush=True)


if __name__ == "__main__":
    main()
