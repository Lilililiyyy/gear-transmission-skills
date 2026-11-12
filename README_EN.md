# Gear Transmission Skills

[简体中文](README.md) · **English**

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/) [![License: AGPL-3.0-only](https://img.shields.io/badge/License-AGPL--3.0--only-blue.svg)](LICENSE) [![GitHub Release](https://img.shields.io/github/v/release/Lilililiyyy/gear-transmission-skills)](https://github.com/Lilililiyyy/gear-transmission-skills/releases)

Agent Skills for gear transmission design, connecting natural-language requests to independent computational kernels for parameter studies, design comparisons, and reproducible calculations.

## Contents

- [Available skills](#available-skills)
- [Installation](#installation)
- [Usage](#usage)
- [Outputs and validation](#outputs-and-validation)
- [Documentation and development](#documentation-and-development)
- [Changelog](#changelog)

## Available skills

<!-- skills:start -->

| Skill                                                                  | Capabilities                                                                                                                     | Version | Independent computational project                                             |
| ---------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ------- | ----------------------------------------------------------------------------- |
| [cylindrical-gear-analysis](skills/cylindrical-gear-analysis/SKILL.md) | Cylindrical LTCA, modification / misalignment comparisons, contact plots, explicit-factor strength and benchmark cases | 0.3.0 | [Cylindrical Gear LTCA](https://github.com/Lilililiyyy/Cylindrical-Gear-LTCA) |
| [bevel-gear-design](skills/bevel-gear-design/SKILL.md) | Zero-offset straight / spiral bevel geometry, forces, nominal strength, staged CAD / TCA / nodal FE-LTCA, and independent reruns | 0.10.0 | [Bevel Gear Design](https://github.com/Lilililiyyy/Bevel-Gear-Design) |
<!-- skills:end -->

## Installation

### Install through an agent

Send this to Codex with local command execution available:

```text
Install Codex skills from this repository:
https://github.com/Lilililiyyy/gear-transmission-skills.git

Clone it to a stable local directory, then run scripts/update-codex-skills.sh --pull.
On Windows PowerShell, use the equivalent .\scripts\update-codex-skills.ps1 --pull.
After installation, run --check through the same entry point.
Retain the complete skill directories under skills/; do not copy only SKILL.md or delete the checkout.
Finally, use the Python interpreter reported by the installer to run the bundled single-condition
example with plots and an independent rerun. Return the report path and result.
```

The script configures the project's `.venv`, installs computational and plotting dependencies, and syncs the complete skill to `~/.agents/skills/`. Python 3.10 / 3.11 and Git are required. `--pull` accepts only fast-forward updates on a clean checkout; `--check` verifies resources, import paths, and kernel hashes without modifying files or solving a case. Retain the checkout for calculations and updates.

### Command-line installation

macOS / Linux or Git Bash:

```bash
git clone https://github.com/Lilililiyyy/gear-transmission-skills.git
cd gear-transmission-skills
bash scripts/update-codex-skills.sh --pull
bash scripts/update-codex-skills.sh --check
```

Windows PowerShell:

```powershell
git clone https://github.com/Lilililiyyy/gear-transmission-skills.git
cd gear-transmission-skills
.\scripts\update-codex-skills.ps1 --pull
.\scripts\update-codex-skills.ps1 --check
```

For later updates, rerun `--pull` followed by `--check` in the same checkout. See the [installation guide](docs/installation_EN.md) for custom directories, existing environments, and other agents.

### Manual installation

From the project root, in a Python 3.10 / 3.11 virtual environment, run:

```bash
python -m pip install -c constraints.txt -e ".[plots]"
```

Then register the complete `skills/cylindrical-gear-analysis/` directory:

| Agent       | Registration                                                        | Explicit invocation                         |
| ----------- | ------------------------------------------------------------------- | ------------------------------------------- |
| Codex       | User-level `~/.agents/skills/`                                      | `$cylindrical-gear-analysis`                |
| Claude Code | User-level `~/.claude/skills/`                                      | `/cylindrical-gear-analysis`                |
| WorkBuddy   | Import a local skill package through skill management and enable it | Mention the skill in natural language       |
| Cursor      | Select the `cursor` target with `npx skills`                        | Mention the skill or use the skill selector |

See the [installation guide](docs/installation_EN.md) for platform commands, the optional Skills CLI, and WorkBuddy integration. Both skill registration and Python environment setup are required.

## Usage

Describe the operating conditions and design objectives directly:

```text
Use cylindrical-gear-analysis for an external spur gear pair: 28/56 teeth, normal module 2.5 mm,
face width 25 mm, pressure angle 20 degrees, E=210 GPa and nu=0.3 for both gears,
and misalignment 12 micrometres. Compare pinion lead crowning of 0, 8, and 12 micrometres
at 20, 80, and 120 N·m using TE and pressure. Use the coupled model with local indentation,
32 phases, and 21 strips; list the remaining defaults.
Save the report, plots, and raw data, then perform an independent rerun.
```

You can also run an example directly in the installed environment:

```bash
python -m gear_transmission run skills/cylindrical-gear-analysis/assets/single_case.json --output output/example --plots
python -m gear_transmission replay output/example --output output/example_replay
```

Output directories must be new. Bundled inputs: [single condition](skills/cylindrical-gear-analysis/assets/single_case.json) · [multi-condition modification comparison](skills/cylindrical-gear-analysis/assets/crown_study.json).

Standard study runs export common-scale contact comparisons to `plots/point_NNN/07_contact_pass_comparison.png`. To redraw an existing run:

```bash
python -m gear_transmission plot output/example --output output/example_replot
```

This validates the saved files and reads their original arrays. It writes new figures and a `plot.json` recording that no solver was executed. The display uses the actual meshing interval; non-meshing root regions and genuine uncovered samples remain distinct. Edge-loaded footprints need not form closed ellipses.

To reproduce the independent project's updated **18/77 helical case and contact figures**:

```bash
python -m gear_transmission case --output output/helical_18_77 --plots
python -m gear_transmission replay output/helical_18_77 --output output/helical_18_77_replay
```

Or ask: “Use cylindrical-gear-analysis to reproduce the updated complete 18/77 case, including the tooth-pass pressure map, three-cycle response, instantaneous pressure band, explicit-factor strength comparison, and an independent rerun.”

Ordinary `run --plots` also produces the updated geometry, cycle, instantaneous-band, and tooth-pass figures from the user's actual inputs. Bundled MASTA data is used only for case comparisons. Strength calculations require explicit rating factors and material stress limits. See the [updated workflow reference](skills/cylindrical-gear-analysis/references/new-workflow.md), currently in Chinese, for plot definitions and scope.

## Outputs and validation

- `report.md`: principal metrics, design trade-offs, and model scope.
- `run.json` and `cases/`: inputs, code identity, raw arrays, and solver diagnostics.
- `replay.json`: input, summary, and complete array comparisons after recomputation.

The complete `case` workflow also exports `workflow_summary.json`, `strength.json`, `ltca/`, and `sensitivity/`. Actual solve directories are listed in `run.json.cases`.

## Documentation and development

- [Installation guide](docs/installation_EN.md): environment setup, agent registration, natural-language installation, and reproducibility checks.
- [Development guide](docs/development_EN.md): directory responsibilities, kernel snapshots, and conventions for adding skills.
- [Input contract](skills/cylindrical-gear-analysis/references/input-contract.md) · [result definitions](skills/cylindrical-gear-analysis/references/results.md), currently in Chinese.

`skills/` contains skills, `packages/` contains kernel snapshots, and `src/gear_transmission/` provides invocation and reproducibility workflows. Computational source and original interfaces are maintained independently.

## License and commercial licensing

This project is licensed under [AGPL-3.0-only](LICENSE). For proprietary or closed-source integration, [contact the author](https://github.com/Lilililiyyy) to discuss a commercial license.

## Changelog

- **2026-11-12 · v0.6.0**

  - 接入圆柱LTCA 0.3.0：保留固定版本响应回归、多工况收敛矩阵及原组件对标的独立状态。
  - 比较请求新增 `study.pressure_metric=center|endpoint`，默认中心压力；所选指标贯穿约束、Pareto、报告、方案比较图、独立复算和只读重绘。
  - 端点指标缺少实际数据时明确拒绝，原中心／端点数组和单位保持不变；快照按来源提交、逐文件哈希和函数行段登记。
  - 快照直接读取固定 Git 提交字节，避免检出换行影响来源哈希。[本地验收](evaluations/ltca-first-round-0.6.0.json)：103项内核测试、50项适配测试及10案响应回归通过；三类工况仍存在离散敏感性，保留未收敛状态。

- **2026-11-12 · v0.5.1**

  - 修复发布校验与README展示格式不一致的问题：保留Release徽章和四列技能表，不再要求已删除的库版本说明区块及状态列。
  - 内核版本继续从快照同步；固定提交与文件哈希保留在快照清单，发布说明继续列出内核提交。
  - 增加四列表格、只读检查及重复同步回归检查；计算内核和历史发布记录不变。

- **2026-11-02 · v0.5.0**

  - Release all skills under one library version; remove per-skill versions. Maintain the version only in the root `pyproject.toml` and read it from distribution metadata at runtime.
  - Generate the integrated kernel version and pinned commit in both skill tables from the existing snapshot manifests. Kernels retain their independent releases and validation.
  - Add metadata consistency checks and release tooling. Successful main-branch CI creates the version tag and bilingual GitHub Release while preserving historical tags and acceptance records.

- **2026-11-02 · v0.4.1**

  - Correct the contact display interval and retain finite centre-cell display support without changing raw pressure or PPTE.
  - Share the upstream contact display code across ordinary plots, benchmark cases and common-scale design comparisons.
  - Add hash-verified, read-only `plot` replotting, separate from solving and numerical replay. See the [0.4.1 acceptance record](evaluations/contact-display-0.4.1.json).

- **2026-11-02 · v0.3.0**

  - Include complete case orchestration, plotting helpers, and public benchmark inputs; record and verify all 47 upstream runtime files.
  - Add `case`, directly invoking the preserved upstream workflow for explicit-factor strength, instantaneous comparisons, tooth-pass pressure maps, and sensitivity cases.
  - Apply updated plotting to ordinary inputs, keeping centre pressure, endpoint pressure, instantaneous bands, and cycle envelopes distinct, with three-cycle display and SI plotting data.
  - Add complete-resource synchronization checks, array and strength reruns, source figure comparisons, and wheel distribution checks.

- **2026-11-01 · v0.2.0**

  - Integrate cylindrical gear kernel `0.2.0`, refreshing 21 files, source hashes, and function mappings from the clean working tree at commit `c63030d`.
  - Include ISO curvature handling for unloaded contact points and extended-contact CSV flags, with regression coverage for both.
  - Align library, skill, and kernel versions; correct the skill entry's license identifier and difference diagnostics for matching infinities, then repeat installation, source-consistency, and independent-rerun checks.

- **2026-10-30 · v0.1.0**

  - First runnable skill, `cylindrical-gear-analysis`: preserved cylindrical gear LTCA kernel, single-condition analysis, modification comparisons, reports, plots, and independent reruns.
  - Cross-platform Codex installation: `--pull` updates and installs; `--check` verifies complete resources and the runtime, while protecting existing same-named skills and local changes.
  - Bilingual READMEs, installation guides, and development notes; detailed computational models remain documented in the independent project.
  - First GitHub publication, passing four CI configurations across Windows / Linux and Python 3.10 / 3.11.
