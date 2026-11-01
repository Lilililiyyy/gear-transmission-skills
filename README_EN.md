# Gear Transmission Skills

[简体中文](README.md) · **English**

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/) [![License: AGPL-3.0-only](https://img.shields.io/badge/License-AGPL--3.0--only-blue.svg)](LICENSE)

Agent Skills for gear transmission design, connecting natural-language requests to independent computational kernels for parameter studies, design comparisons, and reproducible calculations.

The library and integrated LTCA kernel are both **0.2.0**, with the kernel pinned to independent-project commit `c63030d`.

## Contents

- [Available skills](#available-skills)
- [Installation](#installation)
- [Usage](#usage)
- [Outputs and validation](#outputs-and-validation)
- [Documentation and development](#documentation-and-development)
- [Changelog](#changelog)

## Available skills

| Skill | Capabilities | Status | Independent computational project |
|---|---|---|---|
| [cylindrical-gear-analysis](skills/cylindrical-gear-analysis/SKILL.md) | Cylindrical gear LTCA, modification and misalignment analysis, TE / normal stiffness / Hertz pressure, multi-condition design comparisons | Available | [Cylindrical Gear LTCA](https://github.com/Lilililiyyy/Cylindrical-Gear-LTCA) |

The independent project maintains the computational models, complete API, and examples. This library provides agent workflows, input adaptation, reports, and reproducibility checks.

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

| Agent | Registration | Explicit invocation |
|---|---|---|
| Codex | User-level `~/.agents/skills/` | `$cylindrical-gear-analysis` |
| Claude Code | User-level `~/.claude/skills/` | `/cylindrical-gear-analysis` |
| WorkBuddy | Import a local skill package through skill management and enable it | Mention the skill in natural language |
| Cursor | Select the `cursor` target with `npx skills` | Mention the skill or use the skill selector |

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

## Outputs and validation

- `report.md`: principal metrics, design trade-offs, and model scope.
- `run.json` and `cases/`: inputs, code identity, raw arrays, and solver diagnostics.
- `replay.json`: input, summary, and complete array comparisons after recomputation.

There are currently 16 computational workflow checks and 9 installation checks. See the [kernel 0.2.0 synchronization record](evaluations/kernel-sync-0.2.0.json) for the current source comparisons and reruns. The original [computational acceptance record](evaluations/acceptance.json) and [installation record](evaluations/installation.json) retain evidence for their recorded older versions. Numerical reproducibility, discretization convergence, and physical validation are recorded separately; consult the independent project for model accuracy limits.

Installation, actual calculations, and independent reruns across Windows / Linux × Python 3.10 / 3.11 are covered by [continuous integration](https://github.com/Lilililiyyy/gear-transmission-skills/actions/workflows/ci.yml). The initial release results remain in the [CI record](evaluations/ci.json).

## Documentation and development


- [Installation guide](docs/installation_EN.md): environment setup, agent registration, natural-language installation, and reproducibility checks.
- [Development guide](docs/development_EN.md): directory responsibilities, kernel snapshots, and conventions for adding skills.
- [Input contract](skills/cylindrical-gear-analysis/references/input-contract.md) · [result definitions](skills/cylindrical-gear-analysis/references/results.md), currently in Chinese.

`skills/` contains skills, `packages/` contains kernel snapshots, and `src/gear_transmission/` provides invocation and reproducibility workflows. Computational source and original interfaces are maintained independently.

## License and commercial licensing

This project is licensed under [AGPL-3.0-only](LICENSE). For proprietary or closed-source integration, [contact the author](https://github.com/Lilililiyyy) to discuss a commercial license.

## Changelog

- **2026-11-01 · v0.2.0**
  - Integrate cylindrical gear kernel `0.2.0`, refreshing 21 files, source hashes, and function mappings from the clean working tree at commit `c63030d`.
  - Include ISO curvature handling for unloaded contact points and extended-contact CSV flags, with regression coverage for both.
  - Align library, skill, and kernel versions; correct the skill entry's license identifier and difference diagnostics for matching infinities, then repeat installation, source-consistency, and independent-rerun checks.

- **2026-10-30 · v0.1.0**
  - First runnable skill, `cylindrical-gear-analysis`: preserved cylindrical gear LTCA kernel, single-condition analysis, modification comparisons, reports, plots, and independent reruns.
  - Cross-platform Codex installation: `--pull` updates and installs; `--check` verifies complete resources and the runtime, while protecting existing same-named skills and local changes.
  - Bilingual READMEs, installation guides, and development notes; detailed computational models remain documented in the independent project.
  - First GitHub publication, passing four CI configurations across Windows / Linux and Python 3.10 / 3.11.
