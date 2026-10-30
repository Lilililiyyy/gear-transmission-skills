# Gear Transmission Skills

[简体中文](README.md)

Agent Skills for gear transmission design, connecting natural-language requests to independent computational kernels for parameter studies, design comparisons, and reproducible calculations.

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

Existing evidence includes 14 computational workflow checks, 9 installation checks, a benchmark against the original project, and single- and multi-condition reruns. See the [computational acceptance record](evaluations/acceptance.json) and [installation record](evaluations/installation.json). Numerical reproducibility, discretization convergence, and physical validation are recorded separately; consult the independent project for model accuracy limits.

## Documentation and development

- [Installation guide](docs/installation_EN.md): environment setup, agent registration, natural-language installation, and reproducibility checks.
- [Development guide](docs/development_EN.md): directory responsibilities, kernel snapshots, and conventions for adding skills.
- [Input contract](skills/cylindrical-gear-analysis/references/input-contract.md) · [result definitions](skills/cylindrical-gear-analysis/references/results.md), currently in Chinese.

`skills/` contains skills, `packages/` contains kernel snapshots, and `src/gear_transmission/` provides invocation and reproducibility workflows. Computational source and original interfaces are maintained independently.

## Changelog

- **2026-10-30 · v0.1.0**
  - First runnable skill, `cylindrical-gear-analysis`: preserved cylindrical gear LTCA kernel, single-condition analysis, modification comparisons, reports, plots, and independent reruns.
  - Cross-platform Codex installation: `--pull` updates and installs; `--check` verifies complete resources and the runtime, while protecting existing same-named skills and local changes.
  - Bilingual READMEs, installation guides, and development notes; detailed computational models remain documented in the independent project.
