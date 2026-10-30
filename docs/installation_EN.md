# Installation and integration

[简体中文](installation.md)

Installation consists of a computational environment and agent-side skill registration. The instructions assume a complete local checkout, retaining the skill files, resources, and computational packages.

## Requirements

| Component | Requirement |
|---|---|
| Python | Python 3.10 / 3.11 recommended, in an isolated virtual environment |
| Agent | Access to the local project and permission to execute Python commands |
| Node.js | Required only for the optional Skills CLI; the referenced version requires 22.20+ |

## Agent-assisted installation

See the [README](../README_EN.md#install-through-an-agent) for the Codex remote-repository prompt. The script sets up both the computational environment and the complete skill. Other agents can use the general request below.

Send the following request in a session rooted at this project. If using another working directory, include the absolute project path:

```text
Install this project's cylindrical-gear-analysis for the current agent.
Follow this guide to configure a virtual environment, install dependencies, and register
the complete skill directory; reuse an existing working environment.
Run check, run --plots, and replay for the bundled single-condition example in new output directories.
Return the skill registration status, absolute interpreter path, report location, and replay.json.passed.
```

Registration targets and installation conventions are listed below. When an interface import is required, the agent should provide the local package path and remaining actions. Agent-side activation that has not been completed should be reported separately.

## Repository installation scripts

Clone the repository to a stable directory and choose the appropriate entry point. The defaults are `~/.agents/skills/` for skills and the checkout's `.venv` for the runtime. Retain the complete checkout rather than downloading only the skill text.

```bash
git clone https://github.com/Lilililiyyy/gear-transmission-skills.git
cd gear-transmission-skills
bash scripts/update-codex-skills.sh --pull
bash scripts/update-codex-skills.sh --check
```

```powershell
git clone https://github.com/Lilililiyyy/gear-transmission-skills.git
cd gear-transmission-skills
.\scripts\update-codex-skills.ps1 --pull
.\scripts\update-codex-skills.ps1 --check
```

Virtual-environment activation is unnecessary. The script reports the absolute interpreter path and JSON status; a missing compatible Python or dependency installation failure produces a nonzero exit code. If Windows restricts `.ps1` execution, use `python -B scripts/update_codex_skills.py --pull` directly and replace the option with `--check` for verification.

| Option | Behavior |
|---|---|
| No options | Install the current checkout without pulling; also supports a complete ZIP download |
| `--pull` | Run `git pull --ff-only` on a clean checkout, then install or update |
| `--check` | Read-only verification of complete skill resources, the receipt, actual import paths, dependencies, and kernel SHA256; no network access, installation, or solving |
| `--dest PATH` | Custom skill directory; also configurable through `CODEX_SKILLS_DIR` |
| `--runtime PATH` | Virtual-environment directory; pass the same option for later checks |

Set `GEAR_SKILLS_PYTHON` to select the Python interpreter used to create the environment. Only registered skills from this library are updated. The script stops and preserves existing same-named directories with different content or local edits. Previously managed, unmodified resources are synchronized during updates; unrelated skills remain untouched.

Claude Code can use the same installation logic with its user-level target:

```bash
bash scripts/update-codex-skills.sh --dest "$HOME/.claude/skills" --pull
bash scripts/update-codex-skills.sh --dest "$HOME/.claude/skills" --check
```

On Windows, add `--dest "$HOME/.claude/skills"` to the `.ps1` entry point. The destination's `.gear-transmission-skills-install.json` receipt records the checkout and interpreter paths for agents. Rerun installation after moving the checkout or environment. A successful `--check` confirms files and runtime consistency; confirm host discovery and numerical reproducibility as described below.

## Computational environment

Skip this section if the script already configured the environment. For manual installation, run these commands from the project root.

### Windows PowerShell

```powershell
$env:PYTHONUTF8="1"
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -c constraints.txt -e ".[plots]"
```

### macOS / Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -c constraints.txt -e ".[plots]"
```

For computation only, use `-e .` and omit subsequent `--plots` options. Use the same virtual-environment interpreter for later calls. `constraints.txt` records the principal dependency versions from the existing acceptance environment.

## Agent registration

| Agent | User-level target | Project-level target | Invocation |
|---|---|---|---|
| Codex | `~/.agents/skills/cylindrical-gear-analysis/` | `.agents/skills/cylindrical-gear-analysis/` | `$cylindrical-gear-analysis` |
| Claude Code | `~/.claude/skills/cylindrical-gear-analysis/` | `.claude/skills/cylindrical-gear-analysis/` | `/cylindrical-gear-analysis` |
| WorkBuddy | Import a local package through skill management | Use the entry point supported by the client | Mention the skill in natural language |

### Directory installation

Copy or link the complete `skills/cylindrical-gear-analysis/` directory to the selected target. Retain `SKILL.md`, `agents/`, `scripts/`, `references/`, and `assets/`. Check the source of an existing same-named target before updating it.

The Windows example targets Codex. For Claude Code, replace `.agents\skills` with `.claude\skills`:

```powershell
$gearSkillSource=(Resolve-Path '.\skills\cylindrical-gear-analysis').Path
$gearSkillHome=Join-Path $env:USERPROFILE '.agents\skills'
$gearSkillTarget=Join-Path $gearSkillHome 'cylindrical-gear-analysis'
if (Test-Path -LiteralPath $gearSkillTarget) { throw '技能目录已存在 / Skill directory already exists.' }
New-Item -ItemType Directory -Path $gearSkillHome -Force | Out-Null
Copy-Item -LiteralPath $gearSkillSource -Destination $gearSkillTarget -Recurse
```

The macOS / Linux example also targets Codex. For Claude Code, replace `.agents/skills` with `.claude/skills`:

```bash
mkdir -p "$HOME/.agents/skills"
cp -R ./skills/cylindrical-gear-analysis "$HOME/.agents/skills/"
```

Both agents support linked skill directories. Copies must be synchronized after source updates; links require the original project to remain at a stable location. Reopen the session or restart the agent if the new skill has not been discovered.

### Optional Skills CLI

[Vercel Skills CLI](https://github.com/vercel-labs/skills) can replace the directory-copy step. Run one command for the selected agent:

```bash
npx skills add ./ --skill cylindrical-gear-analysis --agent codex --global
npx skills add ./ --skill cylindrical-gear-analysis --agent claude-code --global
npx skills add ./ --skill cylindrical-gear-analysis --agent cursor --global
```

`--global` selects user-level installation; omit it for project-level installation. The installer handles skill files only. Configure the Python computational environment separately.

### WorkBuddy

Use **Add skill → Upload local skill package** in skill management to import the complete skill, then enable it. The agent can prepare a package for the installed client version. Provide the complete project location and virtual-environment interpreter after import.

Package formats and directory conventions depend on the current WorkBuddy client. This library has not been tested in that client; CodeBuddy Code directory conventions are not treated as a universal WorkBuddy standard.

## Verification and invocation

On Windows, use the same interpreter throughout:

```powershell
.\.venv\Scripts\python.exe -m gear_transmission check skills/cylindrical-gear-analysis/assets/single_case.json
.\.venv\Scripts\python.exe -m gear_transmission run skills/cylindrical-gear-analysis/assets/single_case.json --output output/example --plots
.\.venv\Scripts\python.exe -m gear_transmission replay output/example --output output/example_replay
```

On macOS / Linux, use `.venv/bin/python` instead. Output directories must be new.

| Check | Expected result |
|---|---|
| Agent registration | The skill is visible in the current agent and enabled where required |
| Actual calculation | `run.json` contains `status=completed` and `analysis_completed=true`; the report and raw arrays are readable |
| Independent rerun | `replay.json` contains `passed=true` |

You can also request the workflow directly from the agent:

```text
Use cylindrical-gear-analysis with this project's Python virtual environment to run the bundled
single-condition example, save the report and plots, and perform an independent rerun.
```

`check` does not execute the numerical solver. Confirm agent discovery, calculation completion, and numerical reproducibility separately. A working computational program alone does not establish that the agent has loaded the skill.

## Installation conventions

For agent-assisted installation and subsequent maintenance:

- Locate `gear-skills.project.json` and `pyproject.toml` from the working directory or the supplied path. Register only the requested agent.
- Reuse or create an isolated project environment and keep installation, execution, and recomputation on the same Python interpreter. Report missing interpreters as incomplete setup.
- Retain all skill resources and the computational kernel. Do not alter formulas or replace the registered snapshot during installation.
- Execute the calculation and rerun, returning the installation location, interpreter, report, and status. Distinguish deployed files from confirmed agent-side loading.
- Cloud sessions require code and dependencies in their actual execution environment. Local paths and registration do not migrate automatically.

## References and verification scope

Directory and invocation conventions were checked on 2026-10-30.

- [Claude Code Skills](https://code.claude.com/docs/en/skills) · [OpenAI Skills](https://learn.chatgpt.com/docs/build-skills).
- [Official WorkBuddy skill documentation](https://cloud.tencent.com/document/product/1831/134432).
- [Vercel Skills CLI](https://github.com/vercel-labs/skills) · [nature-skills](https://github.com/Yuan1z0825/nature-skills).

Existing evidence covers local Python invocation, numerical reruns, and scripted installation checks through Windows PowerShell and Git Bash; see the [installation record](../evaluations/installation.json). Native Linux / macOS environments and activation in individual clients have not all been tested; the current client state remains authoritative.
