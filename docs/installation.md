# 安装与接入

[English](installation_EN.md)

安装包含计算环境和宿主技能注册两部分。以下操作基于完整的本地项目；技能文件、资源和计算包均须保留。

## 环境要求

| 组件 | 要求 |
|---|---|
| Python | 推荐 3.10／3.11，使用独立虚拟环境 |
| 智能体 | 可读取本地项目并执行 Python 命令 |
| Node.js | 仅使用可选的 Skills CLI 时需要，当前参考版本要求 22.20+ |

## 智能体辅助安装

Codex 的远程仓库安装提示词见 [README](../README.md#通过智能体安装)。脚本同时配置计算环境与完整技能；其他宿主可使用下面的通用请求。

在以本项目为工作目录的会话中发送下面的请求；使用其他工作目录时补充项目绝对路径：

```text
请为当前智能体安装本项目的 cylindrical-gear-analysis。
按本文完成虚拟环境配置、依赖安装和完整技能目录注册；复用已有可用环境。
实际执行随附单工况的 check、run --plots 和 replay，使用新的结果目录。
返回技能注册状态、解释器完整路径、报告位置和 replay.json.passed。
```

安装目标与操作约定见下文。需要界面导入时，应由智能体提供本地包路径和剩余操作；尚未完成宿主启用的状态应单独说明。

## 仓库安装脚本

在稳定目录 clone 本仓库后，任选对应平台入口。默认安装至 `~/.agents/skills/`，运行环境位于仓库 `.venv`；保留完整仓库，不要仅下载 Skill 文本。

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

无需激活虚拟环境。脚本返回完整解释器路径与 JSON 状态；缺少兼容 Python 或依赖安装失败时返回非零退出码。Windows 若限制 `.ps1` 执行，可直接使用 `python -B scripts/update_codex_skills.py --pull`，校验时将参数换为 `--check`。

| 选项 | 行为 |
|---|---|
| 无参数 | 安装当前检出，不拉取远端；可用于完整 ZIP 下载 |
| `--pull` | 干净工作区执行 `git pull --ff-only`，然后安装／更新 |
| `--check` | 只读检查完整技能资源、安装记录、实际导入路径、依赖和内核 SHA256；不联网、不安装、不求解 |
| `--dest PATH` | 自定义技能目录；也可用环境变量 `CODEX_SKILLS_DIR` |
| `--runtime PATH` | 指定虚拟环境目录；后续校验需传入同一参数 |

`GEAR_SKILLS_PYTHON` 可指定创建环境所用的 Python 路径。仅更新本库登记的技能；已有同名目录若来源不同或含本地修改，脚本停止并保留内容。旧版本已登记且未被修改的资源随更新同步；其他技能不受影响。

Claude Code 也可使用相同安装逻辑，将目标设为其用户目录：

```bash
bash scripts/update-codex-skills.sh --dest "$HOME/.claude/skills" --pull
bash scripts/update-codex-skills.sh --dest "$HOME/.claude/skills" --check
```

Windows 同样在 `.ps1` 入口添加 `--dest "$HOME/.claude/skills"`。安装记录保存于目标目录的 `.gear-transmission-skills-install.json`，供智能体读取仓库与解释器位置。仓库或虚拟环境移动后须重新运行安装。`--check` 成功表示文件与运行环境一致，宿主发现和数值复算仍按下文确认。

## 计算环境

脚本已完成环境配置时可跳过本节。手动安装时，在项目根目录执行。

### Windows PowerShell

```powershell
$env:PYTHONUTF8="1"
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -c constraints.txt -e ".[plots]"
```

### macOS／Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -c constraints.txt -e ".[plots]"
```

仅计算时使用 `-e .` 并省略后续 `--plots`。后续调用须使用上述虚拟环境的解释器；`constraints.txt` 记录已有验收环境的主要依赖版本。

## 宿主注册

| 宿主 | 用户级目标 | 项目级目标 | 调用 |
|---|---|---|---|
| Codex | `~/.agents/skills/cylindrical-gear-analysis/` | `.agents/skills/cylindrical-gear-analysis/` | `$cylindrical-gear-analysis` |
| Claude Code | `~/.claude/skills/cylindrical-gear-analysis/` | `.claude/skills/cylindrical-gear-analysis/` | `/cylindrical-gear-analysis` |
| WorkBuddy | 通过技能管理界面导入本地包 | 以客户端支持的入口为准 | 自然语言点名 |

### 目录安装

将完整 `skills/cylindrical-gear-analysis/` 复制或链接到对应目标。保留 `SKILL.md`、`agents/`、`scripts/`、`references/` 和 `assets/`；同名目标已存在时先核对来源。

Windows 示例以 Codex 为目标；Claude Code 将 `.agents\skills` 换为 `.claude\skills`：

```powershell
$gearSkillSource=(Resolve-Path '.\skills\cylindrical-gear-analysis').Path
$gearSkillHome=Join-Path $env:USERPROFILE '.agents\skills'
$gearSkillTarget=Join-Path $gearSkillHome 'cylindrical-gear-analysis'
if (Test-Path -LiteralPath $gearSkillTarget) { throw '技能目录已存在 / Skill directory already exists.' }
New-Item -ItemType Directory -Path $gearSkillHome -Force | Out-Null
Copy-Item -LiteralPath $gearSkillSource -Destination $gearSkillTarget -Recurse
```

macOS／Linux 的 Codex 示例；Claude Code 将 `.agents/skills` 换为 `.claude/skills`：

```bash
mkdir -p "$HOME/.agents/skills"
cp -R ./skills/cylindrical-gear-analysis "$HOME/.agents/skills/"
```

两个宿主均支持链接目录。复制安装需要随源文件更新而同步；链接安装要求原项目位置稳定。新技能尚未被会话发现时，重新打开会话或重启宿主。

### 可选 Skills CLI

[Vercel Skills CLI](https://github.com/vercel-labs/skills) 可替代上述目录复制。选择对应宿主执行一条：

```bash
npx skills add ./ --skill cylindrical-gear-analysis --agent codex --global
npx skills add ./ --skill cylindrical-gear-analysis --agent claude-code --global
npx skills add ./ --skill cylindrical-gear-analysis --agent cursor --global
```

`--global` 表示用户级安装，省略它表示项目级安装。安装器只处理技能文件，Python 计算环境仍须单独配置。

### WorkBuddy

通过技能管理中的 **添加技能 → 上传本地技能包** 导入完整 Skill，并确认已启用。可让智能体依据本机版本准备包；导入后仍须提供完整项目位置和虚拟环境解释器。

WorkBuddy 的包格式与目录约定以当前客户端为准。本库尚未在该客户端实测，不把 CodeBuddy Code 的目录约定作为其统一标准。

## 验证与调用

Windows 下用同一解释器执行：

```powershell
.\.venv\Scripts\python.exe -m gear_transmission check skills/cylindrical-gear-analysis/assets/single_case.json
.\.venv\Scripts\python.exe -m gear_transmission run skills/cylindrical-gear-analysis/assets/single_case.json --output output/example --plots
.\.venv\Scripts\python.exe -m gear_transmission replay output/example --output output/example_replay
```

macOS／Linux 将解释器换成 `.venv/bin/python`。输出目录须为新目录。

| 检查项 | 结果 |
|---|---|
| 宿主注册 | 技能在当前宿主可见，必要时已启用 |
| 实际计算 | `run.json` 中 `status=completed`、`analysis_completed=true`，报告与原数组可读取 |
| 独立复算 | `replay.json` 中 `passed=true` |

在宿主中也可直接请求：

```text
请使用 cylindrical-gear-analysis 和本项目的 Python 虚拟环境，运行随附单工况，
保存报告与图表，并独立复算。
```

`check` 不执行数值求解。宿主发现、计算完成和数值复现分别确认；计算程序可运行不能单独证明宿主已经加载技能。

## 安装约定

用于智能体辅助安装及后续维护：

- 从当前工作目录或用户给定路径定位 `gear-skills.project.json` 与 `pyproject.toml`，仅注册指定宿主。
- 复用或创建项目独立环境，保持安装、运行与复算使用同一 Python；缺少所需解释器时报告未完成项。
- 保留完整技能资源及计算内核，不在安装过程中改写公式或替换登记快照。
- 实际执行计算与复算，返回安装位置、解释器、报告和状态。仅部署文件但尚未确认宿主加载时，明确记录该状态。
- 云端会话需要在实际执行环境安装代码与依赖；本机路径和本地注册不会自动迁移。

## 参考与验证范围

目录与调用方式核对日期：2026-10-30。

- [Claude Code Skills](https://code.claude.com/docs/en/skills) · [OpenAI Skills](https://learn.chatgpt.com/docs/build-skills)。
- [WorkBuddy 官方技能说明](https://cloud.tencent.com/document/product/1831/134432)。
- [Vercel Skills CLI](https://github.com/vercel-labs/skills) · [nature-skills](https://github.com/Yuan1z0825/nature-skills)。

已有记录覆盖本地 Python 调用、数值复算与脚本在 Windows PowerShell／Git Bash 中的安装校验，见[安装验收](../evaluations/installation.json)。Windows／Linux × Python 3.10／3.11 的全新环境安装、实际计算及复算已通过 [CI 验收](../evaluations/ci.json)。macOS 原生环境及各客户端的界面启用尚未逐一实装验收，实际以当前客户端状态为准。
