# Gear Transmission Skills

[中文](README.md) · [English](README_EN.md)

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/) [![License: AGPL-3.0-only](https://img.shields.io/badge/License-AGPL--3.0--only-blue.svg)](LICENSE)

齿轮传动设计的 Agent Skills 库，通过自然语言调用独立计算内核，完成参数分析、方案比较和结果复算。

## 目录

- [可用技能](#可用技能)
- [安装](#安装)
- [使用](#使用)
- [输出与验证](#输出与验证)
- [文档与开发](#文档与开发)
- [更新记录](#更新记录)

## 可用技能

| Skill | 主要能力 | 状态 | 独立计算项目 |
|---|---|---|---|
| [cylindrical-gear-analysis](skills/cylindrical-gear-analysis/SKILL.md) | 圆柱齿轮 LTCA、修形与失准分析、TE／法向刚度／Hertz 压力、多工况方案比较 | 可用 | [Cylindrical Gear LTCA](https://github.com/Lilililiyyy/Cylindrical-Gear-LTCA) |

计算模型、完整接口和案例由独立项目维护。本库提供智能体工作流、输入适配、报告与复算交付。

## 安装

### 通过智能体安装

向可执行本地命令的 Codex 发送：

```text
请从这个仓库安装 Codex skills：
https://github.com/Lilililiyyy/gear-transmission-skills.git

请 clone 到稳定的本地目录，然后运行 scripts/update-codex-skills.sh --pull。
Windows PowerShell 请使用等价入口 .\scripts\update-codex-skills.ps1 --pull。
安装后用同一入口运行 --check 验证。
请保留 skills/ 下的完整技能目录，不要只复制 SKILL.md，也不要删除本地仓库。
最后使用脚本返回的 Python 解释器运行随附单工况（含图表）并独立复算，返回报告路径与结论。
```

脚本配置项目 `.venv`，安装计算与绘图依赖，并将完整 Skill 同步至 `~/.agents/skills/`。需要 Python 3.10／3.11 和 Git；`--pull` 只接受干净工作区的快进更新，`--check` 只读校验资源、导入路径和内核哈希，不执行求解。保留本地仓库用于后续计算与更新。

### 命令行安装

macOS／Linux 或 Git Bash：

```bash
git clone https://github.com/Lilililiyyy/gear-transmission-skills.git
cd gear-transmission-skills
bash scripts/update-codex-skills.sh --pull
bash scripts/update-codex-skills.sh --check
```

Windows PowerShell：

```powershell
git clone https://github.com/Lilililiyyy/gear-transmission-skills.git
cd gear-transmission-skills
.\scripts\update-codex-skills.ps1 --pull
.\scripts\update-codex-skills.ps1 --check
```

后续更新：在同一检出目录重新运行 `--pull`，再运行 `--check`。自定义目录、已有环境与其他宿主见[安装指南](docs/installation.md)。

### 手动安装

在 Python 3.10／3.11 虚拟环境中，于项目根目录执行：

```bash
python -m pip install -c constraints.txt -e ".[plots]"
```

随后注册完整的 `skills/cylindrical-gear-analysis/`：

| 智能体 | 注册方式 | 显式调用 |
|---|---|---|
| Codex | 用户级 `~/.agents/skills/` | `$cylindrical-gear-analysis` |
| Claude Code | 用户级 `~/.claude/skills/` | `/cylindrical-gear-analysis` |
| WorkBuddy | 技能管理界面导入本地技能包，并启用 | 自然语言点名 |
| Cursor | 可通过 `npx skills` 选择 `cursor` 目标 | 自然语言点名或技能选择器 |

平台命令、可选的 Skills CLI 和 WorkBuddy 接入说明见 [安装指南](docs/installation.md)。技能注册与 Python 环境安装均须完成。

## 使用

可以直接描述工况和设计目标：

```text
请使用 cylindrical-gear-analysis 分析外啮合直齿轮：28/56 齿，法向模数 2.5 mm，
齿宽 25 mm，压力角 20°，两轮 E=210 GPa、ν=0.3，失准 12 μm。
比较小轮 0、8、12 μm 齿向鼓形在 20、80、120 N·m 下的 TE 与压力。
采用 coupled 模型与局部压陷，32 相位、21 条带；列出其余默认值。
保存报告、图表与原始数据，并独立复算。
```

也可在已安装环境中直接运行案例：

```bash
python -m gear_transmission run skills/cylindrical-gear-analysis/assets/single_case.json --output output/example --plots
python -m gear_transmission replay output/example --output output/example_replay
```

输出目录须为新目录。随附输入：[单工况](skills/cylindrical-gear-analysis/assets/single_case.json) · [多工况修形比较](skills/cylindrical-gear-analysis/assets/crown_study.json)。

## 输出与验证

- `report.md`：主要指标、方案权衡与模型范围。
- `run.json` 与 `cases/`：输入、代码身份、原始数组和求解诊断。
- `replay.json`：重新求解后的输入、汇总量与全数组比较。

已有 14 项计算工作流检查、9 项安装检查，以及原工程基准对比和单／多工况复算记录，见[计算验收](evaluations/acceptance.json)与[安装验收](evaluations/installation.json)。数值复现、离散收敛和物理验证分别记录；模型精度范围见独立计算项目。

Windows／Linux × Python 3.10／3.11 的全新环境安装、实际计算与独立复算已全部通过 [CI 验收](evaluations/ci.json)。

## 文档与开发

- [源码来源审查（初步）](https://github.com/Lilililiyyy/gear-transmission-skills/blob/docs/unified-readme-2026/docs/source-provenance-review.md)：记录三个仓库的来源线索、内部复用关系与审查边界。

- [安装指南](docs/installation.md)：环境、宿主注册、自然语言安装与复算检查。
- [开发说明](docs/development.md)：目录职责、内核快照与新增 Skill 的约定。
- [输入契约](skills/cylindrical-gear-analysis/references/input-contract.md) · [结果定义](skills/cylindrical-gear-analysis/references/results.md)。

`skills/` 保存技能，`packages/` 保存内核快照，`src/gear_transmission/` 提供调用与复算。计算源码和原接口保持独立维护。

## 许可与商业授权

本项目采用 [AGPL-3.0-only](LICENSE)。如需闭源集成，可[联系作者](https://github.com/Lilililiyyy)获取商业授权。

## 更新记录

- **2026-10-30 · v0.1.0**
  - 首个可运行技能 `cylindrical-gear-analysis`：保留圆柱齿轮 LTCA 内核，支持单工况、修形方案比较、报告、图表与独立复算。
  - 增加跨平台 Codex 安装入口：`--pull` 更新并安装，`--check` 校验完整资源与运行环境，保护已有同名技能和本地修改。
  - 提供中英双语 README、安装与开发说明；详细计算模型链接至独立项目。
  - 首次发布至 GitHub，通过 Windows／Linux 与 Python 3.10／3.11 的四项 CI 验收。
