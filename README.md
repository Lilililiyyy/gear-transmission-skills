# Gear Transmission Skills

[中文](README.md) · [English](README_EN.md)

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/) [![License: AGPL-3.0-only](https://img.shields.io/badge/License-AGPL--3.0--only-blue.svg)](LICENSE)

齿轮传动设计的 Agent Skills 库，通过自然语言调用独立计算内核，完成参数分析、方案比较和结果复算。

当前库版本为 **0.4.1**。圆柱 LTCA 内核 **0.2.1** 固定至 `cb2e1ec`，包含统一的接触图与重绘接口；锥齿轮内核 0.10.0 固定至 `05da3ac`。

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
| [cylindrical-gear-analysis](skills/cylindrical-gear-analysis/SKILL.md) | 圆柱齿轮 LTCA、修形／失准比较、新版接触图、显式因子标准强度复算、完整案例对标 | 可用 | [Cylindrical Gear LTCA](https://github.com/Lilililiyyy/Cylindrical-Gear-LTCA) |
| [bevel-gear-design](skills/bevel-gear-design/SKILL.md) | 零偏置直齿／弧齿的几何、啮合力、名义强度及分阶段CAD／TCA／节点FE-LTCA，结果独立复算 | 可用；接触需可选依赖 | [Bevel Gear Design](https://github.com/Lilililiyyy/Bevel-Gear-Design) |

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

随后注册 `skills/` 下所需的完整技能目录。默认安装器同时注册圆柱与锥齿轮技能；锥齿轮显式调用为 `$bevel-gear-design`：

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

多方案的标准接触比较图位于 `plots/point_NNN/07_contact_pass_comparison.png`。已有结果可按新版显示规则重绘：

```bash
python -m gear_transmission plot output/example --output output/example_replot
```

该命令校验原文件、读取原数组并生成新图，`plot.json` 明确记录未重新求解；原图和原结果保留。显示范围采用两轮实际啮合区，不能将非啮合齿根误判为漏算，或通过补压把偏载斑点画成闭合椭圆。

复现独立项目新版 **18/77 斜齿轮完整案例及接触图**：

```bash
python -m gear_transmission case --output output/helical_18_77 --plots
python -m gear_transmission replay output/helical_18_77 --output output/helical_18_77_replay
```

自然语言也可直接要求：“使用 cylindrical-gear-analysis 复现新版 18/77 完整案例，生成整齿接触压力、三周期响应、同相位压力带及显式因子强度对比，并独立复算。”

普通 `run --plots` 同样输出新版几何、三周期、瞬时压力带与整齿通过程图，采用用户实际输入；随附 MASTA 参照只在完整案例中用于对标。强度使用显式输入的评级因子及材料应力限。图义、目录及适用范围见[新版工作流](skills/cylindrical-gear-analysis/references/new-workflow.md)。

锥齿轮可以从实际几何与强度计算开始：

```text
使用 bevel-gear-design 校核随附17/47弧齿锥齿轮示例，计算几何、啮合力和名义强度。
明确示例中的材料应力限与制造假设，保存图表和原始结果，并独立复算。
```

```bash
python -m gear_transmission.bevel_cli check skills/bevel-gear-design/assets/spiral_17_47_case.json
python -m gear_transmission.bevel_cli run skills/bevel-gear-design/assets/spiral_17_47_case.json --output output/bevel --plots --reference skills/bevel-gear-design/assets/spiral_17_47_reference.json
python -m gear_transmission.bevel_cli replay output/bevel --output output/bevel_replay
```

锥齿轮 `--through geometry/strength/cad/tca/ltca` 默认到strength；CAD与接触另装 `python -m pip install -e ".[bevel-contact,plots]"`。统一LTCA保持节点FE路线；原生分布牵引FE和局部半解析API的适用范围、原始场及复算边界见[锥齿轮运行说明](skills/bevel-gear-design/references/runtime.md)。准双曲面与整机寿命不在当前流程内。

## 输出与验证

- `report.md`：主要指标、方案权衡与模型范围。
- `run.json` 与 `cases/`：输入、代码身份、原始数组和求解诊断。
- `replay.json`：重新求解后的输入、汇总量与全数组比较。

完整 `case` 另存 `workflow_summary.json`、`strength.json`、`ltca/` 与 `sensitivity/`；所有结果以 `run.json.cases` 登记的目录为准。

当前共 38 项工作流与安装检查；接触显示域、源工程一致性、重绘与复算见[0.4.1 验收](evaluations/contact-display-0.4.1.json)。早期的完整资源、图表来源对比与复算见[0.3.0 同步验收](evaluations/full-sync-0.3.0.json)。[0.2.0 内核验收](evaluations/kernel-sync-0.2.0.json)、[首版计算验收](evaluations/acceptance.json)与[安装验收](evaluations/installation.json)保留各自版本记录。数值复现、离散收敛和物理验证分别记录；模型精度范围见独立计算项目。

Windows／Linux × Python 3.10／3.11 的安装、实际计算和独立复算纳入[持续验收](https://github.com/Lilililiyyy/gear-transmission-skills/actions/workflows/ci.yml)，首版结果见 [CI 记录](evaluations/ci.json)。

0.4.0 新增锥齿轮5项行为检查后，当时本地共35项通过。锥齿轮直齿／弧齿强度、同工况参照、安装后及wheel调用与复算、直齿单相位CAD／TCA／节点FE-LTCA及独立请求前向验证见[锥齿轮技能验收](evaluations/bevel-skill-0.1.0.json)。单相位求解及复现通过不代表全周期、离散峰压收敛或物理试验验证。

## 文档与开发


- [安装指南](docs/installation.md)：环境、宿主注册、自然语言安装与复算检查。
- [开发说明](docs/development.md)：目录职责、内核快照与新增 Skill 的约定。
- [输入契约](skills/cylindrical-gear-analysis/references/input-contract.md) · [结果定义](skills/cylindrical-gear-analysis/references/results.md)。

- [锥齿轮输入与单位](skills/bevel-gear-design/references/input-contract.md) · [结果与限制](skills/bevel-gear-design/references/results.md) · [来源映射](docs/bevel-source-map.json)。

`skills/` 保存技能，`packages/` 保存内核快照，`src/gear_transmission/` 提供调用与复算。计算源码和原接口保持独立维护。

## 许可与商业授权

本项目采用 [AGPL-3.0-only](LICENSE)。如需闭源集成，可[联系作者](https://github.com/Lilililiyyy)获取商业授权。

## 更新记录

- **2026-11-02 · v0.4.1**
  - 修正圆柱接触图把完整渐开线范围误作实际啮合范围的问题；中心条带的显示支撑仅延拓到其自身积分边界，原压力与 PPTE 不变。
  - 单图、源项目案例与多方案比较共用源工程 `contact_display`，统一色带、等压线与坐标；标准流程自动生成统一色标的接触比较图。
  - 新增 `plot`：先校验保存结果，再只读原数组重绘；明确区分重绘与重新求解／复算。修复验收见 [0.4.1 记录](evaluations/contact-display-0.4.1.json)。

- **2026-11-02 · v0.4.0**
  - 新增 `bevel-gear-design`，接入0.10.0锥齿轮重建内核、两组原生案例和制造输入，保留逐文件哈希与函数行段来源。
  - 新增 `gear-bevel` 及校验、分阶段计算、绘图、同工况参照比较与实际独立复算；安装器同步两个完整Skill。
  - 直齿／弧齿几何强度、来源及文件篡改检查纳入自动验收。CAD、TCA、LTCA及不同高级接触路线的端到端状态分别记录，不由强度案例通过推断。


- **2026-11-02 · v0.3.0**
  - 补齐独立项目的案例编排、绘图后处理及公开对标输入，完整登记 47 个来源文件并在运行时校验。
  - 新增 `case` 完整案例入口，直接调用原项目新版工作流，支持显式因子强度、同相位接触对标、整齿压力和敏感性分析。
  - 普通工况接入新版绘图，区分中心压力、独立端点压力、瞬时压力带与通过程包络；保留三周期原样展开及 SI 复绘数据。
  - 增加全资源同步／校验脚本、全数组及强度复算、图表来源对比和 wheel 分发验收。

- **2026-11-01 · v0.2.0**
  - 集成独立圆柱齿轮项目 `0.2.0`，按提交 `c63030d` 的干净工作树更新 21 个文件、来源哈希和函数映射。
  - 同步 ISO 未承载接触点的曲率处理与延伸接触 CSV 标记导出，并补充对应回归检查。
  - 统一库、Skill 和内核版本，修正 Skill 入口的许可标记及复算中相同无穷值的差值诊断；重新验证安装、来源一致性与独立复算。

- **2026-10-30 · v0.1.0**
  - 首个可运行技能 `cylindrical-gear-analysis`：保留圆柱齿轮 LTCA 内核，支持单工况、修形方案比较、报告、图表与独立复算。
  - 增加跨平台 Codex 安装入口：`--pull` 更新并安装，`--check` 校验完整资源与运行环境，保护已有同名技能和本地修改。
  - 提供中英双语 README、安装与开发说明；详细计算模型链接至独立项目。
  - 首次发布至 GitHub，通过 Windows／Linux 与 Python 3.10／3.11 的四项 CI 验收。
