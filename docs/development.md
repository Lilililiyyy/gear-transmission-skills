# 开发、内核接入与验证

[English](development_EN.md)

本库维护智能体工作流和计算能力的接入。独立计算项目维护模型、算法、完整接口和案例；当前已接入 [Cylindrical Gear LTCA](https://github.com/Lilililiyyy/Cylindrical-Gear-LTCA)。

## 代码归属

|位置|职责|
|---|---|
|`skills/<skill-name>/`|任务触发条件、执行流程、资源和交付要求|
|`src/gear_transmission/`|显式单位换算、输入检查、函数调用、报告与独立复算|
|`packages/cylindrical-gear-ltca/`|按原字节保留的独立重建计算内核快照|
|`tests/`|输入适配、原数组、失败状态与复算的自动检查|
|`evaluations/`|场景清单、实际验收和独立智能体调用记录|

新增代码说明和诊断采用中文，面向用户的入口文档维护中英版本；原函数、协议字段、单位和索引继续保留。Skill 脚本和工作流调用内核，不重复实现计算公式。

## 新增 Skill

先确定用户任务、必要输入、计算入口和验收条件，再建立 `skills/<skill-name>/SKILL.md`。名称使用小写字母、数字和连字符，描述明确说明何时使用。

按实际需要添加资源：

```text
skills/<skill-name>/
├─ SKILL.md          # 必需：名称、触发描述与必要流程
├─ agents/           # 可选：宿主显示与调用元数据
├─ scripts/          # 可选：可重复执行的辅助入口
├─ references/       # 可选：输入、输出、单位和条件分支
└─ assets/           # 可选：输入案例与输出模板
```

入口保持简洁，条件性细节放进引用文档。依赖、模型选择、缺失输入和失败行为需要明确；安装 Skill 文件不能被描述成同时安装所有计算依赖。

Skill 可先通过显式请求验证，再补充自动触发、缺失信息追问和越界任务的评测。加入 README 技能目录时标明实际状态，不用空目录或占位 Skill 宣称已经接入。

## 更新计算内核

快照身份保存在 [kernel_snapshot.json](../src/gear_transmission/kernel_snapshot.json)，函数、文件哈希及行段对应保存在 [source-map.json](source-map.json)。来源工程含未提交更新时，来源提交号只作背景，实际文件以快照哈希为准。

算法变更先在独立工程完成，再显式更新集成快照、对应记录和验收。保留原公式、数据布局、单位、符号与分支。原工程继续更新时，本库不会自动替换已验收快照，也不在两个位置同时修改算法。

当前集成内核为 `0.2.0`，来源为干净工作树提交 `c63030dc4314a67ba282a3899f7ec8a5f6950fe1`。本次按原字节同步 21 个文件，更新函数映射，并增加 ISO 零载荷内啮合及延伸接触 CSV 回归检查；验收见 [0.2.0 同步记录](../evaluations/kernel-sync-0.2.0.json)。历史验收保持原版本和身份，不用新版本结果改写旧记录。

集成说明见 [内核快照 README](../packages/cylindrical-gear-ltca/README.md)。模型与完整接口优先查阅独立项目；随快照保留的 [模型说明](../packages/cylindrical-gear-ltca/docs/model.md) 和 [接口说明](../packages/cylindrical-gear-ltca/docs/api.md) 对应当前集成版本。

## 本地检查

在完成安装的环境中执行：

```bash
python -m unittest discover -s tests -v
```

已有检查覆盖单位换算、重复输入、材料与几何类型、内啮合斜齿、ISO 设置、零载荷、输出保护、求解失败、原数组和独立复算。修改输入协议或计算调用时，补充与行为有关的检查；普通文档修改核对链接和命令即可。

新增或修改的可执行脚本需要实际运行；复杂 Skill 还应通过独立智能体使用真实请求验证。场景清单见 [scenarios.json](../evaluations/scenarios.json)，列出场景不表示已经执行全部评测。

## 首版验收

实际记录见 [acceptance.json](../evaluations/acceptance.json) 与 [forward-test.json](../evaluations/forward-test.json)：

|检查|已有结果|
|---|---|
|计算工作流自动检查|14 项通过|
|安装自动检查|9 项通过；完整资源、更新、只读校验、修改保护与 Git 快进|
|选定直齿基准与原工程比较|22 组原数组完全一致|
|单工况独立复算|22 组原数组通过|
|3 候选 × 3 工况独立复算|198 组原数组通过|
|独立智能体使用另一组斜齿输入|实际调用、绘图及 22 组数组复算通过|
|wheel 安装后的执行与复算|通过|
|远端 CI|Windows／Linux × Python 3.10／3.11 四项通过，见 [CI 记录](../evaluations/ci.json)|

执行完成、求解退出条件通过、数值复现、空间／相位收敛、原组件响应精度和物理试验验证分别记录。未执行的验证保持未验证。独立项目的模型更新或验证结论，也须核对是否覆盖当前快照和所用工况。

## 组织约定

本项目按职责划分技能入口、工作流适配、独立计算内核快照、测试与验收记录。新增技能应保持入口简洁，将条件性细节放入 references，将可重复执行的辅助逻辑放入 scripts，并将案例输入与输出模板放入 assets。计算算法由对应内核项目维护，集成层负责参数适配、单位转换、调用和结果交付。
