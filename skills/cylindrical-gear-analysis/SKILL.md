---
name: cylindrical-gear-analysis
description: 调用保留的独立计算内核，分析圆柱齿轮单副的修形、给定失准、承载传递误差TE、法向啮合刚度和Hertz压力，比较给定候选与工况并保存可复算结果。适用于直齿或斜齿、内或外啮合的LTCA计算与修形方案比较。锥齿轮、整套轴承箱体系统、标准强度或疲劳寿命不能由本Skill直接给出。
license: MIT
metadata:
  version: "0.1.0"
---

# 圆柱齿轮承载接触分析

使用随本项目保留的 `gear_contact` 重建计算内核完成计算。新增 `gear_transmission` 仅负责输入适配、保存和复算。运行结果必须来自实际调用，不能凭语言模型记忆或样例响应生成数值。

## 输入与方法选择

先读 [输入与单位契约](references/input-contract.md)，按用户任务整理几何、材料、修形、载荷、失准和离散。只追问影响本任务且当前资料中缺失的信息；用户接受默认值时把假设写入请求的 `assumptions`。

单工况使用 `case`；修形候选比较再给 `study.designs` 与 `study.operating_points`。比较任务的载荷、失准和小轮修形分别放在工况及候选中，不能同时在 `case` 指定。压力上限须来自用户明确依据，不自行推断许用值或隐藏评价权重。

默认模型是 `coupled` 固定齿根弹性降阶模型及局部线接触压陷。ISO模式需要明确的 `iso_settings` 和 `local_contact=false`，避免重复计入局部压陷。不要为得到期望曲线而切换模型或修改采样相位。

外啮合两轮螺旋角等值异号，内啮合同号；内齿圈放在 `wheel`，齿数为正整数。扭矩为已选受载面的非负幅值；失准为有效重叠齿宽右端减左端的法向间隙。修形正值表示材料去除。不要根据数值大小猜单位。

## 执行

环境或导入有问题时读 [运行说明](references/runtime.md)。使用安装本项目的Python执行本Skill的 `scripts/run_ltca.py`；脚本及请求采用绝对路径，不依赖当前工作目录。

```text
python <本Skill绝对目录>/scripts/run_ltca.py check <请求.json>
python <本Skill绝对目录>/scripts/run_ltca.py run <请求.json> --output <尚不存在的结果目录> --plots
python <本Skill绝对目录>/scripts/run_ltca.py replay <原结果目录> --output <尚不存在的复算目录>
```

`--plots` 需要绘图依赖；仅计算可以省略。示例输入在 `assets/single_case.json` 与 `assets/crown_study.json`，只有用户要求示例或输入确实一致时才直接运行。

先校验，再实际求解。`check`通过只说明结构与基本参数有效。命令非零退出或 `run.json` 未完成时，读取诊断并说明失败阶段，不把部分输出写成成功报告。修改输入后在新目录运行。

用户要求可复算交付时，必须执行 `replay` 并检查 `replay.json.passed`。复算会校验原文件、源码和数值环境，并重新求解，不读取原数组作为响应。

## 交付与解释

读 [结果定义](references/results.md)，查看 `report.md`、`run.json`、各工况 `summary.json` 与 `diagnostics.json`；有比较任务时还读 `study.json`。向用户给出主要指标、方案权衡、输入／报告／数据的可用路径、求解与复算状态。

TE法向长度与角度TE分别解释；刚度是N/m，不能解释为扭转刚度。峰压取原求解值，显示插值不能改变指标。瞬时Hertz带与周期包络分别标注；未覆盖区域不能补为确定零压。

明确区分执行完成、求解退出条件通过、数值复现、离散收敛、原组件对标和物理试验验证。当前运行未执行的验证保持未验证。给定候选的非支配结果不表示全局最优，也不保证所有指标同时改善。

用户要求标准安全系数、材料寿命、动态噪声、完整二维自由边界接触或整套传动系统结论时，说明当前能力缺口。可交付本Skill覆盖的部分计算，但不能将其冒充所要求的完整结论。
