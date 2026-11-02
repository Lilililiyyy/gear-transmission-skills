---
name: bevel-gear-design
description: 调用独立重建计算内核，完成零偏置直齿、弧齿锥齿轮的几何、啮合力、名义强度、实际齿面与CAD、无载TCA及节点FE承载接触分析，保存输入、原始结果并独立复算。用于锥齿轮参数设计、已给方案校核及接触分析；不适用于准双曲面、圆柱齿轮或未经建模的整机寿命与噪声预测。
license: AGPL-3.0-only
metadata:
  version: "0.1.0"
---

# 锥齿轮设计与分析

使用本库按原字节保留的 `bevel_gear` 重建计算内核。`gear_transmission.bevel_workflow` 是新写调用与复算适配，不在适配层另写接触公式；所有数值结论必须来自实际执行。

## 选择任务与输入

先读 [输入与单位](references/input-contract.md)。确认直齿或圆弧齿线、零偏置相交轴，以及用户要求的终止阶段。已有方案直接校核；参数选型可比较用户给定的候选，记录目标、约束和假设，不把有限候选比较称为全局优化。

只追问当前阶段所缺的关键参数。名义几何不等于制造齿面；实际弧齿齿面需要两轮各自的机床、刀具和齿坯输入，不从螺旋角推造机床设置。材料许用应力、载荷系数和目标安全系数须有明确来源，不自动借用示例值。

采用原生 `schema=bevel_analysis_case_v1`。公开输入统一为m、rad、Pa、N和N·m，转速字段为rpm；保留 `ltca_torque_sign`、工作侧与旋向约定。工程单位换算显式写入输入说明，不根据数值大小猜单位。

## 运行

用完整项目安装环境中的Python执行本Skill的 `scripts/run_bevel.py`；入口、请求与输出使用绝对路径。环境问题读 [运行与模型路线](references/runtime.md)。

```text
python <本Skill绝对目录>/scripts/run_bevel.py check <案例.json> --through strength
python <本Skill绝对目录>/scripts/run_bevel.py run <案例.json> --through strength --output <新结果目录> --plots
python <本Skill绝对目录>/scripts/run_bevel.py replay <原结果目录> --output <新复算目录>
```

`--through` 可选 `geometry`、`strength`、`cad`、`tca`、`ltca`，默认 `strength`。后续阶段执行前置阶段，不能把名义强度输出当作已运行LTCA。`check` 实际校核代数几何及所选评级输入，不执行CAD、接触平衡或有限元。

两组随附案例是演示输入：`assets/straight_16_32_case.json` 与 `assets/spiral_17_47_case.json`。弧齿案例须保留相邻的两份制造JSON；只在用户要求示例或输入一致时直接运行。可用 `--reference <同工况参照.json>` 核对几何与强度，读取 `same_case_inputs` 与逐项比较，不能将参照值作为计算响应。

统一 `ltca` 入口使用原节点FE路线，孔壁固定、准静态、小变形线弹性。新版分布牵引FE与局部半解析通过各自原生API调用，选用前读运行文档；它们没有自动并入此命令的复算协议。不能悄悄改用另一模型来消除报错或改善图形。

命令非零退出或 `run.json.status != completed` 时读取诊断、逐相位状态及残差，说明失败阶段。失败结果保留在原目录，修改输入后使用新目录。可复算交付必须实际执行 `replay` 并检查 `replay.json.passed`。

## 解释与交付

读 [结果定义与限制](references/results.md)，核对 `run.json`、`calculation/workflow.json` 和适用阶段的原始场。向用户交付主要几何／力／应力／安全系数或接触指标，以及输入、报告、图表、数据和复算文件路径。

区分执行完成、各相位求解接受、设计目标、数值复现、网格／相位收敛、外部对标及物理验证。未做的验证保持未验证。标准接触应力不等于LTCA峰压；总TE不等于纯弹性TE；从动轴切线扭转刚度为N·m/rad，不能称为N/m。

局部半解析闭合量不含完整轮体柔度，不能解释为整轮总TE。瞬时受载区与周期最大压力包络分别标注；未解析区域保留未知，显示插值不产生新的平衡解。当前统一强度链是2014名义应力加显式设计限值，不自动给出完整材料寿命评级。准双曲面、真实直齿刀具展成齿根、整套支承、润滑温升、动态噪声或实验精度超出此已接入流程时，明确缺口并交付可完成的部分。
