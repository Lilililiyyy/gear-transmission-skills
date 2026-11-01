# 新版完整案例与绘图

库版本 `0.3.0` 保留独立内核 `0.2.0` 提交 `c63030d` 的完整运行资源：计算模块、案例编排、后处理、案例输入及公开对标文件。公式不在 Skill 适配层重写。

## 复现独立项目的 18/77 案例

使用安装记录中的 Python 和本 Skill 入口脚本：

```text
python <本Skill绝对目录>/scripts/run_ltca.py case --output <新的结果目录> --plots
python <本Skill绝对目录>/scripts/run_ltca.py replay <原结果目录> --output <新的复算目录>
```

`case` 无位置参数时，读取随包分发的 `examples/helical_18_77_case.json`。它直接执行保留的 `run_cylindrical_workflow.py`，包括显式因子 ISO 6336:2006 强度组合链、LTCA、几何／强度对标、同相位同齿对压力带、MASTA 原图对照，以及原案例固定的三个失准／鼓形敏感性方案。几何图调用保留的 `plot_case_geometry`。

完整案例 JSON 沿用独立项目 SI 协议：`schema_version=1`、`analysis`、`loading`、`strength`、`instant`、`case_id`、`label`。`loading` 统一给名义扭矩与 LTCA 目标乘子；不要同时在 `analysis` 或 `iso_settings` 重复指定扭矩。`strength` 需要显式的评级因子及可选材料应力限；不能从 LTCA 峰压反标。

编辑案例后可执行 `check <案例.json>`、`case <案例.json> --output ...` 或 `run <案例.json> --output ...`。随附 MASTA 图只对应绑定的原案例；原入口会拒绝不同几何／相位的同相位图对比。任意齿轮输入应使用常规 LTCA 请求和 `run --plots`，它不会引入 MASTA 数据。自定义参照的原入口与接口完整保留在 `gear_ltca_snapshot` 资源包中，详见独立项目。

## 两种结果目录

|入口|主要结果|
|---|---|
|`run` 常规请求|`cases/` 原数组；`plots/` 新图；多方案各自的图在 `cases/case_NNN/plots/`，同工况曲线比较在 `plots/point_NNN/`|
|`case` 完整案例|`ltca/` 主计算；`strength.json`；`workflow_summary.json`；`sensitivity/`；根目录 `01` 至 `07` 图表和数据|

完整案例的 `--plots` 会执行三个原案例敏感性方案，复算时自动保留相同求解集合，比较四次计算的全数组及强度结果。它们是给定扰动案例，不是自动优化。

## 图义

- `03_ltca_cycle`：同一稳态周期重复三次，末端追加起点；没有连续瞬态求解。
- `04_contact_pattern`：单相位、单齿对的 Hertz 压力带；可用时采用独立端点压力和半宽，保留原中心值及载荷。
- `05_contact_coordinates`：瞬时与周期采样位置；普通与齿顶延伸分支分别标示。
- `07_project_contact_pass`：普通齿廓材料坐标上的通过程压力，重复点取最大值，按新版三角线性插值绘制，凸包外保留未覆盖状态。延伸分支单独保存。
- 完整案例另有 `02_case_benchmarks`、`06_modification_sensitivity`、`07_native_contact_comparison`；MASTA 原位图只复制用于参照，本项目侧全部重新计算。

中心压力、独立端点压力和显示插值最大值分别解释，不能互换。图表与 SI 数据均存入产物哈希；数值复现不等于离散收敛、完整二维接触精度或物理验证。
