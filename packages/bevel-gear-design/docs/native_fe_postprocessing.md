# 原生有限元结果审核与接触压力图源

`bevel_gear.native_results_workflow`读取正式`FE_boundary`工作流结果，独立重算原响应，再审核本态完整材料域并导出压力图源。后处理不建立刚度矩阵、不调用接触QP或数值分解；支持把阶段分开运行以控制内存峰值。

| 阶段 | 输入与结果 |
|---|---|
| `response` | 严格核对summary、phase、输入及最终operator/响应SHA；从原tet和u/F独立重算BTDBu、Wf/WTu、轴矩与应变能。 |
| `geometry` | 恢复同工作流的原模型档案，用OCP三角化已有BRep，保存完整轮体目标及逐原三角AABB测试的保守裁剪。 |
| `domain` | 使用本态实际u/φ和材料参考Rk；逐请求齿及全部其余实际小轮齿重新判域。仅共享原CAD/原FE的静态几何。 |
| `maps` | 每相每个实际承载齿导出P0守恒中心P1瞬时图源，分别保存完整背景、unknown和认证zero。 |
| `cycle` | 完整闭周期、各态全部承载齿及完整材料/外齿证据通过后，构造离散全周期最大压力包络。 |
| `comparison` | 选择明确瞬时相位/物理齿，导出既有四图工具所需的标准FE瞬时与周期两个文件。 |

常规分阶段调用：

```shell
python examples/audit_native_fe_results.py results/distributed_fe_summary.json --output audited --stage response
python examples/audit_native_fe_results.py results/distributed_fe_summary.json --output audited --stage geometry
python examples/audit_native_fe_results.py results/distributed_fe_summary.json --output audited --stage domain
python examples/audit_native_fe_results.py results/distributed_fe_summary.json --output audited --stage maps
python examples/audit_native_fe_results.py results/distributed_fe_summary.json --output audited --stage cycle
python examples/audit_native_fe_results.py results/distributed_fe_summary.json --output audited --stage comparison --instantaneous-phase-index 4 --instantaneous-tooth-offset 0
```

`4`对应随附九相位案例的α=0；其他workflow使用自己的实际相位编号。也可用`--stage all --include-cycle`顺序完成全部阶段；只处理瞬时结果时使用`--stage all`。`--phase-indices`选择原相位子集，周期阶段仍须满足一个主动齿距的均匀升序闭区间并包含重复端点；离散包络不自动代表相位收敛。

`comparison_fe`内生成`full_gear_fe_instantaneous_triangle_pressure.json.gz`和`full_gear_fe_sampled_cycle_maximum_triangle_pressure.json.gz`。它们直接兼容`load_triangle_pressure_map(..., recompute=True)`。与同工况的两份真实半解析数据放在同一目录后，`examples/replot_contact_composite.py --view raw --data ... --output ... --lang both`可重绘原始四面板，各图独立色标。本入口仅生成FE三角数据，不自动生成工程恢复所需P0格源。默认工程视图需要另外提供同状态绑定的`native_pressure_recovery_sources.json.gz`，详见[压力显示说明](contact_pressure_display.md)。`--pressure-levels levels.json`仍可指定三角等值域导出的线性绝对Pa色阶，须从0开始并覆盖原峰值。

输入必须为正式`bevel_distributed_fe_workflow_v1`及其独立phase完成receipt。reader自动从summary取得原phase SHA，分别核对workflow、phase与实际API执行指纹。被拒绝、未完成、缺失或不同来源的状态不能进入后续认证。旧私有诊断状态不通过补字段转换成正式API结果。

response阶段只需要同SHA原mesh、原E/ν和原物理约束。正式齿轮按原blank的孔径及孔壁容差独立选择固定节点，随后与原归档free DOF逐位比较。明确标为`synthetic_FE_fixture_no_CAD`的小块回归可额外提供`--boundary-conditions`的两体plane JSON；此入口不用于真实齿轮，也不从free DOF反推约束。

geometry/domain阶段要求workflow内`actual_manufacturing`、`actual_cad_sha256`、`cad_files`和`generated_model_archives`。每体模型引用包含descriptor、gear BRep和blank BRep各自的安全相对basename及SHA。加载恢复原grid、制造/blank、面组和BRep，构建原开放shell；不重做实体Boolean。mesh原数组、metadata、规范制造、实际CAD字节和当前模型必须一致。

原输入窗口和真正含压力field的物理齿分别记录。空请求的几何身份保留；没有field的物理齿通过完整CAD工作面距离、同次位移界和材料OUT独立排除，生成自身证据，不增加零压力field。完整域按唯一实际小轮齿核对总数。

所有长度用m、角度用rad、力用N、扭矩用N·m、压力用Pa。P0接触面积与P1显示支撑面积分别记录；中心压力及原积分力保留。蓝色表示低压和有证据的零压，红色表示高压；unknown保留独立掩码并由绘图器显示为未认证域。周期最大压力场的面积积分不表示某个瞬时的合力。

原有限φ接受门控的读取与独立材料域审核分开记录。CAD三角化弦差是工程设置，本接口不将其称为连续NURBS的严格Hausdorff界。实际完整域和周期正向验收需要真实已接受正式结果；合成小块及模型档案回读只验证各自接口。

API入口：

```python
from bevel_gear.native_results_workflow import audit_native_fe_results

result = audit_native_fe_results(
    "results/distributed_fe_summary.json", "audited", stage="response"
)
```

核心reader和单态接口保留在`native_results_archive`；材料/图源在`native_results_pressure`，实际模型目标准备在`native_results_target`。CLI需要NumPy、SciPy、Shapely及psutil，CAD域阶段还需要OCP与VTK。`--maximum-rss-gib`和`--minimum-available-gib`控制本后处理进程；资源中止保留`resource_failure.json`及未完成证据，后续认证拒绝该失败目录。
