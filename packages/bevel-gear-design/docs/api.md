# 接口、单位与结果布局

计算模块均为独立实现，函数及数据类的完整参数以源码签名为准；不要求通过 CLI 调用。当前输入和复算命令统一放在项目 README。

## 主要入口

|任务|模块／接口|输入与输出|
|---|---|---|
|统一直齿／弧齿分析|`workflow.load_bevel_case`、`analyze_bevel_case`|案例 JSON → 几何、强度、CAD、TCA、LTCA；`through` 选择终止阶段|
|名义设计与评级|`analyze_design(config)`、`read_config`|显式几何、载荷、两轮材料和评级设置 → `inputs/geometry/forces/rating/diagnostics`|
|独立几何与力|`geometry.calculate_geometry`、`forces.calculate_mesh_forces`|类型化几何与扭矩 → 几何和名义力|
|2014 名义应力|`evaluate_iso2014_nominal`、齿轮对象的 `nominal_strength_2014`|显式材料、齿距偏差、刀尖和载荷 → `stress_pa/coefficients/root_sections/dynamic/contact_lines`|
|几何与 CAD|`manufacturing`、`generated_surface`、`generated_solid`、`straight_cad`|机床／刀具／齿坯或直齿配置 → 有限工作面、实体、STEP/BREP|
|中性 CAD|`cad_import.load_step_model`、`triangulate_imported_model`|STEP → 真实拓扑与诊断、SI 网格|
|整轮分布牵引 FE|`generated_distributed_fe.build_distributed_contact_system`、`solve_distributed_contact_system`、`solve_generated_distributed_fe_ltca`|实际有限NURBS与同一整轮体网格 → 原P0压力、结构响应、几何外迭代、TE及刚度|
|节点 FE 兼容入口|`refined_fe_workflow.run_refined_fe_workflow`|制造输入、材料、节点网格、相位及有符号扭矩 → 摘要和逐状态数组|
|局部半解析|`surface_workflow.run_surface_workflow`|实际工作面、网格／步长和求解器设置 → 局部压力与闭合|
|覆盖审核|`pressure_domain_workflow`、`pressure_domain_batch`|原状态与实际有限齿面 → 已知／未决域、窗口覆盖；不改变原压力|
|原生FE分片与响应审核|`native_pressure_domain.audit_native_partition_geometry`、`audit_native_response_archive`|实际FE节点、交片及原W/u/F档案 → 平均g0、虚功、残差和能量的独立证据；不创建K或求新u|
|周期压力|`contact_cycle`、`contact_cycle_plotting`|原分相状态 → 参考齿面最大值、argmax、原峰与未知掩码|
|原三角压力总图|`contact_composite.build_triangle_pressure_map`、`plot_contact_composite`|原三角压力、实际齿面轮廓及各自覆盖证据 → FE／半解析×瞬时／采样最大四面板|
|工程压力总图|`contact_display.prepare_pressure_comparison`、`plot_pressure_comparison`|原三角压力及FE原P0格 → 半解析在前、各自色标的四图；`view=engineering/raw`区分工程恢复与原始离散场|

## 单位与符号

|物理量|公开接口单位与定义|
|---|---|
|长度、模数、间隙、齿宽、粗糙度|m；不根据数值大小猜单位|
|角度／角速度|rad／rad/s；带 `rpm` 字段的转速仍为 rpm|
|力／扭矩|N／N·m；接触扭矩保留参数轴符号及工作侧|
|弹性模量、压力、应力|Pa；图中可用 GPa／MPa|
|接触柔度|m/N|
|从动轴切线扭转刚度|N·m/rad；不是 N/m 或 T/θ|
|接触／投影面积|m²，分别保留真实面积和参考切平面面积|
|材料密度、运动黏度|kg/m³、m²/s|
|时间、功率|h、W，字段显式标注|

CAD 内部长度约定 mm，`kernel_length_unit_m=0.001`；公开网格顶点为 m，三角形索引从 0 开始。STEP 重导声明 mm；BREP 无内置单位，沿用 mm 约定。开放或无效实体体积为 `null`，不能解释成零体积。

接触案例 driver/follower 表示小轮给定相位／大轮待求相位，不指定功率流。−64 N·m 对应案例的 coast／制动工作侧；相位扫描顺序不是实际转向。名义常幅评级入口另要求正扭矩幅值、正转速及非负持续时间。

`w_f=R_e−(r sinδ₁+z cosδ₁)` 从外端向内端为正，`h=r cosδ₁−z sinδ₁` 沿节锥法向朝齿顶为正；r 为小轮材料坐标的圆柱半径。两者在源数据中为 m，图中为 mm；该子午投影不是等面积展开。

## 案例与结果

统一案例使用 `schema=bevel_analysis_case_v1`。独立瞬时／周期接触配置通过相对 `base_case` 绑定实际制造输入；制造 `request_file` 相对案例文件解析。复制案例时须一并复制其引用输入，不改变数组的相位、齿对与节点对应关系。

`analysis_completed`、各阶段接受状态、几何／强度参照通过、设计目标、网格收敛和物理验证是不同状态。失败、未知和未实现分支保留明确错误或 `null` 及解释，不能以零补齐。

统一评级的 `capacity` 为寿命／表面条件修正后的容量，`allowable` 为容量除以目标安全系数，`safety=capacity/stress`。`meets_targets` 仅说明当前常幅疲劳指标；静载分配设置不执行屈服或完整静强度检查。寿命 `finite/endurance_assumption/out_of_curve` 分别表示有限寿命、曲线末端耐久假设和曲线外未知；AGMA 显式一循环策略须同时读取 `assumption`。

FE 保留每相位的原节点力、压力、三维积分面积、活动集、实际轴力臂、TE 分项、刚度和残差。局部半解析保留参考切平面格、有效掩码、真实曲面点／角点及雅可比，投影面积和真实面积不能互换。具体数组字段和索引由对应状态 JSON／NPZ 的具名字段给出，不按不同相位的数组长度硬拼。

精选结果位于 `examples/results/`。周期 `*.npz` 为无 pickle SI 二维布局，`*.json` 保存模型定义、相位原峰值、几何掩码、来源与数组 SHA256；`data_files` 绑定随附唯一 NPZ。`maximum_known=0` 只表示当前值是下界，不表示确定零压。瞬时压缩 CSV 可用 `gzip.open(..., mode="rt", encoding="utf-8-sig")` 读取，首行保留具名列及单位。

## 重绘与扩展

分布牵引FE结果使用 `contact_representation="distributed_P0_surface_cells"`。`selected_tooth_offsets` 与 `fields` 一一对应；每场的 `pressure_pa`、`forces_n` 和布尔 `valid_mask` 采用二维 `[row, column]` 布局，`cell_area_m2` 为参考面元面积，原力为 `pA`。`geometry.patches` 保存对应实际曲面中心、四角、参数和参考平面。无效掩码外的数组零仅为存储占位，不构成物理零压证明。

显式 `geometry_source="FE_boundary"` 的分支以实际线性FE边界的双三角交片积分平均g0，不能调用NURBS重采样或以中心间隙替代。`audit_native_partition_geometry(system, meshes)` 从全局节点和机械整体R0重建真实面积/质心、P1形函数、g0、W及刚体力臂；`audit_distributed_response` 对此分支先执行该独立几何复算，再核原Wf、Wᵀu、Ku和能量。各物理齿的 `source_material_reference_frames` 保存折回参考材料齿所需Rk，不能用它替换机械R0。

`audit_native_response_archive(state_path, run_root=..., native_system=..., meshes=..., materials=..., fixed_nodes=..., mesh_source_sha256=...)` 读取原生正式JSON/gzip状态及其外置receipt，复核输入指纹、状态/NPZ/W原字节、原压力格列序和固定自由度，默认分块回装Ku而不形成新K。结果只说明这一次响应绑定与平衡的检查范围；完整材料域、逐片非穿透、压力网格及相位收敛另行审核。未求卸载FE网格相切时，`linearization_step_theta_rad` 仅为当次参考角增量，`total_loaded_te_rad=loaded_follower_angle_rad+ratio*driving_angle_rad` 为总受载TE，不从增量推断纯弹性TE。

有限角外迭代若生成了新的交片/W/g0/力臂，最终 `native_system` 必须由该次 operator 档案恢复。`current_operator_manifest_sha256` 与 `operator_fingerprint` 成对绑定结果和真实系统；原 `native_cache_index_sha256` 仅保留初始静态来源身份。正式读取时，结果的 `current_operator_manifest_file` 为运行目录内相对路径，loader 系统保留同一文件的绝对路径，入口核验原文件 SHA。不能把初始缓存的 W 或交片与最终位姿混用。

`pressure_field.reconstruct_distributed_fe_pressure_field(result, tooth_offset=0)` 生成守恒中心P1显示。原P0压力未知量、承载面积和Gauss力矩不改变；缺少完整承载邻域时拒绝外推。节点兼容入口与半解析分别使用自己的重构函数，不能混用压力定义。

`contact_composite.export_triangle_pressure_map(data, "field.json.gz")` 无损保存SI源三角与等值域；`load_triangle_pressure_map(path, recompute=True)` 核验数值SHA并从源标量重建。`examples/replot_contact_composite.py` 只依赖绘图库及随附源表，输出中英总图。默认工程视图额外需要`native_pressure_recovery_sources.json.gz`，采用`contact_recovery.recover_grid_pressure`逐齿保持原参考平面法向合力，`advected_cycle_maximum`仅连接相邻真实承载相位、同一物理齿的平移标量场。输出的原始／显示峰值、实际／显示相位数、压力SI格及恢复误差分别保存。未有该P0源时选择`--view raw`，不能从图片反推。细节见[压力显示说明](contact_pressure_display.md)。

`distributed_fe_results.load_distributed_fe_cycle(summary_file)` 逐相位验证原状态SHA、相位、输入指纹及源身份，再汇总原P0周期响应；`summarize_distributed_fe_cycle(states)` 可直接接收原JSON状态列表。汇总分别保存从动轴无载TE、参考相位修正、弹性TE、总TE和切线刚度，直接核对实际受载角。逐齿压力积分采用各自的参考平面面积；`force_plane_fraction` 不代表共同法向的矢量合力分担。失败相位保留断点，不补零，所选相位完成不代表相位极值收敛。`examples/plot_distributed_cycle.py` 可由原摘要或随附SI汇总重绘双语响应图。

`prepare_generated_ltca(..., refinement_core_radius_m=...)` 将细体网格延伸至显式接触近场厚度，而不是仅加密表面一层。实际深度分辨率仍须检查。`factorization_factory(Kff)` 可接入同一刚度的专业SPD分解器；返回对象须有 `solve`、真实分解完成后的 `positive_definite=True` 与非空 `analysis_diagnostics`，不伪造LU主元。默认SciPy引擎可独立运行。

`contact_cycle_plotting.load_mesh_cycle_pressure_data(metadata_path)` 验证随附数组哈希后恢复周期结果；`plot_mesh_cycle_pressure(data_sets, output, lang="zh")` 直接生成图，不依赖 CAD 或原 FE 求解过程。

求解与绘图分离。局部接触核、有限元准备、约束求解和后处理模块均可直接调用；外部支承或柔度扩展须显式保留单位、边界和能量定义，避免重复计入局部压陷。参考公式见同目录 `model.md`，具体命令、案例图和当前验证限制集中在 README。

## 真实有限元边界的接触积分API

本接口从两侧真实线性三角齿面构造原生`FE_boundary`接触operator，适用于直齿或弧齿锥齿轮已生成的FE齿面。输入几何、姿态和请求网格均参数化；不依赖固定案例、缓存目录或外部软件组件。

### 能力和适用域

两侧齿面须在指定共同平面上构成连续、单值的P1有限元图。程序保留全部真实三角和正面积`cell ∩ tri₁ ∩ tri₂`，用各凸片面积和质心精确积分仿射FE高度、形函数与刚体速度，输出P0压力自由度的稀疏载荷矩阵W、间隙g和力臂ℓ。

同体投影重叠、高度不连续、退化三角、无效刚体旋转、重叠压力格或非法全局索引均明确报错。缺失齿面不会通过凸包、最近点、假节点或重新归一化填补。

这些函数不装配体刚度、不求接触压力，也不把局部积分通过当作整个齿面或整个啮合周期收敛。生产受载分析应将每次实际W/g/ℓ交给FE接触求解器，并在更新参考姿态后重新构造operator和柔度。

### 函数

```python
from bevel_gear.native_fe_operator import (
    compact_fe_boundary,
    build_native_fe_contact_operator,
    build_native_operator_at_pose,
    advance_native_follower_frame,
    check_native_operator_pose,
    save_native_operator,
    load_native_operator,
)
```

- `compact_fe_boundary(mesh, name, *, boundary_triangles=None)`：从含`nodes_m`和`boundary_triangles`的对象或Mapping提取选定实际边界，保留原全局节点ID。应选择真实的单侧齿面；整个封闭体边界通常不是单值图。
- `build_native_fe_contact_operator(compact, patches, frames, reference_follower_angle_rad, *, driving_angle_rad=0., pinion_torque_nm=None, pinion_to_wheel_ratio=None, batch_cells=128, source_identity=None)`：不要求历史缓存的新建入口。
- `build_native_operator_at_pose(base_system, compact, frames, reference_follower_angle_rad, *, batch_cells=128, progress=None)`：保持原请求集合、格和平面，在新姿态重建W/g/ℓ。可用于有限转角外迭代。
- `advance_native_follower_frame(frames, delta_rad)`：绕从动轮实际世界轴更新整体旋转及显式平移；主动轮姿态不变。
- `check_native_operator_pose(system, compact, u, theta_rad, reference_remaining_gaps_m, *, batch_cells=128, mean_tolerance_m=1e-10, maximum_penetration_m=2e-7, progress=None)`：核原交片全部顶点，并在`φ+θ`重新构造真实FE图，核有限姿态平均间隙及顶点非穿透。返回`completed`、`passed`、逐cell残差和覆盖证据；u仍为原局部P1小变形迹，不是大滑动求解。
- `save_native_operator(system, directory, *, actual_mesh_sha256=None, source_identity=None)`与`load_native_operator(directory)`：保存和严格回读typed交片、原W和完整geometry。目录已存在时拒绝覆盖。

更底层的`ProjectedFacetBody`、`ArrayProjectedFacetBody`和`ArrayProjectedFEIntegrator`分别位于同名数学模块。`ArrayProjectedFEIntegrator.build_cells(cells)`接受显式矩形界及整数`tooth_offset/row/column`，也适合不使用规则网格的调用方。

### 输入布局

`compact`为`{"pinion": {k: facet}, "wheel": {k: facet}}`。单一`k=0`时，也可直接令两成员值为facet字典。每个facet包含：

| 字段 | 物理意义与布局 |
|---|---|
| `name` | `pinion`或`wheel` |
| `nodes_m` | 局部材料坐标中的实际节点，N×3，m |
| `triangles` | 局部节点索引，M×3整数；方向确定真实外法向 |
| `global_node_ids` | 长N，非负唯一整数，对应完整体FE自由度 |
| `total_node_count` | 完整体真实节点数，W行数为其三倍 |
| `outward_normals_local` | 可选M×3单位法向，须与实际三角平行 |

`frames`包含两成员的`*_rotation_rows`和`*_axis_world`。R为3×3正向正交旋转，世界点为`X = R x + t`；轴为世界单位向量。可选`*_translation_world_m`、`*_axis_origin_world_m`均为三维米坐标。不得用含义不明的整体非零`translation_m`代替各体平移。

每个patch包含整数`tooth_offset=k`、`grid`、`plane`和可选`requested_mask`。grid为`RectangularContactGrid`或其字典，数组形状为`(y_count, x_count)`，origin为首格中心，长度m。plane含`origin_world_m`、`tangent_rows`（2×3）和`normal_world`，须构成正向正交单位基。`requested_mask`须为同形状布尔数组；缺省请求整个给定网格。

可选`source_material_reference_frames`由调用者提供并原样保留，用于材料参考齿的Rk后处理。机械积分始终使用当前整体R0；两种姿态不可互换。

### 输出、单位和符号

| 字段 | 定义 |
|---|---|
| `records` | 按patch及行列排列的完整cell；含物理索引、g、面积、两体ℓ及原交片惰性视图 |
| `matrices[member]` | CSC稀疏W，`(3*total_node_count, 完整cell数)`，无量纲；全局行顺序为节点xyz |
| `gap_m` | 共同投影射线方向的从动面高度减主动面高度，m |
| `levers_m` | `ℓ = −c·[a×(X−O)]`，m；c为该体世界载荷系数，a为转轴，O为轴原点 |
| `area_plane_m2` | 完整请求cell的投影面积，m²；未知量为`f=p*A_plane`，N |
| `actual_area_m2` | 两体实际FE齿面面积，m² |
| `typed_patch_arrays` | 全部原piece、顶点、FE形函数和投影表的数值数组 |
| `operator_input_fingerprint` | 实际compact节点/三角/全局ID、请求集合、平面和当前姿态的SHA256数值身份；system及geometry一致 |

完整cell的线性剩余间隙为`g + W₁ᵀu₁ + W₂ᵀu₂ − θℓ₂`。系数为两体各自`c=±N/(N·n)`，主动体取负，从动体取正；N为实际FE外法向，n为共同投影法向。不夹取法向分母或FE重心坐标。

`source_valid_mask`始终表示原请求。未完整覆盖的cell保留`partial_records`、实际missing GeoJSON和原piece数据，其面积权重仍为真实片面积／完整cell面积；该cell不进入W。业务层必须拒绝或显式处理partial，不能当作已验证零压力背景。

内部保留`maximum_mapping_error_m`、`maximum_fe_ray_shape_closure_error_m`等兼容字段；本分支指真实FE射线与P1节点映射的闭合误差，不是CAD与FE离散误差，也不能改名为后者。

### 小型可运行示例

该平面算例演示参数化接口和归档，不是齿轮受载分析案例。依赖NumPy、SciPy和Shapely 2.x；可使用项目的`ltca`可选依赖环境。

```python
import tempfile
from pathlib import Path
import numpy as np
from bevel_gear.halfspace_contact import RectangularContactGrid
from bevel_gear.native_fe_operator import (
    compact_fe_boundary, build_native_fe_contact_operator,
    check_native_operator_pose, save_native_operator, load_native_operator,
)

xy = np.array([[-.01, -.01], [.01, -.01], [.01, .01], [-.01, .01]])
tri = np.array([[0, 1, 2], [0, 2, 3]])
compact = {}
for member, height in (("pinion", 0.), ("wheel", 1e-5)):
    mesh = dict(nodes_m=np.column_stack((xy, np.full(4, height))),
                boundary_triangles=tri if member == "pinion" else tri[:, ::-1])
    compact[member] = compact_fe_boundary(mesh, member)
frames = {m + "_rotation_rows": np.eye(3) for m in compact}
frames.update({m + "_axis_world": [0., 1., 0.] for m in compact})
patch = dict(tooth_offset=0,
    grid=RectangularContactGrid(4, 3, .002, .002, (-.003, -.002)),
    plane=dict(origin_world_m=[0., 0., 0.],
               tangent_rows=[[1., 0., 0.], [0., 1., 0.]], normal_world=[0., 0., 1.]))
system = build_native_fe_contact_operator(compact, [patch], frames, 0.)
u = {m: np.zeros((4, 3)) for m in compact}
proof = check_native_operator_pose(system, compact, u, 0.,
    np.array([r["gap_m"] for r in system["records"]]))
assert proof["passed"]
with tempfile.TemporaryDirectory() as tmp:
    save_native_operator(system, Path(tmp) / "operator")
    restored = load_native_operator(Path(tmp) / "operator")
    assert restored["operator_input_fingerprint"] == system["operator_input_fingerprint"]
print(len(system["records"]), proof["passed"])  # 12 True
```

归档不使用pickle，逐文件验证原字节SHA、列索引、geometry和生成姿态。文件路径仅允许归档根目录内的相对路径；绝对路径、父级引用和越界软链接均拒绝。无历史缓存的新调用不要求`native_cache_index_sha256`或外部mesh文件SHA；已有v1归档仍保留原cache及原文件身份的严格读取。该身份兼容性不免除业务层核验compact确实属于准备K时的实际FE网格。

### 原生FE单相位受载接触外迭代

`native_fe_contact.solve_native_fe_contact(prepared, geometry_template, *, compact_facets=None, ...)` 在已有整轮FE系统上直接运行真实边界接触。`geometry_template` 明确给出 `frames`（整体R0及世界轴）、`driving_angle_rad`、`reference_follower_angle_rad`、`pinion_torque_nm`、`pinion_to_wheel_ratio=z1/z2` 和 `patches`；每个patch包含 `tooth_offset`、`grid`、`plane`、布尔 `requested_mask`，以及需要材料投影时的 `source_material_reference_frames`。模板中额外的旧压力、间隙或柔度字段不参与新求解。

`compact_facets=None` 时从prepared的真实完整四面体边界及实际CAD工作面提取；显式传入时仍核对原全局节点、坐标和真实外边界三角形。每次API调用首个QP冷启动，同一相位的后续外步只映射前次承载物理格，不继承旧力或柔度列。每步重建当前φ的W/g0/力臂；出现partial或缺失原请求格直接拒绝，不能通过缩小有效掩码补零。

主要参数为 `maximum_outer_iterations=6`、`outer_step_tolerance_rad=1e-8`、`mean_geometry_tolerance_m=1e-10`、`maximum_piece_penetration_m=2e-7`、`maximum_working_set=2048` 和 `batch_cells=128`。外步、平均几何及片顶点门限可收紧；最终仍须核验真实u下全部原片顶点及实际加载φ重建。接触斑触窗或无效邻格、名义主动轴扭矩残差超限也不会接受。

返回字典的 `accepted` 表示上述单相位门控；`result` 为最终原QP与角度/几何诊断，`system` 为对应最终typed operator，另有 `geometry_outer_history`、`geometry_check`、`input_fingerprint`、`preparation_audit`、`rejection_reasons` 和 `limitations`。若新operator在进入QP前已不完整，`result=None`，不会把前次u/F配给新W。默认移除完整体响应数组但保留其数字身份；归档或实际材料域审核时应显式使用 `include_body_fields=True`。调用方保有prepared及因子，不因本API调用而被关闭。

生产内存结果以 `operator_input_fingerprint` 和 `actual_native_mesh_numeric_identity` 绑定当前算子及两体完整节点/单元/边界；不要求历史cache或文件路径。正式档案额外绑定实际mesh文件SHA、输入JSON、最终operator和原响应文件；最终operator是输出证据，不写回输入指纹形成循环。

参考φ是调用方明确的数值初值。首QP冷启动不表示任意CAD无载窗口都能承受大刚体外步；若迭代路径使有限工作面投影不完整，应重新设计真实覆盖窗口或初值。`total_loaded_te_rad` 为总受载TE，末步δ为当次参考姿态增量；没有求卸载FE网格相切时不提供纯弹性TE。本入口也不认证完整材料零域、压力空间网格或相位极值收敛。

`pinion_to_wheel_ratio` 必须与 `prepared.models` 声明的两轮实际齿数一致，原QP扭矩控制与TE统一使用实际 `z1/z2`。`frames` 与声明主动/从动角的运动学一致性由调用方明确负责；入口核验真实FE投影与位移几何，不据此宣称已反求或验证全部CAD运动学关系。
