# 结果解释与复算

新版图表及完整案例目录见 [新版完整工作流](new-workflow.md)。常规 `run` 使用 `cases/`；完整 `case` 使用 `ltca/` 和可选 `sensitivity/`，实际目录以 `run.json.cases` 为准。`max_endpoint_pressure_pa` 是独立端点最大压力，不能与 `max_pressure_pa` 的中心定义互换。

先读根目录 `run.json`：只有 `status=completed`、`analysis_completed=true`才表示本次所有任务及请求的后处理完成。失败保留 `error`与已有文件，不能将文件存在当作成功。随后逐个查看 `cases`指向的原内核输出。

|字段|定义|
|---|---|
|`te_peak_to_peak_m`|原相位采样下法向接近量的峰峰值，m|
|`angular_te_peak_to_peak_rad`|以小轮虚功一致法向力臂换算的角度TE，rad|
|`max_pressure_pa`|原离散单元最大Hertz压力，Pa|
|`mean_mesh_stiffness_n_per_m`|各相位当前活动接触分支的平均法向切线刚度，N/m|
|`max_force_residual_n`|内核载荷平衡残差最大值，N|
|`max_gap_residual_m`|活动间隙误差／最大穿透的最大值，m|
|`all_converged`|所有相位满足当前求解器退出条件|

展示可将TE乘1e6换为μm、压力乘1e-6换为MPa、刚度乘1e-6换为N/μm。原数组保持SI单位；不要将法向刚度解释为N·m/rad。

`arrays.npz`包含周期数组与变长单元展平数组，不含pickle；`cell_offsets[k]:cell_offsets[k+1]`对应第k相位。`active_mask`、齿对编号、接触位置、原力、压力、间隙均保留。图形不改写原峰压与响应；瞬时Hertz带重建和周期扫掠包络用途不同，包络不能积分为瞬时合力。

`study.json`保持原 `compare_designs` 的目标、逐工况结果和非支配集合。压力约束只来自显式输入；无约束不等于强度安全。修形可能降低峰压同时增大TE波动，报告应解释权衡。

复算命令必须指定新的目录。`replay`先核对原目录所有文件哈希、代码身份和数值环境，再从原请求生成新的求解，逐数组比较。整数索引、偏移、布尔掩码要求一致；浮点容差、最大差及逐数组结论保存在 `replay.json`，并核对每个实际输入、全部汇总指标和方案比较结果。

`replay.json.passed=true`仅说明给定同环境的数值复现。它不意味着网格／连续相位收敛、原组件响应精度或物理试验验证。本次未运行的验证在 `run.json.validation`中保持未验证。
