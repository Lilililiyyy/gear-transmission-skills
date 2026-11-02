# 弧齿锥齿轮整轮有限元接触分析

`examples/run_native_fe_case.py`将实际制造齿面、带孔整轮、TET4 网格、有限齿面相交、多齿接触平衡及结果归档串联起来。同一齿体分解可服务于多个压力网格和主动轮相位，各状态分别求解接触与有限转角外迭代。

随附 `spiral_17_47_relieved_native_fe_case.json` 为 17/47 弧齿副的单相位、三档压力网格案例，采用 30 μm 法向齿顶修缘、−64 N·m 主动轴扭矩、210 GPa 弹性模量和 0.3 泊松比，两轮孔壁固定。输入中的长度、角度、压力和扭矩分别使用 m、rad、Pa 和 N·m。

## 分阶段运行

从项目根目录运行，输出目录应为尚未存在的新目录。

```bash
python -m pip install ".[ltca,plots,render3d]"

# 读取输入，列出制造、相位、压力网格和局部加密计划
python examples/run_native_fe_case.py --through plan --output output/fe_plan

# 生成实际闭实体、网格及接触几何覆盖报告
python examples/run_native_fe_case.py --through templates --output output/fe_geometry
```

`--through`依次支持`plan`、`cad`、`mesh`、`templates`和`solve`。`templates`阶段输出真实几何与完整面元覆盖情况；承载压力由`solve`阶段计算。

`--mesh-dir`可复用此前生成的两份网格 JSON。程序核对制造参数、齿坯、CAD 字节和局部加密区域的一致性；复用网格仍需与本次案例保持一致。

## 计算承载接触

大型算例提供可选 oneMKL PARDISO_64 求解器，当前验证运行库版本为 2025.2。程序使用调用者指定的数值库路径和缓存目录。项目只分发 Python 源码。

以下路径占位符由使用者替换，`output/fe_solve`和缓存目录均须为新目录：

```bash
python examples/run_native_fe_case.py --through solve --output output/fe_solve \
  --sparse-solver pardiso64-ooc --mkl-runtime /path/to/mkl_runtime \
  --ooc-cache-root /path/to/new_ssd_cache \
  --pinion-ooc-core-mib 4096 --wheel-ooc-core-mib 8192 \
  --factor-threads 4 --ooc-maximum-process-rss-gib 40 \
  --ooc-minimum-available-gib 10 --ooc-minimum-free-disk-gib 40
```

该命令采用 Bash 换行符；PowerShell 可合并为一行执行。core 参数约束每体盘外分解的内部工作内存，进程限额覆盖整个 Python 进程。资源门控在数值调用边界检查，实际峰值仍应由运行环境监测。缓存适合放在本地 SSD，计算结果保存在`--output`指定目录。

较小的输入可选择`--sparse-solver scipy`，使用 SciPy SuperLU。实际齿轮的精细三维网格具有较高内存需求，应按所选网格和求解器配置资源。

## 结果与压力图

`levels/<压力网格名称>/distributed_fe_summary.json`是各压力网格的结果入口。每个接受状态保存原压力、总受载传动误差、从动轴切线刚度、齿对分载、实际节点力与位移、原载荷算子以及模型档案。默认使用紧凑数值响应归档；`--include-body-fields`可额外在状态中保留完整节点字段。

完成求解后，可独立运行[接触结果审核与图源导出](native_fe_postprocessing.md)，复核原刚度方程、力和能量，并基于本态位移与受载角审核完整齿面。认证零压域采用蓝色，未认证域保留未知掩码。图源保存原 P0 求解量及守恒 P1 连续显示场。

## 多相位与精度研究

自定义案例可指定主动轮角列表、物理齿窗口、压力网格列表、材料、扭矩及求解门限。完整周期应覆盖一个主动齿距，包含重复端点，并使局部网格加密区域覆盖全部实际相位。区域文件中的相位来源与覆盖范围会被核对。

三档压力网格用于空间精度研究。周期最大压力需要完整周期的逐态接受结果和完整齿面证据，随后还需进行相位加密比较。程序分别保存原求解验收、独立响应审核、材料域覆盖、空间收敛与相位收敛状态。

当前完整工作面相交和有限转角外迭代案例入口适用于弧齿锥齿轮。直齿已有几何、工程根连接、标准强度、理想齿面 TCA 及既有 FE 接口；直齿刀具展成齿根与完整工作面接触链属于后续扩展。
