# 17／47弧齿锥齿轮的九相位整轮FE案例

案例 `examples/spiral_17_47_relieved_native_fe_cycle9_case.json` 给出一个齿距内包含两端点的九个离散姿态，直接用于重建随附有限元采样结果。主动轮扭矩为−64 N·m，材料为E=210 GPa、ν=0.3，两轮采用30 μm齿顶修缘。显式从动角只作为外迭代初值，实际受载角及载荷分配重新求解。

本版案例使用完整压力面元分支。边缘部分面元的接触约束与积分求解列入下一版本；当前入口在QP前拒绝部分面元，不将其删除、缩小请求范围或当作零压力。这里的九个指定姿态不包含该分支，增加其他姿态或改变工况仍须通过原覆盖与接触检查。

## 准备和求解

先读取输入并生成计划；该阶段不创建CAD、网格、刚度矩阵或压力解：

```sh
python examples/run_native_fe_case.py examples/spiral_17_47_relieved_native_fe_cycle9_case.json --output output/cycle9-plan --through plan
```

从制造参数生成双轮实体和体网格：

```sh
python examples/run_native_fe_case.py examples/spiral_17_47_relieved_native_fe_cycle9_case.json --output output/cycle9-mesh --through mesh
```

通过 `--mesh-dir` 精确复用本次网格，并指定独立新缓存目录及实际数值引擎路径：

```sh
python examples/run_native_fe_case.py examples/spiral_17_47_relieved_native_fe_cycle9_case.json --output output/cycle9-fe --through solve --mesh-dir output/cycle9-mesh/meshes --sparse-solver pardiso64-ooc --mkl-runtime <mkl_runtime> --factor-threads 4 --ooc-cache-root <new_cache_directory> --pinion-ooc-core-mib 4096 --wheel-ooc-core-mib 8192 --ooc-maximum-process-rss-gib 40 --ooc-minimum-free-disk-gib 40 --ooc-minimum-available-gib 6
```

数值引擎须满足接口声明的版本和资源要求。两体弹性因子在这些姿态之间复用，各态的真实接触几何和载荷分配独立求解。输出包含逐相位输入指纹、实际网格与算子身份、原响应及最终收据。未计算或未接受的姿态不会补成零压力。接触专用内存策略不输出FE单元应力，标准强度校核独立进行。

网格必须保留原JSON、metadata及绑定的制造、齿坯和CAD文件。CAD／网格引擎版本可能影响重建字节；身份不匹配时会明确拒绝。

## 网格、压力格和工程容差

主体场控制尺寸为1 mm，背景接触轨迹为100 μm。附加50 μm局部轨迹仅由主动角0 rad生成，区域文件如实保留 `source_driver_angles_rad: [0.0]`。显式 `refinement_phase_policy: allow_global_background` 允许其他八态使用该非均匀背景网格；默认 `require_local_sources` 策略仍要求每个所选相位具有局部来源。这些Gmsh场尺寸不保证每个单元都小于对应尺寸，也不构成全周期细化或空间误差证明。

压力格为43×83，格距593.023×50 μm，保持25.5×4.15 mm完整物理窗口。CAD未知查询、制造／齿坯身份、真实FE覆盖、平均几何、力平衡、Ku和接触求解门均保持。

接口默认最大点穿透容差为0.2 μm；显式值须为正且不超过1 μm。本工程案例设置 `solver.maximum_piece_penetration_m = 1e-6`。P0压力采用面元常值，接触平均间隙约束不要求片内所有顶点间隙精确为零；1 μm是法向点间隙验收容差，不是压力误差指标。工程显示许可与原始验收状态分别保存，不改写原记录。

九个显示状态中，最差顶点法向间隙约为-0.573 μm，处于显式1 μm工程容差内。需要瞬时压力细节窗时，可在重绘命令后追加 `--instantaneous-zoom`。

## 索引、图表和结果含义

本案例的 `--phase-index` 取0至8，α=0对应索引4。随附图源的原采样编号 `source_indices = [0,16,32,48,64,80,96,112,128]` 与本案例九个位置按实际主动角依次对应；它们是数据标识，不是本入口的相位位置参数。

瞬时接触图取α=0的参考齿对。原始周期图为九个真实姿态在材料坐标中的采样最大压力。README主图下排另外采用逐齿法向合力守恒恢复，以及9求解相位到129显示相位的随接触中心插值；这是工程显示估计，不新增FE平衡解或声称峰压已收敛。原始视图以`--view raw`单独重绘。两种显示和数值定义见[压力显示说明](contact_pressure_display.md)。总受载TE定义为实际从动角加 `(z1/z2)×主动角`，不是纯弹性TE；从动轴切线刚度是在本态固定几何及承载活动集下的线性化量，SI单位N·m/rad。

仅重绘随附结果，无需CAD或再次求解：

```sh
python examples/replot_contact_composite.py --output output/contact-replot --lang both
python examples/replot_native_response.py --output output/native-response --lang both
```

## English

`examples/spiral_17_47_relieved_native_fe_cycle9_case.json` defines nine discrete positions, including both ends of one tooth pitch, for the 17/47 spiral-bevel pair. The pinion torque is −64 N·m, E = 210 GPa, ν = 0.3, and both members have 30 μm tip relief. The explicit follower angles are numerical initial guesses; the loaded geometry and contact forces are solved at every position.

This release uses complete pressure cells. Contact constraints and integration for partially covered boundary cells are planned for the next version; the current solver rejects them before QP rather than discarding them or treating them as unloaded. The nine specified positions do not require that branch. Other positions and loads must independently satisfy the same coverage and contact checks.

The commands above generate a plan, build the CAD and mesh, and solve using the retained mesh files with explicit numerical-engine resource limits. The bulk size control is 1 mm, the background-track control is 100 μm, and the additional 50 μm local tracks originate only at driver angle zero. Their source list remains `[0.0]`. The explicit `allow_global_background` policy permits the other eight positions to use this nonuniform mesh; it certifies neither fine coverage throughout the cycle nor a spatial error bound. CAD identity, unknown-query, native FE coverage, mean-geometry, equilibrium, Ku and contact checks remain active. Rebuilt CAD or mesh bytes can differ across engine versions; incompatible identities are rejected.

The 43×83 pressure grid has spacings of 593.023 and 50 μm and retains the full 25.5×4.15 mm physical window. The API default point-penetration tolerance remains 0.2 μm; explicit values must be positive and no greater than 1 μm. This engineering case explicitly uses 1 μm. P0 cell-average contact constraints do not require every vertex gap to vanish. This tolerance is not a pressure-error claim. Engineering display admission is recorded separately from the original acceptance status.

The worst vertex normal gap among the nine displayed states is approximately -0.573 μm, within the explicit 1 μm engineering criterion. Append `--instantaneous-zoom` to the pressure replot command for equal-aspect instantaneous details.

The case uses position indices 0–8, with α=0 at index 4. The source-data identifiers `[0,16,32,48,64,80,96,112,128]` map to these positions in the same driver-angle order; they are not the CLI position indices. The raw cycle map is the maximum over nine computed material-coordinate states. The README lower row additionally uses per-tooth normal-force-preserving pressure recovery and contact-center interpolation from nine solved to 129 display phases. These are engineering estimates, with no new FE equilibrium solves or peak-pressure convergence claim. Use `--view raw` for unprocessed results; see [Pressure-display methods](contact_pressure_display.md). Total loaded TE is the actual follower angle plus `(z1/z2)×driver angle`, with no elastic-only split. Follower tangent stiffness is linearized at the fixed geometry and active set and has SI units N·m/rad. The final two commands redraw the supplied figures directly from their source data.
