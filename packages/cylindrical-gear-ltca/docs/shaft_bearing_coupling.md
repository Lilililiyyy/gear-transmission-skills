# 齿轮、阶梯轴与线性轴承准静态耦合

此功能为新写的系统柔度扩展，计算一对平行轴圆柱齿轮及其两根轴。每根轴由实心圆截面阶梯轴段和安装在刚性箱体上的线性轴承构成。齿轮轮体视为刚体，轮齿柔度和局部压陷仍按原 LTCA 模型计算。

## 输入与复算

```bash
python -B -m gear_contact examples/shaft_bearing_case.json --output output/shaft_bearing
python -B -m gear_contact output/shaft_bearing/input.json --output output/shaft_bearing_replay
```

JSON 的 `shaft_system` 包含 `pinion` 和 `wheel` 两根轴。Python 可直接用 `ShaftSystem`、`Shaft`、`ShaftSegment` 与 `LinearBearing` 构造同样的输入，再把 `shaft_system` 传给 `analyze_pair`。用 `solve_shaft(shaft, gear_wrench)` 可单独校核轴，截面载荷按 `Fx, Fy, Fz, Mx, My, Mz` 输入，单位为 N 和 N·m。

- 轴段沿共同的 +z 方向从轴首排列，节点 0 为轴首，节点坐标由 `segments` 的 `length_m` 累加。每段输入 `diameter_m`，可选 `material`。齿轮和轴承位置如处在轴段内部，应在那里分段。
- `gear_node` 是齿宽中心所在的安装截面，`torque_node` 是输入或输出端的角度参考截面。两者为零基节点编号。轴系弹性问题只约束 `torque_node` 的绕 z 转角，由该端承受反扭矩。
- `bearings` 按输入顺序记录，每项给出 `node` 及三向平动刚度 `stiffness_n_per_m=[kx, ky, kz]`。零表示该方向不提供支承。每个径向方向需至少两个支承节点，每根轴至少有一处正轴向刚度。支承为线性弹簧，箱体端位移为零。
- 轴系耦合支持 `independent` 和 `coupled` 齿梁模式。`misalignment_m` 和 `boundary_gap` 仍表示独立给定的几何偏差，已经由轴系算出的位移不应再次输入。ISO/Love、齿顶延伸及外部结构柔度替换接口暂不与本功能组合。

## 力、位移与接触协调

轴段采用两节点 Timoshenko 梁，考虑双向弯曲和剪切，以及轴向拉压、扭转。节点位移自由度为 `ux, uy, uz, θx, θy, θz`。实心圆截面剪切修正系数取 `6(1+ν)/(7+6ν)`。轴承平动弹簧组装进轴刚度，保留齿轮截面六个自由度的柔度。公式来源见 [Timoshenko 梁单元](https://mooseframework.inl.gov/moose/modules/solid_mechanics/C0TimoshenkoBeam.html) 与 [Cowper 剪切系数](https://doi.org/10.1115/1.3625046)。

以小轮齿宽中心为坐标原点。外啮合的大轮中心在 +x 一侧，内啮合的大轮中心在 −x 一侧。轴向坐标与原接触网格一致。按当前受载面约定，`αw` 为工作压力角，`βb` 为小轮带符号基圆螺旋角，作用于大轮的单位法向力为 `(cosβb·sinαw, cosβb·cosαw, −sinβb)`，小轮承受反向力。条带展长和齿宽坐标给出作用点，力矩相对各自齿轮安装截面求取。小轮的绕 z 啮合力矩与输入扭矩平衡，大轮的绕 z 啮合力矩在内外啮合中符号相反。

令 `B` 的每行为单位条带力在两轮安装截面的六分量载荷。对同一组条带力 `f`，用 `Bᵀ f` 得到轴系载荷，用 `B u` 将截面弹性位移投影成法向间隙增量。以 `K` 表示施加支承及角度参考后的轴系刚度，联立方程的新增项为 `C_shaft = B K⁻¹ Bᵀ`。将它加到原齿梁柔度后，用原单边接触约束求 `f`，因此同一相位的轴系变形和载荷分布相互一致。

这里采用小位移线性投影，候选接触线仍按名义几何生成。滚动体载荷、轴承预紧与游隙、轮体及箱体柔度属于后续模型范围。

## 结果与验证

- `diagnostics.json` 的各相位 `elasticity[].shaft_bearing` 保存两个轴的节点位移与转角、齿轮截面载荷、轴承反力、反扭矩及力和力矩平衡残差。
- `arrays.npz` 的 `shaft_normal_displacement_m` 按原 `cell_offsets` 展开，正值表示接触载荷产生的法向间隙增量。初始 `gaps_m` 仍保存原几何间隙。
- 启用轴系时，法向接近量及等效小轮角度 TE 参考两根轴所指定的扭矩端，包含轴系贡献。`mesh_stiffness_n_per_m` 表示该条件下总柔度的法向切线刚度。它随载荷和啮合位置变化。
- 解析验收用非对称简支阶梯轴，比较两段弯曲和剪切挠度、支承沉降、轴向变形、扭转和反力。刚性极限验收同步提高轴与轴承刚度，对照原齿面条带力、活动集合、接近量和切线刚度。另校核内外啮合扭矩符号、螺旋轴向力、平衡残差及导出输入复算。验收为解析与数值验证，尚无整机实测对照。
