# 接触压力显示与原始解

四联图固定顺序为：半解析瞬时、半解析周期最大、FE瞬时、FE周期最大。每幅独立使用从0开始的线性MPa色标，蓝→青→绿→黄→红；颜色不能跨图直接比较，应读取各自色标。横轴w_f为从外端向内端的齿宽位置，纵轴h为沿节锥法向朝齿顶的高度，绘图单位mm；p为法向接触压力。坐标投影不是等面积展开。

## 两种视图

`--view raw`直接在各自原P1三角域内采样压力；不做压力恢复或相位插值。四份`*_triangle_pressure.json.gz`保留原SI标量、三角连接、逐场状态摘要和未知域。随附`comparison_case.json`描述的也是这一原始采样结果。

默认`--view engineering`的半解析上排保持原场，FE下排使用以下独立后处理。原始压力、受载角、TE、刚度及求解接受标志不改变。其目的为低成本展示有限相位案例；后处理场不是新增接触平衡解。

## FE空间压力恢复

在原参考平面面元格上，对非负P0压力p做物理长度核的高斯恢复Gσ，并按每个齿、每个原相位的法向合力F归一化：

```text
F = Σ p_i A
q = Gσ(p)，有效域之外及边缘一格置零
p_rec = q F / (Σ q_i A)
```

A为原参考平面面元面积，单位m²；p和p_rec为Pa，F为N。恢复核使用实际格距，数组布局为`[y,x]`，不对RGB图片模糊。边缘一格零支承用于构造完整中心P1三角域；截边会影响局部形状和中心，位移量一并记录。不能据此声称保留原压力力矩或局部接触互补关系。

本例σ=50 μm，为100 μm背景表面控制尺寸的一半；该尺寸是网格生成控制量，不保证每条实际边长。原压力格距593.023×50 μm。σ是可配置的工程显示参数，不是收敛判据；`--smoothing-length-um`可显式修改。结果记录每个原场的原／恢复峰值、参考平面合力误差及中心漂移。随附案例最大相对合力误差约2.22×10⁻¹⁶。

恢复后在实际材料坐标上使用原格中心构建P1三角场，参考积分权重仍为A/2，逐场核验`∫p_rec dA=F`。完整工作面外不外推，原未解析区仍单独遮罩。显示格或子午投影面积积分不能替代原参考平面积分。

## 接触轨迹随相位插值

设同一物理齿在相邻真实相位的恢复压力为p_i、p_j，合力加权材料坐标中心之差为Δc。对相对相位t∈[0,1]定义：

```text
p(t,x) = (1−t) p_i(x−tΔc) + t p_j(x+(1−t)Δc)
```

先随接触中心移动，再对标量压力作凸组合，最后与所有原始恢复场一起取逐点最大值。这里不借用半解析结果、不拟合椭圆、不截顶或人为匹配两种方法的峰压。只连接同一齿、相邻原相位且两端均承载的场；不跨缺失或零载相位。实现使用非负双线性平移采样，插值包络不会产生高于端点显示场的新峰值。

本例原9个求解相位覆盖一个主动轮齿距，间隔约2.647°；每个区间细分16份，得到129个显示相位，新增FE求解次数为0。此公式假设邻相位间接触区随中心平稳移动，不能恢复未采样的载荷突变或接触拓扑变化，也不提供连续相位最大值的上界。插值场没有重新满足刚度方程、扭矩平衡和互补约束；其积分不能当作某一新增相位的真实载荷。

## 本例的峰压和条带

原FE周期峰值2214.8 MPa位于主动相位约+2.647°、物理齿偏移−1。其邻域出现高压／近零压交错；压力格细于部分FE背景表面离散，存在局部离散振荡，峰压尚无网格收敛证明。单位关系`f_i=p_i A`、扭矩及能量检查通过，只能说明相应平衡和数据一致性，不能证明峰压正确。有限齿体柔度与半空间柔度的物理差别也不能单独解释这一局部峰值。

原周期图的分离条带主要来自9相位采样，另有局部空间离散纹理。加密显示格不能增加真实相位信息。主图分别采用空间压力恢复和随轨迹相位插值，区分处理这两个问题。161×641显示格上的四个峰值依次约为1078、1082、975、1118 MPa；FE显示值不代替原1069.1／2214.8 MPa，也不用于强度验收。

## 复现及接口

```sh
python -m pip install ".[plots]"
python examples/replot_contact_composite.py --output output/contact --lang both
python examples/replot_contact_composite.py --view raw --stem 17_contact_composite_raw --output output/raw --lang both
```

工程显示额外读取`native_pressure_recovery_sources.json.gz`，包含原FE P0格、有效掩码、材料坐标和状态摘要，不依赖CAD、刚度矩阵或厂商运行库。恢复源与原周期数据须有相同相位和状态身份，否则拒绝。输出`data/*_SI.npz`按四图顺序保存SI显示格，`data/*_display.json`保存原峰值、显示峰值、色标、参数及诊断；这两项是显示产物，不是原接触结果。

程序接口为`contact_display.prepare_pressure_comparison`和`plot_pressure_comparison`；数值步骤在`contact_recovery`中，可独立调用。原精确三角等值域导出接口保留在`contact_composite`。

## English

The fixed panel order is semi-analytical instantaneous pressure, semi-analytical cycle maximum, FE instantaneous pressure, and FE cycle maximum. Each panel has an independent linear MPa colorbar, starting at zero blue and ending in red. Read each scale when comparing colors. Face-width coordinate w_f runs from heel to toe; h points toward the tip along the pitch-cone normal. Both plot axes use mm; p is normal pressure. The projection is not area preserving.

`--view raw` samples the original P1 triangle fields without pressure recovery or phase interpolation. The four original SI triangle archives and `comparison_case.json` continue to describe the unprocessed solutions. The default engineering view leaves the semi-analytical fields unchanged and postprocesses the FE lower row. Original pressures, angles, TE, stiffness and acceptance flags remain unchanged.

Spatial recovery applies a Gaussian kernel with physical length σ to the original nonnegative P0 grid, then normalizes each tooth/state to preserve its reference-plane normal resultant `F=Σp_i A`. A one-cell zero boundary ensures complete P1 support; clipping can shift the pressure center, which is recorded. This does not preserve every original moment or contact constraint. The example uses σ=50 μm, half the nominal background surface control size. It is a configurable display choice, not a convergence criterion. Actual grid spacings are 593.023×50 μm. The maximum relative resultant error is approximately 2.22×10⁻¹⁶. Actual material coordinates and reference triangle weights A/2 are retained; projected pixel integrals are distinct from these force integrals.

Temporal interpolation uses the formula above: transport two neighboring fields with their force-weighted contact centers, then form a convex pressure combination. Only the same physical tooth at adjacent solved, loaded phases is connected. Missing or unloaded phases are not bridged. There is no ellipse fit, peak clipping, or transfer of semi-analytical pressures to FE. The envelope includes all original recovered samples and introduces no peak above the endpoint display fields.

Nine solved FE phases, separated by approximately 2.647°, become 129 display phases through 16 subdivisions per interval, with zero new FE solves. This assumes smoothly moving contact between samples. It cannot reconstruct unknown load jumps or topology changes, certify a continuous-phase maximum, or verify equilibrium at interpolated positions. Interpolated pressure integrals are not newly solved physical loads.

The original FE peak, 2214.8 MPa, occurs near pinion phase +2.647° on tooth offset −1, with alternating high and near-zero neighboring values. The pressure grid is finer than parts of the background FE surface mesh; local pressure has not been shown mesh-converged. Cell-force units, torque and energy checks do not establish peak-pressure accuracy. Different FE and half-space compliance scopes alone do not explain this local peak. Separate bands in the raw cycle map mainly arise from sparse phase sampling; spatial pressure oscillations are a different issue. Increasing display resolution alone adds no solved phases.

On the 161×641 display grid the four peaks are approximately 1078, 1082, 975 and 1118 MPa. The FE display estimates do not replace the original instantaneous/cycle peaks of 1069.1/2214.8 MPa and are not intended for strength acceptance.

The commands above reproduce both bilingual views without CAD or a new FE solve. Engineering rendering additionally reads `native_pressure_recovery_sources.json.gz`, containing original FE grid pressures, masks, material coordinates and state identities. The SI raster archive and display JSON explicitly distinguish raw and displayed peaks, settings, phase counts and diagnostics. They are display artifacts, not contact-solver results. The public functions are `contact_display.prepare_pressure_comparison` and `plot_pressure_comparison`; numerical recovery functions reside in `contact_recovery`.
