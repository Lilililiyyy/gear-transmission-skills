# 案例输入与复现

[中文首页](../README.md) · [English README](../README.en.md)

## 输入配置

README 的17/47案例共享宏观参数，名义评级和实际修缘接触分别保存制造输入。名义强度对标不代替修缘齿面的压力验证。

| 用途 | 输入 | 内容 |
|---|---|---|
| 名义弧齿设计 | [spiral_17_47_case.json](../examples/spiral_17_47_case.json) | 宏观几何、材料、载荷、标准强度、名义制造与TCA配置 |
| 名义数值参照 | [spiral_17_47_reference.json](../examples/spiral_17_47_reference.json) | 同量定义的MASTA几何和标准强度保存值 |
| 实际修缘实体 | [spiral_17_47_relieved_contact_design.json](../examples/spiral_17_47_relieved_contact_design.json) | 双轮刀具／机床／齿坯、材料及30 μm齿顶修缘 |
| 局部半解析 | [spiral_17_47_relieved_surface_case.json](../examples/spiral_17_47_relieved_surface_case.json) | 实际修缘输入、129相位及接触求解网格 |
| 整轮FE周期 | [spiral_17_47_relieved_native_fe_cycle9_case.json](../examples/spiral_17_47_relieved_native_fe_cycle9_case.json) | 同一修缘输入、网格、九个姿态、压力格及求解设置 |
| 直齿设计 | [straight_16_32_case.json](../examples/straight_16_32_case.json) | 16/32直齿几何、材料、载荷及强度输入 |

名义对标使用152.4 mm大轮刀盘；修缘接触使用146.304 mm大轮刀盘、2.36 mm刀尖点宽及相应机床设置。两轮的30 μm修缘从归一化有效齿廓位置v=0.65开始，采用C²连续五次材料去除函数；v=0为齿根、v=1为齿顶。设计应力限1500/450 MPa为案例显式输入。

## 分阶段运行

从项目根目录运行，输出写入`output/`。名义几何和标准强度只需基础包，绘图使用`.[plots]`；CAD、TCA及接触安装对应的可选依赖：

```sh
python -m pip install -e ".[ltca,plots,render3d]"
```

名义弧齿对标与直齿强度：

```sh
python examples/run_spiral_workflow.py examples/spiral_17_47_case.json --reference examples/spiral_17_47_reference.json --through strength --output output/strength --plots
python examples/run_spiral_workflow.py examples/straight_16_32_case.json --through strength --output output/straight
```

实际修缘齿面与闭实体，导出STEP/BREP及中英三维图：

```sh
python examples/build_contact_geometry.py examples/spiral_17_47_relieved_contact_design.json --output output/contact_geometry --plots
```

名义制造配置的TCA与实际修缘配置的局部半解析接触分别运行：

```sh
python examples/run_spiral_workflow.py examples/spiral_17_47_case.json --through tca --output output/nominal_tca
python examples/run_surface_contact.py examples/spiral_17_47_relieved_surface_case.json --output output/relieved_halfspace
```

整轮FE可先仅检查输入并生成计划，不创建CAD、网格或刚度矩阵：

```sh
python examples/run_native_fe_case.py examples/spiral_17_47_relieved_native_fe_cycle9_case.json --through plan --output output/fe_plan
```

完整求解的数值引擎、内存、缓存及复用设置见[九相位FE案例](native_cycle9.md)和[FE入口说明](native_fe_case.md)。整轮接触比名义评级和图表重绘需要更多资源。

## 重绘随附结果

安装`.[plots]`即可从SI源重绘压力及周期响应，无需CAD或重新求解接触：

```sh
python examples/replot_contact_composite.py --output output/contact --lang both
python examples/replot_contact_composite.py --view raw --stem 17_contact_composite_raw --output output/contact_raw --lang both
python examples/replot_native_response.py --output output/response --lang both
```

图表输出在`zh/`和`en/`，数值输出在`data/`；绘图生成PNG、SVG、PDF。仓库保留精选PNG及可重绘数值，避免重复存放大型矢量图。

默认压力图的FE下排采用50 μm物理核宽的法向合力守恒恢复，以及九个求解相位到129个显示相位的随轨迹插值。原始视图采用`--view raw`。方法、原峰值、显示峰值、参数与限制见[压力显示说明](contact_pressure_display.md)；瞬时斑的独立Hertz尺度检查见[尺度数据](../examples/results/relieved_contact/contact_scale_reference.json)。

原压力与TE／刚度CSV在[案例结果目录](../examples/results/relieved_contact/)，可独立读取。标准评级应力、LTCA局部压力、总TE与局部半空间闭合量分别定义，不跨模型互换。

## English

The 17/47 example shares macrogeometry across nominal rating and relieved contact analysis, with separate manufacturing inputs. The table links the nominal case and MASTA reference, actual relieved solids, 129-phase local semi-analytical case, nine-phase full-gear FE case, and 16/32 straight-gear design. Nominal strength comparisons do not validate spatial pressure on the relieved geometry.

The nominal wheel cutter diameter is 152.4 mm. The contact configuration uses a 146.304 mm cutter, 2.36 mm blade point width and corresponding machine settings. Both members have 30 μm normal tip relief with a C² quintic law starting at normalized profile position v=0.65, where v=0 is the root and v=1 the tip. The 1500/450 MPa design stress limits are explicit case inputs.

Run the commands above from the project root. The basic package supports nominal geometry and strength; `.[plots]` adds figures; `.[ltca,plots,render3d]` adds optional surface, CAD and contact dependencies. The commands separately compute nominal strength, actual relieved STEP/BREP solids and bilingual views, nominal-configuration TCA, and actual relieved local semi-analytical contact. These are distinct manufacturing configurations.

The FE `--through plan` command checks inputs without building CAD, a mesh or a stiffness matrix. Full solve commands and resource settings are in the [nine-phase case](native_cycle9.md) and [FE runner](native_fe_case.md). Full-gear contact requires more resources than nominal rating or figure reproduction.

The three replot commands reproduce the engineering pressure view, raw pressure view and FE TE/stiffness response without CAD or a new contact solve. Bilingual figures appear in `zh/` and `en/`, numerical outputs in `data/`; PNG, SVG and PDF are generated. The repository retains selected PNGs and numerical sources rather than duplicate large vector figures.

Default FE rendering uses a 50 μm physical recovery length while preserving each original tooth's reference-plane normal resultant, then interpolates along contact tracks from nine solved phases to 129 display phases. The raw view applies neither step. See [Pressure-display methods](contact_pressure_display.md) for raw/displayed peaks, settings and limitations, and [contact-scale data](../examples/results/relieved_contact/contact_scale_reference.json) for an independent Hertz scale check.

Original pressure maps and TE/stiffness CSVs are in the [example results](../examples/results/relieved_contact/). Standard rating stress, local LTCA pressure, total TE and local half-space closure are distinct quantities, not substitutes across models.
