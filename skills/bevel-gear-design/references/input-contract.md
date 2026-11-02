# 输入、单位与阶段

协议沿用 `bevel_analysis_case_v1`，不接受圆柱齿轮Skill的case/study协议。字段完整签名以随内核保留的源码为准；不要静默丢弃无法识别的参数。引用制造文件通过 `load_bevel_case` 相对案例文件展开，结果 `inputs.json` 保存展开后的全部输入。

| 输入块 | 所需内容与阶段 |
|---|---|
| `schema`、`kind` | 固定schema；`straight`或`spiral`。后者要求`tooth_trace=circular_arc`，零螺旋角圆弧齿线仍不同于直齿 |
| `geometry` | 两轮正齿数、大轮外端节径、齿宽、轴交角、压力角、螺旋角、齿线、齿高和齿厚参数；零偏置 |
| `load` | `pinion_torque_n_m`正幅值、`pinion_speed_rpm`正转速；`ltca_torque_sign`为−1或+1，直齿统一入口仅+1 |
| `force_convention` | `flank=drive/coast`及`spiral_hand`，不得将加载符号、扫描方向与真实功率流混同 |
| `materials` | 强度及接触需要两轮材料；E为Pa、ν无量纲；2014名义强度还需要密度、热处理与应力数等显式材料字段 |
| `strength` | 齿距偏差、正反刀尖半径、应用系数、齿宽载荷系数、有效齿宽比例、修形等级、工作侧；具体必需项由原函数签名规定 |
| `stress_limits` | 接触／两轮齿根设计应力限、最低安全系数及`basis`；不得自行推定材料寿命、可靠度或润滑修正 |
| `manufacturing` | 弧齿各自`request`或`request_file`（二者互斥）、`blank`齿坯；名义／制造绑定在统一流程的geometry阶段也会检查 |
| `cad`、`cad_sampling` | 直齿两轮孔径等实体设置；弧齿制造面的采样设置。CAD起需要对应几何依赖 |
| `analysis` | TCA起需要严格递增、无重复的`driver_angles_rad`；保留物理齿索引`tooth_offsets`和载荷扫描幅值 |
| `finite_element`、`loaded_contact` | LTCA网格、节点／通道上限、接触匹配、几何及平衡容差；不能为了成功而暗中放宽门限 |

长度m，角度rad，力N，扭矩N·m，E／压力／应力Pa，密度kg/m³。rpm按字段保留。CAD内部mm，公开网格m、索引从0开始。外端端面模数与中点法向模数不能互换；齿数比与速比方向分别说明。

仅几何任务不应追问不参与计算的材料限值；但统一弧齿案例需要制造绑定。若用户仅有宏观几何，可直接调用 `bevel_gear.geometry.GeometryInput`、`calculate_geometry` 和 `forces.calculate_mesh_forces`，显式记录此路线尚未纳入Skill命令复算协议，不能伪造缺失的制造数据。

示例中的应力限、载荷系数、制造设置和网格仅属于演示。修改名义几何后要重核制造绑定，不能继续套用17/47机床输入。
