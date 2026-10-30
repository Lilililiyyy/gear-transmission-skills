# 请求、单位与默认值

请求是UTF-8 JSON，根字段为 `schema`、`case`、可选 `study` 和 `assumptions`。`schema`固定为 `gear_ltca_request_v1`。未知字段、重复JSON键、非有限数、重复指定同一物理量均报错。

`case`沿用原内核的命名配置。`pair.pinion` 与 `pair.wheel`至少包含 `teeth`、`normal_module_m`、`face_width_m`。常用字段包括 `pressure_angle_rad`、带旋向的 `helix_angle_rad`、`profile_shift`、`internal` 和实际根／顶圆及有效形径。`pair`另可给 `center_distance_m`、`axial_offset_m`。

两轮材料分别用 `material_pinion`、`material_wheel`，含 `youngs_modulus_pa`、`poisson_ratio`。单工况给 `torque_n_m`、`misalignment_m`、`modification_pinion`、`modification_wheel`。参数修形支持 `lead_crown_m`、`profile_crown_m`、`tip_relief_m`等原数据类字段；原测量曲线／曲面／组合修形协议也保留。

离散设置放在 `case.phases`（每周期相位数）和 `case.strips`（每条接触线条带数），均为至少2的整数。`case.phase_offset`以基节周期计，默认0；不要用响应曲线拟合选择这个值。

|输入字段后缀|含义|程序换算|
|---|---|---|
|`_m`|m|保持|
|`_mm`|mm|×1e-3 → `_m`|
|`_um`|μm|×1e-6 → `_m`|
|`_rad`|rad|保持|
|`_deg`|度|×π/180 → `_rad`|
|`_pa`|Pa|保持|
|`_mpa`|MPa|×1e6 → `_pa`|
|`_gpa`|GPa|×1e9 → `_pa`|

例如 `normal_module_mm=2.5` → `normal_module_m=0.0025`；`lead_crown_um=8` → `lead_crown_m=8e-6`。别名只改变单位，不能使内核不支持的字段变为支持。扭矩只用N·m，不能把N·mm直接填入 `torque_n_m`。

默认值由保留内核提供：压力角20°、直齿、零变位、外齿、零轴向错位；中心距按现有几何规则求得；材料E=210 GPa、ν=0.3；无修形与零失准；32相位、21条带，`elastic_mode=coupled`、`local_contact=true`，起采相位0。其他默认字段见实际输出的 `input.json`。这些默认值不表示实际设计已知条件；输入不足时先询问，或在用户接受后列入 `assumptions`。`run.json.defaulted_fields`保留缺省路径。

比较任务给出：

```json
{
  "study": {
    "designs": [
      {"name": "无鼓形", "modification_pinion": {"lead_crown_um": 0}},
      {"name": "8μm鼓形", "modification_pinion": {"lead_crown_um": 8}}
    ],
    "operating_points": [
      {"name": "80Nm", "torque_n_m": 80, "misalignment_um": 12}
    ]
  }
}
```

这个片段须与完整 `schema` 和 `case`组合。比较任务的 `case`不能再给扭矩、失准或小轮修形；大轮修形仍可在 `case`给定。候选名、工况名分别唯一，名称仅作标签。可给显式 `pressure_limit_mpa`／`pressure_limit_pa`，不提供时只按TE和压力做非支配比较。

`elastic_mode`保留 `coupled`、`independent`、`iso`、`iso_love`。ISO相关输入、局部参照深度及齿顶延伸仍按原内核协议；不得擅自补造名义载荷因子。Python外部柔度和间隙回调不能由纯JSON重现，因此本Skill首版入口不接受回调代码。
