# 分布牵引 FE 的复算输入

`examples/run_distributed_fe_contact.py` 从显式制造输入生成真实整轮、体网格和分布 P0 受载接触。网格及输入身份会进入结果指纹；原 P0 压力和用于重绘的实际几何逐相位保存。离散接受、网格研究、完整相位覆盖和完整面域证明分别处理。

## 压力网格方向

案例的 `analysis` 可增加：

```json
{
  "plane_alignment": "gap_hessian",
  "alignment_step_m": 0.0001,
  "alignment_fallback_to_surface_basis": true
}
```

省略时 `plane_alignment="none"`，沿用原曲面投影基。`gap_hessian` 从实际双齿面间隙的两档差分求主曲率方向，仅选择正交网格方向；每个原中心、角点和 Gauss 点仍查询真实有限曲面，不用二次间隙或压力拟合替代。方向保存于 `build_parameters`，同相位的受载外迭代保持不变。

`alignment_fallback_to_surface_basis` 默认 False。显式启用时，主轴差分在有限面外失败的齿仍保留原投影基和全部原候选，原因写入 `state.grid_alignment`。这不是外齿排除或零压力证明。自动主轴不能与显式 `plane_basis_angles_rad` 同时指定。

## 附加局部网格区域

保留原全轨迹细化带后，可在 `mesh` 指定相对案例文件的路径：

```json
{
  "refinement_regions_file": "controls/refinement_regions.json"
}
```

区域文件为 `bevel_refinement_regions_v1`，包含以下必填项：

| 字段 | 含义 |
| --- | --- |
| `units` | 必须是 `"m"`；折线与尺寸均为 SI 长度 |
| `coordinate_system` | 必须是 `"gear_body_local"`；分别使用两轮实际本体坐标 |
| `source_driver_angles_rad` | 有限、不重复的源相位列表，单位 rad |
| `manufacturing` | `pinion`、`wheel` 各自完整标准化的 `request` 和 `blank` |
| `regions` | `pinion`、`wheel` 各自的附加区域列表，可一轮为空 |

每个区域使用以下四个字段：

```python
region = {
    "refinement_polylines_m": actual_body_local_polylines,
    "refined_size_m": 75e-6,
    "refinement_core_radius_m": 300e-6,
    "refinement_radius_m": 1.2e-3,
}
```

折线必须来自真实曲面材料点，至少有两个有限三维点，相邻点不能重合。尺寸须满足 `fine < bulk` 且 `fine <= core < radius`。原全轨迹场和所有附加 Distance/Threshold 场由 Gmsh `Min` 组合，不能以附加场放粗其他区域。上述数值仅说明接口单位；具体尺寸及实际网格分辨率须通过相应算例验证。

制造身份按实际 CAD 模型的标准化数据类结果比较，包含默认值。可从已生成且确有附加区域的实际网格元数据导出身份与区域，避免手抄请求：

```python
import json
from pathlib import Path

meshes = {
    member: json.loads(Path(f"outputs/actual_mesh/{member}_mesh.json").read_text(encoding="utf-8"))
    for member in ("pinion", "wheel")
}
document = {
    "schema": "bevel_refinement_regions_v1",
    "units": "m",
    "coordinate_system": "gear_body_local",
    # 这里仅示意源相位0；应填该组控制区真正生成时的源相位，不能由节点数推断。
    "source_driver_angles_rad": [0.0],
    "manufacturing": {
        member: {
            "request": mesh["metadata"]["actual_manufacturing_request"],
            "blank": mesh["metadata"]["actual_blank_settings"],
        }
        for member, mesh in meshes.items()
    },
    "regions": {
        member: mesh["metadata"].get("refinement_regions", [])
        for member, mesh in meshes.items()
    },
}
if not any(document["regions"].values()):
    raise ValueError("该组网格没有附加区域，应使用标准全轨迹输入")
Path("refinement_regions.json").write_text(
    json.dumps(document, ensure_ascii=False, indent=2), encoding="utf-8"
)
```

生成多物理齿的控制点时，世界坐标必须以该轮**整轮姿态**反变换到本体。不能分别用各齿对的局部参考 frame 把多个实际齿折回同一个参考齿。

## 由实际主轴窄带生成区域

`bevel_gear.contact_refinement.create_contact_refinement` 是新写实现，可直接接入已有的真实有限工作面准备对象。它不需要材料、柔度或压力结果；工作侧符号只用于公共几何参考角约定。

仓库提供无需私有缓存的制造输入例程：

```text
python examples/create_contact_refinement.py examples/spiral_17_47_relieved_contact_design.json --intervals 8 --tooth-offsets -1 0 1 --output outputs/contact_refinement_regions.json
```

依赖可通过 `pip install -e ".[surface]"` 安装。默认实际生成闭实体及有限工作面，在±半主动齿距之间采样9个相位；主轴带长/短半宽为8mm/0.6mm，沿带/跨带采样不大于250µm/100µm，目标局部尺寸75µm、core300µm、影响半径1.2mm。可用 `--angles-rad` 明确给相位，或调整 `--intervals`、`--tooth-offsets`、采样和局部尺寸参数。

每相位每齿对先查询真实间隙主轴，再取真实有限NURBS点。每个有效点会以原UV重新取CAD值，并再次检查有限BRep面。行、列只连接相邻的连续有效采样；遇到无效点立即断线，不能跨空档、不同齿或不同相位连接。所有世界点按该相位的整轮R0反变换到本体，故即使未请求k=0控制齿，也不使用各patch frame折回参考齿。

默认遇到边缘主轴stencil失败会保留原基和真实采样点，完整记录原因；`--strict-alignment` 可选择直接失败。它不删除该齿，也不把失败主轴当作非承载证据。

输出区域JSON符合上述输入协议，另有 `.inputs.json` 和 `.diagnostics.json`。只有所有请求相位均完成后才写出正式区域文件；失败时保留未完成诊断。该窄带只是网格细化控制，不替换原全轨迹背景、不缩小接触候选域，也不保证实际载荷一定在带内。9相位不能冒充65/129相位或连续全过程的加密覆盖；最终每个FE状态仍须检查实际active点、深层单元、压力空间与响应收敛。

所选计算相位必须存在于文件声明的源相位中；只有 `[0.0]` 的区域不能直接用于完整周期运行。源相位列表仅是计划覆盖声明，仍需独立检查实际加载区、深层单元和数值收敛，不能把列表本身当作证明。

## 生成与复用

自动划网时，driver 将区域按轮传给 `prepare_generated_ltca(..., refinement_regions={...})`。直接使用该函数时，可以仅提供需要附加场的成员；默认 None 或空列表不改变原划网调用与诊断数据。

显式 `--mesh-dir` 复用时，driver 要求已有网格的制造身份及 `metadata.refinement_regions` 与区域文件一致。已生成的网格不会再次应用附加场。底层 prepare 函数拒绝“复用网格同时再应用新区域”的含混调用。

区域文件原字节复制为 `refinement_regions_input.json`，实际 SHA 和展开后的数值一起进入 `distributed_fe_inputs.json` 的指纹。之后原外部文件被编辑，也不能替换该次计算已冻结的输入。

`analysis_completed` 只表示所选相位均通过其离散验收。完整一齿距、端点、压力/体网格、几何非线性及完整接触面域仍有独立判据；当前接口不会自动将这些验证标记设为 True。
