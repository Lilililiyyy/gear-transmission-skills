# 环境与计算路线

本Skill随库统一发布，不独立编号；实际内核版本与固定提交见 README 技能表格及 `bevel_snapshot.json`。

保留完整仓库，并在独立Python 3.10／3.11环境安装 `python -m pip install -e ".[plots]"`，或使用库的 `scripts/update-codex-skills.ps1`／`.sh`。安装器会登记两个Skill并校验锥齿轮内核哈希。只复制SKILL.md不够；Skill脚本也不包含Python运行时。

几何和名义强度不需要CAD。统一入口的 `cad`、`tca`、`ltca` 按需安装：

```text
python -m pip install -e ".[bevel-contact,plots]"
```

可选依赖为NumPy、SciPy、cadquery-ocp、Gmsh、VTK、Shapely和psutil，具体范围见项目pyproject.toml。未安装时保留原缺失依赖错误，不退回假齿面或样例数组。大规模FE需先核对节点数、内存和用户的运行约束。

安装后的等价入口为 `python -m gear_transmission.bevel_cli` 或 `gear-bevel`。也可直接使用原函数：

```python
from bevel_gear.workflow import load_bevel_case
from gear_transmission.bevel_workflow import run_case, replay_run
case = load_bevel_case("案例绝对路径.json")
run_case(case, "新输出绝对路径", through="strength")
replay_run("原输出绝对路径", "新复算绝对路径")
```

| 任务路线 | 原内核入口与边界 |
|---|---|
| 统一分阶段流程 | `workflow.analyze_bevel_case`，Skill命令使用此路线；LTCA为节点FE兼容入口，不能称为分布P0牵引解 |
| 宏观设计与其他标准分支 | `api.analyze_design`、`rating`、`rating_b2`、`agma_rating`；读取快照`docs/api.md`及原函数，分支和输入不同于统一2014名义链 |
| 原生分布牵引FE | `generated_distributed_fe`、`native_fe_contact`、`native_fe_operator`；需要真实prepared网格、材料、固定自由度、相位和请求patch，保持当前operator与结果身份 |
| 局部半解析 | `surface_workflow.run_surface_workflow`；矩形面元半空间柔度及局部闭合，不包含整轮体柔度 |
| 原始压力与周期显示 | `contact_display`、`contact_composite`、`distributed_fe_results`；须提供相应原场和覆盖证据，不从图片或其他模型压力补造输入 |

后四类原生API可按用户明确任务新写调用脚本，但没有自动进入本Skill的run/replay协议，也不代表当前发布已验收其全部组合。读取仓库 `packages/bevel-gear-design/docs/api.md`、`model.md` 和对应源码；源码已随本库保留，勿下载不固定版本替换已验收计算。

内核来源及逐文件哈希见 `src/gear_transmission/bevel_snapshot.json`，函数行段映射见 `docs/bevel-source-map.json`。这是公开重建计算内核；适配层属于新写实现，原厂商组件不参与运行。
