# 锥齿轮重建计算内核快照

来自公开 [Bevel Gear Design](https://github.com/Lilililiyyy/Bevel-Gear-Design)，版本0.10.0，干净提交 `05da3ac166a750a6326907af5bd0dba8aa6ab034`。`bevel_gear/`、原接口／模型文档和许可按原字节保留；这是公开重建计算内核，不称为厂商原始源码。

文件来源及SHA256见 `src/gear_transmission/bevel_snapshot.json`，函数与源行段对应见 `docs/bevel-source-map.json`。更新算法先在独立工程完成，再使用 `scripts/sync_bevel_snapshot.py --source <干净来源工程>` 显式同步和验收。

Skill新写适配在 `src/gear_transmission/bevel_workflow.py`，入口与示例在 `skills/bevel-gear-design/`。统一命令的LTCA保持原节点FE路线；高级原生分布牵引与半解析API保留但不自动纳入命令复算。AGPL-3.0-only，保留原LICENSE。
