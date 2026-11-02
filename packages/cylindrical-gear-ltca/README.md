# 保留的圆柱齿轮重建计算内核

本目录保存 `cylindrical-gear-ltca` 的完整运行资源快照：全部 `gear_contact` 模块、案例编排及绘图脚本、案例 JSON、公开 MASTA 对标文件、原安装元数据、依赖约束、许可和模型／接口说明。48 个来源文件保持原字节；未在集成层重写公式或修改分支。不复制已生成的结果和缓存。

当前集成版本为 `0.2.1`，来源提交为 `cb2e1ec8a8d1ac2b6137aa6a9b781b9a797c6a7f`。新增的 `gear_contact/contact_display.py` 由源工程维护，统一实际啮合区、中心单元显示支撑、原数据重绘与着色；案例和 Skill 直接调用。快照身份见 `../../src/gear_transmission/kernel_snapshot.json`，逐文件 SHA256 对应来源工作树原始字节（含换行）。

运行时验证实际导入的内核与全部案例、后处理及对标资源。源码对应表见 `../../docs/source-map.json`。新增 `__init__.py` 仅作为 `gear_ltca_snapshot` 分发标记；它和本 README 属于集成说明，未冒充上游文件。本项目内核属于已有新写重建实现，未标成厂商原始源码。

统一安装使用库根的 `pyproject.toml`，分发 `gear_transmission`、原导入名 `gear_contact` 和完整资源包 `gear_ltca_snapshot`。wheel 为直接运行原案例入口保留工程目录；其中内核源码与公开 `gear_contact` 来自同一份源文件，均核验哈希。不要在同一环境安装另一个占用 `gear_contact` 的发行包。

算法仍在原工程维护，本目录为明确固定版本的集成快照。后续同步必须重新记录工作树身份、文件哈希并通过数值对标；不要同时编辑两处公式。当前未复制原工程大型结果、厂商宿主、私有来源资料或Git目录。
