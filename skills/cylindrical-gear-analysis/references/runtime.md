# 运行环境与安装

本Skill需要可执行本地Python的宿主，以及完整 `gear-transmission-skills` 项目提供的两个包：新写 `gear_transmission` 与保留内核 `gear_contact`。Skill文本本身不提供数值执行环境。

建议Python 3.10／3.11；首版在Python 3.10.11、NumPy 1.26.4、SciPy 1.14.1、Matplotlib 3.7.1下验收。先尝试导入，不在已有可用环境中无故重新安装依赖。

未安装时，在完整项目根执行：

```text
python -m pip install -c constraints.txt -e .
```

需要图表时使用 `-e ".[plots]"`。PowerShell可使用 `.venv/Scripts/python`，其他平台用 `.venv/bin/python`。需要创建环境时用 `python -m venv .venv`，然后始终指定这个Python。

Windows下自动化读取中文终端JSON时，设置 `$env:PYTHONUTF8="1"`，避免系统代码页使管道输出乱码；所有保存的JSON、报告和请求文件始终采用UTF-8。

完整项目根含 `gear-skills.project.json`；从Skill文件位置向上查找这个标记，不假定固定父目录层数。Skill被单独注册时，项目可能不在其父目录，此时使用用户提供的项目位置或已安装环境；不得自动寻找或加载其他厂商组件。

通过仓库安装脚本部署时，技能安装目录的父级保存 `.gear-transmission-skills-install.json`，其中 `root` 和 `python` 记录完整项目与虚拟环境解释器路径。先读取该记录并验证路径可用；缺失时再使用用户提供的位置。可在项目根通过 `scripts/update_codex_skills.py --check` 核验安装，若曾指定自定义 `--dest`／`--runtime`，需传入相同参数。

安装后可在任意目录执行 `python -m gear_transmission`，或用本Skill脚本的绝对路径。运行时核验实际加载的内核逐文件SHA256，若与快照不同则停止并修复环境，不能跳过身份检查。避免与另一个提供 `gear_contact` 的发行包混装。

所有结果写入用户指定的新目录；没有覆盖开关。依赖或求解失败必须报告真实错误；本Skill不下载、调用或修改厂商软件。
