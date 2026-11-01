"""新写资源定位：已登记的独立工程运行资源同时适用于源码安装和 wheel。"""
from pathlib import Path

import gear_ltca_snapshot


def snapshot_root():
    return Path(gear_ltca_snapshot.__file__).resolve().parent
