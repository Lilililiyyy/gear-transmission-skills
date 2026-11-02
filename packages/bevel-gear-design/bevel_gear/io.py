"""设计配置与计算结果的JSON/CSV读写，导出边界处理无穷寿命。"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path


def json_safe(value):
    """无穷寿命写null，具体含义由life.state区分，不输出非标准JSON数值。"""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def read_config(path) -> dict:
    with Path(path).open("r", encoding="utf-8-sig") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError("配置文件顶层必须为对象")
    return value


def write_json(value, path):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(json_safe(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return target


def export_result(result: dict, output_dir) -> list[str]:
    """保存完整输入/结果、两轮强度汇总、载荷与修正系数表。"""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    files = [write_json(result["inputs"], output / "inputs.json"), write_json(result, output / "result.json")]
    for filename, mapping in (("forces.csv", result["forces"]), ("coefficients.csv", result["rating"]["coefficients"]),
                              ("rating.csv", result["rating"])):
        target = output / filename
        with target.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("quantity", "value"))
            for key, value in mapping.items():
                if isinstance(value, (str, int, float, bool)) or value is None:
                    writer.writerow((key, json_safe(value)))
        files.append(target)
    return [str(path.resolve()) for path in files]
