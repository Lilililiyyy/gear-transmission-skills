"""命名JSON配置适配；只构造已声明的数据类，不执行配置中的代码。"""
from __future__ import annotations

from dataclasses import fields

from .geometry import Gear, GearPair
from .elasticity import Material
from .surfaces import SurfaceModification, MeasuredCurve, MeasuredSurface, CompositeModification
from .analysis import analyze_pair


def _construct(cls, value):
    if not isinstance(value, dict):
        raise ValueError(f"{cls.__name__} 配置必须为对象")
    unknown = set(value)-{f.name for f in fields(cls)}
    if unknown:
        raise ValueError(f"{cls.__name__} 含未知字段：{', '.join(sorted(unknown))}")
    return cls(**value)


def modification_from_config(value):
    """参数曲面、测量曲线、测量网格及组合曲面的无损命名配置。"""
    if not isinstance(value, dict):
        raise ValueError("修形配置必须为对象")
    value = dict(value)
    kind = value.pop("kind", None)
    if kind is None:
        kind = ("composite" if "components" in value else
                "surface" if "profile_coordinates" in value else
                "curve" if "coordinates" in value else "parameters")
    if kind == "composite":
        if set(value) != {"components"} or not isinstance(value["components"], (list, tuple)):
            raise ValueError("组合曲面只接受components列表")
        return CompositeModification(tuple(modification_from_config(item) for item in value["components"]))
    classes = {"parameters": SurfaceModification, "curve": MeasuredCurve, "surface": MeasuredSurface}
    if kind not in classes:
        raise ValueError("修形kind必须为parameters、curve、surface或composite")
    return _construct(classes[kind], value)


def analyze_config(config):
    """读取可编辑JSON配置，兼容export_result写出的input.json。"""
    if not isinstance(config, dict) or "pair" not in config or "torque_n_m" not in config:
        raise ValueError("配置必须包含pair和torque_n_m")
    values = dict(config)
    pair_data = dict(values.pop("pair"))
    try:
        pinion = _construct(Gear, pair_data.pop("pinion"))
        wheel = _construct(Gear, pair_data.pop("wheel"))
    except KeyError as error:
        raise ValueError("齿轮副配置必须包含pinion和wheel") from error
    pair = _construct(GearPair, {"pinion": pinion, "wheel": wheel, **pair_data})
    for name in ("material_pinion", "material_wheel"):
        if name in values:
            values[name] = _construct(Material, values[name])
    for name in ("modification_pinion", "modification_wheel"):
        if name in values:
            values[name] = modification_from_config(values[name])
    for name in ("external_compliance_provider", "external_boundary_gap"):
        if values.pop(name, False):
            raise ValueError("含外部回调的配置需要由Python调用方重新提供相应函数，不能仅凭JSON复现")
    allowed = {"torque_n_m", "material_pinion", "material_wheel", "modification_pinion", "modification_wheel",
               "misalignment_m", "phases", "strips", "phase_offset", "elastic_mode", "local_contact", "iso_settings", "extended_contact_settings", "local_reference_depths_m"}
    if set(values)-allowed:
        raise ValueError("含未知计算字段："+", ".join(sorted(set(values)-allowed)))
    return analyze_pair(pair, **values)
