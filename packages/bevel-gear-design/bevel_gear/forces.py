"""新写的锥齿轮中点啮合力分解，力N、扭矩N·m。

每轮局部径向正方向远离轴线，轴向正方向指向共同节锥顶。
切向值给出阻抗扭矩对应的幅值；两轮局部基向量不同，不能把分量
直接相加作全局平衡。drive/coast选择两侧工作齿面，旋向单独指定。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .geometry import GeometryResult


@dataclass(frozen=True)
class MeshForces:
    pinion_torque_n_m: float
    wheel_torque_n_m: float
    pinion_tangential_n: float
    pinion_radial_n: float
    pinion_axial_n: float
    wheel_tangential_n: float
    wheel_radial_n: float
    wheel_axial_n: float
    normal_force_n: float
    flank: str
    spiral_hand: int
    diagnostics: dict


def calculate_mesh_forces(geometry: GeometryResult, pinion_torque_n_m: float, *,
                          flank: str = "drive", spiral_hand: int = 1) -> MeshForces:
    """对称压力角、零偏置、无效率损失的中点局部力。

    扭矩输入为非负幅值；旋向+1/-1及drive/coast共同确定螺旋分量。
    +1定义为drive面上使小轮轴向螺旋分量为正的旋向；转换成加工
    左/右旋标记时必须结合观察方向和装配基向量，不能直接套用。
    本接口不推断转向、驱动件切换或全局轴系坐标；需要反向工况时
    调用方必须明确选择相应工作齿面和基向量。
    """
    if not isinstance(geometry, GeometryResult):
        raise TypeError("啮合力入口需要GeometryResult")
    if not math.isfinite(pinion_torque_n_m) or pinion_torque_n_m < 0:
        raise ValueError("小轮扭矩幅值须为非负有限数，单位N·m")
    if flank not in ("drive", "coast"):
        raise ValueError("工作齿面只能为drive或coast")
    if isinstance(spiral_hand, bool) or spiral_hand not in (-1, 1):
        raise ValueError("旋向标志必须为+1或-1")
    beta = geometry.inputs.spiral_angle_rad
    alpha = geometry.inputs.pressure_angle_rad
    tangent = 2 * pinion_torque_n_m / geometry.pinion.mean_pitch_diameter_m
    pressure = math.tan(alpha) / math.cos(beta)
    spiral = spiral_hand * (1 if flank == "drive" else -1) * math.tan(beta)
    def components(delta, sign):
        radial = tangent * (pressure * math.cos(delta) - sign * spiral * math.sin(delta))
        axial = tangent * (pressure * math.sin(delta) + sign * spiral * math.cos(delta))
        return radial, axial
    radial_p, axial_p = components(geometry.pinion.pitch_cone_angle_rad, 1)
    radial_w, axial_w = components(geometry.wheel.pitch_cone_angle_rad, -1)
    normal = tangent / (math.cos(alpha) * math.cos(beta))
    magnitude_p = math.sqrt(tangent*tangent + radial_p*radial_p + axial_p*axial_p)
    magnitude_w = math.sqrt(tangent*tangent + radial_w*radial_w + axial_w*axial_w)
    return MeshForces(
        pinion_torque_n_m, pinion_torque_n_m * geometry.ratio,
        tangent, radial_p, axial_p, tangent, radial_w, axial_w, normal, flank, spiral_hand,
        {"pinion_normal_magnitude_error_n": magnitude_p - normal,
         "wheel_normal_magnitude_error_n": magnitude_w - normal,
         "efficiency": 1.0,
         "coordinate_convention": "每轮径向远离轴线、轴向指向共同锥顶，切向为幅值"},
    )
