"""第一性原理新写的理想球面渐开线，不含刀具包络或齿根圆角。

在单位球上沿基圆展开切向大圆弧。若基锥半角为δb、球面展开弧长
为s，则δ满足cosδ=cosδb*cos(s)，方位渐开线函数为
s/sinδb-atan2(sin(s),sinδb*cos(s))。这是一条真实解析曲线，
不是把平面渐开线贴到锥体上，也不是实际加工齿面的恢复结果。
局部坐标从共同锥顶向齿轮外端为+Z，齿厚中心方位为0，节点处+X
为径向、+Y为切向。所有长度m、角度rad。
"""
from __future__ import annotations

from dataclasses import dataclass
import math


Vector3 = tuple[float, float, float]


def _finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name}必须为有限数")
    return float(value)


def spherical_involute_function(polar_angle_rad: float, base_cone_angle_rad: float) -> float:
    """球面渐开线方位函数，基圆处取连续极限0。"""
    delta = _finite(polar_angle_rad, "极角")
    base = _finite(base_cone_angle_rad, "基锥角")
    if not 0 < base < math.pi / 2 or not base - 2e-14 <= delta < math.pi / 2:
        raise ValueError("球面渐开线要求0<基锥角<=极角<π/2；基圆以下不是渐开线")
    delta = max(delta, base)
    arc = math.acos(min(1.0, max(-1.0, math.cos(delta) / math.cos(base))))
    return arc / math.sin(base) - math.atan2(math.sin(arc), math.sin(base) * math.cos(arc))


def spherical_point(cone_distance_m: float, polar_angle_rad: float, azimuth_rad: float) -> Vector3:
    """以锥顶为球心的球面坐标，不包含轮系装配变换。"""
    radius = _finite(cone_distance_m, "球面锥距")
    polar = _finite(polar_angle_rad, "极角")
    azimuth = _finite(azimuth_rad, "方位角")
    if radius <= 0 or not 0 <= polar <= math.pi:
        raise ValueError("球面半径必须为正，极角必须位于[0,π]")
    radial = radius * math.sin(polar)
    return radial * math.cos(azimuth), radial * math.sin(azimuth), radius * math.cos(polar)


def local_tooth_trace_angle(cone_distance_m: float, reference_cone_distance_m: float,
                           reference_spiral_angle_rad: float, *, trace: str = "straight",
                           cutter_radius_m: float | None = None) -> float:
    """独立的局部齿线角关系，不等于实际刀具包络齿面。

    constant_lead采用R*sinβ不变；equal_pitch采用tanβ/R不变；
    circular_arc以显式刀盘半径计算圆弧关系。此函数不能据局部齿线角
    直接生成所谓真实弧齿面；本模块的解析齿面仍严格限于直齿。
    """
    radius = _finite(cone_distance_m, "球面锥距")
    reference = _finite(reference_cone_distance_m, "参考锥距")
    beta = _finite(reference_spiral_angle_rad, "参考螺旋角")
    if radius <= 0 or reference <= 0 or abs(beta) >= math.pi / 2:
        raise ValueError("锥距须为正，参考螺旋角绝对值须小于π/2")
    if trace == "straight":
        if beta != 0:
            raise ValueError("直齿线不能指定非零参考螺旋角")
        return 0.0
    if trace == "equal_pitch":
        return math.atan(radius * math.tan(beta) / reference)
    if trace == "constant_lead":
        argument = reference * math.sin(beta) / radius
    elif trace == "circular_arc":
        if cutter_radius_m is None or _finite(cutter_radius_m, "刀盘半径") <= 0:
            raise ValueError("圆弧齿线须显式指定正刀盘半径m")
        argument = (radius*radius + 2*reference*cutter_radius_m*math.sin(beta) - reference*reference) / (2*radius*cutter_radius_m)
    else:
        raise ValueError("齿线只能为straight、constant_lead、equal_pitch或circular_arc")
    if abs(argument) > 1 + 1e-12:
        raise ValueError("该截面局部齿线角超出反正弦实数域")
    return math.asin(min(1.0, max(-1.0, argument)))


@dataclass(frozen=True)
class SphericalInvoluteProfile:
    """理想直齿球面渐开线的一齿双侧齿廓。

    pitch_angular_tooth_thickness_rad是节点圆上的完整齿厚角，等于
    节点端面弧齿厚除以节圆半径。法向和端面压力角仅在直齿时相同；
    此类不以附加扭转冒充实际弧齿刀具生成齿面。
    """
    pitch_cone_angle_rad: float
    pressure_angle_rad: float
    pitch_angular_tooth_thickness_rad: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            _finite(getattr(self, name), name)
        if not 0 < self.pitch_cone_angle_rad < math.pi / 2:
            raise ValueError("理想球面渐开线要求锐角节锥")
        if not 0 < self.pressure_angle_rad < math.pi / 2:
            raise ValueError("节点压力角必须在0和π/2之间")
        if not 0 < self.pitch_angular_tooth_thickness_rad < math.pi:
            raise ValueError("节点完整齿厚角必须在0和π之间")

    @property
    def base_cone_angle_rad(self) -> float:
        return math.asin(math.sin(self.pitch_cone_angle_rad) * math.cos(self.pressure_angle_rad))

    def _side(self, flank: str) -> int:
        if flank not in ("left", "right"):
            raise ValueError("齿面只能为left或right")
        return 1 if flank == "right" else -1

    def half_material_angle(self, polar_angle_rad: float) -> float:
        """正值表示两齿面仍分离；非正值表示达到尖齿或齿面相交。"""
        return (self.pitch_angular_tooth_thickness_rad / 2
                + spherical_involute_function(self.pitch_cone_angle_rad, self.base_cone_angle_rad)
                - spherical_involute_function(polar_angle_rad, self.base_cone_angle_rad))

    def flank_angle(self, polar_angle_rad: float, flank: str = "right") -> float:
        sign = self._side(flank)
        half = self.half_material_angle(polar_angle_rad)
        if half <= 0:
            raise ValueError("该极角处齿厚非正，已达到尖齿或两解析齿面相交")
        return sign * half

    def pressure_angle_at(self, polar_angle_rad: float) -> float:
        """解析单位法向与局部切向的夹角；基圆处为0。"""
        spherical_involute_function(polar_angle_rad, self.base_cone_angle_rad)
        return math.acos(min(1.0, math.sin(self.base_cone_angle_rad) / math.sin(polar_angle_rad)))

    def evaluate(self, cone_distance_m: float, polar_angle_rad: float,
                 flank: str = "right") -> Vector3:
        return spherical_point(cone_distance_m, polar_angle_rad, self.flank_angle(polar_angle_rad, flank))

    def normal(self, polar_angle_rad: float, flank: str = "right") -> Vector3:
        """朝向齿轮材料外部的解析单位法向，基圆处采用单侧极限。

        法向与锥顶径向正交；基圆处曲率奇异仍存在，有限法向并不表示
        可以在该处使用有限曲率的局部接触应力模型。
        """
        sign = self._side(flank)
        angle = self.flank_angle(polar_angle_rad, flank)
        alpha = self.pressure_angle_at(polar_angle_rad)
        cp, sp = math.cos(angle), math.sin(angle)
        cd, sd = math.cos(polar_angle_rad), math.sin(polar_angle_rad)
        # e_delta指向增大极角，e_phi为增大方位角；两者均为单位向量。
        e_delta = cd * cp, cd * sp, -sd
        e_phi = -sp, cp, 0.0
        return tuple(math.sin(alpha) * e_delta[i] + sign * math.cos(alpha) * e_phi[i] for i in range(3))
