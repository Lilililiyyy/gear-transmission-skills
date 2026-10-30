"""新写的渐开线圆柱齿轮几何；内部单位为 m、rad、N。

模型采用理想渐开螺旋面和给定实际根、顶圆。基圆以下齿厚仅作径向
延拓，不生成刀具过渡曲线，不判定根切，也不代替齿侧间隙或干涉检查。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def _involute(angle: float) -> float:
    return math.tan(angle) - angle


def _inverse_involute(value: float) -> float:
    """在物理区间内反求渐开线函数，不引入迭代初值经验式。"""
    if value < 0 or not math.isfinite(value):
        raise ValueError("变位组合不能形成实数工作压力角")
    low, high = 0.0, math.pi / 2 - 1e-9
    for _ in range(80):
        middle = (low + high) / 2
        if _involute(middle) < value:
            low = middle
        else:
            high = middle
    return (low + high) / 2


@dataclass(frozen=True)
class Gear:
    """单轮几何。螺旋角带手性，正值和负值表示相反旋向。

    默认外齿顶高/根高系数为 1/1.25，内齿采用对应的径向相反方向。
    profile_shift 是法向齿条变位系数；内齿参考齿厚的变位项取负号。
    tip/root_diameter_m 定义实体边界，tip/root_form_diameter_m 定义
    可用渐开线的顶侧/根侧边界。梁的固定根仍是实体根圆，接触与
    修形归一化使用完整有效形径范围，不随当前接触区域重新缩放。
    """
    teeth: int
    normal_module_m: float
    face_width_m: float
    pressure_angle_rad: float = math.radians(20)
    helix_angle_rad: float = 0.0
    profile_shift: float = 0.0
    internal: bool = False
    tip_diameter_m: float | None = None
    root_diameter_m: float | None = None
    root_form_diameter_m: float | None = None
    tip_form_diameter_m: float | None = None

    def __post_init__(self) -> None:
        if isinstance(self.teeth, bool) or not isinstance(self.teeth, (int, np.integer)) or self.teeth < 3:
            raise ValueError("齿数必须为不小于3的整数")
        for name in ("normal_module_m", "face_width_m"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} 必须为正有限长度，单位m")
        if not math.isfinite(self.pressure_angle_rad) or not 0 < self.pressure_angle_rad < math.pi / 2:
            raise ValueError("法向压力角必须在0和π/2之间，单位rad")
        if not math.isfinite(self.helix_angle_rad) or abs(self.helix_angle_rad) >= math.pi / 2:
            raise ValueError("螺旋角绝对值必须小于π/2，单位rad")
        if not math.isfinite(self.profile_shift):
            raise ValueError("变位系数必须为有限数")
        for name in ("tip_diameter_m", "root_diameter_m", "root_form_diameter_m", "tip_form_diameter_m"):
            value = getattr(self, name)
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} 必须为正有限长度，单位m")
        inner, outer = sorted((self.tip_radius_m, self.root_radius_m))
        if inner <= 0 or outer <= inner:
            raise ValueError("实际根、顶圆必须为不同的正半径")
        if (self.internal and self.tip_radius_m >= self.root_radius_m) or (
            not self.internal and self.tip_radius_m <= self.root_radius_m
        ):
            raise ValueError("外齿顶圆须大于根圆，内齿顶圆须小于根圆")
        if self.internal and self.tip_radius_m < self.base_radius_m:
            raise ValueError("内齿顶圆不能低于基圆，当前模型只支持完整渐开线内齿廓")
        if outer <= self.base_radius_m:
            raise ValueError("实际齿廓没有可用渐开线区间")
        form_root, form_tip = self.root_form_radius_m, self.tip_form_radius_m
        if not inner <= form_root <= outer or not inner <= form_tip <= outer:
            raise ValueError("根部与顶部形径必须位于实体根圆与顶圆之间")
        if min(form_root, form_tip) < self.base_radius_m:
            raise ValueError("有效渐开线形径不能低于基圆")
        if (self.internal and form_root <= form_tip) or (not self.internal and form_root >= form_tip):
            raise ValueError("根部与顶部形径必须按该齿轮的径向生长方向严格有序")
        if self.tooth_thickness_m(self.tip_radius_m) <= 0:
            raise ValueError("当前参数导致齿顶齿厚非正，请检查变位与实际齿顶直径")

    @property
    def transverse_module_m(self) -> float:
        return self.normal_module_m / math.cos(self.helix_angle_rad)

    @property
    def reference_radius_m(self) -> float:
        return self.transverse_module_m * self.teeth / 2

    @property
    def transverse_pressure_angle_rad(self) -> float:
        return math.atan(math.tan(self.pressure_angle_rad) / math.cos(self.helix_angle_rad))

    @property
    def base_radius_m(self) -> float:
        return self.reference_radius_m * math.cos(self.transverse_pressure_angle_rad)

    @property
    def base_helix_angle_rad(self) -> float:
        return math.atan(math.tan(self.helix_angle_rad) * math.cos(self.transverse_pressure_angle_rad))

    @property
    def tip_radius_m(self) -> float:
        if self.tip_diameter_m is not None:
            return self.tip_diameter_m / 2
        direction = -1 if self.internal else 1
        return self.reference_radius_m + self.normal_module_m * (direction + self.profile_shift)

    @property
    def root_radius_m(self) -> float:
        if self.root_diameter_m is not None:
            return self.root_diameter_m / 2
        direction = 1 if self.internal else -1
        return self.reference_radius_m + self.normal_module_m * (1.25 * direction + self.profile_shift)

    @property
    def root_form_radius_m(self) -> float:
        """根侧有效渐开线形径，与实体根圆及梁的固定边界分开。"""
        if self.root_form_diameter_m is not None:
            return self.root_form_diameter_m / 2
        return self.root_radius_m if self.internal else max(self.root_radius_m, self.base_radius_m)

    @property
    def tip_form_radius_m(self) -> float:
        """顶侧有效渐开线形径，缺省等于实体齿顶圆。"""
        return self.tip_radius_m if self.tip_form_diameter_m is None else self.tip_form_diameter_m / 2

    @property
    def involute_roll_bounds_m(self) -> tuple[float, float]:
        """完整可用渐开线的最小/最大正展长；内齿返回顶部/根部。"""
        inner, outer = sorted((self.tip_form_radius_m, self.root_form_radius_m))
        rb = self.base_radius_m
        return math.sqrt(max(0.0, inner * inner - rb * rb)), math.sqrt(outer * outer - rb * rb)

    def tooth_thickness_m(self, radius):
        """给定半径的端面弧齿厚m，可接受数组。

        基圆以下令展角为0并保持角齿厚，从而给固定齿根梁一个明确的
        近似截面；实际渐开线外的过渡圆角及根部应力集中不在此式中。
        """
        radii = np.asarray(radius, dtype=float)
        inner, outer = sorted((self.tip_radius_m, self.root_radius_m))
        tolerance = 1e-12 * outer
        if np.any(~np.isfinite(radii)) or np.any(radii < inner - tolerance) or np.any(radii > outer + tolerance):
            raise ValueError("齿厚计算半径必须位于实际根圆与顶圆之间，单位m")
        radii = np.clip(radii, inner, outer)
        rb = self.base_radius_m
        tangent = np.sqrt(np.maximum((radii / rb) ** 2 - 1.0, 0.0))
        inv = tangent - np.arctan(tangent)
        direction = -1 if self.internal else 1
        reference_angle_thickness = (
            math.pi + direction * 4 * self.profile_shift * math.tan(self.pressure_angle_rad)
        ) / self.teeth
        angular_thickness = reference_angle_thickness + direction * 2 * (
            _involute(self.transverse_pressure_angle_rad) - inv
        )
        result = radii * angular_thickness
        if np.any(result <= 0):
            raise ValueError("给定半径处渐开线齿厚非正")
        return float(result) if result.ndim == 0 else result


@dataclass(frozen=True)
class GearPair:
    """平行轴啮合副；共同轴向坐标以小轮齿宽中点为零。

    未给中心距时，以变位和无齿侧间隙条件反求工作压力角。显式中心距
    可覆盖它，但当前几何不自动补偿齿侧间隙。只有大轮允许为内齿圈。
    """
    pinion: Gear
    wheel: Gear
    center_distance_m: float | None = None
    axial_offset_m: float = 0.0

    def __post_init__(self) -> None:
        p, w = self.pinion, self.wheel
        if p.internal:
            raise ValueError("小轮必须为外齿，内啮合请设置wheel.internal=True")
        if not math.isclose(p.normal_module_m, w.normal_module_m, rel_tol=1e-10, abs_tol=0):
            raise ValueError("两轮法向模数必须一致")
        if not math.isclose(p.pressure_angle_rad, w.pressure_angle_rad, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("两轮法向压力角必须一致")
        expected = p.helix_angle_rad if w.internal else -p.helix_angle_rad
        if not math.isclose(w.helix_angle_rad, expected, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("外啮合须等值异向螺旋角，内啮合须等值同向螺旋角")
        if w.internal and w.teeth <= p.teeth:
            raise ValueError("内齿圈齿数必须大于小轮齿数")
        if not math.isfinite(self.axial_offset_m):
            raise ValueError("轴向错位必须为有限数，单位m")
        if self.center_distance_m is None:
            if w.internal:
                denominator, shift = w.teeth - p.teeth, w.profile_shift - p.profile_shift
            else:
                denominator, shift = w.teeth + p.teeth, w.profile_shift + p.profile_shift
            inv_working = _involute(p.transverse_pressure_angle_rad) + 2 * shift * math.tan(p.pressure_angle_rad) / denominator
            alpha = _inverse_involute(inv_working)
            base_span = w.base_radius_m + (-p.base_radius_m if w.internal else p.base_radius_m)
            object.__setattr__(self, "center_distance_m", base_span / math.cos(alpha))
        if not math.isfinite(self.center_distance_m) or self.center_distance_m <= 0:
            raise ValueError("中心距必须为正有限长度，单位m")
        base_span = w.base_radius_m + (-p.base_radius_m if w.internal else p.base_radius_m)
        if self.center_distance_m <= base_span:
            raise ValueError("中心距必须大于相应基圆半径之和或之差")
        if self.face_overlap_m <= 0:
            raise ValueError("两轮没有正的轴向齿宽重叠")
        if self.contact_roll_bounds_m[1] <= self.contact_roll_bounds_m[0]:
            raise ValueError("实际根顶圆裁剪后没有可用啮合路径")

    @property
    def base_pitch_m(self) -> float:
        return 2 * math.pi * self.pinion.base_radius_m / self.pinion.teeth

    @property
    def working_pressure_angle_rad(self) -> float:
        base_span = self.wheel.base_radius_m + (-self.pinion.base_radius_m if self.wheel.internal else self.pinion.base_radius_m)
        return math.acos(base_span / self.center_distance_m)

    @property
    def base_helix_angle_rad(self) -> float:
        return self.pinion.base_helix_angle_rad

    @property
    def face_bounds_m(self) -> tuple[float, float]:
        return (
            max(-self.pinion.face_width_m / 2, self.axial_offset_m - self.wheel.face_width_m / 2),
            min(self.pinion.face_width_m / 2, self.axial_offset_m + self.wheel.face_width_m / 2),
        )

    @property
    def face_overlap_m(self) -> float:
        low, high = self.face_bounds_m
        return high - low

    @property
    def contact_roll_bounds_m(self) -> tuple[float, float]:
        """将两轮完整渐开线区间映射到小轮展长并求交。"""
        pl, ph = self.pinion.involute_roll_bounds_m
        wl, wh = self.wheel.involute_roll_bounds_m
        span = self.center_distance_m * math.sin(self.working_pressure_angle_rad)
        if self.wheel.internal:
            return max(pl, wl - span), min(ph, wh - span)
        return max(pl, span - wh), min(ph, span - wl)

    @property
    def transverse_contact_ratio(self) -> float:
        low, high = self.contact_roll_bounds_m
        return (high - low) / self.base_pitch_m

    @property
    def overlap_ratio(self) -> float:
        return self.face_overlap_m * abs(math.tan(self.base_helix_angle_rad)) / self.base_pitch_m

    @property
    def normal_force_lever_arm_m(self) -> float:
        """虚功 F_n·dδ=T_p·dθ_p 对应的法向力臂。"""
        return self.pinion.base_radius_m * math.cos(self.base_helix_angle_rad)

    def normal_force_from_torque(self, torque_n_m: float) -> float:
        if not math.isfinite(torque_n_m) or torque_n_m < 0:
            raise ValueError("扭矩幅值必须为非负有限数，单位N·m")
        return torque_n_m / self.normal_force_lever_arm_m

    def normal_approach_to_angular_te(self, approach_m):
        approach = np.asarray(approach_m, dtype=float)
        if np.any(~np.isfinite(approach)):
            raise ValueError("法向接近量必须为有限数，单位m")
        answer = approach / self.normal_force_lever_arm_m
        return float(answer) if answer.ndim == 0 else answer


@dataclass(frozen=True)
class ContactGrid:
    """原接触线上的中点积分单元，widths_m 是实际线长而非轴向宽度。"""
    roll_pinion_m: np.ndarray
    roll_wheel_m: np.ndarray
    face_coordinate_m: np.ndarray
    widths_m: np.ndarray
    tooth_pair_index: np.ndarray
    profile_factor_pinion: np.ndarray
    profile_factor_wheel: np.ndarray
    lead_factor_pinion: np.ndarray
    lead_factor_wheel: np.ndarray
    curvature_per_m: np.ndarray
    radius_pinion_m: np.ndarray
    radius_wheel_m: np.ndarray

    def __len__(self) -> int:
        return len(self.widths_m)

    @property
    def has_singular_curvature(self) -> bool:
        """基圆尖点的曲率半径为0；几何保留它，局部接触层须拒绝该域。"""
        return bool(np.any(~np.isfinite(self.curvature_per_m)))


def contact_grid(pair: GearPair, phase: float, strips: int = 21) -> ContactGrid:
    """构造一基节周期内的接触线，phase按周期取模。

    每条线先与有效齿宽和两轮实际渐开线范围精确相交，再均匀划分
    strips 个中点单元。半基节相位偏置仅选取相位原点，不改变周期。
    相位沿小轮正展长前进；螺旋角的符号决定接触线轴向斜率。
    在基圆尖点保留0曲率半径和无穷组合曲率，高层局部接触模型必须
    拒绝该域；本函数不以任意微小半径替换真实几何端点。
    """
    if not math.isfinite(phase):
        raise ValueError("相位必须为有限数，单位为基节周期")
    if isinstance(strips, bool) or not isinstance(strips, (int, np.integer)) or strips < 1:
        raise ValueError("每条接触线的积分单元数必须为正整数")
    left, right = pair.face_bounds_m
    low, high = pair.contact_roll_bounds_m
    pitch = pair.base_pitch_m
    beta = pair.base_helix_angle_rad
    slope = math.tan(beta)
    cosine = math.cos(beta)
    origin = low + ((phase % 1.0) + 0.5) * pitch
    offset_low, offset_high = sorted((left * slope, right * slope))
    first = math.ceil((low - offset_high - origin) / pitch)
    last = math.floor((high - offset_low - origin) / pitch)
    face_coordinates, rolls, lengths, indices = [], [], [], []
    for index in range(first, last + 1):
        intercept = origin + index * pitch
        line_left, line_right = left, right
        if abs(slope) > 1e-15:
            clip_low, clip_high = sorted(((low - intercept) / slope, (high - intercept) / slope))
            line_left, line_right = max(left, clip_low), min(right, clip_high)
        elif intercept < low or intercept > high:
            continue
        if line_right <= line_left + 1e-14 * pair.face_overlap_m:
            continue
        dx = (line_right - line_left) / strips
        coordinates = line_left + (np.arange(strips) + 0.5) * dx
        face_coordinates.append(coordinates)
        rolls.append(intercept + coordinates * slope)
        lengths.append(np.full(strips, dx / cosine))
        indices.append(np.full(strips, index, dtype=int))
    if not rolls:
        raise ValueError("当前相位没有理论接触线；该几何存在啮合中断")
    x = np.concatenate(face_coordinates)
    rp = np.concatenate(rolls)
    widths = np.concatenate(lengths)
    ids = np.concatenate(indices)
    span = pair.center_distance_m * math.sin(pair.working_pressure_angle_rad)
    rw = span + rp if pair.wheel.internal else span - rp
    tolerance = 1e-13 * max(pair.pinion.reference_radius_m, pair.wheel.reference_radius_m)
    if np.any(rp < -tolerance) or np.any(rw < -tolerance):
        raise ValueError("接触点展长为负，当前渐开线几何不成立")
    # 只消除边界浮点舍入的负零，不将基圆尖点改成任意有限曲率。
    rp, rw = np.maximum(rp, 0), np.maximum(rw, 0)
    radius_p, radius_w = rp / cosine, rw / cosine
    inverse_p = np.divide(1.0, radius_p, out=np.full_like(radius_p, np.inf), where=radius_p > 0)
    inverse_w = np.divide(1.0, radius_w, out=np.full_like(radius_w, np.inf), where=radius_w > 0)
    curvature = inverse_p + (-1 if pair.wheel.internal else 1) * inverse_w
    if np.any(curvature <= 0):
        raise ValueError("组合法截面曲率必须为正")
    pl, ph = pair.pinion.involute_roll_bounds_m
    wl, wh = pair.wheel.involute_roll_bounds_m
    fp = np.clip((rp - pl) / (ph - pl), 0, 1)
    fw = np.clip(((wh - rw) if pair.wheel.internal else (rw - wl)) / (wh - wl), 0, 1)
    lp = np.clip((x + pair.pinion.face_width_m / 2) / pair.pinion.face_width_m, 0, 1)
    lw = np.clip((x - pair.axial_offset_m + pair.wheel.face_width_m / 2) / pair.wheel.face_width_m, 0, 1)
    return ContactGrid(rp, rw, x, widths, ids, fp, fw, lp, lw, curvature, radius_p, radius_w)
