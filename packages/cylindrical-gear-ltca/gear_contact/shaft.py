"""新写实现：阶梯轴、线性轴承及其向齿面法向的准静态柔度投影。"""
from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, sin, tan

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from .elasticity import Material


@dataclass(frozen=True)
class ShaftSegment:
    """实心圆轴段，长度和直径为 m。"""

    length_m: float
    diameter_m: float
    material: Material = Material()

    def __post_init__(self):
        if any(isinstance(v, bool) or not np.isfinite(v) or v <= 0
               for v in (self.length_m, self.diameter_m)):
            raise ValueError("轴段长度和直径必须为有限正数，单位 m")
        if not isinstance(self.material, Material):
            raise ValueError("轴段材料必须为 Material")


@dataclass(frozen=True)
class LinearBearing:
    """刚性箱体上的三向平动弹簧，刚度顺序为 x、y、z，单位 N/m。"""

    node: int
    stiffness_n_per_m: tuple[float, float, float]

    def __post_init__(self):
        stiffness = np.asarray(self.stiffness_n_per_m, dtype=float)
        if (stiffness.shape != (3,) or np.any(~np.isfinite(stiffness))
                or np.any(stiffness < 0) or not np.any(stiffness > 0)):
            raise ValueError("轴承刚度须为三个非负有限值，至少一向大于零，单位 N/m")
        object.__setattr__(self, "stiffness_n_per_m", tuple(map(float, stiffness)))


@dataclass(frozen=True)
class Shaft:
    """节点从轴首 0 开始沿 +z 编号。扭矩端只约束绕 z 的角位移。"""

    segments: tuple[ShaftSegment, ...]
    bearings: tuple[LinearBearing, ...]
    gear_node: int
    torque_node: int

    def __post_init__(self):
        segments, bearings = tuple(self.segments), tuple(self.bearings)
        if not segments or any(not isinstance(s, ShaftSegment) for s in segments):
            raise ValueError("轴至少需要一个 ShaftSegment")
        if any(not isinstance(b, LinearBearing) for b in bearings):
            raise ValueError("轴承必须为 LinearBearing")
        for node in (self.gear_node, self.torque_node, *(b.node for b in bearings)):
            if isinstance(node, bool) or not isinstance(node, (int, np.integer)) or not 0 <= node <= len(segments):
                raise ValueError("齿轮、轴承和扭矩端均须位于有效的零基轴节点")
        if len({b.node for b in bearings}) != len(bearings):
            raise ValueError("同一轴节点只允许定义一个轴承")
        for axis in (0, 1):
            if sum(b.stiffness_n_per_m[axis] > 0 for b in bearings) < 2:
                raise ValueError("每个径向方向至少需要两个不同节点上的正刚度支承")
        if not any(b.stiffness_n_per_m[2] > 0 for b in bearings):
            raise ValueError("每根轴至少需要一处正轴向刚度用于轴向定位")
        object.__setattr__(self, "segments", segments)
        object.__setattr__(self, "bearings", bearings)


@dataclass(frozen=True)
class ShaftSystem:
    """与小轮和大轮对应的两根平行轴，轴承座固定于刚性箱体。"""

    pinion: Shaft
    wheel: Shaft

    def __post_init__(self):
        if not isinstance(self.pinion, Shaft) or not isinstance(self.wheel, Shaft):
            raise ValueError("轴系必须分别给出小轮轴和大轮轴")


def _element_stiffness(segment):
    """两节点等截面圆梁。每节点顺序为 ux、uy、uz、θx、θy、θz。"""
    length, diameter = segment.length_m, segment.diameter_m
    e, g = segment.material.youngs_modulus_pa, segment.material.shear_modulus_pa
    area, inertia = pi * diameter**2 / 4, pi * diameter**4 / 64
    nu = segment.material.poisson_ratio
    shear_area = area * 6 * (1 + nu) / (7 + 6 * nu)
    phi = 12 * e * inertia / (g * shear_area * length**2)
    l2 = length**2
    # 转角按截面转动定义，x 位移对应 +θy，y 位移对应 −θx。
    bending = e * inertia / (length**3 * (1 + phi)) * np.array([
        [12, 6 * length, -12, 6 * length],
        [6 * length, (4 + phi) * l2, -6 * length, (2 - phi) * l2],
        [-12, -6 * length, 12, -6 * length],
        [6 * length, (2 - phi) * l2, -6 * length, (4 + phi) * l2],
    ])
    matrix = np.zeros((12, 12))
    matrix[np.ix_([0, 4, 6, 10], [0, 4, 6, 10])] = bending
    signs = np.array([1, -1, 1, -1])
    matrix[np.ix_([1, 3, 7, 9], [1, 3, 7, 9])] = bending * np.outer(signs, signs)
    extension = np.array([[1, -1], [-1, 1]]) / length
    matrix[np.ix_([2, 8], [2, 8])] = e * area * extension
    matrix[np.ix_([5, 11], [5, 11])] = g * (2 * inertia) * extension
    return matrix


class _ShaftFlexibility:
    def __init__(self, shaft):
        self.shaft = shaft
        self.positions = np.r_[0., np.cumsum([s.length_m for s in shaft.segments])]
        size = 6 * len(self.positions)
        stiffness = np.zeros((size, size))
        for index, segment in enumerate(shaft.segments):
            indices = slice(6 * index, 6 * index + 12)
            stiffness[indices, indices] += _element_stiffness(segment)
        for bearing in shaft.bearings:
            indices = 6 * bearing.node + np.arange(3)
            stiffness[indices, indices] += bearing.stiffness_n_per_m
        self.gear = slice(6 * shaft.gear_node, 6 * shaft.gear_node + 6)
        free = np.arange(size) != 6 * shaft.torque_node + 5
        reduced = stiffness[np.ix_(free, free)]
        if not np.all(np.isfinite(reduced)) or np.any(reduced.diagonal() <= 0):
            raise ValueError("轴系刚度超出数值范围")
        scale = 1 / np.sqrt(reduced.diagonal())
        unit_loads = np.zeros((size, 6))
        unit_loads[self.gear] = np.eye(6)
        try:
            factor = cho_factor(scale[:, None] * reduced * scale)
            unit_displacements = cho_solve(factor, scale[:, None] * unit_loads[free])
        except np.linalg.LinAlgError as error:
            raise ValueError("轴系刚度非正定，请检查支承和轴段尺寸") from error
        self.unit_displacements = np.zeros((size, 6))
        self.unit_displacements[free] = scale[:, None] * unit_displacements
        gear_flexibility = self.unit_displacements[self.gear]
        self.gear_flexibility = (gear_flexibility + gear_flexibility.T) / 2

    def response(self, wrench):
        displacement = (self.unit_displacements @ wrench).reshape(-1, 6)
        bearings = self.shaft.bearings
        nodes = [bearing.node for bearing in bearings]
        reactions = -np.asarray([b.stiffness_n_per_m for b in bearings]) * displacement[nodes, :3]
        lever = np.zeros_like(reactions)
        lever[:, 2] = self.positions[nodes] - self.positions[self.shaft.gear_node]
        torque_reaction = -float(wrench[5])
        return {
            "node_positions_m": self.positions,
            "translation_m": displacement[:, :3],
            "rotation_rad": displacement[:, 3:],
            "gear_force_n": wrench[:3],
            "gear_moment_n_m": wrench[3:],
            "bearing_nodes": nodes,
            "bearing_reactions_n": reactions,
            "torque_reaction_n_m": torque_reaction,
            "force_balance_n": wrench[:3] + reactions.sum(axis=0),
            "moment_balance_n_m": wrench[3:] + np.cross(lever, reactions).sum(axis=0)
            + [0, 0, torque_reaction],
        }


def solve_shaft(shaft: Shaft, gear_wrench):
    """求单轴静态响应。齿轮节点载荷按 Fx、Fy、Fz、Mx、My、Mz 输入。"""
    if not isinstance(shaft, Shaft):
        raise ValueError("shaft 必须为 Shaft")
    wrench = np.asarray(gear_wrench, dtype=float)
    if wrench.shape != (6,) or not np.all(np.isfinite(wrench)):
        raise ValueError("齿轮截面载荷必须为六个有限数，单位 N 和 N·m")
    return _ShaftFlexibility(shaft).response(wrench)


class _ShaftCoupling:
    def __init__(self, pair, system):
        self.shafts = (_ShaftFlexibility(system.pinion), _ShaftFlexibility(system.wheel))
        alpha, beta = pair.working_pressure_angle_rad, pair.base_helix_angle_rad
        self.transverse_normal = np.array([sin(alpha), cos(alpha), 0.])
        # 接触线展长随 z 以 tanβ 变化，因此法向的轴向分量为 −sinβ。
        self.normal = self.transverse_normal * cos(beta) - [0., 0., sin(beta)]
        self.pitch_point_x = pair.pinion.base_radius_m / cos(alpha)
        self.pitch_roll = pair.pinion.base_radius_m * tan(alpha)
        self.wheel_center = np.array([
            -pair.center_distance_m if pair.wheel.internal else pair.center_distance_m,
            0., pair.axial_offset_m])

    def project(self, grid):
        # 每行是单位条带力产生的截面力和力矩，也是截面位移到间隙的虚功映射。
        position = ((grid.roll_pinion_m - self.pitch_roll)[:, None] * self.transverse_normal
                    + [self.pitch_point_x, 0., 0.])
        position[:, 2] = grid.face_coordinate_m
        normal = np.broadcast_to(self.normal, position.shape)
        pinion = np.column_stack((-normal, np.cross(position, -normal)))
        wheel = np.column_stack((normal, np.cross(position - self.wheel_center, normal)))
        flex = (pinion @ self.shafts[0].gear_flexibility @ pinion.T
                + wheel @ self.shafts[1].gear_flexibility @ wheel.T)
        return (pinion, wheel), (flex + flex.T) / 2

    def diagnostics(self, projections, forces):
        results = {}
        normal_displacement = np.zeros_like(forces)
        for name, shaft, projection in zip(("pinion", "wheel"), self.shafts, projections):
            wrench = projection.T @ forces
            results[name] = shaft.response(wrench)
            normal_displacement += projection @ (shaft.gear_flexibility @ wrench)
        results["normal_displacement_m"] = normal_displacement
        return results
