"""由应变能建立的固定齿根、变截面齿梁结构柔度。

这是独立新写的降阶模型。径向齿梁采用矩形截面，分别积分弯曲、
Timoshenko 剪切和轴向压缩能；不含接触斑局部压陷和轮体/轴承柔度。
齿宽耦合由明确的横向位移梯度剪切能建立，不以经验柔度相加。
内部单位为 m、N、Pa、rad，输出柔度为 m/N。
"""
from __future__ import annotations

from dataclasses import dataclass
from math import cos

import numpy as np
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class Material:
    """各向同性线弹性材料；杨氏模量 Pa，泊松比无量纲。"""

    youngs_modulus_pa: float = 210e9
    poisson_ratio: float = 0.3

    def __post_init__(self):
        if not np.isfinite(self.youngs_modulus_pa) or self.youngs_modulus_pa <= 0:
            raise ValueError("杨氏模量必须为有限正数，单位 Pa")
        if not np.isfinite(self.poisson_ratio) or not -1 < self.poisson_ratio < 0.5:
            raise ValueError("各向同性材料泊松比必须在 (-1, 0.5) 内")

    @property
    def shear_modulus_pa(self):
        """G=E/[2(1+ν)]，单位 Pa。"""
        return self.youngs_modulus_pa / (2 * (1 + self.poisson_ratio))

    @property
    def plane_strain_compliance_per_pa(self):
        """接触半空间的材料系数 (1−ν²)/E，单位 1/Pa。"""
        return (1 - self.poisson_ratio**2) / self.youngs_modulus_pa


def _beam_sections(gear, contact_radii, widths, cosine, material, points):
    """对每个条带以单位法向力积分能量，返回三项对角柔度。"""
    direction = -1.0 if gear.internal else 1.0
    heights = direction * (contact_radii - gear.root_radius_m)
    if np.any(heights <= 0) or not np.all(np.isfinite(heights)):
        raise ValueError("接触点必须位于固定齿根至齿顶的可用齿高内")
    root_to_tip = direction * (gear.tip_radius_m - gear.root_radius_m)
    if np.any(heights > root_to_tip * (1 + 1e-10)):
        raise ValueError("接触点超过齿顶，无法建立齿梁")

    s, ds = _radial_quadrature(gear, heights, points)
    radial = gear.root_radius_m + direction * s
    thickness = np.asarray(gear.tooth_thickness_m(radial), dtype=float)
    thickness = np.broadcast_to(thickness, s.shape)
    if np.any(thickness <= 0) or not np.all(np.isfinite(thickness)):
        raise ValueError("齿梁的每个积分截面都需要有限正齿厚")
    # widths 是接触线积分长度；乘 cosβb 得条带的实际轴向截面宽。
    axial_width = widths * cosine
    area = axial_width[:, None] * thickness
    second_moment = axial_width[:, None] * thickness**3 / 12
    e, g = material.youngs_modulus_pa, material.shear_modulus_pa
    cos_alpha = np.clip(gear.base_radius_m / contact_radii, 0, 1)
    sin_alpha_sq = 1 - cos_alpha**2
    # 法向载荷的横截面分量为 F cosβb，虚功再投影一次，故为 cos²βb。
    transverse_projection = cosine**2
    bending = transverse_projection * cos_alpha**2 * np.sum(
        ds * (heights[:, None] - s)**2 / (e * second_moment), axis=1)
    # 矩形截面的抛物线剪应力分布给出剪切修正系数 5/6。
    shear = transverse_projection * cos_alpha**2 * np.sum(
        ds / ((5 / 6) * g * area), axis=1)
    compression = transverse_projection * sin_alpha_sq * np.sum(ds / (e * area), axis=1)
    return heights, bending, shear, compression


def _radial_quadrature(gear, heights, points):
    """在基圆处拆积分区间，避免将齿根延拓和渐开线厚度跨段积分。"""
    heights = np.atleast_1d(heights)
    direction = -1.0 if gear.internal else 1.0
    split = np.clip(direction * (gear.base_radius_m - gear.root_radius_m), 0, heights)
    starts = np.stack((np.zeros_like(heights), split), axis=1)
    lengths = np.stack((split, heights - split), axis=1)
    nodes, weights = leggauss(points)
    s = starts[:, :, None] + lengths[:, :, None] * (nodes + 1) / 2
    ds = lengths[:, :, None] * weights / 2
    return s.reshape(len(heights), -1), ds.reshape(len(heights), -1)


def _root_shape(s, height):
    """固定端三次位移场，接触点 φ(h)=1；仅在共同齿高内使用。"""
    ratio = s / height
    return ratio**2 * (3 - ratio) / 2


def _one_gear_compliance(gear, contact_radii, widths, faces, groups, cosine,
                         material, mode, points):
    heights, bending, shear, compression = _beam_sections(
        gear, contact_radii, widths, cosine, material, points)
    diagonal = bending + shear + compression
    if np.any(diagonal <= 0) or not np.all(np.isfinite(diagonal)):
        raise ValueError("齿梁能量积分未产生有限正柔度")
    if mode == "independent":
        return np.diag(diagonal), (bending, shear, compression), 0

    stiffness = np.diag(1 / diagonal)
    direction = -1.0 if gear.internal else 1.0
    edges = 0
    for label in np.unique(groups):
        indices = np.flatnonzero(groups == label)
        indices = indices[np.argsort(faces[indices], kind="stable")]
        for first, second in zip(indices[:-1], indices[1:]):
            spacing = faces[second] - faces[first]
            if spacing <= 0:
                raise ValueError("同一对齿的轴向条带坐标必须互异")
            height = min(heights[first], heights[second])
            s, ds = (value[0] for value in _radial_quadrature(gear, [height], points))
            section = np.asarray(gear.tooth_thickness_m(
                gear.root_radius_m + direction * s), dtype=float)
            phi_a = _root_shape(s, heights[first])
            phi_b = _root_shape(s, heights[second])
            # 相邻节点间取轴向线性位移，物理横截面位移为 φ(s)u/cosβb。
            # U_edge=G/(2 Δz cos²βb)∫t(s)[φ_b(s)u_b−φ_a(s)u_a]²ds。
            # 每条边都是非负二次能；与正定局部弹簧能相加后 K 必正定。
            scale = material.shear_modulus_pa / (spacing * cosine**2)
            aa = scale * np.sum(ds * section * phi_a**2)
            bb = scale * np.sum(ds * section * phi_b**2)
            ab = -scale * np.sum(ds * section * phi_a * phi_b)
            stiffness[first, first] += aa
            stiffness[second, second] += bb
            stiffness[first, second] += ab
            stiffness[second, first] += ab
            edges += 1
    # 对称求解后只消除舍入级非对称，不作特征值裁剪或经验稳定项。
    np.linalg.cholesky(stiffness)
    compliance = np.linalg.solve(stiffness, np.eye(len(diagonal)))
    compliance = (compliance + compliance.T) / 2
    return compliance, (bending, shear, compression), edges


def build_structural_compliance(pair, grid, material_pinion=Material(),
                                material_wheel=Material(), mode="coupled", *,
                                quadrature_points=32):
    """建立两轮共同法向结构柔度 C，使 u=C F，C 的单位为 m/N。

    independent：各轴向条带独立的变截面径向齿梁。
    coupled：另加沿齿宽位移梯度的剪切能；自由齿宽边缘为自然边界，
    不同 tooth_pair_index 的齿梁之间不耦合。共同齿高以两条带较小接触
    高度为界，未加载齿顶部分的翘曲没有建模。

    矩形截面的厚度采用弧齿厚近似，齿根固定。基圆以下的几何延拓由
    Gear 明确提供；这里不恢复真实过渡曲线、不计根部圆角应力集中。
    螺旋角只保留横截面力及虚功的 cos²βb 投影，忽略轴向力导致的
    齿/轮体弯曲。内圈齿梁从外侧齿根向内积分，压缩力取压缩能正值。
    该降阶模型不等于三维有限元，也不包含 Hertz 局部接触柔度。
    """
    if mode not in ("independent", "coupled"):
        raise ValueError("弹性模式必须为 independent 或 coupled")
    if not isinstance(quadrature_points, int) or not 4 <= quadrature_points <= 128:
        raise ValueError("高斯积分点数必须为 4 至 128 的整数")
    widths = np.asarray(grid.widths_m, dtype=float)
    faces = np.asarray(grid.face_coordinate_m, dtype=float)
    groups = np.asarray(grid.tooth_pair_index)
    roll_a = np.asarray(grid.roll_pinion_m, dtype=float)
    roll_b = np.asarray(grid.roll_wheel_m, dtype=float)
    arrays = (widths, faces, groups, roll_a, roll_b)
    if widths.ndim != 1 or len(widths) == 0 or any(a.shape != widths.shape for a in arrays):
        raise ValueError("接触网格各数组必须为同长的非空一维数组")
    if not all(np.all(np.isfinite(a)) for a in arrays) or np.any(widths <= 0):
        raise ValueError("接触网格必须有限，条带长度必须为正")
    if np.any(roll_a < 0) or np.any(roll_b < 0):
        raise ValueError("齿面展长必须为非负值")
    cosine = cos(pair.base_helix_angle_rad)
    if not np.isfinite(cosine) or cosine <= 0:
        raise ValueError("基圆螺旋角的余弦必须为有限正数")
    radius_a = np.hypot(pair.pinion.base_radius_m, roll_a)
    radius_b = np.hypot(pair.wheel.base_radius_m, roll_b)
    result_a, parts_a, edges_a = _one_gear_compliance(
        pair.pinion, radius_a, widths, faces, groups, cosine, material_pinion, mode,
        quadrature_points)
    result_b, parts_b, edges_b = _one_gear_compliance(
        pair.wheel, radius_b, widths, faces, groups, cosine, material_wheel, mode,
        quadrature_points)
    compliance = result_a + result_b
    np.linalg.cholesky(compliance)
    diagnostics = {
        "mode": mode,
        "boundary": "两轮齿根固定；齿宽两端为自然边界",
        "local_contact_included": False,
        "normal_projection": "横截面梁柔度乘 cos²βb；未计轴向力导致的轮体弯曲",
        "section_model": "变弧齿厚矩形截面；基圆以下由几何模块给出径向延拓",
        "coupling_energy": "G/(2Δz cos²βb)∫t(s)(φ_b u_b−φ_a u_a)²ds",
        "coupled_edges_per_gear": [edges_a, edges_b],
        "quadrature_points_per_segment": quadrature_points,
        "quadrature_segments": "在基圆处分开齿根延拓段与渐开线段",
        "pinion_bending_m_per_n": parts_a[0].tolist(),
        "pinion_shear_m_per_n": parts_a[1].tolist(),
        "pinion_compression_m_per_n": parts_a[2].tolist(),
        "wheel_bending_m_per_n": parts_b[0].tolist(),
        "wheel_shear_m_per_n": parts_b[1].tolist(),
        "wheel_compression_m_per_n": parts_b[2].tolist(),
    }
    return compliance, diagnostics
