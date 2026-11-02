"""理想直齿球面渐开线的真实解析齿面网格及OBJ/JSON导出。

网格只含两侧有效渐开线。低于基圆的区域、实际根圆角、齿槽底面和
轮体实体均不填造；完整圆周是同一解析齿面的严格周期阵列。
宽向使用以共同锥顶为球心的球面端界；顶/根侧边界由宏观锥面与
这些球面的解析交线确定。球面端界与平面端盖是不同的建模约定。
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

from .geometry import GeometryResult
from .spherical_involute import SphericalInvoluteProfile, Vector3


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _cross(a, b):
    return a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]


def _sample_count(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 2:
        raise ValueError(f"{name}必须为不小于2的整数")


@dataclass(frozen=True)
class ToothSurfaceMesh:
    """三角网格，法向与三角形均朝材料外；不代表闭合实体。"""
    vertices_m: tuple[Vector3, ...]
    normals: tuple[Vector3, ...]
    triangles: tuple[tuple[int, int, int], ...]
    uv: tuple[tuple[float, float], ...]
    flank: tuple[str, ...]
    tooth_indices: tuple[int, ...]
    region_labels: tuple[str, ...]
    metadata: dict


@dataclass(frozen=True)
class SphericalToothSurface:
    """有效球面渐开线齿面；局部+Z沿轴从锥顶指向齿轮外端。

    齿厚中心位于XZ平面，节点+X是径向、+Y是切向；left位于负Y侧。
    root/tip reference heights沿节锥法线定义，锥度决定各圆锥面的斜率。
    """
    profile: SphericalInvoluteProfile
    teeth: int
    inner_cone_distance_m: float
    outer_cone_distance_m: float
    reference_cone_distance_m: float
    mean_addendum_m: float
    mean_dedendum_m: float
    addendum_taper_angle_rad: float
    dedendum_taper_angle_rad: float
    member: str = "pinion"
    requested_root_fillet_radius_m: float | None = None

    def __post_init__(self):
        if not isinstance(self.profile, SphericalInvoluteProfile):
            raise TypeError("齿廓必须为SphericalInvoluteProfile")
        if isinstance(self.teeth, bool) or not isinstance(self.teeth, int) or self.teeth < 3:
            raise ValueError("齿数必须为不小于3的整数")
        lengths = (self.inner_cone_distance_m, self.outer_cone_distance_m,
                   self.reference_cone_distance_m, self.mean_addendum_m, self.mean_dedendum_m)
        if any(not math.isfinite(v) or v <= 0 for v in lengths):
            raise ValueError("球面锥距及节点齿高必须为正有限长度m")
        if not self.inner_cone_distance_m < self.reference_cone_distance_m < self.outer_cone_distance_m:
            raise ValueError("参考球面锥距必须严格位于内、外球面端界之间")
        for angle in (self.addendum_taper_angle_rad, self.dedendum_taper_angle_rad):
            if not math.isfinite(angle) or not 0 <= angle < math.pi / 2:
                raise ValueError("齿高锥度必须在[0,π/2)内")
        if self.requested_root_fillet_radius_m is not None and (
            not math.isfinite(self.requested_root_fillet_radius_m) or self.requested_root_fillet_radius_m <= 0
        ):
            raise ValueError("给定的根圆角请求必须为正有限长度；真实过渡仍未实现")
        if self.profile.pitch_angular_tooth_thickness_rad >= 2 * math.pi / self.teeth:
            raise ValueError("节点完整齿厚角必须小于一个齿距角")
        for radius in (self.inner_cone_distance_m, self.reference_cone_distance_m, self.outer_cone_distance_m):
            self.profile_bounds(radius)

    def _boundary_polar_angle(self, radius: float, height: float, taper: float, sign: int) -> float:
        # 宏观锥面的截线为 h=A+B*s；s是节锥母线坐标，不是点的球面半径。
        # 联立 s²+h²=radius²，正根保持对应的节锥外向分支。
        b = math.tan(taper)
        a = height - self.reference_cone_distance_m * b
        discriminant = (1 + b*b) * radius*radius - a*a
        if discriminant <= 0:
            raise ValueError("宏观齿顶/根锥面与该球面没有有效交线")
        s = (-a*b + math.sqrt(discriminant)) / (1 + b*b)
        h = a + b*s
        if s <= 0 or h <= 0:
            raise ValueError("宏观齿顶/根边界落在非正母线坐标或非正齿高处")
        polar = self.profile.pitch_cone_angle_rad + sign * math.atan2(h, s)
        if not 0 < polar < math.pi / 2:
            raise ValueError("球面截交边界离开有效锐角锥体")
        return polar

    def physical_profile_bounds(self, cone_distance_m: float) -> tuple[float, float]:
        radius = float(cone_distance_m)
        if not math.isfinite(radius) or not self.inner_cone_distance_m <= radius <= self.outer_cone_distance_m:
            raise ValueError("球面锥距必须在指定内、外球面端界内，单位m")
        root = self._boundary_polar_angle(radius, self.mean_dedendum_m, self.dedendum_taper_angle_rad, -1)
        tip = self._boundary_polar_angle(radius, self.mean_addendum_m, self.addendum_taper_angle_rad, 1)
        return root, tip

    def profile_bounds(self, cone_distance_m: float) -> tuple[float, float]:
        """仅返回有效渐开线范围；基圆以下实体根部不补线或补圆角。"""
        root, tip = self.physical_profile_bounds(cone_distance_m)
        lower = max(root, self.profile.base_cone_angle_rad)
        if tip <= lower:
            raise ValueError("该截面不存在正长度有效球面渐开线")
        if self.profile.half_material_angle(tip) <= 0:
            raise ValueError("齿顶边界超出尖齿位置，不能生成相交的有效齿面")
        if 2 * self.profile.half_material_angle(lower) >= 2 * math.pi / self.teeth:
            raise ValueError("根侧有效渐开线齿厚已占满齿距，当前参数不形成分离齿槽")
        return lower, tip

    def evaluate(self, cone_distance_m: float, polar_angle_rad: float, flank: str = "right") -> Vector3:
        lower, upper = self.profile_bounds(cone_distance_m)
        if not math.isfinite(polar_angle_rad) or not lower - 2e-14 <= polar_angle_rad <= upper + 2e-14:
            raise ValueError("极角超出该截面的有效齿面边界；根过渡区未实现")
        return self.profile.evaluate(cone_distance_m, min(upper, max(lower, polar_angle_rad)), flank)

    def mesh(self, width_count: int = 17, profile_count: int = 41) -> ToothSurfaceMesh:
        """生成单齿两面。预览建议9×21；一般显示17×41；细查25×65。

        采样数不是精度保证；需要距离/曲率分析时应逐级加密并核验收敛。
        """
        _sample_count(width_count, "宽向采样数")
        _sample_count(profile_count, "齿廓采样数")
        if 2 * width_count * profile_count > 2_000_000:
            raise ValueError("单齿网格超过两百万点，请显式降低采样规模")
        vertices, normals, parameters, flanks, indices, regions, triangles = [], [], [], [], [], [], []
        missing_root = False
        for flank in ("left", "right"):
            start = len(vertices)
            for i in range(width_count):
                u = i / (width_count - 1)
                radius = self.inner_cone_distance_m + u * (self.outer_cone_distance_m - self.inner_cone_distance_m)
                # 固定端点可避免插值最后一位浮点舍入造成伪越界。
                radius = min(self.outer_cone_distance_m, max(self.inner_cone_distance_m, radius))
                root, _ = self.physical_profile_bounds(radius)
                low, high = self.profile_bounds(radius)
                missing_root = missing_root or root < self.profile.base_cone_angle_rad
                for j in range(profile_count):
                    v = j / (profile_count - 1)
                    polar = low + v * (high - low)
                    vertices.append(self.profile.evaluate(radius, polar, flank))
                    normals.append(self.profile.normal(polar, flank))
                    parameters.append((u, v))
                    flanks.append(flank)
                    indices.append(0)
                    regions.append("base_limit" if abs(polar - self.profile.base_cone_angle_rad) < 2e-13 else "spherical_involute")
            for i in range(width_count - 1):
                for j in range(profile_count - 1):
                    a = start + i * profile_count + j
                    b, c, d = a + profile_count, a + 1, a + profile_count + 1
                    for tri in ((a, b, c), (c, b, d)):
                        normal = _cross(_sub(vertices[tri[1]], vertices[tri[0]]), _sub(vertices[tri[2]], vertices[tri[0]]))
                        if _dot(normal, normals[tri[0]]) < 0:
                            tri = tri[0], tri[2], tri[1]
                        triangles.append(tri)
        metadata = {
            "units": "m", "angle_units": "rad", "member": self.member,
            "tooth_count": self.teeth, "represented_tooth_indices": [0],
            "width_count": width_count, "profile_count": profile_count,
            "coordinate_frame": "锥顶为原点；+Z从锥顶指向外端；齿厚中心在XZ平面，节点+X径向、+Y切向",
            "surface_model": "第一性原理理想直齿球面渐开线；两侧有效齿面",
            "width_boundary_model": "以共同锥顶为球心的内外球面端界，非平面端盖",
            "profile_boundary_model": "宏观齿顶/齿根圆锥面与球面截交；基圆以下不生成",
            "root_below_base_missing": missing_root,
            "requested_root_fillet_radius_m": self.requested_root_fillet_radius_m,
            "root_fillet_implemented": False, "closed_solid": False,
            "base_curvature_singular": any(region == "base_limit" for region in regions),
            "pitch_cone_angle_rad": self.profile.pitch_cone_angle_rad,
            "base_cone_angle_rad": self.profile.base_cone_angle_rad,
            "pressure_angle_rad": self.profile.pressure_angle_rad,
            "warnings": ["网格不含齿根圆角、槽底、轮体或实际刀具包络", "解析几何及采样不等于接触/承载求解或制造认证"],
        }
        return ToothSurfaceMesh(tuple(vertices), tuple(normals), tuple(triangles), tuple(parameters),
                                tuple(flanks), tuple(indices), tuple(regions), metadata)

    def mesh_full_gear(self, width_count: int = 17, profile_count: int = 41) -> ToothSurfaceMesh:
        return periodic_gear_mesh(self.mesh(width_count, profile_count), self.teeth)


def build_spherical_tooth_surface(geometry: GeometryResult, member: str = "pinion", *,
                                 include_root_transition: bool = False) -> SphericalToothSurface:
    """从已算宏观几何构造理想直齿齿面；绝不将弧齿扭转近似作真实齿面。"""
    if not isinstance(geometry, GeometryResult):
        raise TypeError("齿面入口需要GeometryResult")
    if member not in ("pinion", "wheel"):
        raise ValueError("轮成员只能为pinion或wheel")
    if geometry.inputs.spiral_angle_rad != 0 or geometry.diagnostics.get("tooth_trace") != "straight":
        raise NotImplementedError("实际弧齿/零度圆弧的刀具包络齿面尚未实现，不能用扭转的直齿面代替")
    if include_root_transition:
        raise NotImplementedError("真实齿根圆角和基圆以下过渡尚未实现，本接口仅提供有效球面渐开线")
    side = getattr(geometry, member)
    angular_thickness = side.mean_transverse_tooth_thickness_m / (side.mean_pitch_diameter_m / 2)
    profile = SphericalInvoluteProfile(side.pitch_cone_angle_rad, geometry.inputs.pressure_angle_rad, angular_thickness)
    return SphericalToothSurface(
        profile, side.teeth, geometry.inner_cone_distance_m, geometry.outer_cone_distance_m,
        geometry.mean_cone_distance_m, side.mean_addendum_m, side.mean_dedendum_m,
        side.addendum_taper_angle_rad, side.dedendum_taper_angle_rad, member,
        side.root_fillet_radius_m,
    )


def _rotate_z(vector: Vector3, angle: float) -> Vector3:
    cosine, sine = math.cos(angle), math.sin(angle)
    return cosine * vector[0] - sine * vector[1], sine * vector[0] + cosine * vector[1], vector[2]


def periodic_gear_mesh(single_tooth: ToothSurfaceMesh, tooth_count: int) -> ToothSurfaceMesh:
    """精确旋转同一解析齿面作周期阵列，不增加未实现的根部或实体。"""
    if isinstance(tooth_count, bool) or not isinstance(tooth_count, int) or tooth_count < 3:
        raise ValueError("完整圆周齿数必须为不小于3的整数")
    if set(single_tooth.tooth_indices) != {0}:
        raise ValueError("周期阵列入口必须是单齿网格，不能重复阵列完整齿轮")
    if tooth_count != single_tooth.metadata.get("tooth_count"):
        raise ValueError("阵列齿数必须与解析齿厚/齿距参数一致")
    count = len(single_tooth.vertices_m)
    if count * tooth_count > 2_000_000:
        raise ValueError("完整圆周网格超过两百万点，请减小采样数")
    vertices, normals, triangles, uv, flanks, teeth, regions = [], [], [], [], [], [], []
    for tooth in range(tooth_count):
        angle = 2 * math.pi * tooth / tooth_count
        vertices.extend(_rotate_z(p, angle) for p in single_tooth.vertices_m)
        normals.extend(_rotate_z(n, angle) for n in single_tooth.normals)
        triangles.extend(tuple(index + tooth * count for index in tri) for tri in single_tooth.triangles)
        uv.extend(single_tooth.uv)
        flanks.extend(single_tooth.flank)
        teeth.extend([tooth] * count)
        regions.extend(single_tooth.region_labels)
    metadata = dict(single_tooth.metadata)
    metadata.update(represented_tooth_indices=list(range(tooth_count)), periodic_array=True, closed_solid=False)
    return ToothSurfaceMesh(tuple(vertices), tuple(normals), tuple(triangles), tuple(uv), tuple(flanks),
                            tuple(teeth), tuple(regions), metadata)


def export_surface_obj(mesh: ToothSurfaceMesh, path) -> Path:
    """仅用标准库导出三角齿面及解析顶点法向，组名保留齿号与左右侧。"""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# 单位m；理想球面渐开线有效齿面；非闭合实体\n")
        stream.write("# 不含真实根圆角、槽底和轮体\n")
        for p in mesh.vertices_m:
            stream.write("v " + " ".join(format(v, ".17g") for v in p) + "\n")
        for n in mesh.normals:
            stream.write("vn " + " ".join(format(v, ".17g") for v in n) + "\n")
        group = None
        for triangle in mesh.triangles:
            index = triangle[0]
            current = mesh.tooth_indices[index], mesh.flank[index]
            if current != group:
                stream.write(f"g tooth_{current[0]}_{current[1]}\n")
                group = current
            stream.write("f " + " ".join(f"{i+1}//{i+1}" for i in triangle) + "\n")
    return target


def export_surface_json(mesh: ToothSurfaceMesh, path) -> Path:
    """保存完整离散点、法向、三角形、参数域、齿号与模型边界。"""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": mesh.metadata, "vertices_m": mesh.vertices_m, "normals": mesh.normals,
               "triangles": mesh.triangles, "uv": mesh.uv, "flank": mesh.flank,
               "tooth_indices": mesh.tooth_indices, "region_labels": mesh.region_labels}
    target.write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":")), encoding="utf-8")
    return target
