"""实际三角离散曲面的最近距离；完整BVH候选，不以顶点距离代替整片距离。

本模块按需调用VTK原生三角相交、点面距离及有限边段距离。数值结果
属于给定三角网。名义弦差可另扣作工程余量，不能因此称为连续CAD的
严格距离证明；若要后者，调用方必须另外给出真实Hausdorff误差上界。
"""
from __future__ import annotations

import math


def _vtk():
    try:
        import vtk
    except ImportError as error:
        raise ImportError('整片三角距离需要可选VTK依赖，请安装render3d扩展') from error
    return vtk


def _triangle(value):
    import numpy as np
    result = np.asarray(value, dtype=float)
    if result.shape != (3, 3) or not np.isfinite(result).all():
        raise ValueError('三角形须为有限3×3坐标，单位m')
    if float(np.linalg.norm(np.cross(result[1]-result[0], result[2]-result[0]))) <= 1e-30:
        raise ValueError('不能用退化三角形建立整片距离证据')
    return result


class _TrianglePair:
    """复用原生工作区；单对象不跨线程共享。"""
    def __init__(self):
        vtk = _vtk()
        self.vtk = vtk
        self.first = vtk.vtkTriangle()
        self.second = vtk.vtkTriangle()
        self.first.GetPoints().SetDataTypeToDouble()
        self.second.GetPoints().SetDataTypeToDouble()
        self.closest = [0., 0., 0.]
        self.parametric = [0., 0., 0.]
        self.weights = [0., 0., 0.]
        self.subid = vtk.reference(0)
        self.distance2 = vtk.reference(0.)
        self.segment_first = [0., 0., 0.]
        self.segment_second = [0., 0., 0.]
        self.t = vtk.reference(0.)
        self.s = vtk.reference(0.)

    def set_source(self, triangle):
        import numpy as np
        self.source_array = np.asarray(triangle, dtype=float)
        self.source = self.source_array.tolist()
        self.source_edge_length = float(np.linalg.norm(self.source_array-np.roll(self.source_array, 1, axis=0), axis=1).max())

    def distance2_to(self, triangle):
        import numpy as np
        # 原生相交代码的固定内部容差会把SI微小FE片误判为相交。
        # 双方共同平移并按最大边长归一化；距离平方最后还原到m²。
        scale = max(self.source_edge_length,
            float(np.linalg.norm(triangle-np.roll(triangle, 1, axis=0), axis=1).max()))
        if not math.isfinite(scale) or scale <= 0:
            raise ArithmeticError('有限三角对未获得有效归一化长度')
        origin = self.source_array[0]
        source = ((self.source_array-origin)/scale).tolist()
        target = ((triangle-origin)/scale).tolist()
        if self.vtk.vtkTriangle.TrianglesIntersect(*source, *target):
            return 0.
        for i, point in enumerate(source):
            self.first.GetPoints().SetPoint(i, point)
        for i, point in enumerate(target):
            self.second.GetPoints().SetPoint(i, point)
        # 另逐条检验有限边穿过对面内部；不能让相交快速分支的漏判
        # 变成仅顶点/边边正距离，从而高估分离下界。
        for cell, points in ((self.first, target), (self.second, source)):
            for i in range(3):
                if cell.IntersectWithLine(points[i], points[(i+1) % 3], 1e-12,
                        self.t, self.segment_first, self.parametric, self.subid):
                    return 0.
        best = math.inf
        for cell, points in ((self.first, target), (self.second, source)):
            for point in points:
                status = cell.EvaluatePosition(point, self.closest, self.subid,
                    self.parametric, self.distance2, self.weights)
                distance2 = float(self.distance2)
                if status < 0 or not math.isfinite(distance2) or distance2 < -1e-25:
                    raise ArithmeticError('VTK点面最近距离未获得有效有限结果')
                best = min(best, max(0., distance2))
        # 非相交三角形的最近点还可能同时位于两条边内部，必须检查九组有限边段。
        for a in range(3):
            for b in range(3):
                value = self.vtk.vtkLine.DistanceBetweenLineSegments(
                    source[a], source[(a+1) % 3], target[b], target[(b+1) % 3],
                    self.segment_first, self.segment_second, self.t, self.s)
                if not math.isfinite(value) or value < -1e-25:
                    raise ArithmeticError('VTK有限边段距离未获得有效有限结果')
                best = min(best, max(0., value))
        return best*scale*scale


def triangle_triangle_distance_m(first, second):
    """两完整有限三角形的数值最近距离m，含穿越相交与内部边边最近点。"""
    pair = _TrianglePair()
    pair.set_source(_triangle(first))
    return math.sqrt(pair.distance2_to(_triangle(second)))


class TriangleSurfaceDistance:
    """固定真实三角网的BVH距离引擎；同一实例不可并发跨线程调用。

    首先由源片重心求一个真实距离上界D。任一可能更近的目标三角形
    都必须与源片AABB向六向各扩展D的盒相交，因此完整查询该盒中的
    所有BVH单元后逐片求最近距离；不限制候选数，不构造凸包补孔洞。
    """
    def __init__(self, triangles_m, *, numerical_tolerance_m=1e-10, verify_bvh=False):
        import numpy as np
        vtk = _vtk()
        from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        triangles = np.asarray(triangles_m, dtype=float)
        if triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not len(triangles) or not np.isfinite(triangles).all():
            raise ValueError('对面须为非空有限N×3×3真实三角网，单位m')
        if np.any(np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]), axis=1) <= 1e-30):
            raise ValueError('对面含退化三角形，不能静默丢弃后认证最近距离')
        if (isinstance(numerical_tolerance_m, bool) or not isinstance(numerical_tolerance_m, (int, float)) or
                not math.isfinite(numerical_tolerance_m) or numerical_tolerance_m <= 0):
            raise ValueError('数值距离余量须为正有限m')
        self.triangles = triangles.copy()
        self.np = np
        self.vtk = vtk
        self.tolerance = numerical_tolerance_m
        self.verify_bvh = bool(verify_bvh)
        points = vtk.vtkPoints()
        points.SetData(numpy_to_vtk(triangles.reshape(-1, 3), deep=True, array_type=vtk.VTK_DOUBLE))
        cells = vtk.vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, 3*len(triangles)+1, 3, dtype=np.int64), deep=True),
                      numpy_to_vtkIdTypeArray(np.arange(3*len(triangles), dtype=np.int64), deep=True))
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(points)
        self.polydata.SetPolys(cells)
        self.locator = vtk.vtkStaticCellLocator()
        self.locator.SetDataSet(self.polydata)
        self.locator.BuildLocator()
        self.ids = vtk.vtkIdList()
        self.closest = [0., 0., 0.]
        self.cellid, self.subid = vtk.reference(0), vtk.reference(0)
        self.distance2 = vtk.reference(0.)
        self.pair = _TrianglePair()
        self.bounds_min = triangles.min(axis=1)
        self.bounds_max = triangles.max(axis=1)

    def query(self, source_triangle_m):
        """返回整片数值距离及完整候选证据，不把该值自动称为连续CAD下界。"""
        triangle = _triangle(source_triangle_m)
        np = self.np
        self.locator.FindClosestPoint(triangle.mean(axis=0), self.closest, self.cellid, self.subid, self.distance2)
        value = float(self.distance2)
        if not math.isfinite(value) or value < 0:
            raise ArithmeticError('源片重心的真实网格距离未获得有效上界')
        upper = math.sqrt(value)
        radius = upper+self.tolerance
        low, high = triangle.min(axis=0)-radius, triangle.max(axis=0)+radius
        bounds = tuple(float(x) for axis in range(3) for x in (low[axis], high[axis]))
        self.locator.FindCellsWithinBounds(bounds, self.ids)
        bucket_candidates = {int(self.ids.GetId(i)) for i in range(self.ids.GetNumberOfIds())}
        candidate_indices = np.asarray(sorted(bucket_candidates), dtype=np.int64)
        # 桶查询包含盒外单元；按真实单元AABB过滤，不改变完整候选集合。
        active = np.all(self.bounds_max[candidate_indices] >= low, axis=1) & np.all(self.bounds_min[candidate_indices] <= high, axis=1)
        candidates = set(candidate_indices[active].tolist())
        if int(self.cellid) not in candidates:
            raise ArithmeticError('BVH候选未包含构造上界的最近单元，拒绝给出不完整下界')
        if self.verify_bvh:
            required = set(np.flatnonzero(np.all(self.bounds_max >= low, axis=1) & np.all(self.bounds_min <= high, axis=1)).tolist())
            if not required.issubset(candidates):
                raise ArithmeticError('BVH未覆盖所有相交单元AABB候选')
        self.pair.set_source(triangle)
        best = upper*upper
        tested = 0
        for index in sorted(candidates):
            best = min(best, self.pair.distance2_to(self.triangles[index]))
            tested += 1
            if best <= self.tolerance*self.tolerance:
                best = 0.
                break
        distance = math.sqrt(max(0., best))
        return dict(distance_m=distance, centroid_distance_upper_bound_m=upper,
                    candidate_count=len(candidates), tested_candidate_count=tested,
                    bvh_bucket_candidate_count=len(bucket_candidates), primitive_pair_scale_normalized=True,
                    complete_bvh_candidates=True, all_aabb_candidates_independently_checked=self.verify_bvh,
                    numerical_tolerance_m=self.tolerance,
                    distance_scope='给定实际三角网的完整有限三角形距离，非仅顶点距离；不填补孔洞')


def engineering_clearance_margin(distance_m, relative_motion_m, *, source_allowance_m,
        target_allowance_m, numerical_tolerance_m=1e-10):
    """双面弦差参数只作显式工程余量；没有独立全域上界时不宣称连续CAD严格证明。"""
    values = (distance_m, relative_motion_m, source_allowance_m, target_allowance_m, numerical_tolerance_m)
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('距离、相对运动和两面余量须为非负有限m')
    lower = max(0., distance_m-source_allowance_m-target_allowance_m-numerical_tolerance_m)
    return dict(distance_lower_under_declared_allowances_m=lower,
                relative_motion_upper_bound_m=relative_motion_m,
                remaining_margin_m=lower-relative_motion_m,
                separated_under_declared_allowances=lower > relative_motion_m,
                source_allowance_m=source_allowance_m, target_allowance_m=target_allowance_m,
                allowance_evidence='调用方声明的工程弦差余量，可能仅为网格设置；并非本函数独立验证的连续曲面全域上界',
                strict_continuous_cad_proof=False)
