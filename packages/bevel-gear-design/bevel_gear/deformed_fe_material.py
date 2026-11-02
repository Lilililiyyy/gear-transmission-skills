"""真实线性四面体变形闭边界的材料域审核；不重建刚度或求解位移。

本模块以给定有限元三角边界为材料域。CAD点通过原网格最近边界迹
携带同源P1位移；CAD与FE迹之间的几何偏差另作完整三角形上界。
名义CAD弦差属于显式工程余量，不能冒称连续NURBS严格Hausdorff界。
"""
from __future__ import annotations

import math
import hashlib


def closed_boundary_evidence(nodes_m, boundary_triangles):
    """以共享节点索引逐条核验闭边界；不以重复坐标拼接或补洞。"""
    import numpy as np
    from .distributed_pressure_domain import numeric_array_identity
    nodes = np.asarray(nodes_m, dtype=float)
    triangles = np.asarray(boundary_triangles)
    if (nodes.ndim != 2 or nodes.shape[1:] != (3,) or not np.isfinite(nodes).all()
            or triangles.ndim != 2 or triangles.shape[1:] != (3,) or not len(triangles)
            or not np.issubdtype(triangles.dtype, np.integer)):
        raise ValueError('闭边界需要有限N×3节点和非空整数M×3共享节点索引')
    if triangles.min() < 0 or triangles.max() >= len(nodes):
        raise ValueError('闭边界三角形引用了不存在的节点')
    canonical = np.sort(triangles, axis=1)
    if np.any(canonical[:, 1:] == canonical[:, :-1]):
        raise ValueError('闭边界含重复节点退化三角形')
    unique = np.unique(canonical, axis=0)
    if len(unique) != len(triangles):
        raise ValueError('闭边界含重复三角形')
    xyz = nodes[triangles]
    if np.any(np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1) <= 1e-30):
        raise ValueError('闭边界含几何退化三角形')
    edges = np.sort(np.concatenate((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]])), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    if np.any(counts != 2):
        raise ValueError('材料边界未闭合或非流形：每条共享索引边必须恰好属于两片')
    return dict(passed=True, triangle_count=len(triangles), edge_count=len(counts),
        shared_index_edge_incidence_exactly_two=True, duplicate_triangles=0,
        nodes_identity=numeric_array_identity(nodes),
        boundary_identity=numeric_array_identity(triangles, integer=True),
        definition='原共享节点索引闭合核验；不焊接、不封孔、不构造凸包')


def deformation_jacobian_evidence(nodes_m, tetrahedra, displacements_m, *, chunk_size=20000):
    """分块复核原四面体到实际变形四面体的Jacobian比值，无K或分解。"""
    import numpy as np
    nodes, displacement = np.asarray(nodes_m), np.asarray(displacements_m)
    tets = np.asarray(tetrahedra)
    if nodes.shape != displacement.shape or tets.ndim != 2 or tets.shape[1:] != (4,):
        raise ValueError('变形Jacobian需要同节点实际u和四节点四面体')
    low, high = math.inf, -math.inf
    for first in range(0, len(tets), chunk_size):
        ids = tets[first:first+chunk_size]
        x = nodes[ids]; moved = x+displacement[ids]
        original = np.linalg.det(x[:, 1:]-x[:, :1])
        current = np.linalg.det(moved[:, 1:]-moved[:, :1])
        if np.any(original == 0) or not np.isfinite(original).all() or not np.isfinite(current).all():
            raise ValueError('原或变形四面体Jacobian退化/非有限')
        ratio = current/original
        low, high = min(low, float(ratio.min())), max(high, float(ratio.max()))
    return dict(passed=low > 0, minimum_determinant_ratio=low, maximum_determinant_ratio=high,
        tetrahedron_count=len(tets), factorization_performed=False,
        global_self_intersection_proof=False,
        definition='同源u下每个原四面体有向体积比为正；不因此声称任意全局自交严格定理')


def tetrahedral_boundary_evidence(tetrahedra, boundary_triangles, *, chunk_size=100000):
    """核验给定闭边界恰为原四面体单次出现的全部面，使用紧凑64位键。

    不建立百万项Python字典；每个面三个排序节点索引打包为一个64位
    整数，完整排序计数。当前实现最多支持2^21个节点索引范围，超限
    明确拒绝而非截断。内部面须恰有两个单元共享。
    """
    import numpy as np
    from .distributed_pressure_domain import numeric_array_identity
    tets = np.asarray(tetrahedra); boundary = np.asarray(boundary_triangles)
    if (tets.ndim != 2 or tets.shape[1:] != (4,) or not len(tets) or
            boundary.ndim != 2 or boundary.shape[1:] != (3,) or not len(boundary)
            or not np.issubdtype(tets.dtype, np.integer) or not np.issubdtype(boundary.dtype, np.integer)
            or tets.min() < 0 or boundary.min() < 0):
        raise ValueError('四面体边界核验需要非空非负整数单元与有限三角面索引')
    maximum = max(int(tets.max()), int(boundary.max())); bits = max(1, maximum.bit_length())
    if bits*3 > 63:
        raise ValueError('紧凑完整面计数当前支持节点索引小于2^21，不截断更大索引')
    def packed(faces):
        ordered = np.sort(faces, axis=1).astype(np.uint64)
        return (ordered[:, 0] << (2*bits)) | (ordered[:, 1] << bits) | ordered[:, 2]
    keys = np.empty(4*len(tets), dtype=np.uint64)
    table = np.asarray([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    for first in range(0, len(tets), chunk_size):
        value = tets[first:first+chunk_size]
        keys[4*first:4*(first+len(value))] = packed(value[:, table].reshape(-1, 3))
    keys.sort()
    starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1])+1]
    counts = np.diff(np.r_[starts, len(keys)])
    if np.any((counts != 1) & (counts != 2)):
        raise ValueError('原四面体面邻接并非一个或两个单元，不能当作有效材料网格')
    actual = keys[starts[counts == 1]]
    supplied = np.sort(packed(boundary))
    passed = bool(np.array_equal(actual, supplied))
    return dict(passed=passed, tetrahedron_count=len(tets), actual_boundary_triangle_count=len(actual),
        supplied_boundary_triangle_count=len(boundary), internal_face_count=int((counts == 2).sum()),
        all_tetrahedron_faces_counted=True, exact_shared_node_indices=True, packed_index_bits=bits,
        tetrahedra_identity=numeric_array_identity(tets, integer=True),
        boundary_identity=numeric_array_identity(boundary, integer=True),
        definition='全部4N原单元面排序计数；单次面集合与提供边界逐整数完全一致')


class ClosedTriangleMaterial:
    """完整有限FE闭边界的整片距离和材料分类；按需加载VTK。"""
    def __init__(self, nodes_m, boundary_triangles, *, tolerance_m=1e-10):
        import numpy as np
        from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        from .triangle_surface_distance import TriangleSurfaceDistance, _vtk
        self.closure = closed_boundary_evidence(nodes_m, boundary_triangles)
        self.nodes = np.asarray(nodes_m, dtype=float)
        self.indices = np.asarray(boundary_triangles, dtype=np.int64)
        self.tolerance = tolerance_m
        self.distance = TriangleSurfaceDistance(self.nodes[self.indices], numerical_tolerance_m=tolerance_m)
        vtk = _vtk()
        points = vtk.vtkPoints()
        points.SetData(numpy_to_vtk(self.nodes, deep=True, array_type=vtk.VTK_DOUBLE))
        cells = vtk.vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, 3*len(self.indices)+1, 3, dtype=np.int64), deep=True),
                      numpy_to_vtkIdTypeArray(self.indices.reshape(-1), deep=True))
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(points); self.polydata.SetPolys(cells)
        self.classifier = vtk.vtkSelectEnclosedPoints()
        # 此容差为相对包围盒尺度；整片距离的SI容差另行扣除。
        self.classifier.SetTolerance(1e-10)
        self.classifier.Initialize(self.polydata)

    def classify(self, point_m):
        """返回有限闭体IN/OUT；边界邻域须由调用方用距离排除。"""
        return 'IN' if self.classifier.IsInsideSurface(point_m) else 'OUT'

    def triangle_evidence(self, triangle_m):
        """整片正距离结合单点OUT，证明连通三角片全在有限FE材料外。"""
        import numpy as np
        triangle = np.asarray(triangle_m, dtype=float)
        distance = self.distance.query(triangle)
        state = self.classify(triangle.mean(axis=0)) if distance['distance_m'] > self.tolerance else 'ON_OR_INTERSECTING'
        return dict(**distance, material_state=state, outside=state == 'OUT',
            closed_shared_index_boundary=True)


def _nearest(engine, point):
    engine.locator.FindClosestPoint(point, engine.closest, engine.cellid, engine.subid, engine.distance2)
    value = float(engine.distance2)
    if not math.isfinite(value) or value < 0:
        raise ArithmeticError('原FE迹最近点未得到有效距离')
    return math.sqrt(value), int(engine.cellid)


def _point_triangle_distance(engine, point, triangle):
    cell = engine.pair.second
    for i, value in enumerate(triangle):
        cell.GetPoints().SetPoint(i, value)
    pair = engine.pair
    status = cell.EvaluatePosition(point, pair.closest, pair.subid, pair.parametric, pair.distance2, pair.weights)
    value = float(pair.distance2)
    if status < 0 or not math.isfinite(value) or value < -1e-25:
        raise ArithmeticError('凸距离上界的有限三角最近点计算失败')
    return math.sqrt(max(0., value))


def _split(triangle):
    import numpy as np
    a, b, c = triangle
    ab, bc, ca = (a+b)*.5, (b+c)*.5, (c+a)*.5
    return tuple(np.asarray(points) for points in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))


def canonical_midpoint_leaf(vertices_world_m, root_material_m, rotation_rows, depth, *, tolerance_m=1e-10):
    """恢复四分树整数重心坐标材料叶片，供跨phase静态几何缓存。"""
    import numpy as np
    root, world, rotation = (np.asarray(x, dtype=float) for x in (root_material_m, vertices_world_m, rotation_rows))
    if root.shape != (3, 3) or world.shape != (3, 3) or rotation.shape != (3, 3) or not 0 <= depth <= 30:
        raise ValueError('原CAD三角形、世界叶片、旋转与四分树深度不合法')
    basis = (root[1:]-root[0]).T
    xy = (world@rotation-root[0])@basis@np.linalg.inv(basis.T@basis)
    bary = np.column_stack((1.-xy.sum(axis=1), xy)); scale = 2**depth
    integral = np.rint(bary*scale)
    if np.any(integral < 0) or np.any(integral.sum(axis=1) != scale):
        raise ValueError('静态缓存叶片不能恢复为原四分树整数重心坐标')
    canonical = (integral/scale)@root
    if float(np.linalg.norm(canonical-world@rotation, axis=1).max()) > tolerance_m:
        raise ValueError('静态缓存叶片与当前实际CAD分片不一致')
    return canonical


def projected_surface_cover_upper_bound(triangle_m, target, known_distance_upper_m):
    """实际FE有限片沿源法向的构造性距离上界；必须完整覆盖源三角。

    对每个有限FE片Q，在源三角平面内取其真实正投影。交叠多边形中
    x沿源法向到Q的交点距离是仿射函数的绝对值，最大值位于交叠顶点。
    只有这些有限投影片的GEOS并集完整包含源三角才认可；不buffer、
    不补孔洞、不用采样面积代替集合覆盖。该构造给出一个真实可选Q点，
    因而是最近距离的上界，而非近似最近点或下界。
    """
    import numpy as np
    try:
        from shapely.geometry import Polygon
        from shapely.ops import unary_union
    except ImportError:
        return dict(passed=False, reason='缺少可选Shapely/GEOS完整平面覆盖引擎')
    triangle = np.asarray(triangle_m, dtype=float)
    origin = triangle[0]; first = triangle[1]-origin
    first /= np.linalg.norm(first)
    normal = np.cross(first, triangle[2]-origin); normal /= np.linalg.norm(normal)
    basis = np.asarray([first, np.cross(normal, first)])
    source = Polygon((triangle-origin)@basis.T)
    if not source.is_valid or source.area <= 0:
        return dict(passed=False, reason='源三角投影无效')
    expansion = float(known_distance_upper_m)+target.tolerance
    low, high = triangle.min(axis=0)-expansion, triangle.max(axis=0)+expansion
    bounds = tuple(float(x) for axis in range(3) for x in (low[axis], high[axis]))
    target.locator.FindCellsWithinBounds(bounds, target.ids)
    indices = sorted({int(target.ids.GetId(i)) for i in range(target.ids.GetNumberOfIds())})
    indices = [i for i in indices if bool(np.all(target.bounds_max[i] >= low) and np.all(target.bounds_min[i] <= high))]
    pieces = []
    for index in indices:
        value = target.triangles[index]
        n = np.cross(value[1]-value[0], value[2]-value[0]); denominator = float(n@normal)
        if abs(denominator) <= float(np.linalg.norm(n))*1e-6:
            continue
        projected = Polygon((value-origin)@basis.T)
        if not projected.is_valid or projected.area <= 0: continue
        overlap = source.intersection(projected)
        if overlap.is_empty or overlap.area <= 0: continue
        if overlap.geom_type != 'Polygon' or overlap.interiors:
            return dict(passed=False, reason='两个真实凸三角投影交集未得到无孔凸多边形')
        xy = np.asarray(overlap.exterior.coords)
        actual = origin+xy@basis
        upper = float(np.abs((value[0]-actual)@n/denominator).max())+target.tolerance
        # 覆盖并集用原共享节点FE片的完整投影，避免各片独立裁到T边界
        # 后引入浮点构造缝；距离上界仍由其与T的真实交叠多边形计算。
        pieces.append((upper, projected, index, denominator/float(np.linalg.norm(n))))
    if not pieces:
        return dict(passed=False, reason='未找到有限FE面片真实投影覆盖')
    thresholds = sorted({p[0] for p in pieces})
    def covered(threshold):
        # 不容忍面积小的孔洞；covers检验整个源集合（含边界）。
        return unary_union([p[1] for p in pieces if p[0] <= threshold]).covers(source)
    if not covered(thresholds[-1]):
        return dict(passed=False, reason='真实有限FE片投影并集未完整覆盖源三角，保留原保守上界',
            candidate_count=len(indices), projected_piece_count=len(pieces))
    left, right = 0, len(thresholds)-1
    while left < right:
        middle = (left+right)//2
        if covered(thresholds[middle]): right = middle
        else: left = middle+1
    upper = thresholds[left]
    used = [dict(target_triangle_index=i, distance_upper_bound_m=u, normal_direction_cosine=c)
        for u, _, i, c in pieces if u <= upper]
    return dict(passed=True, distance_upper_bound_m=upper, full_source_triangle_set_covered=True,
        candidate_count=len(indices), projected_piece_count=len(pieces), selected_pieces=used,
        buffer_or_convex_hull_applied=False,
        definition='实际有限FE投影片完整覆盖源T；每交叠多边形上真实法向射线长度的顶点仿射上界')


def triangle_to_surface_upper_bound(triangle_m, target, *, accuracy_m=5e-7, maximum_depth=7,
        maximum_leaves=1024):
    """完整源三角→目标有限三角网距离的上界，包含穿越FE片边的情况。

    对任一固定目标三角Q，dist(x,Q)为凸函数，故源三角T内的最大值
    不超过三顶点dist(v,Q)最大值。选择真实目标片即可得到有效上界；
    中点四分只收紧上界，达到资源上限也不会以采样最大值冒充上界。
    """
    import heapq
    import numpy as np
    triangle = np.asarray(triangle_m, dtype=float)
    if triangle.shape != (3, 3) or not np.isfinite(triangle).all():
        raise ValueError('几何桥接需要有限3×3三角形m')
    if accuracy_m <= 0 or maximum_depth < 0 or maximum_leaves < 1:
        raise ValueError('几何桥接精度及资源上限不合法')
    serial = 0; evaluated = 0; observed = 0.
    def bound(value):
        nonlocal observed, evaluated
        queries = [_nearest(target, p) for p in (*value, value.mean(axis=0))]
        observed = max(observed, *(q[0] for q in queries))
        candidates = sorted({q[1] for q in queries})
        upper = min(max(_point_triangle_distance(target, p, target.triangles[k]) for p in value) for k in candidates)
        evaluated += 1
        return upper+target.tolerance
    heap = [(-bound(triangle), serial, 0, triangle)]
    # 先尝试真实FE片的构造性完整投影覆盖；跨片长CAD三角往往可在
    # 一次GEOS覆盖中给出紧界，避免先生成上千个四分叶。失败才走原树。
    early_cover = None
    if -heap[0][0]-observed > accuracy_m:
        early_cover = projected_surface_cover_upper_bound(triangle, target, -heap[0][0])
        if early_cover['passed']:
            upper = min(-heap[0][0], early_cover['distance_upper_bound_m'])
            if upper < observed-target.tolerance*4:
                raise ArithmeticError('原片完整FE投影覆盖界低于真实观察距离')
            return dict(distance_upper_bound_m=upper, observed_distance_lower_bound_m=observed,
                remaining_bound_gap_m=max(0., upper-observed), terminal_leaf_count=1,
                evaluated_region_count=evaluated, maximum_used_depth=0,
                projected_FE_cover=early_cover, valid_upper_even_if_not_tight=True,
                target_complete_surface=True,
                formula='真实有限FE片沿源法向构造距离；GEOS投影并集完整覆盖T；交叠顶点仿射上界')
    while -heap[0][0]-observed > accuracy_m and len(heap)+3 <= maximum_leaves:
        item = heapq.heappop(heap)
        if item[2] >= maximum_depth:
            heapq.heappush(heap, item); break
        for child in _split(item[3]):
            serial += 1
            heapq.heappush(heap, (-bound(child), serial, item[2]+1, child))
    upper = -heap[0][0]
    cover = None
    if upper-observed > accuracy_m:
        cover = projected_surface_cover_upper_bound(triangle, target, upper)
        if cover['passed']:
            upper = min(upper, cover['distance_upper_bound_m'])
            if upper < observed-target.tolerance*4:
                raise ArithmeticError('完整FE投影覆盖上界低于已观察实际距离，拒绝不一致几何证据')
    return dict(distance_upper_bound_m=upper, observed_distance_lower_bound_m=observed,
        remaining_bound_gap_m=max(0., upper-observed), terminal_leaf_count=len(heap),
        evaluated_region_count=evaluated, maximum_used_depth=max(item[2] for item in heap),
        projected_FE_cover=cover,
        valid_upper_even_if_not_tight=True, target_complete_surface=True,
        formula='min_Q max_顶点 dist(v,Q)；dist(x,Q)凸；完整中点四分覆盖源三角')


def complete_trace_carriers(triangle_m, target, *, source_allowance_m=0., proven_distance_upper_m=None):
    """可能承载整个源片最近点迹的所有原FE面片；完整BVH盒候选。"""
    import numpy as np
    triangle = np.asarray(triangle_m, dtype=float)
    distance, closest = _nearest(target, triangle.mean(axis=0))
    radius = float(np.linalg.norm(triangle-triangle.mean(axis=0), axis=1).max())
    # dist(x,FE)的1-Lipschitz性给出每个源点的最近点距离上界。
    upper = distance+radius if proven_distance_upper_m is None else float(proven_distance_upper_m)
    if not math.isfinite(upper) or upper < distance-target.tolerance:
        raise ValueError('最近迹候选的已证明整片距离上界不合法')
    expansion = upper+2*source_allowance_m+2*target.tolerance
    low, high = triangle.min(axis=0)-expansion, triangle.max(axis=0)+expansion
    bounds = tuple(float(x) for axis in range(3) for x in (low[axis], high[axis]))
    target.locator.FindCellsWithinBounds(bounds, target.ids)
    indices = sorted({int(target.ids.GetId(i)) for i in range(target.ids.GetNumberOfIds())})
    # VTK桶查询可返回盒外单元；逐单元真实AABB过滤只移除不相交项。
    indices = [i for i in indices if bool(np.all(target.bounds_max[i] >= low)
        and np.all(target.bounds_min[i] <= high))]
    if closest not in indices:
        raise ArithmeticError('FE最近迹候选遗漏了构造上界的原面片')
    return dict(indices=indices, complete_bvh_candidates=True, expansion_m=expansion,
        query_low_m=low.tolist(), query_high_m=high.tolist(),
        source_allowance_m=source_allowance_m, nearest_centroid_triangle=closest)


def clip_triangle_to_box_with_weights(triangle_m, low_m, high_m):
    """真实有限FE三角与完整投影候选盒求交，保留原P1重心权重。"""
    import numpy as np
    triangle, low, high = (np.asarray(x, dtype=float) for x in (triangle_m, low_m, high_m))
    if triangle.shape != (3, 3) or low.shape != (3,) or high.shape != (3,) or np.any(low >= high):
        raise ValueError('FE三角裁剪需要有限三角与严格正尺寸盒')
    polygon = [(triangle[i], np.eye(3)[i]) for i in range(3)]
    for axis in range(3):
        for side, limit in ((1., low[axis]), (-1., high[axis])):
            if not polygon: return []
            output = []; previous = polygon[-1]; old_distance = side*(previous[0][axis]-limit)
            for current in polygon:
                distance = side*(current[0][axis]-limit)
                if (old_distance >= 0.) != (distance >= 0.):
                    t = old_distance/(old_distance-distance)
                    output.append((previous[0]+t*(current[0]-previous[0]), previous[1]+t*(current[1]-previous[1])))
                if distance >= 0.: output.append(current)
                previous, old_distance = current, distance
            polygon = output
    # 只删除逐字节相同的相邻重复点，不以几何容差焊接不同材料点。
    result = []
    for point, weights in polygon:
        if not result or not np.array_equal(point, result[-1][0]): result.append((point, weights))
    if len(result) > 1 and np.array_equal(result[0][0], result[-1][0]): result.pop()
    return result


class DeformedFEMaterialAudit:
    """同源真实u的有限FE全材料远域审核上下文，无刚度装配与分解。

    target_cad_context为DistributedPressureDomainAudit。必须有完整整轮
    CAD候选盒；所有潜在接触目标CAD片都逐片核对CAD→原FE迹上界。
    """
    def __init__(self, meshes, result, *, response_evidence, target_cad_context,
            bridge_accuracy_m=5e-7, bridge_maximum_depth=7, bridge_maximum_leaves=1024,
            static_geometry_cache=None):
        import numpy as np
        from .distributed_pressure_domain import (distributed_body_identity, distributed_response_fingerprint,
            numeric_array_identity)
        from .triangle_surface_distance import TriangleSurfaceDistance
        from .pressure_background import _digest
        from pathlib import Path
        from . import triangle_surface_distance as distance_module
        fingerprint = distributed_response_fingerprint(result)
        if response_evidence.get('passed') is not True or response_evidence.get('full_response_fingerprint') != fingerprint:
            raise ValueError('变形材料域必须绑定同一次通过W/u/F/Ku核验的P0响应')
        if 'body_fields' not in result:
            raise ValueError('变形FE材料域需要原实际位移数组，不能只有全局范数')
        identity = distributed_body_identity(result)
        self.context = target_cad_context; self.result = result; self.response = response_evidence
        self.tolerance = target_cad_context.tolerance
        self.bridge_options = dict(accuracy_m=bridge_accuracy_m, maximum_depth=bridge_maximum_depth,
            maximum_leaves=bridge_maximum_leaves)
        # 缓存只承载静态原网格/CAD几何，不能保存其他相位的u、距离或IN/OUT。
        geometry_key = _digest(dict(meshes={m: dict(nodes=numeric_array_identity(mesh.nodes_m),
            tetrahedra=numeric_array_identity(mesh.tetrahedra, integer=True),
            boundary=numeric_array_identity(mesh.boundary_triangles, integer=True)) for m, mesh in meshes.items()},
            CAD=target_cad_context.source_geometry_identity, target=target_cad_context.target_identity,
            source_chord_setting_m=target_cad_context.source_deflection,
            target_chord_setting_m=target_cad_context.target_deflection,
            complete_target_box=target_cad_context.complete_box,
            implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            distance_engine_source_sha256=hashlib.sha256(Path(distance_module.__file__).read_bytes()).hexdigest(),
            bridge_options=self.bridge_options, tolerance_m=self.tolerance))
        if static_geometry_cache is None:
            static_geometry_cache = {}
        if not isinstance(static_geometry_cache, dict):
            raise TypeError('静态几何缓存须为调用方明确持有的字典')
        cache = static_geometry_cache.setdefault(geometry_key, dict(original={}, closures={}, target_bridge={}, source_trace={},
            tetrahedral_boundary={}))
        self.static_cache = cache; self.static_geometry_key = geometry_key
        self.original = {}; self.meshes = meshes; self.u = {}; self.closures = {}; jacobians = {}
        for member, mesh in meshes.items():
            check = response_evidence['body_checks'][member]
            if numeric_array_identity(mesh.nodes_m) != check['mesh_nodes_identity']:
                raise ValueError('变形材料边界的原节点与真实响应不一致：'+member)
            mechanical = check.get('mechanical_recheck', {})
            if (check.get('mesh_tetrahedra_identity') or mechanical.get('mesh_tetrahedra_identity')) != numeric_array_identity(mesh.tetrahedra, integer=True):
                raise ValueError('变形材料边界的原四面体与Ku证据不一致：'+member)
            self.u[member] = np.asarray(result['body_fields'][member]['displacements_m'], dtype=float)
            jacobians[member] = deformation_jacobian_evidence(mesh.nodes_m, mesh.tetrahedra, self.u[member])
            if not jacobians[member]['passed']:
                raise ValueError('实际位移使原四面体翻转，不能认证材料零域：'+member)
            if member not in cache['closures']:
                cache['closures'][member] = closed_boundary_evidence(mesh.nodes_m, mesh.boundary_triangles)
            if member not in cache['tetrahedral_boundary']:
                cache['tetrahedral_boundary'][member] = tetrahedral_boundary_evidence(mesh.tetrahedra, mesh.boundary_triangles)
            if not cache['tetrahedral_boundary'][member]['passed']:
                raise ValueError('提供的闭边界并非原四面体材料完整边界：'+member)
            self.closures[member] = cache['closures'][member]
            if member not in cache['original']:
                cache['original'][member] = TriangleSurfaceDistance(np.asarray(mesh.nodes_m)[np.asarray(mesh.boundary_triangles)],
                    numerical_tolerance_m=self.tolerance)
            self.original[member] = cache['original'][member]
        frames = self.context._loaded_frames(result)
        self.rotations = {m: np.asarray(frames[m+'_rotation_rows']) for m in meshes}
        wheel = meshes['wheel']
        self.target_material = ClosedTriangleMaterial(np.asarray(wheel.nodes_m)+self.u['wheel'], wheel.boundary_triangles,
            tolerance_m=self.tolerance)
        pin = meshes['pinion']
        self.deformed_source = (np.asarray(pin.nodes_m)+self.u['pinion'])[np.asarray(pin.boundary_triangles)]
        self.source_distance_cache = {}; self.target_bridge_cache = cache['target_bridge']
        self.clipped_distance_cache = {}
        self.motion = sum(item['maximum_displacement_norm_m'] for item in identity.values())
        self.binding = dict(full_response_fingerprint=fingerprint, body_identity=identity,
            original_boundary_closure=self.closures, deformed_wheel_closure=self.target_material.closure,
            actual_deformation_jacobians=jacobians,
            actual_tetrahedral_boundary=cache['tetrahedral_boundary'],
            static_geometry_cache_identity=geometry_key,
            loaded_rotations={m: r.tolist() for m, r in self.rotations.items()},
            finite_FE_material_model=True, strict_continuous_cad_proof=False,
            cad_chord_allowance_evidence='名义CAD三角化弦差参数；本模块未声称独立连续NURBS全域Hausdorff界')

    def _target_bridge(self, source_wheel_local):
        """只裁剪不可能接触的CAD目标；所有可能近片均有原FE迹上界。"""
        import numpy as np
        target = self.context.target; box = self.context.complete_box
        if target is None:
            return dict(passed=False, reason='缺少完整目标CAD三角候选网')
        expansion = self.motion+self.context.source_deflection+self.context.target_deflection+self.tolerance
        low, high = source_wheel_local.min(axis=0)-expansion, source_wheel_local.max(axis=0)+expansion
        if box is not None and (np.any(low < box[0]) or np.any(high > box[1])):
            return dict(passed=False, reason='潜在接触CAD候选盒超出已证明完整的裁剪盒')
        bounds = tuple(float(x) for axis in range(3) for x in (low[axis], high[axis]))
        target.locator.FindCellsWithinBounds(bounds, target.ids)
        indices = sorted({int(target.ids.GetId(i)) for i in range(target.ids.GetNumberOfIds())})
        indices = [i for i in indices if bool(np.all(target.bounds_max[i] >= low)
            and np.all(target.bounds_min[i] <= high))]
        worst = 0.; worst_index = None
        for index in indices:
            if index not in self.target_bridge_cache:
                self.target_bridge_cache[index] = triangle_to_surface_upper_bound(target.triangles[index],
                    self.original['wheel'], **self.bridge_options)
            upper = self.target_bridge_cache[index]['distance_upper_bound_m']
            if upper > worst:
                worst, worst_index = upper, index
        return dict(passed=True, distance_upper_bound_m=worst+self.context.target_deflection,
            candidate_count=len(indices), candidate_indices=indices, worst_candidate_index=worst_index,
            complete_bvh_candidates=True, query_low_m=low.tolist(), query_high_m=high.tolist(),
            potential_contact_expansion_m=expansion, target_nominal_CAD_chord_allowance_m=self.context.target_deflection)

    def triangle_evidence(self, vertices_world_m, *, canonical_source_local_m=None):
        """判断一个原CAD叶片，失败保留未知，不把材料IN认作接触。"""
        import numpy as np
        triangle = np.asarray(vertices_world_m, dtype=float)
        local_source = triangle@self.rotations['pinion']
        input_roundoff = 0.
        if canonical_source_local_m is not None:
            canonical = np.asarray(canonical_source_local_m, dtype=float)
            if canonical.shape != (3, 3):
                raise ValueError('静态源叶片须为同物理齿局部3×3坐标')
            input_roundoff = float(np.linalg.norm(canonical-local_source, axis=1).max())
            if input_roundoff > self.tolerance:
                raise ValueError('静态源叶片与当前相位同物理叶片不匹配')
            # 对齐误差显式加入上下界，不通过截断坐标制造不真实缓存命中。
            local_source = canonical
        local_target = triangle@self.rotations['wheel']
        source_key = local_source.astype('<f8', copy=False).tobytes()
        if source_key not in self.static_cache['source_trace']:
            bridge = triangle_to_surface_upper_bound(local_source, self.original['pinion'], **self.bridge_options)
            self.static_cache['source_trace'][source_key] = dict(carriers=complete_trace_carriers(local_source,
                self.original['pinion'], source_allowance_m=self.context.source_deflection+self.tolerance,
                proven_distance_upper_m=bridge['distance_upper_bound_m']), bridge=bridge)
        entry = self.static_cache['source_trace'][source_key]
        carriers, bridge = entry['carriers'], entry['bridge']
        source_bound = bridge['distance_upper_bound_m']+self.context.source_deflection+input_roundoff
        target = self._target_bridge(local_target)
        if not target['passed']:
            return dict(passed=False, source_bridge=bridge, source_carriers=carriers, target_bridge=target)
        lower = math.inf; states = []; min_index = None
        for index in carriers['indices']:
            if index not in self.source_distance_cache:
                facet = self.deformed_source[index]@self.rotations['pinion'].T@self.rotations['wheel']
                self.source_distance_cache[index] = self.target_material.triangle_evidence(facet)
            value = self.source_distance_cache[index]
            if not value['outside'] or value['distance_m'] <= source_bound+target['distance_upper_bound_m']+self.tolerance:
                # 整个FE载体片可能碰到相邻已解接触区；真实最近迹只在完整盒内。
                # 逐原片精确裁剪并用原P1权重施加u，不能把最近迹改成插值贴图。
                key = (source_key, index)
                if key not in self.clipped_distance_cache:
                    polygon = clip_triangle_to_box_with_weights(self.original['pinion'].triangles[index],
                        carriers['query_low_m'], carriers['query_high_m'])
                    pieces = []
                    if len(polygon) >= 3:
                        weights = np.asarray([item[1] for item in polygon])
                        moved = weights@self.deformed_source[index]@self.rotations['pinion'].T@self.rotations['wheel']
                        for i in range(1, len(moved)-1):
                            fan = moved[[0, i, i+1]]
                            if float(np.linalg.norm(np.cross(fan[1]-fan[0], fan[2]-fan[0]))) <= 1e-30:
                                continue
                            pieces.append(self.target_material.triangle_evidence(fan))
                    # 若裁剪片退化或数值失败，回退整原片保守值，绝不跳过它。
                    if not polygon:
                        self.clipped_distance_cache[key] = dict(empty_trace_intersection=True, outside=True,
                            exact_convex_halfspace_clipping=True, original_boundary_triangle_index=index)
                    elif pieces:
                        self.clipped_distance_cache[key] = dict(distance_m=min(p['distance_m'] for p in pieces),
                            outside=all(p['outside'] for p in pieces), clipped_to_complete_trace_box=True,
                            clipped_polygon_material_m=[p.tolist() for p, _ in polygon],
                            clipped_polygon_P1_weights=[w.tolist() for _, w in polygon],
                            pieces=pieces, original_boundary_triangle_index=index)
                    else:
                        self.clipped_distance_cache[key] = dict(value, clipping_degenerate_kept_full_triangle=True)
                value = self.clipped_distance_cache[key]
            states.append(dict(index=index, **value))
            if value.get('empty_trace_intersection'):
                continue
            if value['distance_m'] < lower:
                lower, min_index = value['distance_m'], index
        if not math.isfinite(lower):
            raise ArithmeticError('完整最近迹候选裁剪后无任何材料点，拒绝虚假正间隙')
        remaining = lower-source_bound-target['distance_upper_bound_m']-self.tolerance
        outside = all(item['outside'] for item in states)
        return dict(passed=outside and remaining > 0., deformed_FE_distance_m=lower,
            distance_lower_under_declared_trace_allowances_m=max(0., remaining), remaining_margin_m=remaining,
            source_CAD_to_FE_trace_upper_m=source_bound, target_CAD_to_FE_trace_upper_m=target['distance_upper_bound_m'],
            source_bridge=bridge, target_bridge=target, source_carriers=carriers,
            source_material_triangle_sha256=hashlib.sha256(source_key).hexdigest(),
            controlling_source_FE_triangle=min_index, source_FE_material_evidence=states,
            all_complete_carriers_OUT=outside, actual_u_applied=True, relative_motion_upper_bound_m=0.,
            canonical_source_alignment_allowance_m=input_roundoff,
            strict_continuous_cad_proof=False)

    def refine_certificate(self, parent_certificate, *, maximum_depth=5, maximum_regions=500000,
            parent_sha256=None, progress=None):
        """只细化父证书未知片；近域、远域与未决完整覆盖原工作面。

        已知父叶原样保留。所有新增零域都依赖本次真实变形闭材料及
        完整几何桥接；失败叶仍未知。该方法不改变原压力、u/F或父证书。
        """
        from copy import deepcopy
        from pathlib import Path
        import hashlib
        import numpy as np
        from .distributed_pressure_domain import (distributed_pressure_fingerprint,
            distributed_pressure_domain_background)
        from .pressure_domain_audit import _grid_cover, _periodic_pinion_triangles
        from .pressure_background import _digest
        if isinstance(maximum_depth, bool) or not isinstance(maximum_depth, int) or not 0 <= maximum_depth <= 12:
            raise ValueError('材料域细化深度必须为0至12整数')
        parent_identity = _digest(parent_certificate)
        offset = parent_certificate['tooth_offset']
        before = distributed_pressure_fingerprint(self.result, tooth_offset=offset)
        source_before = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        distributed_pressure_domain_background(self.result, self.context.background, tooth_offset=offset,
            certificate=parent_certificate, pinion_request=self.context.models['pinion'].request)
        parent_depth = parent_certificate['regions']['maximum_depth']
        maximum_depth = max(parent_depth, maximum_depth)
        patch = next(p for p in self.result['geometry']['patches'] if p['tooth_offset'] == offset)
        if offset in self.result['selected_tooth_offsets']:
            field = self.result['fields'][self.result['selected_tooth_offsets'].index(offset)]
            near = _grid_cover(patch, field, padding_m=self.context.source_deflection)
        else:
            near = lambda triangle: False
        parent_rows = parent_certificate['regions']['regions']
        physical_roots, _, _ = _periodic_pinion_triangles(self.context.triangles, offset,
            self.context.models['pinion'].request.tooth_count, self.context.models['wheel'].request.tooth_count)
        def canonical_leaf(value, source_index, depth):
            # 以已验证的原CAD四分树整数重心坐标恢复同一材料坐标叶片。
            # 跨phase缓存不对世界坐标随意舍入；变换舍入差仍显式扣除。
            return canonical_midpoint_leaf(value, physical_roots[source_index], self.rotations['pinion'], depth,
                tolerance_m=self.tolerance)
        rows = []; evaluated = 0; refined = 0; zero = 0; solved = 0
        def visit(parent_index, value, depth):
            nonlocal evaluated, refined, zero, solved
            original = parent_rows[parent_index]
            row = dict(original, vertices_world_m=value.tolist(), centre_world_m=value.mean(axis=0).tolist(),
                area_m2=float(np.linalg.norm(np.cross(value[1]-value[0], value[2]-value[0]))*.5), depth=depth,
                parent_region_index=parent_index, state='unresolved')
            evaluated += 1
            if near(value):
                row.update(state='original_solved_domain', near_domain_rechecked=True)
                solved += 1
            else:
                proof = self.triangle_evidence(value,
                    canonical_source_local_m=canonical_leaf(value, original['source_triangle_index'], depth))
                # 完整FE面片材料证据按真实面片ID去重保存，叶片保留全部候选ID。
                proof.pop('source_FE_material_evidence', None)
                row['deformed_material_evidence'] = proof
                if proof['passed']:
                    row.update(state='excluded_by_distance_bound',
                        distance_lower_bound_m=proof['distance_lower_under_declared_trace_allowances_m'],
                        relative_motion_upper_bound_m=0.,
                        exclusion_model='实际变形FE闭材料＋完整双面CAD/FE迹上界；非原全局位移范数法')
                    zero += 1
                elif depth < maximum_depth and len(rows)+evaluated+4 < maximum_regions:
                    refined += 1
                    for child in _split(value): visit(parent_index, child, depth+1)
                    return
            rows.append(row)
        for index, row in enumerate(parent_rows):
            if row['state'] == 'unresolved':
                visit(index, np.asarray(row['vertices_world_m'], dtype=float), row['depth'])
            else:
                rows.append(deepcopy(row))
            if progress and index % 500 == 0:
                progress(dict(stage='deformed_material_refinement', tooth_offset=offset,
                    original_rows_processed=index+1, original_row_count=len(parent_rows),
                    evaluated_regions=evaluated, retained_regions=len(rows), refined_regions=refined,
                    newly_excluded=zero, newly_near_solved=solved,
                    source_FE_cache_count=len(self.source_distance_cache), target_bridge_cache_count=len(self.target_bridge_cache)))
        unknown = [row for row in rows if row['state'] == 'unresolved']
        regions = dict(parent_certificate['regions'], regions=rows, maximum_depth=maximum_depth,
            passed=not unknown, unresolved_region_count=len(unknown), terminal_region_count=len(rows),
            evaluated_region_count=parent_certificate['regions']['evaluated_region_count']+evaluated,
            unresolved_source_triangles=sorted({row['source_triangle_index'] for row in unknown}),
            area_m2={name: math.fsum(row['area_m2'] for row in rows if row['state'] == name) for name in
                ('original_solved_domain', 'excluded_by_distance_bound', 'unresolved')})
        child = deepcopy(parent_certificate)
        child.update(regions=regions, unresolved_region_count=len(unknown), passed=not unknown,
            complete_working_flank=not unknown,
            evidence_kind='原P0有效自由度已解近域；父远域；同次实际变形FE完整闭材料及双面CAD→FE迹上界新增远域',
            deformed_material_binding=deepcopy(self.binding),
            deformed_source_FE_evidence={**parent_certificate.get('deformed_source_FE_evidence', {}),
                **{str(k): v for k, v in self.source_distance_cache.items()}},
            deformed_clipped_trace_evidence={**parent_certificate.get('deformed_clipped_trace_evidence', {}),
                **{hashlib.sha256(k[0]).hexdigest()+':'+str(k[1]): v for k, v in self.clipped_distance_cache.items()}},
            target_CAD_FE_trace_evidence={**parent_certificate.get('target_CAD_FE_trace_evidence', {}),
                **{str(k): v for k, v in self.target_bridge_cache.items()}},
            strict_continuous_cad_proof=False)
        if 'whole_material_components' in child:
            child['parent_whole_material_components'] = dict(components=child.pop('whole_material_components'),
                region_indices_scope='父证书原regions行号；不得当作细分子证书的新行号',
                parent_certificate_sha256=parent_sha256)
        for key in ('full_working_flank_covered', 'omitted_regions_excluded', 'unresolved_regions_excluded', 'far_field_motion_bound_passed'):
            child['checks'][key] = not unknown
        after = distributed_pressure_fingerprint(self.result, tooth_offset=offset)
        source_after = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        if before != after or source_before != source_after or _digest(parent_certificate) != parent_identity:
            raise ArithmeticError('变形材料审核期间原态、父证书或执行源码发生改变')
        child['refinement_identity'] = dict(parent_certificate_sha256=parent_sha256,
            parent_semantic_fingerprint=parent_identity, parent_unmodified=True,
            state_fingerprint_before=before, state_fingerprint_after=after,
            source_sha256_before=source_before, source_sha256_after=source_after, unchanged=True,
            parent_original_row_count=len(parent_rows), newly_evaluated_regions=evaluated,
            unchanged_parent_known_regions=sum(row['state'] != 'unresolved' for row in parent_rows),
            near_domain_definition_unchanged=True, original_pressure_modified=False)
        child['limitations'] = ['原P0弱平均约束及独立Gauss门控不等于连续点态非穿透证明',
            '实际变形材料域属于原线性四面体边界；CAD/FE最近迹模型和名义CAD弦差明确为工程几何桥接',
            '几何桥接逐完整三角保守上界，非仅采样最大误差；没有独立连续NURBS全域Hausdorff定理',
            '本域审核不证明压力空间网格或相位收敛']
        distributed_pressure_domain_background(self.result, self.context.background, tooth_offset=offset,
            certificate=child, pinion_request=self.context.models['pinion'].request)
        return child
