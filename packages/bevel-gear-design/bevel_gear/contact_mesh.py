"""新写FE边界面、解析直齿侧面筛选和真实三角形最近点接触离散。

坐标m、面积m²，三角形为0起始FE节点索引。只在调用时加载NumPy。
本模块不生成柔度或载荷；第一面单重心积分是单向离散，不是双向mortar认证。
"""
from __future__ import annotations
import math


def _numpy():
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError("接触面离散需要可选NumPy依赖") from error
    return np


def _finite_nonnegative(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name}须为{'正' if positive else '非负'}有限实数")
    return float(value)


def _indices(value, width, nodes, name, np):
    raw = np.asarray(value)
    if raw.size == 0:
        return np.empty((0, width), dtype=np.int64)
    if raw.ndim != 2 or raw.shape[1] != width or not np.issubdtype(raw.dtype, np.integer):
        raise ValueError(f"{name}须为整数索引的N×{width}数组")
    result = raw.astype(np.int64, copy=False)
    if result.min() < 0 or result.max() >= nodes or np.any(np.diff(np.sort(result, axis=1), axis=1) == 0):
        raise ValueError(f"{name}含越界或重复节点索引")
    return result


def boundary_facets(mesh) -> dict:
    """以四面体单邻接确定边界，用面/单元重心校正朝外，保留原面实体标签。

    triangles是FE全局节点索引；nodes_m保留原全节点坐标。未知标签为0，
    boundary_indices稳定标识本次全部边界面的序号。非流形面和退化单元明确拒绝。
    """
    from .elasticity import TetraMesh
    if not isinstance(mesh, TetraMesh):
        raise TypeError("mesh须为TetraMesh")
    np = _numpy()
    nodes = np.asarray(mesh.nodes_m, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) < 4 or not np.isfinite(nodes).all():
        raise ValueError("FE节点须为至少4个有限三维坐标m")
    tetrahedra = _indices(mesh.tetrahedra, 4, len(nodes), "四面体", np)
    if not len(tetrahedra):
        raise ValueError("至少需要一个四面体")
    cells = nodes[tetrahedra]
    edges = cells[:, 1:] - cells[:, :1]
    six_volume = np.abs(np.einsum("ij,ij->i", np.cross(edges[:, 0], edges[:, 1]), edges[:, 2]))
    scale = np.max(np.linalg.norm(edges, axis=2), axis=1)
    if np.any(six_volume <= 64*np.finfo(float).eps*scale**3):
        raise ValueError("四面体退化或过度扁平，不能确定可靠边界法向")
    if len({tuple(sorted(map(int, row))) for row in tetrahedra}) != len(tetrahedra):
        raise ValueError("存在重复四面体")
    adjacency = {}
    for cell_index, row in enumerate(tetrahedra):
        for local in ((1,2,3),(0,2,3),(0,1,3),(0,1,2)):
            triangle = tuple(int(row[i]) for i in local)
            adjacency.setdefault(tuple(sorted(triangle)), []).append((triangle, cell_index))
    if any(len(items) > 2 for items in adjacency.values()):
        raise ValueError("存在多于两个四面体共享同一三角面的非流形连接")
    given = _indices(mesh.boundary_triangles, 3, len(nodes), "边界三角形", np)
    raw_tags = np.asarray(mesh.boundary_tags)
    if raw_tags.size and (raw_tags.ndim != 1 or len(raw_tags) != len(given) or not np.issubdtype(raw_tags.dtype, np.integer)):
        raise ValueError("边界标签须为与给定三角形等长的整数数组")
    tag_map = {}
    for index, row in enumerate(given):
        key = tuple(sorted(map(int, row)))
        tag = int(raw_tags[index]) if raw_tags.size else 0
        if key in tag_map and tag_map[key] != tag:
            raise ValueError("同一边界面有冲突标签")
        tag_map[key] = tag
    triangles, centroids, normals, areas, tags = [], [], [], [], []
    for key in sorted(adjacency):
        items = adjacency[key]
        if len(items) != 1:
            continue
        triangle, cell_index = items[0]
        points = nodes[list(triangle)]
        centroid = points.mean(axis=0)
        normal = np.cross(points[1]-points[0], points[2]-points[0])
        magnitude = float(np.linalg.norm(normal))
        if np.dot(normal, centroid-cells[cell_index].mean(axis=0)) < 0:
            triangle = (triangle[0], triangle[2], triangle[1])
            normal = -normal
        triangles.append(triangle)
        centroids.append(centroid)
        normals.append(normal/magnitude)
        areas.append(.5*magnitude)
        tags.append(tag_map.get(key, 0))
    return {"triangles":np.asarray(triangles,dtype=np.int64).reshape(-1,3),
        "centroids_m":np.asarray(centroids,dtype=float).reshape(-1,3),
        "normals":np.asarray(normals,dtype=float).reshape(-1,3),"areas_m2":np.asarray(areas,dtype=float),
        "face_tags":np.asarray(tags,dtype=np.int64),"nodes_m":nodes,
        "boundary_indices":np.arange(len(triangles),dtype=np.int64),
        "metadata":{"units":"m","normal_definition":"由四面体单元重心指向外部",
                    "boundary_definition":"四面体单邻接三角面"}}


def select_straight_flank_facets(mesh, surface, *, flank="right", tooth_index=0,
                                 matching_tolerance_m=3e-6) -> dict:
    """在实际FE节点上逐顶点匹配解析方位角，选有效球面渐开线侧面。

    不调用旧球面端界profile_bounds裁切FE节点；实体CAD的轴向端界由FE网格
    本身保留。只要求基锥以上、正齿厚的解析侧面方位匹配。
    tooth_index=None选择完整周期中的全部齿，单个三角面三顶点必须归属于同一齿。
    """
    from .tooth_surface import SphericalToothSurface
    if not isinstance(surface, SphericalToothSurface):
        raise TypeError("surface须为SphericalToothSurface")
    if flank not in ("left","right"):
        raise ValueError("flank须为left或right")
    tolerance = _finite_nonnegative(matching_tolerance_m,"matching_tolerance_m",positive=True)
    if tooth_index is not None and (isinstance(tooth_index,bool) or not isinstance(tooth_index,int) or not 0 <= tooth_index < surface.teeth):
        raise ValueError("tooth_index须为齿数范围内的整数或None")
    np = _numpy()
    facets = boundary_facets(mesh)
    nodes = facets["nodes_m"]
    radial = np.hypot(nodes[:,0],nodes[:,1])
    polar = np.arctan2(radial,nodes[:,2])
    azimuth = np.arctan2(nodes[:,1],nodes[:,0])
    pitch = 2*math.pi/surface.teeth
    indices = (np.floor(azimuth/pitch+.5).astype(np.int64) % surface.teeth) if tooth_index is None else np.full(len(nodes),tooth_index,dtype=np.int64)
    valid = np.zeros(len(nodes),dtype=bool)
    distances = np.full(len(nodes),math.inf)
    for index, delta in enumerate(polar):
        if radial[index] <= 0 or not surface.profile.base_cone_angle_rad-2e-14 <= delta < math.pi/2:
            continue
        try:
            expected = surface.profile.flank_angle(max(float(delta),surface.profile.base_cone_angle_rad),flank)
        except ValueError:
            continue
        error = math.atan2(math.sin(float(azimuth[index])-indices[index]*pitch-expected),
                           math.cos(float(azimuth[index])-indices[index]*pitch-expected))
        distances[index] = 2*radial[index]*abs(math.sin(error*.5))
        valid[index] = distances[index] <= tolerance
    triangles = facets["triangles"]
    mask = np.all(valid[triangles],axis=1)
    if tooth_index is None:
        mask &= np.all(indices[triangles] == indices[triangles[:,0,None]],axis=1)
    result = {key:(value if key in ("nodes_m","metadata") else value[mask]) for key,value in facets.items()}
    result["tooth_indices"] = indices[triangles[mask,0]]
    result["maximum_vertex_matching_distance_m"] = np.max(distances[triangles[mask]],axis=1) if np.any(mask) else np.empty(0)
    result["metadata"] = {**facets["metadata"],"flank":flank,"tooth_index":tooth_index,
        "matching_tolerance_m":tolerance,"selection":"逐FE顶点的解析齿面方位弦距；不裁切球面宽端界"}
    return result


def _rotation(value,np):
    rotation = np.asarray(value,dtype=float)
    if rotation.shape != (3,3) or not np.isfinite(rotation).all() or not np.allclose(rotation.T@rotation,np.eye(3),rtol=0,atol=1e-10) or not math.isclose(float(np.linalg.det(rotation)),1.,abs_tol=1e-10):
        raise ValueError("旋转须为按行存储的有限3×3正交矩阵且行列式为+1")
    return rotation


def _facet_arrays(facets,np):
    nodes = np.asarray(facets["nodes_m"],dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or not np.isfinite(nodes).all():
        raise ValueError("facets.nodes_m须为有限N×3坐标")
    triangles = _indices(facets["triangles"],3,len(nodes),"接触三角形",np)
    count = len(triangles)
    centroids = np.asarray(facets["centroids_m"],dtype=float)
    normals = np.asarray(facets["normals"],dtype=float)
    areas = np.asarray(facets["areas_m2"],dtype=float)
    if centroids.shape != (count,3) or normals.shape != (count,3) or areas.shape != (count,) or not all(np.isfinite(v).all() for v in (centroids,normals,areas)):
        raise ValueError("接触面重心、法向或面积布局无效")
    lengths = np.linalg.norm(normals,axis=1)
    if np.any(lengths <= 0) or np.any(areas <= 0):
        raise ValueError("接触面法向与面积须非退化")
    vertices = nodes[triangles]
    cross = np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0])
    edge_scale = np.maximum(np.sum((vertices[:,1]-vertices[:,0])**2,axis=1),np.sum((vertices[:,2]-vertices[:,0])**2,axis=1))
    if np.any(np.linalg.norm(cross,axis=1) <= 32*np.finfo(float).eps*edge_scale):
        raise ValueError("接触三角形退化或过度扁平")
    if not np.allclose(vertices.mean(axis=1),centroids,rtol=1e-10,atol=1e-12):
        raise ValueError("提供的面重心与实际FE节点不一致")
    return triangles,vertices,centroids,normals/lengths[:,None],areas


def _closest_on_triangles(point,vertices,np):
    """真实点到三角形：内部投影与三条有限线段共同取最小。

    point为(3,)时共享同一点；为(N,3)时逐行对应N个三角形，复用
    完全相同的几何公式，避免逐点构造大量NumPy小数组。
    """
    a,b,c = vertices[:,0],vertices[:,1],vertices[:,2]
    ab,ac,ap = b-a,c-a,point-a
    dot = lambda u,v:np.einsum("ij,ij->i",u,v)
    normal = np.cross(ab,ac)
    denominator = dot(normal,normal)
    v = dot(np.cross(ap,ac),normal)/denominator
    w = dot(np.cross(ab,ap),normal)/denominator
    bary = np.stack((1-v-w,v,w),axis=1)
    plane = a+v[:,None]*ab+w[:,None]*ac
    points = [plane]
    weights = [bary]
    for first,second,edge_id in ((a,b,0),(b,c,1),(c,a,2)):
        edge = second-first
        t = np.clip(dot(point-first,edge)/dot(edge,edge),0.,1.)
        points.append(first+t[:,None]*edge)
        if edge_id == 0:weights.append(np.stack((1-t,t,np.zeros_like(t)),axis=1))
        elif edge_id == 1:weights.append(np.stack((np.zeros_like(t),1-t,t),axis=1))
        else:weights.append(np.stack((t,np.zeros_like(t),1-t),axis=1))
    points = np.stack(points,axis=1)
    weights = np.stack(weights,axis=1)
    query=np.asarray(point)
    if query.ndim==2:
        if query.shape[1:]!=(3,) or len(query) not in (1,len(vertices)):
            raise ValueError('逐行最近点查询须与三角形数量相同或仅含一个共享点')
        query=query[:,None,:]
    squared = np.sum((points-query)**2,axis=2)
    squared[np.any(bary < 0,axis=1),0] = math.inf
    candidate = np.argmin(squared,axis=1)
    rows = np.arange(len(vertices))
    return points[rows,candidate],weights[rows,candidate],squared[rows,candidate]


def pair_facet_contacts(first, second, first_rotation_rows, second_rotation_rows, *,
                        maximum_distance_m=2e-4, maximum_normal_angle_rad=.25) -> list:
    """每个第一面重心对所有法向合格第二三角形作真实最近点搜索。

    world=R@local，共同原点；平移须由调用者预先在节点坐标中处理。
    signed_gap保留正负；两外法向相反才配对。返回节点权重与力方向互为功共轭，
    不保证一一配对或已形成非冗余接触约束，未执行实体侵入的体积检查。
    """
    distance_limit = _finite_nonnegative(maximum_distance_m,"maximum_distance_m")
    angle_limit = _finite_nonnegative(maximum_normal_angle_rad,"maximum_normal_angle_rad")
    if angle_limit >= math.pi/2:
        raise ValueError("最大法向偏角须小于π/2，以保留相反的外法向条件")
    np = _numpy()
    r1,r2 = _rotation(first_rotation_rows,np),_rotation(second_rotation_rows,np)
    tri1,verts1,cent1,norm1,areas1 = _facet_arrays(first,np)
    tri2,verts2,cent2,norm2,areas2 = _facet_arrays(second,np)
    if not len(tri1) or not len(tri2):
        return []
    centers = cent1@r1.T
    vertices = verts2@r2.T
    normals1,normals2 = norm1@r1.T,norm2@r2.T
    contacts = []
    for index,point in enumerate(centers):
        cosines = np.clip((-normals2)@normals1[index],-1.,1.)
        eligible = np.flatnonzero(cosines >= math.cos(angle_limit)-4e-15)
        if not len(eligible):
            continue
        nearest,weights,squared = _closest_on_triangles(point,vertices[eligible],np)
        chosen = int(np.argmin(squared))
        distance = math.sqrt(max(0.,float(squared[chosen])))
        if distance > distance_limit:
            continue
        second_index = int(eligible[chosen])
        normal = normals1[index]-normals2[second_index]
        normal /= np.linalg.norm(normal)
        bary = weights[chosen]
        zeros = int(np.count_nonzero(np.abs(bary) <= 1e-10))
        feature = "vertex" if zeros >= 2 else "edge" if zeros == 1 else "face"
        contacts.append({"first_facet_index":index,"second_facet_index":second_index,
            "first_node_indices":tuple(map(int,tri1[index])),"first_weights":(1/3,1/3,1/3),
            "second_node_indices":tuple(map(int,tri2[second_index])),"second_weights":tuple(map(float,bary)),
            "first_point_world_m":tuple(map(float,point)),"second_point_world_m":tuple(map(float,nearest[chosen])),
            "normal_world":tuple(map(float,normal)),"first_force_normal_local":tuple(map(float,-r1.T@normal)),
            "second_force_normal_local":tuple(map(float,r2.T@normal)),"area_m2":float(areas1[index]),
            "distance_m":distance,"signed_gap_m":float(np.dot(nearest[chosen]-point,normal)),
            "normal_angle_rad":math.acos(float(cosines[second_index])),"second_closest_feature":feature})
    return contacts
