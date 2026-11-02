"""真实有限FE投影分片的几何/响应审核；不调用NURBS，不创建K或求压力。

协议名称只定位分支。几何来自原全局节点、整体机械R和原有限交片，
必须逐片复算，不能仅凭integration_method声称算法或压力通过。
"""
from __future__ import annotations

import hashlib
import math

NATIVE_GEOMETRY_SOURCE = 'FE_boundary'
NATIVE_INTEGRATION_METHOD = 'P0_cell_intersection_with_both_actual_FE_triangles_centroid_exact'
NATIVE_METHOD_ALIASES = ('projected_FE_partition_exact',)
NATIVE_OPERATOR_IDENTITY_FIELDS = ('current_operator_manifest_sha256', 'operator_fingerprint')
NATIVE_IDENTITY_FIELDS = ('geometry_source', 'integration_method', 'native_cache_index_sha256',
    'source_input_fingerprint', 'source_window_input_fingerprint', 'actual_native_mesh_sha256',
    'operator_input_fingerprint', 'actual_native_mesh_numeric_identity',
    *NATIVE_OPERATOR_IDENTITY_FIELDS)


def _sha_identity(value, name):
    if (not isinstance(value, str) or len(value) != 64 or
            any(c not in '0123456789abcdef' for c in value.lower())):
        raise ValueError('原生FE缺少明确的SHA256身份：'+name)


def _native_numeric_mesh_identity(value):
    """生产内存接口绑定实际完整网格数组，不要求历史缓存或文件来源。"""
    if not isinstance(value, dict) or set(value) != {'pinion', 'wheel'}:
        raise ValueError('生产原生FE须绑定两体完整网格数字身份')
    for member in value.values():
        if not isinstance(member, dict) or set(member) != {'nodes_m', 'tetrahedra', 'boundary_triangles'}:
            raise ValueError('生产网格身份须含nodes_m/tetrahedra/boundary_triangles')
        for name, width, dtype in (('nodes_m', 3, '<f8'), ('tetrahedra', 4, '<i8'), ('boundary_triangles', 3, '<i8')):
            item = member[name]
            if not isinstance(item, dict): raise ValueError('生产网格数字身份须为具名字典：'+name)
            shape = item.get('shape', [])
            if item.get('dtype') != dtype or len(shape) != 2 or shape[1] != width or isinstance(shape[0], bool) or not isinstance(shape[0], int) or shape[0] < 1:
                raise ValueError('生产网格身份的dtype/shape不完整：'+name)
            _sha_identity(item.get('sha256'), name)
    return value


def _native_operator_identity(value):
    """新 outer 的最终 operator 单独绑定；原静态 cache SHA 不能替代它。"""
    present = {name: value[name] for name in NATIVE_OPERATOR_IDENTITY_FIELDS if name in value}
    if present or 'current_operator_manifest_file' in value:
        if len(present) != len(NATIVE_OPERATOR_IDENTITY_FIELDS):
            raise ValueError('当前原生 operator 须同时提供 manifest SHA 和 operator_fingerprint')
        for name, digest in present.items():
            if (not isinstance(digest, str) or len(digest) != 64 or
                    any(c not in '0123456789abcdef' for c in digest.lower())):
                raise ValueError('当前原生 operator 身份须为明确 SHA256：'+name)
    return present


def native_geometry_mode(result):
    """识别显式原生分支，标签冲突/缺失不回退NURBS。"""
    geometry = result.get('geometry', {})
    top, inner = result.get('geometry_source'), geometry.get('geometry_source')
    if top is not None and inner is not None and top != inner:
        raise ValueError('result与geometry的geometry_source不一致')
    source = top if top is not None else inner
    if source == NATIVE_GEOMETRY_SOURCE:
        if top != source or inner != source:
            raise ValueError('原生FE结果须在result和geometry都明确geometry_source=FE_boundary')
        method = result.get('integration_method')
        if method != geometry.get('integration_method'):
            raise ValueError('原生FE result/geometry的integration_method不一致')
        if method not in (NATIVE_INTEGRATION_METHOD, *NATIVE_METHOD_ALIASES):
            raise NotImplementedError('未支持的原生FE积分协议：'+str(method))
        if result.get('operator_input_fingerprint') is not None and result['operator_input_fingerprint'] != geometry.get('operator_input_fingerprint'):
            raise ValueError('生产原生result/geometry的operator_input_fingerprint不同')
        return True
    if result.get('native_cache_index_sha256') is not None or result.get('operator_input_fingerprint') is not None or _native_operator_identity(result):
        raise ValueError('具有native缓存身份却未声明FE_boundary，不能按NURBS读取')
    if source not in (None, 'NURBS', 'CAD_NURBS'):
        raise NotImplementedError('未支持的压力几何来源：'+str(source))
    return False


def native_result_identity(result):
    """新字段存在时才进入新payload；历史缺省NURBS身份保持原字节算法。"""
    active = native_geometry_mode(result)
    value = {name: result[name] for name in NATIVE_IDENTITY_FIELDS if name in result}
    if active:
        _native_operator_identity(result)
        production = 'operator_input_fingerprint' in value
        if production:
            _sha_identity(value['operator_input_fingerprint'], 'operator_input_fingerprint')
            _native_numeric_mesh_identity(value.get('actual_native_mesh_numeric_identity'))
        for key in ('native_cache_index_sha256', 'source_input_fingerprint'):
            if not production or key in value: _sha_identity(value.get(key), key)
        meshes = value.get('actual_native_mesh_sha256')
        if not production or meshes is not None:
            if not isinstance(meshes, dict) or set(meshes) != {'pinion', 'wheel'}:
                raise ValueError('原生FE结果须绑定两体actual_native_mesh_sha256')
            for digest in meshes.values(): _sha_identity(digest, 'actual_native_mesh_sha256')
    return value


def native_geometry_definition_fingerprint(geometry):
    """只比较实际参考几何，不把求解完成标记当成几何改变。"""
    from .pressure_background import _digest
    names=('tooth_offset','grid','valid_mask','frames','plane_origin_world_m','plane_tangent_rows',
        'plane_outward_normal_world','gaps_m','first_points_world_m','second_points_world_m',
        'first_corner_points_world_m','second_corner_points_world_m','corner_valid_mask','native_partial_cell_mask')
    return _digest(dict(geometry_source=geometry.get('geometry_source'),integration_method=geometry.get('integration_method'),
        reference={key:geometry.get(key) for key in ('driving_angle_rad','pinion_torque_nm','pinion_to_wheel_ratio','reference_follower_angle_rad')},
        patches=[{key:p[key] for key in names if key in p} for p in geometry.get('patches',[])]))


def _polygon_geometry(vertices):
    import numpy as np
    points = np.asarray(vertices, dtype=float)
    if points.ndim != 2 or points.shape[1:] != (2,) or len(points) < 3 or not np.isfinite(points).all():
        raise ValueError('原生交片顶点须为至少3个有限二维SI点')
    local = points-points[0]; following = np.roll(local, -1, axis=0)
    cross = local[:, 0]*following[:, 1]-local[:, 1]*following[:, 0]
    signed = float(cross.sum())
    if not math.isfinite(signed) or signed == 0:
        raise ValueError('原生交片无正面积，不能归一化或补齐')
    centroid = points[0]+np.sum((local+following)*cross[:, None], axis=0)/(3*signed)
    return points, abs(signed)*.5, centroid


def _finite_vector(value,shape,name):
    import numpy as np
    array=np.asarray(value,dtype=float)
    if array.shape!=shape or not np.isfinite(array).all():raise ValueError(name+'须为有限指定形状数组')
    return array


def _finite_scalar(value,name):
    if isinstance(value,bool):raise ValueError(name+'不能用布尔值替代SI数值')
    number=float(value)
    if not math.isfinite(number):raise ValueError(name+'须为有限数')
    return number


def _body_boundary_keys(mesh):
    import numpy as np
    nodes = np.asarray(mesh.nodes_m, dtype=float); triangles = np.asarray(mesh.boundary_triangles)
    if nodes.ndim != 2 or nodes.shape[1:] != (3,) or not np.isfinite(nodes).all():
        raise ValueError('原生审核需要真实全体节点N×3，单位m')
    if triangles.ndim != 2 or triangles.shape[1:] != (3,) or not len(triangles) or not np.issubdtype(triangles.dtype, np.integer):
        raise ValueError('原生审核需要实际全体边界三角全局索引')
    if triangles.min() < 0 or triangles.max() >= len(nodes):
        raise ValueError('原生边界三角引用越界节点')
    bits = max(1, (len(nodes)-1).bit_length())
    if bits*3 > 63:
        raise NotImplementedError('当前紧凑边界身份核验支持小于2^21个节点')
    def packed(value):
        ordered = np.sort(value, axis=-1).astype(np.uint64)
        return (ordered[..., 0] << np.uint64(2*bits)) | (ordered[..., 1] << np.uint64(bits)) | ordered[..., 2]
    return nodes, np.sort(packed(triangles)), packed


def audit_native_partition_geometry(system, meshes, *, result=None, mapping_tolerance_m=2e-10,
        coverage_relative_tolerance=2e-12, barycentric_tolerance=5e-11,
        include_recomputed_values=False, progress=None):
    """从同次实际节点/R0重建交片平均g0、原W列与刚体lever；不求u。

    records可以是按批懒加载Sequence。每次仅保留一个cell及其GEOS交片，
    不把全部quadrature展开成巨型JSON。源Rk只允许作为后处理身份保留，
    机械映射只读patch.frames中的整体R0。返回(g0,物理格keys,几何证据)。
    """
    from .pressure_background import _digest, _canonical
    from .distributed_pressure_domain import numeric_array_identity
    import json
    if system.get('geometry_source') != NATIVE_GEOMETRY_SOURCE:
        raise ValueError('原生FE审核需要显式system.geometry_source=FE_boundary')
    method = system.get('integration_method')
    if method not in (NATIVE_INTEGRATION_METHOD, *NATIVE_METHOD_ALIASES):
        raise NotImplementedError('未支持的原生FE积分协议：'+str(method))
    geometry, records = system.get('geometry'), system.get('records')
    if not isinstance(geometry, dict) or records is None or not len(records):
        raise ValueError('原生FE审核缺少真实geometry或完整records')
    if geometry.get('geometry_source') != NATIVE_GEOMETRY_SOURCE or geometry.get('integration_method') != method:
        raise ValueError('原生geometry头须与system的来源/积分协议一致')
    if system.get('operator_input_fingerprint') is not None:
        _sha_identity(system['operator_input_fingerprint'], 'operator_input_fingerprint')
        if system['operator_input_fingerprint'] != geometry.get('operator_input_fingerprint'):
            raise ValueError('生产原生system/geometry的operator_input_fingerprint不同')
    try:
        import numpy as np
        from shapely.geometry import Polygon, box
        from shapely.ops import unary_union
    except ImportError as error:
        raise ImportError('真实原生FE分片审核需要可选NumPy/GEOS依赖') from error
    if result is not None:
        if not native_geometry_mode(result):raise ValueError('原生系统不能审核NURBS结果')
        native_result_identity(result)
        if result.get('operator_input_fingerprint') is not None:
            if result['operator_input_fingerprint'] != system.get('operator_input_fingerprint'):
                raise ValueError('生产原生结果与当前真实records的operator_input_fingerprint不同')
        elif result['native_cache_index_sha256'] != system.get('native_cache_index_sha256'):
            raise ValueError('原生结果与真实records的cache index SHA不一致')
        if result['integration_method'] != method:
            raise ValueError('原生结果与真实records的积分协议不一致')
        if native_geometry_definition_fingerprint(result['geometry']) != native_geometry_definition_fingerprint(geometry):
            raise ValueError('原生结果与真实records使用的geometry不一致')
        if _native_operator_identity(result) != _native_operator_identity(system):
            raise ValueError('原生结果与真实交片/W的当前 operator 身份不同，不能复用 outer0 系统')
    operator_identity = _native_operator_identity(system)
    members = ('pinion', 'wheel')
    if set(meshes) != set(members):raise ValueError('原生齿轮审核须提供pinion/wheel两体真实网格')
    bodies = {m: _body_boundary_keys(meshes[m]) for m in members}
    if result is not None and 'actual_native_mesh_numeric_identity' in result:
        numeric_meshes = {m: dict(nodes_m=numeric_array_identity(meshes[m].nodes_m),
            tetrahedra=numeric_array_identity(meshes[m].tetrahedra, integer=True),
            boundary_triangles=numeric_array_identity(meshes[m].boundary_triangles, integer=True)) for m in members}
        if numeric_meshes != result['actual_native_mesh_numeric_identity']:
            raise ValueError('生产原生结果与实际读取的完整节点/四面体/边界数字身份不同')
    patches = {p['tooth_offset']: p for p in geometry.get('patches', [])}
    if not patches:raise ValueError('原生几何缺少真实物理齿patch')
    matrices = system.get('matrices')
    if matrices is not None:
        if set(matrices) != set(members):raise ValueError('原生W须完整提供两体')
        matrices = {m: w.tocsc(copy=True) for m, w in matrices.items()}
        for m, w in matrices.items():
            w.sum_duplicates(); w.sort_indices()
            if not np.isfinite(w.data).all():raise ValueError('原生W含非有限数')
            if w.shape != (3*len(bodies[m][0]), len(records)):
                raise ValueError('原生W尺寸与真实全部节点/原cell列数不符')
    expected = [(k, int(j), int(i)) for k, p in patches.items()
        for j, i in zip(*np.nonzero(np.asarray(p['valid_mask'], dtype=bool)))]
    if len(expected) != len(records):raise ValueError('原生有效mask与完整records列数不符')
    if not isinstance(include_recomputed_values, bool):raise ValueError('重建数值返回开关须为布尔值')
    keys, gaps = [], []; mean_levers = {m: [] for m in members}; digest = hashlib.sha256(); piece_count = 0
    maximum = dict(height_m=0., gap_m=0., shape_weight=0., ray_shape_closure_m=0.,
        coefficient=0., lever_m=0., W_column_relative_error=0., coverage_relative=0.)
    mechanical_frames = {}; world_nodes = {}; boundary_cache = {}; per_patch = {}
    for k, patch in patches.items():
        if patch.get('geometry_source') != NATIVE_GEOMETRY_SOURCE or patch.get('integration_method') != method:
            raise ValueError('原生patch必须与系统的几何/积分协议一致')
        frames = patch['frames']
        for m in members:
            rotation = np.asarray(frames[m+'_rotation_rows'], dtype=float)
            translation = np.asarray(frames.get(m+'_translation_world_m', [0., 0., 0.]), dtype=float)
            if translation.shape!=(3,) or not np.isfinite(translation).all():raise ValueError('两体独立刚体平移须为有限3向量')
            if (rotation.shape != (3, 3) or not np.isfinite(rotation).all() or
                    not np.allclose(rotation.T@rotation, np.eye(3), atol=1e-11, rtol=0) or np.linalg.det(rotation) <= 0):
                raise ValueError('原生机械R0须为真实正向刚体旋转')
            if m in mechanical_frames:
                if not np.allclose(rotation, mechanical_frames[m][0], atol=1e-12, rtol=0) or not np.array_equal(translation, mechanical_frames[m][1]):
                    raise ValueError('不同k的机械R必须为同一整体R0，不能使用source Rk')
            else:
                mechanical_frames[m] = (rotation, translation)
                world_nodes[m] = bodies[m][0]@rotation.T+translation
    for column, record in enumerate(records):
        key = (record['tooth_offset'], record['row'], record['column'])
        if key != expected[column]:raise ValueError('原生records物理k/j/i与mask/W列顺序不一致')
        keys.append(key);patch = patches[key[0]]
        grid = patch['grid'].as_dict() if hasattr(patch['grid'], 'as_dict') else patch['grid']
        dx, dy = float(grid['spacing_x_m']), float(grid['spacing_y_m'])
        centre = np.asarray(grid['origin_xy_m'])+np.asarray([key[2]*dx, key[1]*dy])
        cell = box(centre[0]-dx/2, centre[1]-dy/2, centre[0]+dx/2, centre[1]+dy/2)
        area = dx*dy
        if not math.isclose(float(record['area_plane_m2']), area, rel_tol=coverage_relative_tolerance, abs_tol=1e-24):
            raise ValueError('原生cell面积与原矩形格不一致')
        origin = _finite_vector(patch['plane_origin_world_m'],(3,),'投影原点'); basis = _finite_vector(patch['plane_tangent_rows'],(2,3),'投影切基')
        normal = _finite_vector(patch['plane_outward_normal_world'],(3,),'投影法向')
        if (origin.shape != (3,) or basis.shape != (2, 3) or normal.shape != (3,) or
                not np.allclose(basis@basis.T, np.eye(2), atol=1e-11, rtol=0) or
                not np.allclose(np.cross(*basis), normal, atol=1e-11, rtol=0)):
            raise ValueError('原生共同平面须为正向正交单位坐标')
        parts = []; gap = 0.; levers = {m: 0. for m in members}; W = {m: {} for m in members}
        record_digest = hashlib.sha256()
        for sample in record['quadrature']:
            vertices, piece_area, xy = _polygon_geometry(sample['piece_vertices_xy_m'])
            polygon = Polygon(vertices)
            if not polygon.is_valid or polygon.interiors:raise ValueError('原生交片必须是有效无孔有限多边形')
            if abs(piece_area-polygon.area) > coverage_relative_tolerance*area:
                raise ValueError('原生交片鞋带面积与真实有限多边形不一致')
            if (not math.isclose(piece_area, float(sample['piece_area_plane_m2']), rel_tol=coverage_relative_tolerance, abs_tol=1e-24) or
                    not math.isclose(float(sample['integration_weight']), piece_area/area, rel_tol=coverage_relative_tolerance, abs_tol=1e-16)):
                raise ValueError('原生交片面积/完整cell权重不符，不能重归一化partial')
            if np.linalg.norm(xy-_finite_vector(sample['target_xy_m'],(2,),'真实片质心')) > mapping_tolerance_m:
                raise ValueError('原生积分点不是真实交片质心')
            parts.append(polygon); weight = piece_area/area; heights = {}; sample_identity = {}
            for m, sign in zip(members, (-1., 1.)):
                mapping = sample['mappings'][m]; ids = np.asarray(mapping['nodes'])
                if ids.shape != (3,) or not np.issubdtype(ids.dtype, np.integer) or len(np.unique(ids)) != 3:
                    raise ValueError('原生射线须保留3个真实整数全局节点，不虚构或夹取形函数')
                nodes, boundary, packed = bodies[m]
                if ids.min() < 0 or ids.max() >= len(nodes):raise ValueError('原生射线全局节点越界')
                node_key = tuple(sorted(map(int, ids)))
                if node_key not in boundary_cache.setdefault(m, {}):
                    code = packed(ids); position = int(np.searchsorted(boundary, code))
                    if position >= len(boundary) or boundary[position] != code:
                        raise ValueError('原生射线三节点不是原真实FE边界三角')
                    boundary_cache[m][node_key] = True
                world = world_nodes[m][ids]; projected = (world-origin)@basis.T
                matrix = (projected[1:]-projected[0]).T
                try:ab = np.linalg.solve(matrix, xy-projected[0]); vertex_ab = np.linalg.solve(matrix, (vertices-projected[0]).T).T
                except np.linalg.LinAlgError as error:raise ValueError('原生边界三角投影退化') from error
                shapes = np.r_[1.-ab.sum(), ab]; vertex_shapes = np.column_stack((1.-vertex_ab.sum(axis=1), vertex_ab))
                if min(shapes.min(), vertex_shapes.min()) < -barycentric_tolerance or max(shapes.max(), vertex_shapes.max()) > 1.+barycentric_tolerance:
                    raise ValueError('原生交片或质心超出真实FE三角，不能nearest clamp')
                height = float(shapes@((world-origin)@normal)); mapped = shapes@world
                ray = origin+xy@basis+height*normal; closure = float(np.linalg.norm(mapped-ray))
                maximum['ray_shape_closure_m'] = max(maximum['ray_shape_closure_m'], closure)
                maximum['height_m'] = max(maximum['height_m'], abs(height-_finite_scalar(mapping['height_m'],'原生高度')))
                maximum['shape_weight'] = max(maximum['shape_weight'], float(np.abs(shapes-_finite_vector(mapping['weights'],(3,),'原P1权重')).max()))
                if np.linalg.norm(mapped-_finite_vector(mapping['mapped_world'],(3,),'原世界射线点')) > mapping_tolerance_m:
                    raise ValueError('原生实际全局节点/R0与记录world点不一致')
                face_normal = np.cross(world[1]-world[0], world[2]-world[0]);length=float(np.linalg.norm(face_normal))
                edge_scale=float(np.linalg.norm(world-np.roll(world,1,axis=0),axis=1).max())
                if length<=64*np.finfo(float).eps*edge_scale**2:raise ValueError('实际FE边界三角退化')
                face_normal /= length
                supplied = _finite_vector(mapping['normal_world'],(3,),'原FE法向')
                if abs(float(supplied@supplied)-1.) > 1e-10 or abs(abs(float(face_normal@supplied))-1.) > 1e-10:
                    raise ValueError('原生法向未平行真实FE有限面')
                if abs(float(supplied@normal)-_finite_scalar(mapping['normal_cosine'],'原法向投影余弦')) > 1e-10:
                    raise ValueError('原生法向投影分母不符')
                cosine = float(face_normal@normal)
                if abs(cosine) < 1e-6:raise ValueError('原生实际FE图对投影方向退化')
                coefficient_world = sign*face_normal/cosine
                rotation = mechanical_frames[m][0]; coefficient = rotation.T@coefficient_world
                coefficient_error = float(np.linalg.norm(coefficient-_finite_vector(sample['coefficients_local'][m],(3,),'原工作共轭系数')))
                maximum['coefficient'] = max(maximum['coefficient'], coefficient_error)
                axis = _finite_vector(patch['frames'][m+'_axis_world'],(3,),'刚体轴');axis_origin = _finite_vector(patch['frames'].get(m+'_axis_origin_world_m', [0., 0., 0.]),(3,),'刚体轴原点')
                if abs(float(axis@axis)-1.) > 1e-10:raise ValueError('原生刚体旋转轴须为单位向量')
                lever = -float(coefficient_world@np.cross(axis, mapped-axis_origin))
                maximum['lever_m'] = max(maximum['lever_m'], abs(lever-_finite_scalar(sample['levers_m'][m],'原lever m')))
                levers[m] += weight*lever; heights[m] = height
                if matrices is not None:
                    for node, shape in zip(ids, shapes):
                        for direction in range(3):
                            dof = 3*int(node)+direction
                            W[m][dof] = W[m].get(dof, 0.)+weight*float(shape)*float(coefficient[direction])
                sample_identity[m] = dict(nodes=ids.tolist(),weights=shapes.tolist(),height_m=height,
                    coefficient_local=coefficient.tolist(),lever_m=lever)
            local_gap = heights['wheel']-heights['pinion'];gap += weight*local_gap
            maximum['gap_m'] = max(maximum['gap_m'], abs(local_gap-_finite_scalar(sample['gap_m'],'原交片gap m')))
            record_digest.update(json.dumps(_canonical(dict(vertices=vertices,area=piece_area,centroid=xy,bodies=sample_identity)),
                sort_keys=True,separators=(',',':')).encode('utf-8'));record_digest.update(b'\n');piece_count += 1
        if not parts:raise ValueError('有效原生cell没有实际正面积交片')
        union = unary_union(parts);summed = math.fsum(p.area for p in parts)
        discrepancy = max(float(cell.difference(union).area),float(union.difference(cell).area),
            max(0.,summed-float(union.area)),abs(summed-area))/area
        maximum['coverage_relative'] = max(maximum['coverage_relative'], discrepancy)
        if discrepancy > coverage_relative_tolerance:raise ValueError('原生cell真实交片不满格或重复覆盖，不可作为原P0自由度')
        maximum['gap_m'] = max(maximum['gap_m'], abs(gap-_finite_scalar(record['gap_m'],'原cell平均gap m')))
        if 'gaps_m' in patch:
            maximum['gap_m'] = max(maximum['gap_m'], abs(gap-_finite_scalar(patch['gaps_m'][key[1]][key[2]],'原patch平均gap m')))
        for m in members:
            mean_levers[m].append(levers[m])
            maximum['lever_m'] = max(maximum['lever_m'], abs(levers[m]-_finite_scalar(record['levers_m'][m],'原cell lever m')))
            if matrices is not None:
                matrix = matrices[m];first,last = matrix.indptr[column:column+2]
                indices = np.asarray(sorted(W[m]),dtype=np.int64);values = np.asarray([W[m][int(i)] for i in indices])
                stored_ids,stored = matrix.indices[first:last],matrix.data[first:last]
                all_ids = np.union1d(indices,stored_ids);expected_column=np.zeros(len(all_ids));actual=np.zeros(len(all_ids))
                expected_column[np.searchsorted(all_ids,indices)] = values;actual[np.searchsorted(all_ids,stored_ids)] = stored
                error = float(np.linalg.norm(expected_column-actual)/max(float(np.linalg.norm(actual)),1.))
                maximum['W_column_relative_error'] = max(maximum['W_column_relative_error'],error)
        gaps.append(gap);digest.update(json.dumps([*map(int,key),record_digest.hexdigest()],separators=(',',':')).encode('ascii'));digest.update(b'\n')
        per_patch[key[0]] = per_patch.get(key[0],0)+1
        if progress is not None and (column+1)%128==0:
            progress(dict(stage='native_FE_partition_geometry_recomputed',cells=column+1,pieces=piece_count))
    passed = (maximum['height_m'] <= mapping_tolerance_m and maximum['gap_m'] <= mapping_tolerance_m
        and maximum['shape_weight'] <= 1e-10 and maximum['ray_shape_closure_m'] <= mapping_tolerance_m
        and maximum['coefficient'] <= 1e-10 and maximum['lever_m'] <= mapping_tolerance_m
        and maximum['W_column_relative_error'] <= 1e-10)
    proof = dict(schema='native_FE_partition_geometry_audit_v1',passed=passed,geometry_only=True,geometry_source=NATIVE_GEOMETRY_SOURCE,
        integration_method=method,integration_method_alias_used=method in NATIVE_METHOD_ALIASES,cell_count=len(gaps),piece_count=piece_count,
        cells_by_tooth_offset=per_patch,maximum_reconstruction_errors=maximum,record_geometry_fingerprint=digest.hexdigest(),
        mean_gap_identity=numeric_array_identity(np.asarray(gaps)),body_nodes_identity={m:numeric_array_identity(bodies[m][0]) for m in members},
        native_cache_index_sha256=system.get('native_cache_index_sha256'),W_reconstructed_from_original_pieces=matrices is not None,
        source_Rk_used_for_mechanics=False,NURBS_gap_queried=False,CAD_error_set_from_FE_closure=False,
        maximum_native_FE_ray_shape_closure_error_m=maximum['ray_shape_closure_m'],
        loaded_finite_phi_geometry_checked=False,physical_pressure_verified=False,
        stiffness_assembled=False,factorization_performed=False,pressure_solved=False,
        definition='真实全局FE节点＋机械整体R0；完整交片面积/质心、原三角P1形函数、g0、W和lever重建；非NURBS')
    if include_recomputed_values:
        # 只供同次响应审核计算Wᵀu/力矩；证据保存数组身份后可移除完整数值。
        proof['recomputed_lever_arms_m'] = {m: np.asarray(v) for m, v in mean_levers.items()}
    if operator_identity:
        proof.update(operator_identity)
    if system.get('operator_input_fingerprint') is not None:
        _sha_identity(system['operator_input_fingerprint'], 'operator_input_fingerprint')
        proof['operator_input_fingerprint'] = system['operator_input_fingerprint']
    return np.asarray(gaps),np.asarray(keys,dtype=np.int64),proof


def audit_native_response_archive(state_path, *, run_root, native_system, meshes, materials,
        fixed_nodes, mesh_source_sha256, recheck_body=True, progress=None):
    """只读原生FE档案并重建平均g0/W/u响应，不走NURBS，不创建K或求u。

state可为JSON或gzip JSON，须为源业务接受且身份稳定的原状态。调用方从
同一实际网格档案读取meshes时提供原字节SHA；本入口另绑定全部节点/单元
数字身份、NPZ/W原字节、固定自由度及分片几何。SHA字段本身不证明物理。
默认分块回装Ku；关闭时只复核原保存残差并在证据中明确降低核查范围。
"""
    from pathlib import Path
    from types import SimpleNamespace
    import gzip
    import json
    import numpy as np
    from scipy import sparse
    from .refined_fe_postprocessing import _relative_file
    from .distributed_pressure_domain import (audit_distributed_response, distributed_body_identity,
        matrix_free_body_response, numeric_array_identity)
    state_path, root = Path(state_path).resolve(), Path(run_root).resolve()
    members = ('pinion', 'wheel')
    if not isinstance(recheck_body, bool):raise ValueError('原生Ku独立回验开关须为布尔值')
    if any(set(value) != set(members) for value in (meshes, materials, fixed_nodes, mesh_source_sha256)):
        raise ValueError('原生响应须明确提供两体网格、材料、固定节点及原网格档案SHA')
    def sha(path):
        value = hashlib.sha256()
        with Path(path).open('rb') as stream:
            while block := stream.read(4*1024*1024):value.update(block)
        return value.hexdigest()
    state_sha = sha(state_path)
    with (gzip.open(state_path, 'rt', encoding='utf-8-sig') if state_path.suffix == '.gz'
          else state_path.open('r', encoding='utf-8-sig')) as stream:
        state = json.load(stream)
    if state.get('accepted') is not True or state.get('source_identity_unchanged') is not True:
        raise ValueError('原生正式响应档案只接受源业务接受且身份稳定的状态；候选须另作诊断')
    result = dict(state['result'])
    if not native_geometry_mode(result):raise ValueError('本入口只接受明确FE_boundary原生结果')
    native_result_identity(result)
    operator_file = None
    if _native_operator_identity(result):
        if _native_operator_identity(result) != _native_operator_identity(native_system):
            raise ValueError('原生最终响应与调用方最终 operator 身份不同')
        if any(not isinstance(value.get('current_operator_manifest_file'), str) or
                not value['current_operator_manifest_file'] for value in (result, native_system)):
            raise ValueError('原生最终响应和 loader 须明确最终 operator manifest 文件路径')
        operator_file = _relative_file(root, result['current_operator_manifest_file']).resolve()
        loader_file = Path(native_system['current_operator_manifest_file'])
        if (not loader_file.is_absolute() or loader_file.resolve() != operator_file or
                sha(operator_file) != result['current_operator_manifest_sha256']):
            raise ValueError('最终 operator manifest 原文件路径或字节 SHA 不一致')
    inputs_file = root/'inputs.json'
    inputs = json.loads(inputs_file.read_text(encoding='utf-8'))
    inputs_fingerprint = hashlib.sha256(json.dumps(inputs, ensure_ascii=False, sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()
    if inputs_fingerprint != state['input_fingerprint']:
        raise ValueError('原生实际输入JSON与状态fingerprint不一致')
    if result.get('source_input_fingerprint') != inputs_fingerprint:
        raise ValueError('原生结果没有绑定本档案实际输入指纹')
    if 'native_cache_index_sha256' in result and inputs.get('native_cache_index_sha256') != result['native_cache_index_sha256']:
        raise ValueError('原生实际输入与压力结果的分片缓存SHA不同')
    expected_mesh_sha = result.get('actual_native_mesh_sha256', native_system.get('actual_native_mesh_sha256'))
    if expected_mesh_sha is None or mesh_source_sha256 != expected_mesh_sha:
        raise ValueError('原生结果与调用方实际读取的两体网格档案SHA不一致')
    if 'native_cache_index_sha256' in result and state.get('native_cache_index_sha256') != result['native_cache_index_sha256']:
        raise ValueError('原生状态与内部压力结果的分片缓存SHA不一致')
    final_file = root/'final_state_receipt.json'
    final = json.loads(final_file.read_text(encoding='utf-8'))
    # 完成后的外置receipt绑定状态字节，避免把状态的自哈希写进自身。
    if final.get('source_state_sha256') != state_sha or final.get('input_fingerprint') != state['input_fingerprint']:
        raise ValueError('原生最终receipt没有绑定本状态原字节及输入身份')
    if final.get('accepted') is not True or final.get('source_identity_unchanged') is not True:
        raise ValueError('原生最终receipt未确认源业务接受及身份稳定')
    if _relative_file(root, final['state_file']).resolve() != state_path:
        raise ValueError('原生最终receipt指向另一状态文件')
    receipt_path = root/'response_receipt.json'
    if (_relative_file(root, final['response_receipt_file']).resolve() != receipt_path or
            final.get('response_receipt_sha256') != sha(receipt_path)):
        raise ValueError('原生最终receipt绑定的原W/u/F响应receipt已改变')
    receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
    if receipt != state.get('response_receipt') or receipt.get('passed') is not True or receipt.get('input_fingerprint') != state['input_fingerprint']:
        raise ValueError('原生内嵌/外置响应receipt不同或不属于本输入')
    force = np.asarray(result['equilibrium']['forces_n'], dtype=float)
    keys = np.asarray([(f['tooth_offset'], int(j), int(i)) for f in result['fields']
        for j, i in zip(*np.nonzero(np.asarray(f['valid_mask'], dtype=bool)))], dtype=np.int64)
    body_fields = {}; systems = {}; body_evidence = {}; original_files = {}; matrices = {}
    def sparse_identity(matrix):
        value = matrix.tocsc(copy=True); value.sum_duplicates(); value.sort_indices()
        return dict(shape=list(value.shape), data=numeric_array_identity(value.data),
            indices=numeric_array_identity(value.indices, integer=True), indptr=numeric_array_identity(value.indptr, integer=True))
    for member in members:
        item = receipt['bodies'][member]
        response_path = _relative_file(root, item['response_file'])
        W_path = _relative_file(root, item['load_matrix_file'])
        if sha(response_path) != item['response_sha256'] or sha(W_path) != item['load_matrix_sha256']:
            raise ValueError('原生W或位移/载荷外存档案SHA改变：'+member)
        with np.load(response_path, allow_pickle=False) as stream:
            data = {name:stream[name] for name in stream.files}
        nodes = np.asarray(meshes[member].nodes_m, dtype=float)
        if (not np.array_equal(data['contact_cell_keys'], keys) or not np.array_equal(data['force_plane_n'], force)
                or data['displacements_m'].shape != nodes.shape or data['nodal_forces_n'].shape != nodes.shape):
            raise ValueError('原生原W列/P0力/物理格索引/节点布局不一致：'+member)
        matrix = sparse.load_npz(W_path)
        if matrix.shape != (3*len(nodes), len(force)) or matrix.shape != tuple(item['load_matrix_shape']):
            raise ValueError('原生W维度与完整体网格/压力格不一致')
        if sparse_identity(matrix) != sparse_identity(native_system['matrices'][member]):
            raise ValueError('原生档案W与独立分片重建系统不同：'+member)
        matrices[member] = matrix
        fixed = np.asarray(fixed_nodes[member])
        if fixed.ndim != 1 or not np.issubdtype(fixed.dtype, np.integer) or not len(fixed) or fixed.min() < 0 or fixed.max() >= len(nodes) or len(np.unique(fixed)) != len(fixed):
            raise ValueError('原生固定节点须为完整实际网格内的唯一整数索引')
        fixed_dofs = (3*fixed[:, None]+np.arange(3)).ravel()
        free = np.setdiff1d(np.arange(3*len(nodes)), fixed_dofs)
        if not len(free) or not np.array_equal(data['free_dofs'], free):
            raise ValueError('原生响应保存自由度与实际固定节点不一致')
        if recheck_body:
            mechanical = matrix_free_body_response(meshes[member], materials[member],
                data['displacements_m'], data['nodal_forces_n'], fixed)
        else:
            relative = float(np.linalg.norm(data['ku_minus_f_n'].reshape(-1)[free])/
                max(float(np.linalg.norm(data['nodal_forces_n'].reshape(-1)[free])), 1.))
            mechanical = dict(passed=relative <= 1e-8, relative_free_residual=relative,
                maximum_fixed_displacement_m=float(np.abs(data['displacements_m'].reshape(-1)[fixed_dofs]).max()),
                strain_energy_j=.5*float(np.sum(data['displacements_m']*(data['nodal_forces_n']+data['ku_minus_f_n']))),
                external_half_work_j=.5*float(np.sum(data['displacements_m']*data['nodal_forces_n'])),
                displacement_identity=numeric_array_identity(data['displacements_m']),
                nodal_force_identity=numeric_array_identity(data['nodal_forces_n']),
                mesh_nodes_identity=numeric_array_identity(nodes),
                mesh_tetrahedra_identity=numeric_array_identity(meshes[member].tetrahedra, integer=True))
        body_evidence[member] = mechanical
        saved = result['body_response_checks'][member]
        body_fields[member] = dict(displacements_m=data['displacements_m'], nodal_forces_n=data['nodal_forces_n'],
            **{key:saved[key] for key in ('strain_energy_j','external_work_j','maximum_free_residual_n','relative_free_residual')})
        systems[member] = SimpleNamespace(mesh=meshes[member], fixed_dofs=fixed_dofs, free_dofs=free)
        original_files[member] = dict(response_sha256=item['response_sha256'], load_matrix_sha256=item['load_matrix_sha256'])
        if progress:progress(dict(stage='native_archived_body_bound', member=member,
            relative_Ku_residual=mechanical['relative_free_residual'], independent_Ku_recheck=recheck_body))
    result['body_fields'] = body_fields
    result['distributed_body_identity'] = distributed_body_identity(result)
    system = dict(native_system, matrices=matrices)
    proof = audit_distributed_response(SimpleNamespace(systems=systems), system, result,
        body_response_evidence=body_evidence)
    proof.update(schema='native_FE_archive_response_audit_v1', source_state_sha256=state_sha,
        final_state_receipt_sha256=sha(final_file), response_receipt_sha256=sha(receipt_path),
        actual_native_mesh_sha256=dict(mesh_source_sha256), archived_files_identity=original_files,
        input_document_sha256=sha(inputs_file),
        independent_Ku_recheck=recheck_body, original_state_unchanged=sha(state_path) == state_sha,
        audit_source_sha256=sha(__file__), NURBS_gap_queried=False)
    if operator_file is not None:
        if sha(operator_file) != result['current_operator_manifest_sha256']:
            raise ArithmeticError('原生响应审核期间最终 operator manifest 原字节改变')
        proof.update(_native_operator_identity(result))
    if not proof['original_state_unchanged']:raise ArithmeticError('原生响应审核期间原状态字节改变')
    return result, proof
