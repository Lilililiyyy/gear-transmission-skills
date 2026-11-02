"""完整闭边界的远场材料域审核；距离本体和实际材料分类由显式oracle提供。

每个源三角先给出覆盖完整对面边界、已扣数值和几何余量的整片距离下界。
只有该下界大于相对位移上界的片才进入连通分量。精确共享顶点的片形成
无边界穿越的连通集，因此每个分量一次真实材料分类可以推广至该分量。
正无符号距离本身不能证明材料外；失败/缺失证据统一保留未决。
"""
from __future__ import annotations

import math
import operator


def _nonnegative(value, name):
    if isinstance(value, bool):
        raise ValueError(name+'须为非负有限数，单位m')
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(name+'须为非负有限数，单位m') from error
    if not math.isfinite(value) or value < 0:
        raise ValueError(name+'须为非负有限数，单位m')
    return value


def _triangles(values):
    try:
        result = [tuple(tuple(float(v) for v in point) for point in triangle) for triangle in values]
    except (TypeError, ValueError) as error:
        raise ValueError('源三角须为非空N×3×3坐标，单位m') from error
    if not result or any(len(t) != 3 or any(len(p) != 3 or any(not math.isfinite(v) for v in p) for p in t) for t in result):
        raise ValueError('源三角须为非空有限N×3×3坐标，单位m')
    for triangle in result:
        a = [triangle[1][i]-triangle[0][i] for i in range(3)]
        b = [triangle[2][i]-triangle[0][i] for i in range(3)]
        cross = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
        if not math.isfinite(sum(v*v for v in cross)) or sum(v*v for v in cross) <= 1e-60:
            raise ValueError('源三角退化，不能默默删除后声称完整覆盖')
    return result


def exact_vertex_components(triangles_m, selected_indices=None):
    """按精确相同坐标顶点连接，返回原三角索引分组；不用距离容差拼接。

只连接selected_indices中的三角。未获边界排除的中间片不能作为连接桥。
不相交的近邻顶点即使只差浮点末位，也保持为不同分量。
"""
    triangles = _triangles(triangles_m)
    requested = list(range(len(triangles))) if selected_indices is None else list(selected_indices)
    try:
        selected = [operator.index(index) for index in requested]
    except (TypeError, ValueError) as error:
        raise ValueError('连通片索引须为不重复且存在的整数') from error
    if (len(set(selected)) != len(selected) or any(isinstance(i, bool) for i in requested)
            or any(not 0 <= i < len(triangles) for i in selected)):
        raise ValueError('连通片索引须为不重复且存在的整数')
    parent = {i: i for i in selected}
    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    vertices = {}
    for index in selected:
        for point in triangles[index]:
            if point in vertices:
                parent[find(index)] = find(vertices[point])
            else:
                vertices[point] = index
    groups = {}
    for index in selected:
        groups.setdefault(find(index), []).append(index)
    return sorted((sorted(group) for group in groups.values()), key=lambda group: group[0])


def audit_triangle_material_domain(triangles_m, relative_motion_upper_m, *, distance_lower_oracle,
        material_classifier, boundary_sha256, boundary_is_closed):
    """构建三角远域的完整边界距离＋实际材料分类证据，不重新求解接触。

triangles_m为N×3×3源片坐标m；relative_motion_upper_m给每片相对位移范数
上界m，包含调用方需要的刚体余项及弹性运动。distance_lower_oracle(i,t)
须返回done、complete_target_search和distance_lower_bound_m；下界必须
针对完整对面闭边界，并已扣双面离散误差与查询余量。没有完成或完整候选
证据时，该片不会被认证。

material_classifier(point_m)必须查询同一加载闭体，返回state为outside、
inside、on或unknown及可选原生证据。仅outside可排除；不以法向或无符号
点距猜测材料内外。boundary_sha256绑定两个oracle共同的实际几何。

返回每片与每分量的原证据。此函数不替代距离oracle、CAD闭体有效性、
连续曲面误差界或接触力学验证；调用方须保留上述绑定和实际计算来源。
"""
    triangles = _triangles(triangles_m)
    if boundary_is_closed is not True:
        raise ValueError('连通材料传播需要实际完整闭边界证据，不能使用单个开放工作面')
    if (not isinstance(boundary_sha256, str) or len(boundary_sha256) != 64
            or any(c not in '0123456789abcdef' for c in boundary_sha256)):
        raise ValueError('须绑定两个oracle共同的实际闭体SHA256')
    motions = [_nonnegative(v, '相对运动上界') for v in relative_motion_upper_m]
    if len(motions) != len(triangles):
        raise ValueError('每个源三角必须有且只有一个相对运动上界')
    rows, selected = [], []
    for index, (triangle, motion) in enumerate(zip(triangles, motions)):
        query = distance_lower_oracle(index, triangle)
        if not isinstance(query, dict):
            raise ValueError('整片距离oracle须返回结构化完成状态及完整候选证据')
        lower = None
        if query.get('done') is True and query.get('complete_target_search') is True:
            lower = _nonnegative(query.get('distance_lower_bound_m'), '整片距离下界')
        separated = lower is not None and lower > motion
        rows.append(dict(triangle_index=index, distance_lower_bound_m=lower,
            relative_motion_upper_bound_m=motion, whole_boundary_separated=bool(separated),
            distance_evidence=dict(query), state='unresolved'))
        if separated:
            selected.append(index)
    components = []
    for number, members in enumerate(exact_vertex_components(triangles, selected)):
        representative = min(members, key=lambda i: (rows[i]['distance_lower_bound_m']-motions[i], i))
        point = tuple(sum(vertex[j] for vertex in triangles[representative])/3. for j in range(3))
        classification = material_classifier(point)
        if not isinstance(classification, dict) or classification.get('state') not in ('outside', 'inside', 'on', 'unknown'):
            raise ValueError('材料分类须显式返回outside/inside/on/unknown，不能只提供无符号距离')
        outside = classification['state'] == 'outside'
        components.append(dict(component_index=number, triangle_indices=members,
            representative_triangle_index=representative, representative_point_m=point,
            classification=dict(classification), outside=outside, exact_shared_vertices_only=True,
            all_triangles_whole_boundary_separated=True))
        for index in members:
            rows[index]['component_index'] = number
            rows[index]['state'] = 'outside_with_clearance' if outside else 'unresolved'
    unknown = [row['triangle_index'] for row in rows if row['state'] == 'unresolved']
    return dict(schema='triangle_material_domain_v1', boundary_sha256=boundary_sha256,
        boundary_is_closed=True, triangle_count=len(triangles), rows=rows, components=components,
        unresolved_triangle_indices=unknown, passed=not unknown,
        source_pressure_modified=False, strict_continuous_cad_proof=False,
        scope='调用方已绑定的完整闭边界整片下界与实际材料分类；压力和离散收敛另行验证')
