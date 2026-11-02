"""由完整四面体边界提取原生FE接触面，保留整轮全局节点编号。

边界拓扑按全部四面体面键核验；工作面再由实际有限CAD面选择。
本模块不组装K、不分解、不计算压力，也不证明网格或接触域收敛。
"""
from __future__ import annotations


def verified_tetra_boundary(mesh, *, chunk_size=100000, progress=None):
    """核对Gmsh边界恰为所有tet的单邻接面，并恢复一致朝外法向。

    节点单位m；返回三角形使用原整轮全局节点编号。为采用无碰撞
    uint64面键，要求节点数小于2**21。缺失、重复、非流形或退化输入
    显式拒绝，不以声明的标签代替真实四面体邻接关系。
    """
    import numpy as np
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError('边界审核分块数须为正整数')
    nodes = np.asarray(mesh.nodes_m, dtype=float)
    def indices(value):
        raw = np.asarray(value)
        if raw.dtype.kind not in 'iu' or (raw.dtype.kind == 'u' and raw.size and raw.max() > np.iinfo(np.int64).max):
            raise ValueError('四面体、边界和标签须为整数数组，不能截断实数索引')
        return np.asarray(raw,dtype=np.int64)
    tets = indices(mesh.tetrahedra)
    given = indices(mesh.boundary_triangles)
    tags = indices(mesh.boundary_tags)
    if (nodes.ndim != 2 or nodes.shape[1] != 3 or not len(nodes) or
            not np.isfinite(nodes).all() or len(nodes) >= 2**21):
        raise ValueError('完整边界审核要求有限N×3节点，且N小于2**21')
    if tets.ndim != 2 or tets.shape[1] != 4 or not len(tets):
        raise ValueError('完整边界审核需要非空四节点四面体')
    if given.ndim != 2 or given.shape[1] != 3 or not len(given) or tags.shape != (len(given),):
        raise ValueError('需要原Gmsh完整边界三角形及一一对应标签')
    if tets.min() < 0 or tets.max() >= len(nodes) or given.min() < 0 or given.max() >= len(nodes):
        raise ValueError('原四面体或边界节点索引越界')

    def encode(tri):
        values = np.sort(tri, axis=1).astype(np.uint64)
        return (values[:,0] << np.uint64(42)) | (values[:,1] << np.uint64(21)) | values[:,2]

    keys = encode(given); order = np.argsort(keys); keys = keys[order]
    if np.any(keys[1:] == keys[:-1]):
        raise ValueError('原边界包含重复三角形')
    triangles = given[order].copy(); tags = tags[order].copy()
    matches = np.zeros(len(triangles), dtype=np.int32)
    owners = np.full(len(triangles), -1, dtype=np.int64)
    all_keys = np.empty(4*len(tets), dtype=np.uint64)
    faces = ((1,2,3),(0,2,3),(0,1,3),(0,1,2))
    for start in range(0,len(tets),chunk_size):
        block = tets[start:start+chunk_size]; points = nodes[block]
        edges = points[:,1:] - points[:,:1]
        volume6 = np.abs(np.einsum('ij,ij->i', np.cross(edges[:,0],edges[:,1]),edges[:,2]))
        scale = np.max(np.linalg.norm(edges,axis=2),axis=1)
        if np.any(volume6 <= 64*np.finfo(float).eps*scale**3):
            raise ValueError('原四面体退化，不能核验实际外边界')
        for side,face in enumerate(faces):
            face_keys = encode(block[:,face])
            all_keys[side*len(tets)+start:side*len(tets)+start+len(block)] = face_keys
            found = np.searchsorted(keys,face_keys)
            positions = np.flatnonzero(found < len(keys))
            positions = positions[keys[found[positions]] == face_keys[positions]]
            hits = found[positions]
            np.add.at(matches,hits,1); owners[hits] = start+positions
        if progress and start % (5*chunk_size) == 0:
            progress(dict(stage='native_boundary_keys',processed_tetrahedra=start+len(block),total_tetrahedra=len(tets)))
    all_keys.sort()
    begins = np.r_[0,np.flatnonzero(all_keys[1:] != all_keys[:-1])+1]
    counts = np.diff(np.r_[begins,len(all_keys)])
    if np.any(counts > 2):
        raise ValueError('原tet面存在非流形多邻接')
    actual = all_keys[begins[counts == 1]]
    if not np.array_equal(actual,keys) or np.any(matches != 1):
        raise ValueError('原边界不等于全部四面体的真实单邻接边界')
    del all_keys,begins,counts,actual
    xyz = nodes[triangles]; centroids = xyz.mean(axis=1)
    normals = np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
    magnitude = np.linalg.norm(normals,axis=1)
    outward = centroids-nodes[tets[owners]].mean(axis=1)
    flip = np.einsum('ij,ij->i',normals,outward) < 0
    triangles[flip] = triangles[flip][:,[0,2,1]]; normals[flip] *= -1
    if np.any(magnitude <= 0):
        raise ValueError('真实外边界存在退化三角形')
    return dict(nodes_m=nodes,triangles=triangles,centroids_m=centroids,normals=normals/magnitude[:,None],
        areas_m2=.5*magnitude,face_tags=tags,boundary_indices=np.arange(len(triangles),dtype=np.int64),
        metadata=dict(units='m',normal_definition='由唯一相邻四面体重心指向外部',
            boundary_definition='全部tet面键核验的完整单邻接三角面',
            source_Gmsh_boundary_completely_verified=True,large_adjacency_dictionary_used=False))


def extract_native_flank_compacts(prepared, tooth_offsets, *, progress=None):
    """从已准备的整轮FE提取各实际物理齿；不复制完整体网格到每个齿。

    返回(compact_facets, audit)。compact按pinion/wheel和有符号齿对k
    索引，记录局部表面节点、局部三角形及原整轮global_node_ids。
    小轮物理齿为−k，大轮为+k；面选择及外法向遵循实际CAD工作面。
    整体运动R0由求解时显式提供，材料显示Rk由几何模板另行保留。
    """
    from copy import deepcopy
    from types import SimpleNamespace
    import numpy as np
    from .generated_distributed_fe import _working_facets
    offsets = list(tooth_offsets)
    if (not offsets or len(set(offsets)) != len(offsets) or
            any(isinstance(k,bool) or not isinstance(k,int) for k in offsets)):
        raise ValueError('原生接触齿对偏移须为非空、不重复整数序列')
    if set(prepared.models) != {'pinion','wheel'} or set(prepared.systems) != {'pinion','wheel'}:
        raise ValueError('需要pinion/wheel两轮完整FE准备对象')
    compact = {}; audits = {}
    for member in ('pinion','wheel'):
        mesh = prepared.systems[member].mesh
        boundary = verified_tetra_boundary(mesh,progress=progress)
        local = SimpleNamespace(models=prepared.models,systems={member:prepared.systems[member]},
            _facet_cache={('distributed_boundary',member):boundary})
        pieces = {}; selections = []
        for k in offsets:
            facet = _working_facets(local,member,k)
            ids = np.unique(facet['triangles'])
            if not len(ids):
                raise ValueError('实际CAD工作面没有匹配的FE三角形：'+member+' k='+str(k))
            pieces[k] = dict(name=member,nodes_m=np.array(facet['nodes_m'][ids],copy=True),
                triangles=np.searchsorted(ids,facet['triangles']),global_node_ids=ids,
                total_node_count=len(mesh.nodes_m),outward_normals_local=np.array(facet['normals'],copy=True),
                source_selection_metadata=deepcopy(facet['metadata']),
                source_boundary_indices=np.array(facet['boundary_indices'],copy=True))
            selections.append(dict(tooth_offset=k,physical_tooth_index=int(facet['metadata']['tooth_index']),
                compact_nodes=len(ids),triangles=len(facet['triangles']),
                maximum_vertex_CAD_selection_distance_m=float(np.max(facet['maximum_vertex_matching_distance_m']))))
            if progress:
                progress(dict(stage='native_flank_compacted',member=member,**selections[-1]))
            for key in list(local._facet_cache):
                if key[0] == 'distributed_working_flank':
                    local._facet_cache.pop(key)
        compact[member] = pieces
        audits[member] = dict(source_nodes=len(mesh.nodes_m),source_tetrahedra=len(mesh.tetrahedra),
            verified_boundary_faces=len(boundary['triangles']),selections=selections)
        del local,boundary
    return compact,dict(geometry_only=True,K_assembled=False,pressure_solved=False,bodies=audits)
