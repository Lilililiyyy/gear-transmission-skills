"""全量真实FE投影图单值性和连续性先验的批量实现，不裁剪facet。

保留冻结ProjectedFacetBody的数据、法向、旋转、索引及阈值；把逐face/
逐pair的GEOS交集与高度连续性检查改为数组运算。mapping/barycentric/
ray直接继承原实现。新初始化不冒称执行过原逐pair循环，独立接受对照。
"""
from __future__ import annotations

import math
from .projected_fe_partition import ProjectedFacetBody,_vector,_validated_plane


class ArrayProjectedFacetBody(ProjectedFacetBody):
    """完整facet集合的批量单值/连续先验；无FE体/K/压力计算。"""
    def __init__(self,data,plane,*,minimum_normal_cosine=1e-6,overlap_relative_tolerance=256*2.220446049250313e-16,barycentric_tolerance=5e-11):
        import numpy as np
        import shapely
        from shapely.strtree import STRtree
        for value in (minimum_normal_cosine,overlap_relative_tolerance,barycentric_tolerance):
            if isinstance(value,bool) or not math.isfinite(value) or value<=0:raise ValueError('真实图审核门限须为正有限值')
        self.name=data['name'];self.origin,self.basis,self.normal=_validated_plane(plane)
        self.barycentric_tolerance=barycentric_tolerance
        self.nodes=np.asarray(data['nodes_m'],dtype=float);raw=np.asarray(data['triangles'])
        if self.nodes.ndim!=2 or self.nodes.shape[1]!=3 or not len(self.nodes) or not np.isfinite(self.nodes).all():raise ValueError('真实facet节点须为有限N×3米坐标')
        if raw.ndim!=2 or raw.shape[1]!=3 or not len(raw) or raw.dtype.kind not in 'iu':raise ValueError('真实facet须为非空整数M×3节点索引')
        self.triangles=raw.astype(int)
        if self.triangles.min()<0 or self.triangles.max()>=len(self.nodes):raise ValueError('真实facet节点索引越界')
        self.rotation=np.asarray(data.get('rotation_rows',np.eye(3)),dtype=float)
        if self.rotation.shape!=(3,3) or not np.isfinite(self.rotation).all() or not np.allclose(self.rotation.T@self.rotation,np.eye(3),rtol=0,atol=1e-11) or np.linalg.det(self.rotation)<=0:raise ValueError('R须为三维正交正向旋转')
        self.translation=_vector(data.get('translation_world_m',[0.,0.,0.]),3,'刚体平移')
        self.axis=_vector(data.get('axis_world',[0.,0.,1.]),3,'旋转轴')
        if abs(float(self.axis@self.axis)-1.)>1e-11:raise ValueError('旋转轴须为单位向量，角度单位rad')
        self.axis_origin=_vector(data.get('axis_origin_world_m',[0.,0.,0.]),3,'旋转轴原点')
        ids=np.asarray(data.get('global_node_ids',np.arange(len(self.nodes))))
        if ids.shape!=(len(self.nodes),) or ids.dtype.kind not in 'iu' or ids.min()<0 or len(np.unique(ids))!=len(ids):raise ValueError('全局节点号须为非负唯一整数，保持真实FE节点对应')
        self.node_ids=ids.astype(int);self.total_node_count=data.get('total_node_count',int(ids.max())+1)
        if isinstance(self.total_node_count,bool) or not isinstance(self.total_node_count,int) or self.total_node_count<=ids.max():raise ValueError('真实体总节点数不足以包含facet全局节点号')
        self.world_nodes=self.nodes@self.rotation.T+self.translation
        relative=self.world_nodes-self.origin;self.xy_nodes=relative@self.basis.T;self.heights=relative@self.normal
        xyz=self.world_nodes[self.triangles];xy=self.xy_nodes[self.triangles]
        cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);norm=np.linalg.norm(cross,axis=1)
        edge_scale=np.max(np.linalg.norm(xyz-np.roll(xyz,1,axis=1),axis=2),axis=1)
        if np.any(norm<=64*np.finfo(float).eps*edge_scale**2):raise ValueError('真实FE面三角退化')
        normals=cross/norm[:,None]
        if 'outward_normals_local' in data:
            supplied=np.asarray(data['outward_normals_local'],dtype=float)@self.rotation.T
            if supplied.shape!=normals.shape or not np.isfinite(supplied).all():raise ValueError('外法向布局须与真实facet一致')
            lengths=np.linalg.norm(supplied,axis=1)
            if np.any(np.abs(lengths-1)>1e-10):raise ValueError('外法向须为真实FE面的单位法向')
            alignment=np.einsum('ij,ij->i',supplied,normals)
            if np.any(np.abs(alignment)<1-1e-10):raise ValueError('传入法向不平行真实线性FE面，不能使用CAD法向代替')
            normals*=np.where(alignment<0,-1.,1.)[:,None]
        self.normals=normals;self.cosines=normals@self.normal
        if np.any(np.abs(self.cosines)<minimum_normal_cosine):raise ValueError('FE面相对投影射线退化，不夹取法向分母')
        if np.any(self.cosines*self.cosines[0]<=0):raise ValueError('同体facet投影折叠或外法向方向不一致')
        self.polygons=shapely.polygons(xy)
        areas=np.asarray(shapely.area(self.polygons))
        if np.any(~shapely.is_valid(self.polygons)) or np.any(areas<=0):raise ValueError('真实FE投影三角无效或面积退化')
        self.tree=STRtree(self.polygons)
        self.projected_origins=xy[:,0]
        self.projected_inverse=np.linalg.inv(np.stack((xy[:,1]-xy[:,0],xy[:,2]-xy[:,0]),axis=2))
        self.face_edge_lengths=np.linalg.norm(xyz-np.roll(xyz,1,axis=1),axis=2)
        pairs=self.tree.query(self.polygons,predicate='intersects');first,second=pairs[:,pairs[1]>pairs[0]]
        intersection=shapely.intersection(self.polygons[first],self.polygons[second])
        intersection_area=np.asarray(shapely.area(intersection));limit=overlap_relative_tolerance*np.minimum(areas[first],areas[second])
        if np.any(intersection_area>limit):raise ValueError('同体FE三角投影存在正面积重叠，非单值图')
        self.roundoff_overlap_area_m2=float(intersection_area.sum())
        # get_coordinates会取复合几何中所有分量的真实坐标，比原冻结循环
        # 跳过Multi*/GeometryCollection更强；该差别以计数单独报告。
        coordinates,owner=shapely.get_coordinates(intersection,return_index=True)
        scale_global=max(float(np.max(np.abs(relative))),1e-300)
        if len(coordinates):
            face_first=first[owner];face_second=second[owner]
            def height(facets):
                ab=np.einsum('pij,pj->pi',self.projected_inverse[facets],coordinates-self.projected_origins[facets])
                bary=np.column_stack((1-ab.sum(axis=1),ab))
                return np.einsum('pi,pi->p',bary,self.heights[self.triangles[facets]])
            difference=np.abs(height(face_first)-height(face_second))
            scales=np.maximum(np.maximum(edge_scale[face_first],edge_scale[face_second]),scale_global)
            allowed=512*np.finfo(float).eps*scales
            if np.any(difference>allowed):raise ValueError('同体投影共享边/点对应不同FE高度，非连续单值图')
            maximum_difference=float(difference.max());maximum_ratio=float(np.max(difference/allowed))
        else:maximum_difference=0.;maximum_ratio=0.
        types=np.asarray(shapely.get_type_id(intersection));additional=int(np.sum((types>=4)&~shapely.is_empty(intersection)))
        self.projection_audit=dict(facet_count=len(self.polygons),queried_neighbor_pairs=len(first),
            tolerated_roundoff_overlap_area_m2=self.roundoff_overlap_area_m2,
            minimum_absolute_normal_cosine=float(np.min(np.abs(self.cosines))),
            original_FE_linear_facets_used=True,closest_point_clamping_used=False,
            precondition_backend='batched_full_facet_STRtree_pair_intersections_v1',
            all_original_facets_retained=True,checked_shared_intersection_coordinates=len(coordinates),
            additional_collection_intersections_checked=additional,
            collection_check_is_stronger_than_frozen_skip=True,
            maximum_height_continuity_error_m=maximum_difference,
            maximum_continuity_error_over_original_limit=maximum_ratio,
            area_threshold_same_as_frozen=True,height_threshold_same_as_frozen=True,
            K_assembled=False,factorizations=0)
