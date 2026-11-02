"""真实双FE交片的批量精确P0积分与typed数组，不求接触QP。

数学与冻结projected_fe_partition一致：完整P0矩形的每单位平面合力，
在所有真实正面积凸片上积分P1位移迹。部分覆盖不归一化、不补零。
原单值/连续图先验仍由冻结ProjectedFacetBody执行；本模块只批量化
交片、重心坐标、g/W/力臂与固定仿射片顶点检查，不宣称全周期已验证。
"""
from __future__ import annotations

from numbers import Integral
import math
import time


def _point_data(body,faces,xy,*,barycentric_tolerance):
    """给定真实facet与同投影xy，批量求原P1形函数，不作最近点或夹取。"""
    import numpy as np
    faces=np.asarray(faces,dtype=np.int64);xy=np.asarray(xy,dtype=float)
    if xy.shape!=(len(faces),2):raise ValueError('facet/xy批量布局不一致')
    difference=xy-body.projected_origins[faces]
    ab=np.einsum('pij,pj->pi',body.projected_inverse[faces],difference)
    bary=np.column_stack((1-ab.sum(axis=1),ab))
    if len(bary) and (bary.min() < -barycentric_tolerance or bary.max()>1+barycentric_tolerance):
        raise ArithmeticError('真实交片点越过原FE三角，不夹取重心坐标')
    nodes=body.triangles[faces]
    height=np.einsum('pi,pi->p',bary,body.heights[nodes])
    world=body.origin+xy@body.basis+height[:,None]*body.normal
    via=np.einsum('pi,pij->pj',bary,body.world_nodes[nodes])
    return dict(global_node_ids=body.node_ids[nodes],shape_weights=bary,height_m=height,
        mapped_world_m=world,closure_error_m=np.linalg.norm(via-world,axis=1))


class ArrayProjectedFEIntegrator:
    """复用已验证真实单值FE图；每个交片的数据按数组集中存储。"""
    def __init__(self,first,second,plane,*,minimum_normal_cosine=1e-6,
            overlap_relative_tolerance=256*2.220446049250313e-16,barycentric_tolerance=5e-11):
        import numpy as np
        from .projected_fe_body_arrays import ArrayProjectedFacetBody
        origin=np.asarray(plane['origin_world_m'],dtype=float);basis=np.asarray(plane['tangent_rows'],dtype=float);normal=np.asarray(plane['normal_world'],dtype=float)
        if origin.shape!=(3,) or basis.shape!=(2,3) or normal.shape!=(3,) or not all(np.isfinite(x).all() for x in (origin,basis,normal)):
            raise ValueError('共同平面须为有限三维原点、2×3切基及三维法向')
        if not np.allclose(basis@basis.T,np.eye(2),rtol=0,atol=1e-11) or not np.allclose(np.cross(*basis),normal,rtol=0,atol=1e-11) or abs(float(normal@normal)-1)>1e-11:
            raise ValueError('共同平面须为正向正交单位基')
        for x in (minimum_normal_cosine,overlap_relative_tolerance,barycentric_tolerance):
            if isinstance(x,bool) or not math.isfinite(x) or x<=0:raise ValueError('投影审核门限须为正有限值')
        if first['name']==second['name']:raise ValueError('两体名称不能相同')
        options=dict(minimum_normal_cosine=minimum_normal_cosine,overlap_relative_tolerance=overlap_relative_tolerance,barycentric_tolerance=barycentric_tolerance)
        self.bodies=[ArrayProjectedFacetBody(data,(origin,basis,normal),**options) for data in (first,second)]
        self.barycentric_tolerance=barycentric_tolerance;self.overlap_relative_tolerance=overlap_relative_tolerance

    @classmethod
    def from_verified_bodies(cls,bodies,*,barycentric_tolerance=5e-11,overlap_relative_tolerance=256*2.220446049250313e-16):
        """仅复用冻结类已经验过的图，不接受绕过先验的任意假对象。"""
        import numpy as np
        from .projected_fe_partition import ProjectedFacetBody
        if any(isinstance(x,bool) or not math.isfinite(x) or x<=0 for x in (barycentric_tolerance,overlap_relative_tolerance)):
            raise ValueError('重用图的审核门限须为正有限值')
        if len(bodies)!=2 or any(not isinstance(b,ProjectedFacetBody) for b in bodies) or bodies[0].name==bodies[1].name:
            raise ValueError('须提供两侧已验证ProjectedFacetBody')
        for field in ('origin','basis','normal'):
            if not np.array_equal(getattr(bodies[0],field),getattr(bodies[1],field)):raise ValueError('两侧真实图不共用相同投影平面')
        result=cls.__new__(cls);result.bodies=list(bodies)
        result.barycentric_tolerance=barycentric_tolerance;result.overlap_relative_tolerance=overlap_relative_tolerance
        return result

    def build_cells(self,cells,*,coverage_relative_tolerance=2e-12,maximum_pieces=1000000):
        """保留所有正面积交片；返回稀疏W、轻量header、typed片/顶点数组。

        不直接接入原QP。cell_statuses的缺失GeoJSON与partial记录为真实
        不完整区域；只有完整cell进入W，partial的原交片仍保留在arrays。
        """
        import numpy as np
        import shapely
        from scipy import sparse
        from shapely.strtree import STRtree
        started=time.perf_counter()
        if isinstance(coverage_relative_tolerance,bool) or not math.isfinite(coverage_relative_tolerance) or coverage_relative_tolerance<=0:
            raise ValueError('覆盖相对容差须为正有限值')
        if isinstance(maximum_pieces,bool) or not isinstance(maximum_pieces,int) or maximum_pieces<1:raise ValueError('交片上限须为正整数')
        requests=list(cells);count=len(requests);keys=[];bounds=[]
        for index,cell in enumerate(requests):
            key=tuple(cell.get(k,d) for k,d in (('tooth_offset',0),('row',0),('column',index)))
            value=np.asarray(cell['bounds_xy_m'],dtype=float)
            if any(isinstance(k,bool) or not isinstance(k,Integral) for k in key):raise ValueError('物理齿/行列须为真实整数索引')
            if value.shape!=(4,) or not np.isfinite(value).all() or value[2]<=value[0] or value[3]<=value[1]:raise ValueError('压力矩形须为有限正面积')
            keys.append(tuple(int(k) for k in key));bounds.append(value)
        if len(set(keys))!=count:raise ValueError('压力cell物理索引重复')
        bounds=np.asarray(bounds,dtype=float).reshape((-1,4))
        polygons=shapely.box(bounds[:,0],bounds[:,1],bounds[:,2],bounds[:,3]);areas=np.asarray(shapely.area(polygons))
        cell_tree=STRtree(polygons);pairs=cell_tree.query(polygons,predicate='intersects')
        if pairs.size:
            first,second=pairs[:,pairs[1]>pairs[0]]
            overlap=shapely.area(shapely.intersection(polygons[first],polygons[second]))
            if np.any(overlap>self.overlap_relative_tolerance*np.minimum(areas[first],areas[second])):raise ValueError('压力cell存在正面积重叠')
        bodies=self.bodies;body_polygons=[np.asarray(b.polygons,dtype=object) for b in bodies]
        pairs=bodies[0].tree.query(polygons,predicate='intersects')
        cell1,face1=pairs;one=shapely.intersection(polygons[cell1],body_polygons[0][face1]);positive=shapely.area(one)>0
        one=one[positive];cell1=cell1[positive];face1=face1[positive]
        pairs=bodies[1].tree.query(one,predicate='intersects');q1,face2=pairs
        two=shapely.intersection(one[q1],body_polygons[1][face2])
        pieces,owner=shapely.get_parts(two,return_index=True);positive=shapely.area(pieces)>0
        pieces=pieces[positive];owner=owner[positive]
        if len(pieces)>maximum_pieces:raise ValueError('真实正面积交片数超过显式上限')
        if len(pieces) and (np.any(shapely.get_type_id(pieces)!=3) or np.any(shapely.get_num_interior_rings(pieces))):
            raise ArithmeticError('三凸域正面积交集不是无孔Polygon，不能删片或用单质心近似')
        cell_ids=cell1[q1[owner]];facet_ids=[face1[q1[owner]],face2[owner]]
        # 按请求cell稳定分组；不合并不同的真实facet交片。
        order=np.argsort(cell_ids,kind='stable');pieces=pieces[order];cell_ids=cell_ids[order]
        facet_ids=[v[order] for v in facet_ids]
        piece_areas=np.asarray(shapely.area(pieces));xy=np.asarray(shapely.get_coordinates(shapely.centroid(pieces))).reshape((-1,2))
        weights=piece_areas/areas[cell_ids]
        body_data={};maximum_closure=np.zeros(count)
        for body,facets,sign in zip(bodies,facet_ids,(-1.,1.)):
            data=_point_data(body,facets,xy,barycentric_tolerance=self.barycentric_tolerance)
            coefficient_world=sign*body.normals[facets]/body.cosines[facets,None]
            coefficient_local=coefficient_world@body.rotation
            lever=-np.einsum('pi,pi->p',coefficient_world,np.cross(body.axis,data['mapped_world_m']-body.axis_origin))
            data.update(facet_indices=facets,coefficients_local=coefficient_local,coefficients_world=coefficient_world,
                levers_m=lever,actual_area_piece_m2=piece_areas/np.abs(body.cosines[facets]))
            body_data[body.name]=data
            np.maximum.at(maximum_closure,cell_ids,data['closure_error_m'])
        gap=body_data[bodies[1].name]['height_m']-body_data[bodies[0].name]['height_m']
        mean_gap=np.bincount(cell_ids,weights=weights*gap,minlength=count)
        mean_levers={b.name:np.bincount(cell_ids,weights=weights*body_data[b.name]['levers_m'],minlength=count) for b in bodies}
        actual_areas={b.name:np.bincount(cell_ids,weights=body_data[b.name]['actual_area_piece_m2'],minlength=count) for b in bodies}
        piece_counts=np.bincount(cell_ids,minlength=count);offsets=np.r_[0,np.cumsum(piece_counts)]
        statuses=[];complete=np.zeros(count,dtype=bool)
        for index,polygon in enumerate(polygons):
            selected=pieces[offsets[index]:offsets[index+1]]
            covered=shapely.union_all(selected);missing=polygon.difference(covered);excess=covered.difference(polygon)
            summed=float(piece_areas[offsets[index]:offsets[index+1]].sum());overlap=max(0.,summed-float(covered.area))
            if overlap>coverage_relative_tolerance*areas[index]:raise ArithmeticError('同cell真实交片有超限面积重叠')
            complete[index]=missing.area<=coverage_relative_tolerance*areas[index] and excess.area<=coverage_relative_tolerance*areas[index] and abs(summed-areas[index])<=coverage_relative_tolerance*areas[index]
            statuses.append(dict(requested_cell_index=index,tooth_offset=keys[index][0],row=keys[index][1],column=keys[index][2],
                status='complete' if complete[index] else 'partial',area_plane_m2=float(areas[index]),covered_area_plane_m2=float(covered.area),
                summed_piece_area_plane_m2=summed,missing_area_plane_m2=float(missing.area),overlap_area_plane_m2=overlap,
                missing_geometry_geojson=missing.__geo_interface__,piece_count=int(piece_counts[index]),
                integration_weight_sum=float(weights[offsets[index]:offsets[index+1]].sum()),record_index=None))
        cell_to_column=np.full(count,-1,dtype=np.int64);cell_to_column[complete]=np.arange(int(complete.sum()))
        piece_columns=cell_to_column[cell_ids];usable=piece_columns>=0;matrices={}
        for body in bodies:
            data=body_data[body.name]
            rows=(3*data['global_node_ids'][usable,:,None]+np.arange(3)).ravel()
            columns=np.repeat(piece_columns[usable],9)
            values=(weights[usable,None,None]*data['shape_weights'][usable,:,None]*data['coefficients_local'][usable,None,:]).ravel()
            matrices[body.name]=sparse.coo_matrix((values,(rows,columns)),shape=(3*body.total_node_count,int(complete.sum()))).tocsc()
        # 与冻结逐片实现一样保留完整cell四角；批量射线只取实际相交面，
        # 在共享边上任一相邻面给出相同P1位置，不外推或夹取。
        cb=bounds[complete];corner_xy=np.stack((cb[:,[0,1]],cb[:,[2,1]],cb[:,[2,3]],cb[:,[0,3]]),axis=1).reshape((-1,2))
        corner_query=bodies[0].tree.query(shapely.points(corner_xy),predicate='intersects')
        hit,first_hit=np.unique(corner_query[0],return_index=True)
        if not np.array_equal(hit,np.arange(len(corner_xy))):raise ArithmeticError('完整cell角点缺真实FE覆盖，不补虚假射线')
        corner_world=_point_data(bodies[0],corner_query[1,first_hit],corner_xy,barycentric_tolerance=self.barycentric_tolerance)['mapped_world_m'].reshape((-1,4,3))
        records=[];partial_records=[]
        for index,status in enumerate(statuses):
            record=dict(tooth_offset=keys[index][0],row=keys[index][1],column=keys[index][2],requested_cell_index=index,
                gap_m=float(mean_gap[index]),area_plane_m2=float(areas[index]),
                actual_area_m2={b.name:float(actual_areas[b.name][index]) for b in bodies},
                levers_m={b.name:float(mean_levers[b.name][index]) for b in bodies},
                geometric_levers_m={b.name:float(mean_levers[b.name][index]) for b in bodies},
                maximum_mapping_error_m=float(maximum_closure[index]),coverage=status)
            if complete[index]:
                status['record_index']=len(records);record['first_corner_points_world_m']=corner_world[len(records)]
                records.append(record)
            else:record['interpretation']='不满格原积分诊断；权重未归一化，不进入QP';partial_records.append(record)
        # 原始每片外环去掉重复闭合点，保留ragged偏移；所有正面积片均保留。
        coordinates=np.asarray(shapely.get_coordinates(pieces)).reshape((-1,2));vertex_counts=np.asarray(shapely.get_num_coordinates(pieces),dtype=np.int64)
        end=np.cumsum(vertex_counts);keep=np.ones(len(coordinates),dtype=bool)
        if len(end):keep[end-1]=False
        vertices=coordinates[keep];vertex_offsets=np.r_[0,np.cumsum(vertex_counts-1)]
        return dict(schema='projected_FE_partition_typed_arrays_v1',records=records,partial_records=partial_records,
            matrices=matrices,cell_statuses=statuses,requested_cell_keys=np.asarray(keys,dtype=np.int64).reshape((-1,3)),
            cell_bounds_xy_m=bounds,cell_areas_m2=areas,cell_complete_mask=complete,cell_to_column=cell_to_column,
            piece_cell_indices=cell_ids,piece_facet_indices={b.name:ids for b,ids in zip(bodies,facet_ids)},
            piece_centroids_xy_m=xy,piece_areas_m2=piece_areas,piece_integration_weights=weights,piece_gap_m=gap,
            complete_cell_first_corner_points_world_m=corner_world,
            piece_vertices_xy_m=vertices,piece_vertex_offsets=vertex_offsets,body_data=body_data,
            body_projection_tables={b.name:dict(projected_origins=b.projected_origins,projected_inverse=b.projected_inverse,
                local_triangles=b.triangles,local_node_heights_m=b.heights,global_node_ids=b.node_ids,
                rotation_rows=b.rotation,axis_world=b.axis,axis_origin_world_m=b.axis_origin,
                plane_origin_world_m=b.origin,plane_basis_rows=b.basis,plane_normal_world=b.normal) for b in bodies},
            projection_audits={b.name:b.projection_audit for b in bodies},piece_count=len(pieces),requested_cell_count=count,
            complete_cell_count=int(complete.sum()),partial_cell_count=int((~complete).sum()),
            original_NURBS_gap_used=False,closest_point_clamping_used=False,partial_cell_renormalization_used=False,
            geometry_source='FE_boundary',integration_method='P0_cell_intersection_with_both_actual_FE_triangles_centroid_exact',
            quadrature_object_dictionaries_created=False,QP_connected=False,K_assembled=False,factorizations=0,
            barycentric_tolerance=self.barycentric_tolerance,piece_vertex_shape_validation_performed=False,
            seconds=time.perf_counter()-started)


def evaluate_array_piece_vertices(system,displacements_local_m,theta_rad,*,follower_member='wheel',return_vertex_values=False):
    """向量化核查固定P1交片的所有顶点；不作加载后有限旋转重建。"""
    import numpy as np
    if isinstance(theta_rad,bool) or not math.isfinite(theta_rad):raise ValueError('线性化角须为有限rad')
    members=set(system['matrices'])
    if set(displacements_local_m)!=members or follower_member not in members:raise ValueError('位移体成员或从动体与原W不符')
    displacements={}
    for m,w in system['matrices'].items():
        u=np.asarray(displacements_local_m[m],dtype=float)
        if u.size!=w.shape[0] or not np.isfinite(u).all():raise ValueError('位移须为原真实全局节点布局')
        displacements[m]=u.reshape((-1,3))
    counts=np.diff(system['piece_vertex_offsets']);piece_of_vertex=np.repeat(np.arange(system['piece_count']),counts)
    xy=system['piece_vertices_xy_m'];heights={};opening=np.zeros(len(xy));lever=None
    for m in members:
        table=system['body_projection_tables'][m];data=system['body_data'][m]
        faces=data['facet_indices'][piece_of_vertex]
        ab=np.einsum('vij,vj->vi',table['projected_inverse'][faces],xy-table['projected_origins'][faces])
        shapes=np.column_stack((1-ab.sum(axis=1),ab));nodes=table['local_triangles'][faces]
        if len(shapes) and (shapes.min() < -system['barycentric_tolerance'] or shapes.max()>1+system['barycentric_tolerance']):
            raise ArithmeticError('真实交片顶点越过原FE三角，不夹取形函数')
        heights[m]=np.einsum('vi,vi->v',shapes,table['local_node_heights_m'][nodes])
        global_ids=table['global_node_ids'][nodes];coef=data['coefficients_local'][piece_of_vertex]
        opening+=np.einsum('vi,vij,vj->v',shapes,displacements[m][global_ids],coef)
        if m==follower_member:
            world=table['plane_origin_world_m']+xy@table['plane_basis_rows']+heights[m][:,None]*table['plane_normal_world']
            lever=-np.einsum('vi,vi->v',data['coefficients_world'][piece_of_vertex],np.cross(table['axis_world'],world-table['axis_origin_world_m']))
    ordered=list(system['body_data']);gap=heights[ordered[1]]-heights[ordered[0]]+opening-theta_rad*lever
    vertex_cell_ids=system['piece_cell_indices'][piece_of_vertex];minimum=np.full(system['requested_cell_count'],np.inf)
    np.minimum.at(minimum,vertex_cell_ids,gap)
    means=sum(np.asarray(w.T@displacements[m].ravel()).ravel() for m,w in system['matrices'].items())
    rows=[]
    for index,record in enumerate(system['records']):
        request=record['requested_cell_index'];value=record['gap_m']+float(means[index])-theta_rad*record['levers_m'][follower_member]
        rows.append(dict(record_index=index,tooth_offset=record['tooth_offset'],row=record['row'],column=record['column'],
            mean_remaining_gap_m=value,minimum_remaining_gap_m=float(minimum[request]),
            checked_piece_vertices=int(counts[system['piece_cell_indices']==request].sum())))
    result=dict(minimum_remaining_gap_m=min((r['minimum_remaining_gap_m'] for r in rows),default=None),cell_checks=rows,
        fixed_piecewise_affine_extrema_evaluated=True,loaded_finite_rotation_rebuilt=False,whole_flank_nonpenetration_certified=False,
        piece_vertex_shape_validation_performed=True,
        partial_cells_excluded_from_QP_but_preserved=True)
    if return_vertex_values:result.update(vertex_remaining_gaps_m=gap,piece_of_vertex=piece_of_vertex)
    return result


def typed_array_payload(system):
    """导出普通NPZ可保存的原数组；不含pickle、GEOS对象或逐q字典。"""
    import numpy as np
    result={k:v for k,v in system.items() if isinstance(v,np.ndarray)}
    for group in ('body_data','body_projection_tables'):
        for member,fields in system[group].items():
            for name,value in fields.items():result[group+'__'+member+'__'+name]=np.asarray(value)
    return result
