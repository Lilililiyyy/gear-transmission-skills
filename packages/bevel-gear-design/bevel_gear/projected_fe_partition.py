"""真实线性FE边界的投影分片精确P0积分，不创建或分解FE体刚度。

两侧边界须在共同平面上为单值P1图。cell∩tri1∩tri2是凸片，片内
FE高度、形函数、刚体速度均为仿射函数，面积与质心足以精确积分。
本分支明确使用真实FE线性边界gap，既不调用NURBS也不最近点夹取。
"""
from __future__ import annotations

import math
from numbers import Integral


def _vector(value,length,name):
    import numpy as np
    result=np.asarray(value,dtype=float)
    if result.shape!=(length,) or not np.isfinite(result).all():raise ValueError(name+'须为有限同长向量')
    return result


def _polygons(geometry):
    """只取真实二维正面积片，线/点交集的面积严格为零。"""
    if geometry.is_empty:return
    if geometry.geom_type=='Polygon':
        if geometry.area>0:yield geometry
    elif hasattr(geometry,'geoms'):
        for member in geometry.geoms:yield from _polygons(member)


def _validated_plane(plane):
    """公共构造器统一核验正向正交投影平面，输入坐标单位为米。"""
    import numpy as np
    if isinstance(plane,dict):plane=(plane['origin_world_m'],plane['tangent_rows'],plane['normal_world'])
    if len(plane)!=3:raise ValueError('共同投影平面须包含原点、两个切基与法向')
    origin=_vector(plane[0],3,'投影原点');basis=np.asarray(plane[1],dtype=float);normal=_vector(plane[2],3,'投影法向')
    if basis.shape!=(2,3) or not np.isfinite(basis).all() or abs(float(normal@normal)-1)>1e-11 or not np.allclose(basis@basis.T,np.eye(2),rtol=0,atol=1e-11) or not np.allclose(np.cross(*basis),normal,rtol=0,atol=1e-11):
        raise ValueError('共同投影basis/n须构成正向正交单位坐标系')
    return origin,basis,normal


class ProjectedFacetBody:
    """保存真实facet、全局节点号及刚体姿态；无自由度/形函数虚构。"""
    def __init__(self,data,plane,*,minimum_normal_cosine=1e-6,overlap_relative_tolerance=256*2.220446049250313e-16,barycentric_tolerance=5e-11):
        import numpy as np
        from shapely.geometry import Polygon
        from shapely.strtree import STRtree
        for value in (minimum_normal_cosine,overlap_relative_tolerance,barycentric_tolerance):
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:raise ValueError('真实图审核门限须为正有限值')
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
        self.polygons=[Polygon(tri) for tri in xy]
        if any(not p.is_valid or p.area<=0 for p in self.polygons):raise ValueError('真实FE投影三角无效或面积退化')
        self.tree=STRtree(self.polygons)
        self.projected_origins=xy[:,0]
        self.projected_inverse=np.linalg.inv(np.stack((xy[:,1]-xy[:,0],xy[:,2]-xy[:,0]),axis=2))
        self.face_edge_lengths=np.linalg.norm(xyz-np.roll(xyz,1,axis=1),axis=2)
        self.roundoff_overlap_area_m2=0.;pairs=0
        # 同体投影须是单值图；既不以凸包填域，也不把重叠面选最近的一片。
        for first,p in enumerate(self.polygons):
            for second in self.tree.query(p,predicate='intersects'):
                second=int(second)
                if second<=first:continue
                intersection=p.intersection(self.polygons[second]);pairs+=1
                limit=overlap_relative_tolerance*min(p.area,self.polygons[second].area)
                if intersection.area>limit:raise ValueError('同体FE三角投影存在正面积重叠，非单值图')
                self.roundoff_overlap_area_m2+=intersection.area
                if intersection.is_empty:continue
                # 共享边/点也必须具有一致FE高度；不接受投影边界上的双层面。
                if intersection.geom_type=='Polygon':points=list(intersection.exterior.coords)
                elif intersection.geom_type=='LineString':points=list(intersection.coords)
                elif intersection.geom_type=='Point':points=[intersection.coords[0]]
                else:continue
                h1=self.barycentric(first,points)@self.heights[self.triangles[first]]
                h2=self.barycentric(second,points)@self.heights[self.triangles[second]]
                scale=max(float(edge_scale[first]),float(edge_scale[second]),float(np.max(np.abs(relative))),1e-300)
                if np.max(np.abs(h1-h2))>512*np.finfo(float).eps*scale:raise ValueError('同体投影共享边/点对应不同FE高度，非连续单值图')
        self.projection_audit=dict(facet_count=len(self.polygons),queried_neighbor_pairs=pairs,
            tolerated_roundoff_overlap_area_m2=self.roundoff_overlap_area_m2,
            minimum_absolute_normal_cosine=float(np.min(np.abs(self.cosines))),
            original_FE_linear_facets_used=True,closest_point_clamping_used=False)

    def barycentric(self,index,xy):
        import numpy as np
        points=np.asarray(xy,dtype=float);ab=(points-self.projected_origins[index])@self.projected_inverse[index].T
        return np.column_stack((1-ab.sum(axis=1),ab)) if points.ndim==2 else np.r_[1-ab.sum(),ab]

    def mapping(self,index,xy):
        import numpy as np
        bary=self.barycentric(index,xy)
        if bary.min() < -self.barycentric_tolerance or bary.max()>1+self.barycentric_tolerance:raise ArithmeticError('凸片质心射线未落在其真实FE三角内；不夹取重心坐标')
        nodes=self.triangles[index];height=float(bary@self.heights[nodes])
        world=self.origin+np.asarray(xy)@self.basis+self.normal*height
        via_nodes=bary@self.world_nodes[nodes]
        return dict(nodes=self.node_ids[nodes],weights=bary,mapped_world=world,height_m=height,
            distance_m=float(np.linalg.norm(via_nodes-world)),face_edge_lengths_m=self.face_edge_lengths[index],
            facet_index=int(index),normal_world=self.normals[index],normal_cosine=float(self.cosines[index]))

    def ray(self,xy):
        from shapely.geometry import Point
        point=Point(xy)
        for index in self.tree.query(point,predicate='intersects'):
            index=int(index)
            if self.polygons[index].covers(point):return self.mapping(index,xy)
        raise ArithmeticError('完整cell角点无真实FE投影覆盖，不补虚假射线点')


def build_projected_fe_contact_system(first,second,plane,cells,*,coverage_relative_tolerance=2e-12,
        overlap_relative_tolerance=256*2.220446049250313e-16,minimum_normal_cosine=1e-6,
        barycentric_tolerance=5e-11,maximum_pieces=1000000,include_vertex_validation=True):
    """输出兼容分布接触的records/稀疏W；不满格只列partial，不进入QP。

    first/second：name、nodes_m、triangles、rotation_rows、axis_world；可给
    全局节点号/真实体节点总数、平移和轴原点。plane为origin_world_m、
    tangent_rows(2×3)、normal_world。每个cell给bounds_xy_m(xmin,ymin,
    xmax,ymax)及tooth_offset/row/column。f=p*完整cell面积，故片权始终
    为真实片面积/cell面积，绝不将缺失面积重新归一化。
    """
    import numpy as np
    from scipy import sparse
    from shapely.geometry import box
    from shapely.ops import unary_union
    from shapely.strtree import STRtree
    for value in (coverage_relative_tolerance,overlap_relative_tolerance,minimum_normal_cosine,barycentric_tolerance):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:raise ValueError('几何审核容差须为正有限显式数')
    if isinstance(maximum_pieces,bool) or not isinstance(maximum_pieces,int) or maximum_pieces<1:raise ValueError('分片上限须为正整数')
    if not isinstance(include_vertex_validation,bool):raise ValueError('分片顶点审核数据开关须为布尔值')
    origin=_vector(plane['origin_world_m'],3,'投影原点');normal=_vector(plane['normal_world'],3,'投影法向')
    basis=np.asarray(plane['tangent_rows'],dtype=float)
    if basis.shape!=(2,3) or not np.isfinite(basis).all() or abs(float(normal@normal)-1)>1e-11 or not np.allclose(basis@basis.T,np.eye(2),rtol=0,atol=1e-11) or not np.allclose(np.cross(basis[0],basis[1]),normal,rtol=0,atol=1e-11):raise ValueError('共同投影basis/n须构成正向正交单位坐标系')
    if first['name']==second['name']:raise ValueError('两体名称须不同')
    requests=list(cells);cell_polygons=[];cell_keys=set()
    for index,cell in enumerate(requests):
        key=tuple(cell.get(name,default) for name,default in (('tooth_offset',0),('row',0),('column',index)))
        if any(isinstance(v,bool) or not isinstance(v,Integral) for v in key) or key in cell_keys:raise ValueError('压力cell的物理齿/行/列索引须为唯一整数，不静默取整')
        cell_keys.add(key);bounds=_vector(cell['bounds_xy_m'],4,'cell矩形界')
        if bounds[2]<=bounds[0] or bounds[3]<=bounds[1]:raise ValueError('cell矩形面积须为正')
        cell_polygons.append(box(*bounds))
    cell_tree=STRtree(cell_polygons)
    for i,p in enumerate(cell_polygons):
        for j in cell_tree.query(p,predicate='intersects'):
            j=int(j)
            if j>i and p.intersection(cell_polygons[j]).area>overlap_relative_tolerance*min(p.area,cell_polygons[j].area):raise ValueError('压力cell存在正面积重叠，不能重复定义同一物理域的P0未知量')
    options=dict(minimum_normal_cosine=minimum_normal_cosine,overlap_relative_tolerance=overlap_relative_tolerance,barycentric_tolerance=barycentric_tolerance)
    bodies=[ProjectedFacetBody(data,(origin,basis,normal),**options) for data in (first,second)]
    rows={b.name:[] for b in bodies};cols={b.name:[] for b in bodies};values={b.name:[] for b in bodies}
    records=[];statuses=[];partial_records=[];total_pieces=0
    for cell_index,cell in enumerate(requests):
        bounds=_vector(cell['bounds_xy_m'],4,'cell矩形界')
        if bounds[2]<=bounds[0] or bounds[3]<=bounds[1]:raise ValueError('cell矩形面积须为正')
        polygon=box(*bounds);area=float(polygon.area);pieces=[];samples=[]
        gap=0.;levers={b.name:0. for b in bodies};actual={b.name:0. for b in bodies}
        for i in bodies[0].tree.query(polygon,predicate='intersects'):
            i=int(i);first_piece=polygon.intersection(bodies[0].polygons[i])
            if first_piece.area==0:continue
            for j in bodies[1].tree.query(first_piece,predicate='intersects'):
                j=int(j)
                for piece in _polygons(first_piece.intersection(bodies[1].polygons[j])):
                    total_pieces+=1
                    if total_pieces>maximum_pieces:raise ValueError('真实分片数超过显式上限')
                    if len(piece.interiors):raise ArithmeticError('三凸域交集产生孔洞，不能仅用一个质心积分')
                    xy=np.asarray(piece.centroid.coords[0]);weight=float(piece.area/area)
                    mappings={b.name:b.mapping(index,xy) for b,index in zip(bodies,(i,j))}
                    coefficients={};sample_levers={}
                    for b,sign in zip(bodies,(-1.,1.)):
                        mapped=mappings[b.name];world_coef=sign*mapped['normal_world']/mapped['normal_cosine']
                        coefficients[b.name]=b.rotation.T@world_coef
                        sample_levers[b.name]=-float(world_coef@np.cross(b.axis,mapped['mapped_world']-b.axis_origin))
                        levers[b.name]+=weight*sample_levers[b.name]
                        actual[b.name]+=float(piece.area)/abs(mapped['normal_cosine'])
                    local_gap=mappings[bodies[1].name]['height_m']-mappings[bodies[0].name]['height_m']
                    gap+=weight*local_gap;pieces.append(piece)
                    sample=dict(gap_m=local_gap,integration_weight=weight,target_xy_m=xy,
                        piece_area_plane_m2=float(piece.area),piece_vertices_xy_m=np.asarray(piece.exterior.coords[:-1]),
                        mappings=mappings,coefficients_local=coefficients,levers_m=sample_levers,
                        geometric_levers_m=sample_levers.copy(),integration_rule='真实FE投影凸片面积/质心，仿射场精确')
                    if include_vertex_validation:
                        vertex_weights={};vertex_heights={};vertex_levers={}
                        for b,index in zip(bodies,(i,j)):
                            shapes=b.barycentric(index,sample['piece_vertices_xy_m'])
                            if shapes.min() < -barycentric_tolerance or shapes.max()>1+barycentric_tolerance:raise ArithmeticError('真实分片顶点越过原FE三角，不夹取形函数')
                            local_nodes=b.triangles[index];world=shapes@b.world_nodes[local_nodes]
                            vertex_weights[b.name]=shapes;vertex_heights[b.name]=shapes@b.heights[local_nodes]
                            coef_world=b.rotation@coefficients[b.name]
                            vertex_levers[b.name]=-np.cross(b.axis,world-b.axis_origin)@coef_world
                        sample.update(piece_vertex_shape_weights=vertex_weights,
                            piece_vertex_gap_m=vertex_heights[bodies[1].name]-vertex_heights[bodies[0].name],
                            piece_vertex_levers_m=vertex_levers)
                    samples.append(sample)
        covered=unary_union(pieces) if pieces else polygon.difference(polygon)
        summed=sum(float(piece.area) for piece in pieces);covered_area=float(covered.area)
        overlap=max(0.,summed-covered_area);missing=polygon.difference(covered);excess=covered.difference(polygon)
        if overlap>coverage_relative_tolerance*area:raise ArithmeticError('cell分片存在超审核门限的重叠，不重复施加载荷')
        complete=bool(missing.area<=coverage_relative_tolerance*area and excess.area<=coverage_relative_tolerance*area and abs(summed-area)<=coverage_relative_tolerance*area)
        audit=dict(requested_cell_index=cell_index,tooth_offset=int(cell.get('tooth_offset',0)),row=int(cell.get('row',0)),column=int(cell.get('column',cell_index)),
            status='complete' if complete else 'partial',area_plane_m2=area,covered_area_plane_m2=covered_area,
            summed_piece_area_plane_m2=summed,missing_area_plane_m2=float(missing.area),overlap_area_plane_m2=overlap,
            missing_geometry_geojson=missing.__geo_interface__,piece_count=len(samples),integration_weight_sum=sum(q['integration_weight'] for q in samples),
            record_index=len(records) if complete else None)
        statuses.append(audit)
        record=dict(tooth_offset=audit['tooth_offset'],row=audit['row'],column=audit['column'],gap_m=gap,
            area_plane_m2=area,actual_area_m2=actual,levers_m=levers,geometric_levers_m=levers.copy(),quadrature=samples,
            maximum_mapping_error_m=max((q['mappings'][b.name]['distance_m'] for q in samples for b in bodies),default=0.),
            geometry_branch='真实线性FE面投影分片；非NURBS精确gap',coverage=audit)
        if not complete:
            record['interpretation']='不满格诊断积分；权重未归一化，不能作为已接受QP变量'
            partial_records.append(record);continue
        # 仅满格单元进入同协议的W/records；四角也从真实FE射线取得。
        corners=[(bounds[0],bounds[1]),(bounds[2],bounds[1]),(bounds[2],bounds[3]),(bounds[0],bounds[3])]
        record['first_corner_points_world_m']=[bodies[0].ray(xy)['mapped_world'].tolist() for xy in corners]
        column=len(records)
        for sample in samples:
            for b in bodies:
                mapping=sample['mappings'][b.name];coef=sample['coefficients_local'][b.name]
                for node,shape in zip(mapping['nodes'],mapping['weights']):
                    for direction in range(3):
                        rows[b.name].append(3*int(node)+direction);cols[b.name].append(column)
                        values[b.name].append(sample['integration_weight']*float(shape)*float(coef[direction]))
        records.append(record)
    matrices={b.name:sparse.coo_matrix((values[b.name],(rows[b.name],cols[b.name])),shape=(3*b.total_node_count,len(records))).tocsc() for b in bodies}
    return dict(records=records,matrices=matrices,cell_statuses=statuses,partial_records=partial_records,
        geometry_branch='projected_linear_FE_partition_v1',pressure_basis='P0投影矩形合力；真实FE三角交片精确积分',
        original_NURBS_gap_used=False,closest_point_clamping_used=False,partial_cell_renormalization_used=False,
        piece_vertex_validation_data_included=include_vertex_validation,
        projection_audits={b.name:b.projection_audit for b in bodies},requested_cell_count=len(statuses),
        complete_cell_count=len(records),partial_cell_count=len(partial_records),piece_count=total_pieces,
        tolerances=dict(coverage_relative=coverage_relative_tolerance,overlap_relative=overlap_relative_tolerance,
            minimum_normal_cosine=minimum_normal_cosine,barycentric=barycentric_tolerance),
        units=dict(gap='m',force='N',pressure='Pa',area='m^2',lever='m',angle='rad'),
        factorization_performed=False)


def evaluate_partition_remaining_gaps(system,displacements_local_m,theta_rad,*,follower_member='wheel'):
    """核查固定分片上的线性化剩余gap，极小值在凸片顶点，非仅看平均。

    不执行加载刚体几何重建，也不替代外层有限旋转验收；审核范围仅为
    已纳入QP的满格压力cell，不能当成整个有限齿面的完整零压证明。
    """
    import numpy as np
    if not system.get('piece_vertex_validation_data_included'):raise ValueError('本system未保留真实交片顶点数据，不能宣称处处仿射gap审核')
    if isinstance(theta_rad,bool) or not isinstance(theta_rad,(int,float)) or not math.isfinite(theta_rad):raise ValueError('线性化受载角须为有限rad')
    if set(displacements_local_m)!=set(system['matrices']) or follower_member not in system['matrices']:raise ValueError('位移体成员须与原W一致，且明确从动体')
    u={}
    for member,w in system['matrices'].items():
        values=np.asarray(displacements_local_m[member],dtype=float)
        if values.size!=w.shape[0] or not np.isfinite(values).all():raise ValueError('真实节点位移须与原W全局自由度布局一致')
        u[member]=values.reshape(-1,3)
    means=sum(np.asarray(w.T@u[m].ravel()).ravel() for m,w in system['matrices'].items())
    rows=[];minimum=math.inf
    for index,record in enumerate(system['records']):
        local_min=math.inf;location=None;count=0
        for sample in record['quadrature']:
            gap=np.array(sample['piece_vertex_gap_m'],copy=True)
            for member in u:
                mapping=sample['mappings'][member]
                gap+=sample['piece_vertex_shape_weights'][member]@u[member][mapping['nodes']]@sample['coefficients_local'][member]
            gap-=theta_rad*sample['piece_vertex_levers_m'][follower_member]
            j=int(np.argmin(gap));count+=len(gap)
            if gap[j]<local_min:local_min=float(gap[j]);location=sample['piece_vertices_xy_m'][j].tolist()
        mean=record['gap_m']+float(means[index])-theta_rad*record['levers_m'][follower_member]
        rows.append(dict(record_index=index,row=record['row'],column=record['column'],tooth_offset=record['tooth_offset'],
            mean_remaining_gap_m=mean,minimum_remaining_gap_m=local_min,minimum_location_xy_m=location,checked_piece_vertices=count))
        minimum=min(minimum,local_min)
    return dict(minimum_remaining_gap_m=None if not rows else minimum,cell_checks=rows,
        fixed_piecewise_affine_extrema_evaluated=True,loaded_finite_rotation_rebuilt=False,
        whole_flank_nonpenetration_certified=False,checked_domain='仅已纳入QP的完整压力cell')
