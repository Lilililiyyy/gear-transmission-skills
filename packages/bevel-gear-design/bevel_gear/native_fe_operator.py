"""参数化原生FE接触operator：同rest/格/平面在新phi重建及typed归档。

只构造几何与W，不准备K或求压力。legacy record/q是原typed数组的惰性
视图，供冻结原QP及独立g0/W审核读取；不复制旧C，不调用NURBS gap。
"""
from __future__ import annotations

from collections.abc import Mapping,Sequence
from copy import deepcopy
from dataclasses import asdict,is_dataclass
from pathlib import Path,PureWindowsPath
from numbers import Integral
import gzip
import hashlib
import json
import math
import time

MEMBERS=('pinion','wheel')
METHOD='P0_cell_intersection_with_both_actual_FE_triangles_centroid_exact'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def safe(value):
    import numpy as np
    if is_dataclass(value):return safe(asdict(value))
    if isinstance(value,np.ndarray):return safe(value.tolist())
    if isinstance(value,np.generic):return safe(value.item())
    if isinstance(value,Mapping):return {str(k):safe(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [safe(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):return None
    return value


def fingerprint(value):return hashlib.sha256(json.dumps(safe(value),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf8')).hexdigest()


def _input_identity(value):
    """对实际数值输入取身份；大数组按dtype/shape/原字节SHA绑定。"""
    import numpy as np
    if is_dataclass(value):return _input_identity(asdict(value))
    if isinstance(value,np.ndarray):
        a=np.ascontiguousarray(value)
        if a.dtype.kind=='O':raise ValueError('operator输入数组不能使用对象dtype')
        return dict(dtype=a.dtype.str,shape=list(a.shape),sha256=hashlib.sha256(a.tobytes()).hexdigest())
    if isinstance(value,Mapping):return {str(k):_input_identity(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [_input_identity(v) for v in value]
    return safe(value)


def _archive_path(directory,name):
    """只接受归档根内相对普通文件；拒绝绝对路径、父级和越界软链接。"""
    if not isinstance(name,str) or not name or '\x00' in name:raise ValueError('归档文件名须为非空相对路径')
    normalized=name.replace('\\','/');part=Path(normalized)
    if part.is_absolute() or PureWindowsPath(name).drive or '..' in part.parts:raise ValueError('归档文件路径不得为绝对路径或引用父级')
    candidate=(directory/part).resolve()
    if directory not in candidate.parents or not candidate.is_file():raise ValueError('归档文件缺失或路径越过归档根目录')
    return candidate


def compact_fe_boundary(mesh,name,*,boundary_triangles=None):
    """从全局FE网格选择真实边界三角；返回局部facet及原全局节点ID。

    mesh可为含nodes_m/boundary_triangles的对象或Mapping。可显式给单侧
    齿面三角；整封闭边界通常不是共同平面上的单值图，会被后续拒绝。
    """
    import numpy as np
    get=(lambda key:mesh[key]) if isinstance(mesh,Mapping) else (lambda key:getattr(mesh,key))
    nodes=np.asarray(get('nodes_m'),dtype=float);tri=np.asarray(get('boundary_triangles') if boundary_triangles is None else boundary_triangles)
    if nodes.ndim!=2 or nodes.shape[1:]!=(3,) or not len(nodes) or not np.isfinite(nodes).all():raise ValueError('全局FE节点须为非空有限N×3米坐标')
    if tri.ndim!=2 or tri.shape[1:]!=(3,) or not len(tri) or tri.dtype.kind not in 'iu' or tri.min()<0 or tri.max()>=len(nodes):raise ValueError('真实边界三角须为全局节点范围内非空整数M×3')
    if not isinstance(name,str) or not name:raise ValueError('真实FE体须有非空名称')
    ids,inverse=np.unique(tri,return_inverse=True)
    return dict(name=name,nodes_m=nodes[ids].copy(),triangles=inverse.reshape(-1,3),global_node_ids=ids,total_node_count=len(nodes))


def _body_at_frames(data,frames,member):
    result=dict(data,rotation_rows=frames[member+'_rotation_rows'],axis_world=frames[member+'_axis_world'])
    for suffix in ('translation_world_m','axis_origin_world_m'):
        if member+'_'+suffix in frames:result[suffix]=frames[member+'_'+suffix]
    return result


def build_native_fe_contact_operator(compact,patches,frames,reference_follower_angle_rad,*,
        driving_angle_rad=0.,pinion_torque_nm=None,pinion_to_wheel_ratio=None,batch_cells=128,source_identity=None):
    """新建显式FE_boundary operator，无cache/私有文件/体刚度依赖。

    compact为pinion/wheel各自的{k:真实facet字典}，也可各给单个k=0
    facet字典。patch含tooth_offset、grid、plane及requested_mask；未给
    mask时请求整个给定格。源请求partial仍保留，不能将其当零压力。
    """
    import numpy as np
    from .halfspace_contact import RectangularContactGrid
    if set(compact)!=set(MEMBERS):raise ValueError('原生齿轮接触须提供pinion/wheel两体实际facet')
    compact={m:({0:compact[m]} if 'nodes_m' in compact[m] else compact[m]) for m in MEMBERS}
    prepared=[];seen=set()
    for original in patches:
        p=deepcopy(original);k=p.get('tooth_offset',0)
        if isinstance(k,bool) or not isinstance(k,Integral) or k in seen:raise ValueError('patch物理齿偏移须为唯一整数')
        seen.add(k);p['tooth_offset']=int(k);grid=p['grid'];grid=RectangularContactGrid(**grid) if isinstance(grid,dict) else grid;p['grid']=grid
        mask=np.asarray(p.pop('requested_mask',p.get('source_valid_mask',p.get('valid_mask',np.ones(grid.shape,dtype=bool)))))
        if mask.shape!=grid.shape or mask.dtype.kind!='b':raise ValueError('请求mask须为与grid同形状的布尔数组')
        p['source_valid_mask']=mask.copy();p['valid_mask']=mask.copy()
        if 'plane' in p:
            plane=p.pop('plane');p.update(plane_origin_world_m=plane['origin_world_m'],plane_tangent_rows=plane['tangent_rows'],plane_outward_normal_world=plane['normal_world'])
        p.setdefault('source_material_reference_frames',deepcopy(frames));p['frames']=deepcopy(frames);prepared.append(p)
    if not prepared:raise ValueError('须提供至少一个实际接触patch')
    geometry=dict(patches=prepared,driving_angle_rad=driving_angle_rad,reference_follower_angle_rad=reference_follower_angle_rad)
    for key,value in (('driving_angle_rad',driving_angle_rad),('pinion_torque_nm',pinion_torque_nm),('pinion_to_wheel_ratio',pinion_to_wheel_ratio)):
        if value is not None:
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError(key+'须为有限SI数值')
            geometry[key]=float(value)
    base=dict(geometry=geometry,operator_source_identity=source_identity or {})
    return build_native_operator_at_pose(base,compact,frames,reference_follower_angle_rad,batch_cells=batch_cells)


def _write_json(path,value):
    path=Path(path);raw=(json.dumps(safe(value),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode('utf8')
    if path.suffix=='.gz':
        with path.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as packed:packed.write(raw)
    else:path.write_bytes(raw)
    return dict(file=path.name,sha256=sha(path))


def _read_json(path):
    path=Path(path)
    if path.suffix=='.gz':
        with gzip.open(path,'rt',encoding='utf8') as stream:return json.load(stream)
    return json.loads(path.read_text(encoding='utf8'))


def _restore_geometry(data):
    import numpy as np
    from .halfspace_contact import RectangularContactGrid
    result=dict(data);result['patches']=[]
    bool_fields={'valid_mask','source_valid_mask','native_partial_cell_mask','native_center_ray_valid_mask','corner_valid_mask'}
    array_fields={'gaps_m','wheel_lever_arms_m','pinion_lever_arms_m','wheel_normal_moment_arms_m','pinion_normal_moment_arms_m',
        'surface_slopes','opposing_normal_angles_rad','first_surface_area_jacobian','second_surface_area_jacobian',
        'first_points_world_m','second_points_world_m','first_outward_normals_world','second_outward_normals_world',
        'first_corner_points_world_m','second_corner_points_world_m'}
    for original in data['patches']:
        patch=dict(original);patch['grid']=RectangularContactGrid(**patch['grid']) if isinstance(patch['grid'],dict) else patch['grid']
        for field in bool_fields:
            if field in patch:patch[field]=np.asarray(patch[field],dtype=bool)
        for field in array_fields:
            if field in patch:patch[field]=np.asarray(patch[field],dtype=float)
        result['patches'].append(patch)
    return result


class TypedQPoint(Mapping):
    """单个原交片的懒视图；只有独立审核请求时才生成顶点形函数。"""
    KEYS=('gap_m','integration_weight','target_xy_m','piece_area_plane_m2','piece_vertices_xy_m','mappings',
          'coefficients_local','levers_m','geometric_levers_m','integration_rule','piece_vertex_shape_weights','piece_vertex_gap_m','piece_vertex_levers_m')
    def __init__(self,arrays,index):self.a=arrays;self.i=index;self.cache={}
    def __iter__(self):return iter(self.KEYS)
    def __len__(self):return len(self.KEYS)
    def __getitem__(self,key):
        import numpy as np
        a=self.a;i=self.i
        if key in self.cache:return self.cache[key]
        if key=='gap_m':return float(a['piece_gap_m'][i])
        if key=='integration_weight':return float(a['piece_integration_weights'][i])
        if key=='target_xy_m':return a['piece_centroids_xy_m'][i]
        if key=='piece_area_plane_m2':return float(a['piece_areas_m2'][i])
        if key=='piece_vertices_xy_m':return a['piece_vertices_xy_m'][a['piece_vertex_offsets'][i]:a['piece_vertex_offsets'][i+1]]
        if key=='integration_rule':return '真实FE投影凸片面积/质心，仿射场精确；typed数组视图'
        if key in ('coefficients_local','levers_m','geometric_levers_m'):
            name='levers_m' if key=='geometric_levers_m' else key
            value={m:a['body_data'][m][name][i] for m in MEMBERS};self.cache[key]=value;return value
        if key=='mappings':
            value={}
            for m in MEMBERS:
                d=a['body_data'][m];t=a['body_projection_tables'][m];face=int(d['facet_indices'][i])
                value[m]=dict(nodes=d['global_node_ids'][i],weights=d['shape_weights'][i],mapped_world=d['mapped_world_m'][i],
                    height_m=float(d['height_m'][i]),distance_m=float(d['closure_error_m'][i]),facet_index=face,
                    face_edge_lengths_m=t['face_edge_lengths_m'][face],normal_world=t['normals_world'][face],normal_cosine=float(t['normal_cosines'][face]))
            self.cache[key]=value;return value
        if key.startswith('piece_vertex_'):
            xy=self['piece_vertices_xy_m'];shapes={};heights={};levers={}
            for m in MEMBERS:
                t=a['body_projection_tables'][m];d=a['body_data'][m];face=int(d['facet_indices'][i])
                ab=(xy-t['projected_origins'][face])@t['projected_inverse'][face].T;w=np.column_stack((1-ab.sum(axis=1),ab))
                if w.min() < -a['barycentric_tolerance'] or w.max()>1+a['barycentric_tolerance']:raise ArithmeticError('typed原交片顶点越过真实FE面')
                local=t['local_triangles'][face];shapes[m]=w;heights[m]=w@t['local_node_heights_m'][local]
                world=t['plane_origin_world_m']+xy@t['plane_basis_rows']+heights[m][:,None]*t['plane_normal_world']
                levers[m]=-np.cross(t['axis_world'],world-t['axis_origin_world_m'])@d['coefficients_world'][i]
            self.cache.update(piece_vertex_shape_weights=shapes,piece_vertex_gap_m=heights['wheel']-heights['pinion'],piece_vertex_levers_m=levers)
            if key in self.cache:return self.cache[key]
        raise KeyError(key)


class TypedQSequence(Sequence):
    def __init__(self,arrays,start,end):self.a=arrays;self.start=int(start);self.end=int(end)
    def __len__(self):return self.end-self.start
    def __getitem__(self,index):
        if isinstance(index,slice):return [self[i] for i in range(*index.indices(len(self)))]
        if index<0:index+=len(self)
        if not 0<=index<len(self):raise IndexError(index)
        return TypedQPoint(self.a,self.start+index)
    def __iter__(self):return (TypedQPoint(self.a,i) for i in range(self.start,self.end))


class TypedRecord(Mapping):
    def __init__(self,header,arrays,offsets):self.header=header;self.a=arrays;self.offsets=offsets
    def __iter__(self):yield from self.header;yield 'quadrature'
    def __len__(self):return len(self.header)+1
    def __getitem__(self,key):
        if key=='quadrature':
            ci=self.header['requested_cell_index'];return TypedQSequence(self.a,self.offsets[ci],self.offsets[ci+1])
        if key=='maximum_fe_ray_shape_closure_error_m':return self.header['maximum_mapping_error_m']
        return self.header[key]


class TypedRecords(Sequence):
    def __init__(self,patch_arrays):
        import numpy as np
        self.items=[]
        for a in patch_arrays:
            if a is None:continue
            offsets=np.r_[0,np.cumsum(np.bincount(a['piece_cell_indices'],minlength=a['requested_cell_count']))]
            self.items.extend(TypedRecord(h,a,offsets) for h in a['records'])
    def __len__(self):return len(self.items)
    def __iter__(self):return iter(self.items)
    def __getitem__(self,index):return self.items[index]


def _add_projection_tables(arrays,bodies):
    for b in bodies:
        arrays['body_projection_tables'][b.name].update(world_nodes_m=b.world_nodes,xy_nodes_m=b.xy_nodes,
            face_edge_lengths_m=b.face_edge_lengths,normals_world=b.normals,normal_cosines=b.cosines)


def _ray_many(body,xy):
    import numpy as np
    import shapely
    from .projected_fe_partition_arrays import _point_data
    pairs=body.tree.query(shapely.points(xy),predicate='intersects');hit,first=np.unique(pairs[0],return_index=True)
    if not np.array_equal(hit,np.arange(len(xy))):raise ArithmeticError('完整cell中心/角点缺真实FE射线，不能外推')
    faces=pairs[1,first];out=_point_data(body,faces,xy,barycentric_tolerance=5e-11)
    out.update(normal_world=body.normals[faces],normal_cosine=body.cosines[faces]);return out


def _merge_batches(batches,bodies):
    import numpy as np
    from scipy import sparse
    first=batches[0];out={k:v for k,v in first.items() if k not in ('records','partial_records','matrices','cell_statuses') and not isinstance(v,np.ndarray)}
    arrays=[k for k,v in first.items() if isinstance(v,np.ndarray)];records=[];partial=[];statuses=[];cell_offset=0;column_offset=0;vertex_offset=0
    for b in batches:
        for h in b['records']+b['partial_records']:
            h['requested_cell_index']+=cell_offset;h['coverage']['requested_cell_index']+=cell_offset
            if h['coverage']['record_index'] is not None:h['coverage']['record_index']+=column_offset
        records.extend(b['records']);partial.extend(b['partial_records']);statuses.extend(b['cell_statuses'])
        cell_offset+=b['requested_cell_count'];column_offset+=b['complete_cell_count']
    for key in arrays:
        if key in ('piece_cell_indices','cell_to_column','piece_vertex_offsets'):continue
        out[key]=np.concatenate([b[key] for b in batches],axis=0)
    offsets=[];cell_offset=0;column_offset=0;vertex_offset=0;ids=[];columns=[]
    for b in batches:
        ids.append(b['piece_cell_indices']+cell_offset);c=b['cell_to_column'].copy();c[c>=0]+=column_offset;columns.append(c)
        offsets.append(b['piece_vertex_offsets'][:-1]+vertex_offset)
        cell_offset+=b['requested_cell_count'];column_offset+=b['complete_cell_count'];vertex_offset+=len(b['piece_vertices_xy_m'])
    out['piece_cell_indices']=np.concatenate(ids);out['cell_to_column']=np.concatenate(columns)
    out['piece_vertex_offsets']=np.r_[np.concatenate(offsets),vertex_offset]
    out['body_data']={m:{k:np.concatenate([b['body_data'][m][k] for b in batches]) for k in first['body_data'][m]} for m in MEMBERS}
    out['piece_facet_indices']={m:out['body_data'][m]['facet_indices'] for m in MEMBERS}
    out['body_projection_tables']=first['body_projection_tables'];_add_projection_tables(out,bodies)
    out.update(records=records,partial_records=partial,cell_statuses=statuses,
        matrices={m:sparse.hstack([b['matrices'][m] for b in batches],format='csc') for m in MEMBERS},
        requested_cell_count=cell_offset,complete_cell_count=column_offset,partial_cell_count=cell_offset-column_offset,
        piece_count=sum(b['piece_count'] for b in batches))
    return out


def _patch_geometry(source,frames,arrays,bodies):
    import numpy as np
    p=deepcopy(source);p['frames']=deepcopy(frames);grid=p['grid'];shape=grid.shape
    p['source_valid_mask']=np.asarray(source.get('source_valid_mask',source['valid_mask']),dtype=bool).copy()
    for name in ('valid_mask','native_center_ray_valid_mask'):p[name]=np.zeros(shape,dtype=bool)
    p['native_partial_cell_mask']=p['source_valid_mask'].copy();p['corner_valid_mask']=np.zeros((shape[0]+1,shape[1]+1),dtype=bool)
    scalars=('gaps_m','wheel_lever_arms_m','pinion_lever_arms_m','wheel_normal_moment_arms_m','pinion_normal_moment_arms_m',
        'surface_slopes','opposing_normal_angles_rad','first_surface_area_jacobian','second_surface_area_jacobian')
    for name in scalars:p[name]=np.full(shape,np.nan)
    for prefix in ('first','second'):
        p[prefix+'_points_world_m']=np.full(shape+(3,),np.nan);p[prefix+'_outward_normals_world']=np.full(shape+(3,),np.nan)
        p[prefix+'_corner_points_world_m']=np.full((shape[0]+1,shape[1]+1,3),np.nan)
    p.update(geometry_source='FE_boundary',integration_method=METHOD,native_loaded_geometry_verified=False,
        native_gap_definition='当前phi下真实双FE交片平均；中心/四角为同一实际FE射线',native_operator_rebuilt_at_current_pose=True)
    if arrays is None:return p
    rows=np.asarray([h['row'] for h in arrays['records']],dtype=int);cols=np.asarray([h['column'] for h in arrays['records']],dtype=int)
    if not len(rows):return p
    xy=np.column_stack((grid.origin_xy_m[0]+cols*grid.spacing_x_m,grid.origin_xy_m[1]+rows*grid.spacing_y_m))
    sampled={b.name:_ray_many(b,xy) for b in bodies};p['valid_mask'][rows,cols]=True;p['native_center_ray_valid_mask'][rows,cols]=True;p['native_partial_cell_mask'][rows,cols]=False
    p['gaps_m'][rows,cols]=[h['gap_m'] for h in arrays['records']]
    normal=np.asarray(p['plane_outward_normal_world']);basis=np.asarray(p['plane_tangent_rows'])
    for body,prefix,sign in zip(bodies,('first','second'),(1.,-1.)):
        m=body.name;d=sampled[m];p[prefix+'_points_world_m'][rows,cols]=d['mapped_world_m'];p[prefix+'_outward_normals_world'][rows,cols]=d['normal_world']
        p[prefix+'_surface_area_jacobian'][rows,cols]=[h['actual_area_m2'][m]/h['area_plane_m2'] for h in arrays['records']]
        p[m+'_lever_arms_m'][rows,cols]=[h['levers_m'][m] for h in arrays['records']]
        moment=sign*np.cross(body.axis,arrays['body_data'][m]['mapped_world_m']-body.axis_origin)@normal
        average=np.bincount(arrays['piece_cell_indices'],weights=arrays['piece_integration_weights']*moment,minlength=arrays['requested_cell_count'])
        p[m+'_normal_moment_arms_m'][rows,cols]=[average[h['requested_cell_index']] for h in arrays['records']]
    p['surface_slopes'][rows,cols]=np.maximum.reduce([np.linalg.norm(d['normal_world']@basis.T,axis=1)/np.abs(d['normal_cosine']) for d in sampled.values()])
    p['opposing_normal_angles_rad'][rows,cols]=np.arccos(np.clip(-np.einsum('ij,ij->i',sampled['pinion']['normal_world'],sampled['wheel']['normal_world']),-1,1))
    for j,i in zip(rows,cols):p['corner_valid_mask'][j,i]=p['corner_valid_mask'][j,i+1]=p['corner_valid_mask'][j+1,i]=p['corner_valid_mask'][j+1,i+1]=True
    jj,ii=np.nonzero(p['corner_valid_mask']);qxy=np.column_stack((grid.origin_xy_m[0]+(ii-.5)*grid.spacing_x_m,grid.origin_xy_m[1]+(jj-.5)*grid.spacing_y_m))
    for body,prefix in zip(bodies,('first','second')):p[prefix+'_corner_points_world_m'][jj,ii]=_ray_many(body,qxy)['mapped_world_m']
    return p


def build_native_operator_at_pose(base_system,compact,frames,reference_follower_angle_rad,*,batch_cells=128,progress=None):
    """同源请求集合/格/平面，按新整体R0重建g/W/ell；绝不取旧C。"""
    import numpy as np
    from scipy import sparse
    from .projected_fe_body_arrays import ArrayProjectedFacetBody
    from .projected_fe_partition_arrays import ArrayProjectedFEIntegrator
    from .projected_fe_partition import _validated_plane
    if isinstance(reference_follower_angle_rad,bool) or not isinstance(reference_follower_angle_rad,(int,float)) or not math.isfinite(reference_follower_angle_rad) or isinstance(batch_cells,bool) or not isinstance(batch_cells,int) or batch_cells<1:raise ValueError('新参考phi须为有限rad、批量须为正整数')
    if set(compact)!=set(MEMBERS):raise ValueError('原生operator须提供pinion/wheel实际facet')
    if not np.allclose(frames.get('translation_m',[0,0,0]),0,rtol=0,atol=1e-16):raise ValueError('请用各体translation_world_m明确平移，不能使用歧义的整体translation_m')
    geometry={k:deepcopy(v) for k,v in base_system['geometry'].items() if k!='patches'};geometry['patches']=[]
    geometry.update(reference_follower_angle_rad=float(reference_follower_angle_rad),geometry_source='FE_boundary',integration_method=METHOD,native_operator_rebuilt_at_current_pose=True)
    frames=deepcopy(frames);frames['native_reference_follower_angle_rad']=float(reference_follower_angle_rad)
    patch_arrays=[];matrices={m:[] for m in MEMBERS};started=time.perf_counter()
    options=dict(minimum_normal_cosine=1e-6,overlap_relative_tolerance=256*np.finfo(float).eps,barycentric_tolerance=5e-11)
    seen=set()
    for source in base_system['geometry']['patches']:
        k=source['tooth_offset'];grid=source['grid'];mask=np.asarray(source['source_valid_mask'] if 'source_valid_mask' in source else source['valid_mask']);cells=[]
        if isinstance(k,bool) or not isinstance(k,Integral) or k in seen:raise ValueError('patch物理齿偏移须为唯一整数')
        seen.add(k)
        if mask.shape!=grid.shape or mask.dtype.kind!='b':raise ValueError('源请求mask须为与grid同布局的布尔数组')
        if any(k not in compact[m] or compact[m][k].get('name')!=m for m in MEMBERS):raise ValueError('实际facet物理齿偏移或体名称与patch不符')
        for j,i in zip(*np.nonzero(mask)):
            x=grid.origin_xy_m[0]+i*grid.spacing_x_m;y=grid.origin_xy_m[1]+j*grid.spacing_y_m
            cells.append(dict(tooth_offset=k,row=int(j),column=int(i),bounds_xy_m=[x-grid.spacing_x_m/2,y-grid.spacing_y_m/2,x+grid.spacing_x_m/2,y+grid.spacing_y_m/2]))
        if not cells:
            patch_arrays.append(None);geometry['patches'].append(_patch_geometry(source,frames,None,[]))
            for m in MEMBERS:matrices[m].append(sparse.csc_matrix((3*compact[m][k]['total_node_count'],0)))
            continue
        plane=_validated_plane((source['plane_origin_world_m'],source['plane_tangent_rows'],source['plane_outward_normal_world']))
        bodies=[ArrayProjectedFacetBody(_body_at_frames(compact[m][k],frames,m),plane,**options) for m in MEMBERS]
        integrator=ArrayProjectedFEIntegrator.from_verified_bodies(bodies);batches=[]
        for start in range(0,len(cells),batch_cells):batches.append(integrator.build_cells(cells[start:start+batch_cells]))
        arrays=_merge_batches(batches,bodies);patch_arrays.append(arrays)
        geometry['patches'].append(_patch_geometry(source,frames,arrays,bodies))
        for m in MEMBERS:matrices[m].append(arrays['matrices'][m])
        if progress:progress(dict(stage='native_outer_operator_patch_built',tooth_offset=k,requested_cells=len(cells),complete_cells=arrays['complete_cell_count'],partial_cells=arrays['partial_cell_count'],pieces=arrays['piece_count']))
    if not seen:raise ValueError('至少需要一个实际接触patch')
    requests=[{key:p[key] for key in ('tooth_offset','grid','source_valid_mask','plane_origin_world_m','plane_tangent_rows','plane_outward_normal_world')} for p in geometry['patches']]
    input_fp=fingerprint(_input_identity(dict(compact=compact,requested_patches=requests,frames=frames,
        driving_angle_rad=geometry['driving_angle_rad'],reference_follower_angle_rad=reference_follower_angle_rad)))
    geometry['operator_input_fingerprint']=input_fp;records=TypedRecords(patch_arrays)
    result=dict(records=records,matrices={m:sparse.hstack(matrices[m],format='csc') for m in MEMBERS},geometry=geometry,frames=frames,
        typed_patch_arrays=patch_arrays,geometry_source='FE_boundary',integration_method=METHOD,
        operator_input_fingerprint=input_fp,operator_source_identity=base_system.get('operator_source_identity',{}),
        maximum_cad_fe_mapping_error_m=max((r['maximum_mapping_error_m'] for r in records),default=0.),
        pressure_basis='真实P1 FE面双三角交片精确P0；当前phi下独立新W，不复用旧C',
        ray_backend='actual_FE_triangle_intersection_typed_arrays',quadrature_order='exact_affine_partition_centroids',
        piece_vertex_validation_data_included=True,operator_generation_seconds=time.perf_counter()-started,
        original_NURBS_gap_used=False,halfspace_compliance_used=False,old_compliance_reused=False)
    if 'native_cache_index_sha256' in base_system:result['native_cache_index_sha256']=base_system['native_cache_index_sha256']
    return result


def advance_native_follower_frame(frames,delta_rad):
    """从动轮整体姿态绕其真实世界轴有限转动，增量单位rad。

    轴原点固定；若给各体平移，按同一刚体旋转更新从动轮平移。pinion
    姿态及来源参考齿Rk不变。此操作不把增量命名为总弹性TE。
    """
    import numpy as np
    if isinstance(delta_rad,bool) or not isinstance(delta_rad,(int,float)) or not math.isfinite(delta_rad):raise ValueError('从动轮有限转角须为有限rad')
    out=deepcopy(frames);axis=np.asarray(frames['wheel_axis_world'],dtype=float)
    if axis.shape!=(3,) or not np.isfinite(axis).all() or abs(float(axis@axis)-1)>1e-11:raise ValueError('从动轮旋转轴须为真实三维单位向量')
    x,y,z=axis;cross=np.asarray([[0,-z,y],[z,0,-x],[-y,x,0]])
    rotation=np.eye(3)+math.sin(delta_rad)*cross+(1-math.cos(delta_rad))*(cross@cross)
    out['wheel_rotation_rows']=rotation@np.asarray(frames['wheel_rotation_rows'],dtype=float)
    origin=np.asarray(frames.get('wheel_axis_origin_world_m',[0,0,0]),dtype=float)
    if 'wheel_translation_world_m' in frames or np.any(origin):
        translation=np.asarray(frames.get('wheel_translation_world_m',[0,0,0]),dtype=float)
        if origin.shape!=(3,) or translation.shape!=(3,) or not np.isfinite(origin).all() or not np.isfinite(translation).all():raise ValueError('轴原点及刚体平移须为有限三维米坐标')
        out['wheel_translation_world_m']=origin+rotation@(translation-origin)
    if 'native_reference_follower_angle_rad' in frames:out['native_reference_follower_angle_rad']=frames['native_reference_follower_angle_rad']+delta_rad
    return out


def check_native_operator_pose(system,compact,displacements_local_m,theta_rad,reference_remaining_gaps_m,*,
        batch_cells=128,mean_tolerance_m=1e-10,maximum_penetration_m=2e-7,progress=None):
    """原typed片全部顶点与实际有限phi重建审核；不创建K或重求压力。

    u为原全局节点局部小变形P1迹，长度m。先核原g+Wᵀu-theta*ell，
    再在phi+theta重建真实双FE交片并核g_loaded+W_loadedᵀu。检查覆盖
    全部原请求；partial不会因有效mask缩小而消失。不声称大滑动分析。
    """
    import numpy as np
    from .projected_fe_partition_arrays import evaluate_array_piece_vertices
    if isinstance(batch_cells,bool) or not isinstance(batch_cells,int) or batch_cells<1:raise ValueError('审核批量须为正整数')
    for v in (mean_tolerance_m,maximum_penetration_m):
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<=0:raise ValueError('平均误差及穿入门限须为正有限米值')
    if 'typed_patch_arrays' not in system:raise ValueError('原生pose审核需要原typed交片数据，不能以中心采样代替')
    remaining=np.asarray(reference_remaining_gaps_m,dtype=float)
    if remaining.shape!=(len(system['records']),) or not np.isfinite(remaining).all():raise ValueError('原QP剩余gap须对应全部完整原cell')
    if set(displacements_local_m)!=set(MEMBERS):raise ValueError('位移须提供pinion/wheel原全局节点数组')
    u={}
    for m,w in system['matrices'].items():
        value=np.asarray(displacements_local_m[m],dtype=float)
        if value.size!=w.shape[0] or not np.isfinite(value).all():raise ValueError('位移须与真实全局FE节点布局一致')
        u[m]=value.reshape(-1,3)
    loaded_frames=advance_native_follower_frame(system['frames'],theta_rad);started=time.perf_counter()
    rows=[];oldmin=math.inf;newmin=math.inf;old_error=0.;new_error=0.;q_W_error=0.;old_vertices=0;new_vertices=0;partial=0;cursor=0;patches=[]
    try:
        initial_partial=sum(a['partial_cell_count'] for a in system['typed_patch_arrays'] if a is not None)
        for a in system['typed_patch_arrays']:
            if a is None:continue
            checked=evaluate_array_piece_vertices(a,u,theta_rad);qopening=np.zeros(a['piece_count'])
            for m in MEMBERS:
                d=a['body_data'][m]
                qopening+=np.einsum('pi,pij,pj->p',d['shape_weights'],u[m][d['global_node_ids']],d['coefficients_local'])
            qmeans=np.bincount(a['piece_cell_indices'],weights=a['piece_integration_weights']*qopening,minlength=a['requested_cell_count'])
            for record,item in zip(a['records'],checked['cell_checks']):
                original=system['records'][cursor];key=(record['tooth_offset'],record['row'],record['column'])
                if key!=(original['tooth_offset'],original['row'],original['column']):raise ValueError('typed原片与W物理cell列顺序不同')
                mean=record['gap_m']+qmeans[record['requested_cell_index']]-theta_rad*record['levers_m']['wheel']
                q_W_error=max(q_W_error,abs(mean-item['mean_remaining_gap_m']));old_error=max(old_error,abs(mean-remaining[cursor]))
                oldmin=min(oldmin,item['minimum_remaining_gap_m']);old_vertices+=item['checked_piece_vertices']
                rows.append(dict(tooth_offset=key[0],row=key[1],column=key[2],original_record_index=cursor,
                    old_mean_remaining_m=float(mean),old_minimum_piece_vertex_gap_m=item['minimum_remaining_gap_m'],reference_remaining_gap_m=float(remaining[cursor])))
                cursor+=1
        if cursor!=len(system['records']):raise ValueError('未核全部原完整候选cell')
        compact={m:({0:compact[m]} if 'nodes_m' in compact[m] else compact[m]) for m in MEMBERS}
        loaded=build_native_operator_at_pose(system,compact,loaded_frames,
            system['geometry']['reference_follower_angle_rad']+theta_rad,batch_cells=batch_cells,progress=progress)
        lookup={(r['tooth_offset'],r['row'],r['column']):r for r in rows}
        for patch,a in zip(loaded['geometry']['patches'],loaded['typed_patch_arrays']):
            if a is None:continue
            partial+=a['partial_cell_count'];checked=evaluate_array_piece_vertices(a,u,0.)
            for record in a['partial_records']:
                key=(record['tooth_offset'],record['row'],record['column'])
                if key in lookup:lookup[key].update(loaded_complete=False,loaded_missing_area_m2=record['coverage']['missing_area_plane_m2'])
            for record,item in zip(a['records'],checked['cell_checks']):
                key=(record['tooth_offset'],record['row'],record['column'])
                if key not in lookup:continue
                row=lookup[key];error=abs(item['mean_remaining_gap_m']-row['reference_remaining_gap_m']);new_error=max(new_error,error)
                newmin=min(newmin,item['minimum_remaining_gap_m']);new_vertices+=item['checked_piece_vertices']
                row.update(loaded_complete=True,loaded_mean_gap_m=item['mean_remaining_gap_m'],loaded_vs_old_predicted_mean_error_m=error,
                    loaded_minimum_piece_vertex_gap_m=item['minimum_remaining_gap_m'])
            patches.append(dict(tooth_offset=patch['tooth_offset'],original_native_cells=int(np.asarray(patch['source_valid_mask']).sum()),
                loaded_complete_cells=a['complete_cell_count'],loaded_partial_cells=a['partial_cell_count']))
        old_ok=bool(initial_partial==0 and math.isfinite(oldmin) and oldmin>=-maximum_penetration_m and old_error<=mean_tolerance_m and q_W_error<=mean_tolerance_m)
        new_ok=bool(partial==0 and math.isfinite(newmin) and newmin>=-maximum_penetration_m and new_error<=mean_tolerance_m and all(r.get('loaded_complete') for r in rows))
        return dict(completed=True,passed=old_ok and new_ok,initial_piece_vertex_check=dict(passed=old_ok,minimum_remaining_gap_m=None if not math.isfinite(oldmin) else oldmin,
            maximum_mean_vs_QP_error_m=old_error,checked_vertices=old_vertices,maximum_original_q_mean_vs_W_mean_error_m=q_W_error,partial_cells=initial_partial),
            loaded_FE_phi_geometry_check=dict(passed=new_ok,minimum_remaining_gap_m=None if not math.isfinite(newmin) else newmin,
                maximum_cell_mean_error_m=new_error,partial_cells=partial,checked_vertices=new_vertices,finite_phi_geometry_rebuilt=True,deformed_large_sliding_geometry_used=False),
            delta_theta_rad=theta_rad,loaded_whole_body_frames=loaded_frames,cell_checks=rows,patch_checks=patches,
            tolerances=dict(mean_error_m=mean_tolerance_m,penetration_m=maximum_penetration_m),whole_flank_nonpenetration_certified=False,
            source_NURBS_gap_reused=False,implementation_revision='native_FE_typed_operator',seconds=time.perf_counter()-started)
    except Exception as error:
        return dict(completed=False,passed=False,failure=dict(type=type(error).__name__,message=str(error)),cell_checks=rows,patch_checks=patches,
            original_candidate_retained=True,implementation_revision='native_FE_typed_operator',seconds=time.perf_counter()-started)


def _manifest_fingerprint(manifest):
    value=dict(generation_pose=manifest['whole_body_R0_frames'],actual_mesh_sha256=manifest['actual_mesh_sha256'],
        geometry_fingerprint=manifest['geometry_fingerprint'],column_keys_sha256=manifest['column_keys_sha256'],
        patches=manifest['patches'],matrices=manifest['matrices'],original_native_cache_index_sha256=manifest.get('native_cache_index_sha256'))
    if 'operator_input_fingerprint' in manifest:value['operator_input_fingerprint']=manifest['operator_input_fingerprint']
    return fingerprint(value)


def save_native_operator(system,directory,*,actual_mesh_sha256=None,source_identity=None):
    """普通typed NPZ+JSON+W归档；manifest不嵌自身SHA，返回独立身份。"""
    import numpy as np
    from scipy import sparse
    from .projected_fe_partition_arrays import typed_array_payload
    if 'typed_patch_arrays' not in system:raise ValueError('当前operator没有typed交片；必须先显式重建并重新求解，不能替换已求解的W')
    directory=Path(directory).resolve()
    if directory.exists():raise FileExistsError('当前operator归档目录已存在，不能覆盖')
    directory.mkdir(parents=True);patches=[]
    for i,a in enumerate(system['typed_patch_arrays']):
        if a is None:patches.append(None);continue
        payload=typed_array_payload(a);file=directory/f'patch_{i:02d}_arrays.npz'
        with file.open('wb') as stream:np.savez_compressed(stream,**payload)
        meta={k:v for k,v in a.items() if k not in ('body_data','body_projection_tables','matrices','piece_facet_indices') and not isinstance(v,np.ndarray)}
        record=_write_json(directory/f'patch_{i:02d}_metadata.json.gz',meta)
        patches.append(dict(arrays_file=file.name,arrays_sha256=sha(file),metadata=record))
    matrices={}
    for m,w in system['matrices'].items():
        file=directory/(m+'_W.npz');sparse.save_npz(file,w,compressed=True);matrices[m]=dict(file=file.name,sha256=sha(file),shape=list(w.shape))
    geometry=_write_json(directory/'geometry.json.gz',system['geometry'])
    keys=np.asarray([(r['tooth_offset'],r['row'],r['column']) for r in system['records']],dtype=np.int64).reshape(-1,3)
    column_sha=hashlib.sha256(keys.tobytes()).hexdigest();geometry_fp=fingerprint(system['geometry'])
    manifest=dict(schema='native_FE_typed_operator_manifest_v2',geometry_source='FE_boundary',integration_method=METHOD,
        actual_mesh_sha256=actual_mesh_sha256,generation_driving_angle_rad=system['geometry']['driving_angle_rad'],
        generation_reference_follower_angle_rad=system['geometry']['reference_follower_angle_rad'],whole_body_R0_frames=system['frames'],
        geometry=geometry,geometry_fingerprint=geometry_fp,column_keys_sha256=column_sha,column_keys_shape=list(keys.shape),
        patches=patches,matrices=matrices,source_identity=source_identity or {},units=dict(length='m',angle='rad',force='N',W='N/N'),
        original_NURBS_gap_used=False,old_compliance_reused=False)
    for name in ('native_cache_index_sha256','operator_input_fingerprint'):
        if name in system:manifest[name]=system[name]
    operator_fp=_manifest_fingerprint(manifest);manifest['operator_fingerprint']=operator_fp
    manifest_file=directory/'manifest.json';_write_json(manifest_file,manifest)
    identity=dict(current_operator_manifest_file=str(manifest_file),current_operator_manifest_sha256=sha(manifest_file),operator_fingerprint=operator_fp)
    system.update(identity);return identity


def load_native_operator(directory):
    """从原manifest字节及全部外存SHA恢复实际g/W/原片视图，不补outer0。"""
    import numpy as np
    from scipy import sparse
    directory=Path(directory).resolve();path=_archive_path(directory,'manifest.json');manifest=_read_json(path)
    if manifest.get('schema') not in ('native_FE_typed_operator_manifest_v1','native_FE_typed_operator_manifest_v2'):raise ValueError('未知原生typed operator归档schema')
    if manifest.get('geometry_source')!='FE_boundary' or manifest.get('integration_method')!=METHOD:raise ValueError('归档须明确真实FE边界积分协议')
    geometry_file=_archive_path(directory,manifest['geometry']['file'])
    if sha(geometry_file)!=manifest['geometry']['sha256']:raise ValueError('operator geometry原字节SHA不符')
    geometry=_restore_geometry(_read_json(geometry_file));patch_arrays=[]
    for entry in manifest['patches']:
        if entry is None:patch_arrays.append(None);continue
        array_file=_archive_path(directory,entry['arrays_file']);meta_file=_archive_path(directory,entry['metadata']['file'])
        if sha(array_file)!=entry['arrays_sha256'] or sha(meta_file)!=entry['metadata']['sha256']:raise ValueError('operator typed数组/metadata SHA不符')
        a=_read_json(meta_file);a['body_data']={};a['body_projection_tables']={}
        with np.load(array_file,allow_pickle=False) as archive:
            for key in archive.files:
                if '__' in key:
                    group,member,name=key.split('__',2);a[group].setdefault(member,{})[name]=archive[key].copy()
                else:a[key]=archive[key].copy()
        a['piece_facet_indices']={m:a['body_data'][m]['facet_indices'] for m in MEMBERS};patch_arrays.append(a)
    matrices={}
    for m,entry in manifest['matrices'].items():
        file=_archive_path(directory,entry['file'])
        if sha(file)!=entry['sha256']:raise ValueError('当前operator W原字节SHA不符')
        matrices[m]=sparse.load_npz(file)
        if list(matrices[m].shape)!=entry['shape']:raise ValueError('当前W布局不符')
    records=TypedRecords(patch_arrays);keys=np.asarray([(r['tooth_offset'],r['row'],r['column']) for r in records],dtype=np.int64).reshape(-1,3)
    if fingerprint(geometry)!=manifest['geometry_fingerprint'] or hashlib.sha256(keys.tobytes()).hexdigest()!=manifest['column_keys_sha256']:
        raise ValueError('当前geometry/列索引数学身份不符')
    expected_fp=_manifest_fingerprint(manifest)
    if expected_fp!=manifest['operator_fingerprint']:raise ValueError('当前operator指纹与manifest内容不符')
    if geometry['driving_angle_rad']!=manifest['generation_driving_angle_rad'] or geometry['reference_follower_angle_rad']!=manifest['generation_reference_follower_angle_rad']:
        raise ValueError('当前operator生成角度与geometry不符')
    if manifest['whole_body_R0_frames'].get('native_reference_follower_angle_rad')!=geometry['reference_follower_angle_rad']:
        raise ValueError('当前operator整体R0参考phi不符')
    offset=0
    for patch,a in zip(geometry['patches'],patch_arrays):
        if fingerprint(patch['frames'])!=fingerprint(manifest['whole_body_R0_frames']):raise ValueError('当前patch未使用同一整体R0姿态')
        count=0 if a is None else a['complete_cell_count']
        if a is not None:a['matrices']={m:w[:,offset:offset+count] for m,w in matrices.items()}
        offset+=count
    if offset!=len(records) or any(w.shape[1]!=len(records) for w in matrices.values()):raise ValueError('当前operator交片与W列数不符')
    if manifest.get('operator_input_fingerprint')!=geometry.get('operator_input_fingerprint'):raise ValueError('当前operator输入指纹与geometry不符')
    result=dict(records=records,matrices=matrices,geometry=geometry,frames=manifest['whole_body_R0_frames'],typed_patch_arrays=patch_arrays,
        geometry_source='FE_boundary',integration_method=METHOD,
        current_operator_manifest_file=str(path),current_operator_manifest_sha256=sha(path),operator_fingerprint=manifest['operator_fingerprint'],
        actual_native_mesh_sha256=manifest['actual_mesh_sha256'],
        maximum_cad_fe_mapping_error_m=max((r['maximum_mapping_error_m'] for r in records),default=0.),
        pressure_basis='当前phi下原typed FE交片P0归档',ray_backend='actual_FE_triangle_intersection_typed_arrays',
        quadrature_order='exact_affine_partition_centroids',piece_vertex_validation_data_included=True,
        original_NURBS_gap_used=False,halfspace_compliance_used=False,old_compliance_reused=False)
    for name in ('native_cache_index_sha256','operator_input_fingerprint'):
        if name in manifest:result[name]=manifest[name]
    return result
