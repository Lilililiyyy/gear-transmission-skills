"""原生FE案例的纯几何／网格准备；不装配K、不求压力。

新写业务接线，复用现有CAD、TCA、有限NURBS／OCP查询及Gmsh内核。
固定φ／plane／R0／Rk都是显式几何参数，旧压力、gap、W、C不作为输入。
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import fields
import gc
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

MEMBERS=('pinion','wheel')
FRAME_FIELDS={'pinion_rotation_rows','wheel_rotation_rows','pinion_axis_world','wheel_axis_world',
    'translation_m','pinion_translation_world_m','wheel_translation_world_m',
    'pinion_axis_origin_world_m','wheel_axis_origin_world_m','native_reference_follower_angle_rad'}


def file_sha256(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block:=stream.read(2**20):value.update(block)
    return value.hexdigest()


def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def json_equal(first,second):
    """只规范tuple/list的JSON表示，不放宽制造数值身份。"""
    return json.loads(json.dumps(first,allow_nan=False))==json.loads(json.dumps(second,allow_nan=False))


def write_json(path,value):
    """流式输出，避免asdict深复制整个百万节点网格。"""
    def convert(item):
        for name in ('tolist','item'):
            if callable(getattr(item,name,None)):return getattr(item,name)()
        raise TypeError('不能保存的案例数值类型：'+type(item).__name__)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf8') as stream:
        json.dump(value,stream,ensure_ascii=False,allow_nan=False,separators=(',',':'),default=convert);stream.write('\n')


def _real(value,name,*,positive=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or positive and value<=0:
        raise ValueError(name+'须为有限'+('正' if positive else '')+'实数')
    return float(value)


def _integer(value,name,minimum=0):
    if isinstance(value,bool) or not isinstance(value,int) or value<minimum:raise ValueError(name+'须为范围内整数')
    return value


def _relative(parent,name):
    if not isinstance(name,str) or not name or Path(name).anchor or '..' in Path(name).parts:raise ValueError('案例引用须为目录内相对路径')
    parent=Path(parent).resolve();path=(parent/name).resolve()
    if not path.is_relative_to(parent):raise ValueError('案例引用越出其目录')
    return path


def _vector(value,length,name):
    if not isinstance(value,(list,tuple)) or len(value)!=length:raise ValueError(name+'维数不符')
    return [_real(x,name) for x in value]


def validate_frames(frames):
    """只接收整体运动／材料参考几何，不接受旧物理结果字段。"""
    if not isinstance(frames,dict) or set(frames)-FRAME_FIELDS:raise ValueError('frames含未知或非几何字段')
    for member in MEMBERS:
        matrix=frames.get(member+'_rotation_rows')
        if not isinstance(matrix,(list,tuple)) or len(matrix)!=3:raise ValueError('需要三维旋转矩阵')
        R=[_vector(row,3,'旋转矩阵') for row in matrix]
        if any(abs(sum(a*b for a,b in zip(R[i],R[j]))-(1 if i==j else 0))>2e-12 for i in range(3) for j in range(3)):
            raise ValueError('R0／Rk须为正交旋转')
        determinant=sum(R[0][j]*(R[1][(j+1)%3]*R[2][(j+2)%3]-R[1][(j+2)%3]*R[2][(j+1)%3]) for j in range(3))
        if abs(determinant-1)>2e-12:raise ValueError('旋转矩阵行列式须为1')
        axis=_vector(frames.get(member+'_axis_world'),3,'世界轴')
        if abs(sum(x*x for x in axis)-1)>2e-12:raise ValueError('世界轴须为单位向量')
    for name in FRAME_FIELDS-set(('pinion_rotation_rows','wheel_rotation_rows','pinion_axis_world','wheel_axis_world')):
        if name not in frames:continue
        if name=='native_reference_follower_angle_rad':_real(frames[name],name)
        else:
            vector=_vector(frames[name],3,name)
            if any(abs(x)>1e-14 for x in vector):raise ValueError('当前圆锥CAD准备仅支持共顶点轴系，不支持偏置／准双曲面')


def validate_grid(grid):
    required={'x_count','y_count','spacing_x_m','spacing_y_m','origin_xy_m'}
    if not isinstance(grid,dict) or set(grid)!=required:raise ValueError('压力格仅允许明确的维数／间距／原点')
    nx=_integer(grid['x_count'],'x_count',3);ny=_integer(grid['y_count'],'y_count',3)
    if nx*ny>500000:raise ValueError('单片压力格超过五十万格上限')
    for name in ('spacing_x_m','spacing_y_m'):_real(grid[name],name,positive=True)
    _vector(grid['origin_xy_m'],2,'压力格原点')


def validate_plane(plane):
    if not isinstance(plane,dict) or set(plane)!={'origin_world_m','tangent_rows','normal_world'}:raise ValueError('plane只允许原点／切基／法向')
    _vector(plane['origin_world_m'],3,'plane原点')
    if not isinstance(plane['tangent_rows'],(list,tuple)) or len(plane['tangent_rows'])!=2:raise ValueError('plane须有两个切向基')
    rows=[_vector(v,3,'plane基') for v in plane['tangent_rows']]+[_vector(plane['normal_world'],3,'plane法向')]
    if any(abs(sum(a*b for a,b in zip(rows[i],rows[j]))-(1 if i==j else 0))>2e-12 for i in range(3) for j in range(3)):
        raise ValueError('plane切基／法向须正交归一')


def load_case(path,*,phase_indices=None):
    """读取明确相位列表；返回原文档、base及路径绑定，不创建输出或引擎。"""
    path=Path(path).resolve();case=read_json(path)
    allowed={'schema','name','description','base_case','mesh','phases','pressure_grids','solver','memory_policy','refinement_phase_policy'}
    if not isinstance(case,dict) or set(case)-allowed or case.get('schema')!='bevel_native_fe_case_v1':raise ValueError('需要bevel_native_fe_case_v1及明确字段')
    if case.get('memory_policy','full') not in ('full','contact_only'):raise ValueError('memory_policy只允许full或contact_only')
    if case.get('refinement_phase_policy','require_local_sources') not in ('require_local_sources','allow_global_background'):
        raise ValueError('refinement_phase_policy只允许require_local_sources或allow_global_background')
    base_file=_relative(path.parent,case['base_case']);base=read_json(base_file)
    if base.get('schema')!='bevel_contact_manufacturing_case_v1' or base.get('kind')!='spiral':raise ValueError('当前制造案例准备需要明确弧齿锥齿轮输入')
    phases=case.get('phases')
    if not isinstance(phases,list) or not phases:raise ValueError('phases须为非空实际相位列表')
    angles=[];offsets=None
    for phase in phases:
        if not isinstance(phase,dict) or set(phase)!={'driving_angle_rad','reference_follower_angle_rad','frames','patches'}:raise ValueError('每相仅允许明确角／frames／patches，不接收旧gap／压力／C')
        angles.append(_real(phase['driving_angle_rad'],'主动角rad'));phi=_real(phase['reference_follower_angle_rad'],'参考从动角rad')
        validate_frames(phase['frames'])
        if abs(phase['frames'].get('native_reference_follower_angle_rad',phi)-phi)>1e-12:raise ValueError('frames参考φ与相位φ不同')
        patches=phase['patches']
        if not isinstance(patches,list) or not patches:raise ValueError('每相须有明确物理齿patch')
        current=[]
        for p in patches:
            if not isinstance(p,dict) or set(p)!={'tooth_offset','grid','plane','source_material_reference_frames'}:raise ValueError('patch只允许纯几何定义，mask必须fresh查询')
            if isinstance(p['tooth_offset'],bool) or not isinstance(p['tooth_offset'],int):raise ValueError('物理齿偏移须为整数')
            current.append(p['tooth_offset']);validate_grid(p['grid']);validate_plane(p['plane']);validate_frames(p['source_material_reference_frames'])
        if len(current)!=len(set(current)) or 0 not in current:raise ValueError('物理齿偏移不重复且须含参考齿0')
        if offsets is not None and current!=offsets:raise ValueError('同案例各相须保持同一明确物理齿窗')
        offsets=current
    if any(b<=a for a,b in zip(angles,angles[1:])):raise ValueError('实际相位角须严格递增，不允许重复／占位相位')
    selected=list(range(len(phases))) if phase_indices is None else list(phase_indices)
    if not selected or len(selected)!=len(set(selected)) or any(isinstance(i,bool) or not isinstance(i,int) or not 0<=i<len(phases) for i in selected):raise ValueError('选择相位须为原列表中不重复索引')
    levels=case.get('pressure_grids');names=[]
    if not isinstance(levels,list) or not levels:raise ValueError('pressure_grids须为非空列表')
    for level in levels:
        if not isinstance(level,dict) or set(level)!={'name','spacing_m'}:raise ValueError('每档仅给名称与两个SI格距')
        name=level['name']
        if not isinstance(name,str) or not name or any(not (x.isascii() and (x.isalnum() or x in '_-')) for x in name):raise ValueError('压力档名须为安全独立名称')
        values=_vector(level['spacing_m'],2,'压力格距m')
        if any(x<=0 for x in values):raise ValueError('压力格距须为正')
        names.append(name)
    if len(names)!=len(set(names)):raise ValueError('压力档名不能重复')
    mesh=case.get('mesh');mesh_allowed={'element_size_m','contact_refined_size_m','refinement_core_radius_m','refinement_radius_m','maximum_nodes',
        'refinement_phase_intervals','refinement_half_widths_m','track_spacing_m','sample_spacing_m','refinement_regions_file'}
    if not isinstance(mesh,dict) or set(mesh)!=mesh_allowed:raise ValueError('mesh须明确给出主体／背景／区域及轨迹控制')
    for key in ('element_size_m','contact_refined_size_m','refinement_core_radius_m','refinement_radius_m','track_spacing_m','sample_spacing_m'):_real(mesh[key],key,positive=True)
    if not mesh['contact_refined_size_m']<mesh['element_size_m'] or not mesh['contact_refined_size_m']<=mesh['refinement_core_radius_m']<mesh['refinement_radius_m']:raise ValueError('背景场须fine<bulk且fine≤core<radius')
    if any(x<=0 for x in _vector(mesh['refinement_half_widths_m'],2,'轨迹半宽m')):raise ValueError('轨迹半宽须为正')
    _integer(mesh['maximum_nodes'],'maximum_nodes',4);_integer(mesh['refinement_phase_intervals'],'refinement_phase_intervals',1)
    regions_file=_relative(path.parent,mesh['refinement_regions_file'])
    return dict(case=case,base_case=base,case_file=path,base_file=base_file,regions_file=regions_file,
        selected_phase_indices=selected,tooth_offsets=offsets,case_sha256=file_sha256(path),base_sha256=file_sha256(base_file),regions_sha256=file_sha256(regions_file))



def load_case_refinement_regions(loaded,loader):
    """显式允许背景网格时，只区分局部轨迹来源与求解相位，不伪造来源。

    默认仍要求所选每相都有声明的局部来源。允许背景网格仅改变案例入场策略；
    CAD未知、制造/blank身份及每相原生FE覆盖门均不因此豁免。
    """
    case=loaded['case'];policy=case.get('refinement_phase_policy','require_local_sources')
    if policy not in ('require_local_sources','allow_global_background'):raise ValueError('未知的局部加密相位使用策略')
    selected=loaded['selected_phase_indices'];angles=[case['phases'][i]['driving_angle_rad'] for i in selected]
    document=read_json(loaded['regions_file'])
    validation_angles=angles if policy=='require_local_sources' else document.get('source_driver_angles_rad',[])
    regions,digest,raw=loader(loaded['regions_file'],validation_angles,case['mesh']['element_size_m'])
    if not json_equal(regions['manufacturing'],normalized_manufacturing(loaded['base_case'])):raise ValueError('局部区域的制造/齿坯与本案例不一致，背景策略不能绕过身份门')
    sources=regions['source_driver_angles_rad'];rows=[]
    for index,angle in zip(selected,angles):
        local=any(math.isclose(angle,value,rel_tol=0,abs_tol=1e-12) for value in sources)
        rows.append(dict(phase_index=index,driving_angle_rad=angle,local_source_declared=local,
            mesh_use='declared_local_source' if local else 'global_background_without_this_phase_local_source'))
    scope=dict(policy=policy,source_driver_angles_rad=list(sources),phases=rows,
        phases_using_global_background=[r['phase_index'] for r in rows if not r['local_source_declared']],
        background_track_refined_size_m=case['mesh']['contact_refined_size_m'],bulk_element_size_m=case['mesh']['element_size_m'],
        local_region_refined_sizes_m={m:[r['refined_size_m'] for r in regions['regions'][m]] for m in MEMBERS},
        all_selected_phases_have_declared_local_sources=all(r['local_source_declared'] for r in rows),
        whole_cycle_local_fine_coverage_certified=False,spatial_error_certified=False,
        original_region_source_angles_preserved=True,original_API_pre_QP_FE_coverage_required=True,
        note='尺寸为Gmsh场控制，不保证每个单元均小于该尺寸；局部来源相位也不是实际受载覆盖或误差证明。')
    return regions,digest,raw,scope


def repartition_grid(original,dx,dy):
    """中心不动，以最小整数格覆盖原cell外边界；不缩放／搬运旧mask。"""
    validate_grid(original);dx=_real(dx,'dx',positive=True);dy=_real(dy,'dy',positive=True)
    nx,ny=original['x_count'],original['y_count'];sx,sy=original['spacing_x_m'],original['spacing_y_m']
    origin=original['origin_xy_m'];centre=[origin[0]+(nx-1)*sx/2,origin[1]+(ny-1)*sy/2]
    shape=[math.ceil(nx*sx/dx-1e-12),math.ceil(ny*sy/dy-1e-12)]
    result=dict(x_count=shape[0],y_count=shape[1],spacing_x_m=dx,spacing_y_m=dy,
        origin_xy_m=[centre[0]-(shape[0]-1)*dx/2,centre[1]-(shape[1]-1)*dy/2])
    validate_grid(result)
    return result,dict(original_cell_window_m=[nx*sx,ny*sy],new_cell_window_m=[shape[0]*dx,shape[1]*dy],
        minimum_symmetric_padding_each_side_m=[max(0.,(shape[0]*dx-nx*sx)/2),max(0.,(shape[1]*dy-ny*sy)/2)])


def case_plan(loaded):
    """只列尺寸与相位作用域；不把实际相位列表当连续周期证书。"""
    case=loaded['case'];rows=[]
    for level in case['pressure_grids']:
        for i in loaded['selected_phase_indices']:
            phase=case['phases'][i]
            rows.append(dict(level=level['name'],phase_index=i,driving_angle_rad=phase['driving_angle_rad'],
                patches=[dict(tooth_offset=p['tooth_offset'],grid=repartition_grid(p['grid'],*level['spacing_m'])[0]) for p in phase['patches']]))
    return dict(schema='bevel_native_FE_case_plan_v1',selected_phase_indices=loaded['selected_phase_indices'],rows=rows,
        mesh_K_pressure_created=False,geometry_source='FE_boundary',full_cycle_certified=False,
        note='phases为实际离散列表；多压力grid不是多啮合相位，φ为数值初值')


def build_case_models(loaded):
    """由公开制造输入重建真实双闭实体，不读取既有BRep或受载状态。"""
    from .generated_solid import build_generated_gear
    base=loaded['base_case'];sampling=base.get('cad_sampling',{})
    return {m:build_generated_gear(base['manufacturing'][m]['request'],base['manufacturing'][m]['blank'],**sampling) for m in MEMBERS}


def normalized_manufacturing(base_case):
    """按现有输入类型展开默认值；已显式给出的raw值不得被另一个制造身份替换。"""
    from dataclasses import asdict
    from .generated_surface import _request
    from .generated_solid import GeneratedBlankSettings
    return {m:dict(request=asdict(_request(base_case['manufacturing'][m]['request'])),
        blank=asdict(GeneratedBlankSettings(**base_case['manufacturing'][m]['blank']))) for m in MEMBERS}


def generate_case_meshes(loaded,models,regions,output,*,cad_brep_files,progress=None):
    """只调用CAD／TCA／Gmsh及JSON保存，始终不prepare弹性或创建K。"""
    from .cad_import import ImportedCadModel,export_model_step
    from .generated_tca import prepare_generated_tca
    from .generated_nodal_ltca import generated_refinement_tracks
    from .elasticity import mesh_step_tetra
    from dataclasses import asdict
    case=loaded['case'];controls=case['mesh'];output=Path(output)
    if set(cad_brep_files)!=set(MEMBERS):raise ValueError('网格必须绑定本次双轮真实BRep文件')
    cad_sha={m:file_sha256(cad_brep_files[m]) for m in MEMBERS}
    output.mkdir(parents=True,exist_ok=False)
    tca=prepare_generated_tca(models['pinion'].request,models['wheel'].request,cad_models=models)
    tracks={m:[] for m in MEMBERS};z=models['pinion'].request.tooth_count;n=controls['refinement_phase_intervals']
    for index in range(n+1):
        angle=-math.pi/z+2*math.pi/z*index/n
        result=generated_refinement_tracks(tca,driving_angle_rad=angle,tooth_offsets=loaded['tooth_offsets'],
            half_widths_m=controls['refinement_half_widths_m'],track_spacing_m=controls['track_spacing_m'],sample_spacing_m=controls['sample_spacing_m'])
        for m in MEMBERS:tracks[m].extend(result['tracks'][m])
    write_json(output/'refinement_tracks.json',dict(tracks=tracks,geometry_only=True))
    meshes={};rows={}
    for m in MEMBERS:
        model=models[m]
        if not json_equal(regions['manufacturing'][m],dict(request=asdict(model.request),blank=asdict(model.settings))):raise ValueError('加密区域与当前标准化制造／齿坯身份不同')
        step=output/(m+'.step');export_model_step(ImportedCadModel(model.shape,model.diagnostics),step)
        if progress:progress(dict(stage='mesh_started',member=m))
        mesh=mesh_step_tetra(step,element_size_m=controls['element_size_m'],maximum_nodes=controls['maximum_nodes'],
            refinement_polylines_m=tracks[m],refined_size_m=controls['contact_refined_size_m'],
            refinement_core_radius_m=controls['refinement_core_radius_m'],refinement_radius_m=controls['refinement_radius_m'],
            refinement_regions=regions['regions'][m])
        mesh.metadata.update(actual_manufacturing_request=asdict(model.request),actual_blank_settings=asdict(model.settings),
            actual_relief_cad_sha256=cad_sha[m])
        target=output/(m+'_mesh.json');write_json(target,{f.name:getattr(mesh,f.name) for f in fields(mesh)})
        rows[m]=dict(file=target.name,sha256=file_sha256(target),nodes=len(mesh.nodes_m),tetrahedra=len(mesh.tetrahedra))
        meshes[m]=mesh
        if progress:progress(dict(stage='mesh_completed',member=m,**rows[m]))
    write_json(output/'mesh_manifest.json',dict(schema='bevel_native_FE_case_mesh_v1',meshes=rows,K_assembled=False,pressure_solved=False,
        case_sha256=loaded['case_sha256'],regions_sha256=loaded['regions_sha256'],actual_cad_sha256=cad_sha))
    return meshes


def empty_geometry_patch_evidence(evidence,grid):
    """仅明确全查询的有限面外／法向不配对允许空请求；不签发零压力证明。"""
    validate_grid(grid)
    if (evidence.get('completed') is not True or evidence.get('all_center_and_corner_positions_requeried') is not True or
            evidence.get('requested_cells')!=0 or evidence.get('total_cells')!=grid['x_count']*grid['y_count']):
        raise ArithmeticError('空patch缺少全部中心／角点几何查询证据')
    allowed={'proven_outside_positive_NURBS_control_hull','opposing_normal_limit'}
    allowed.update(m+':'+s for m in MEMBERS for s in ('outside_trimmed_face','OCP_proven_outside_trimmed_face'))
    totals={'centers':grid['x_count']*grid['y_count'],'corners':(grid['x_count']+1)*(grid['y_count']+1)}
    for name,total in totals.items():
        row=evidence.get(name,{})
        counts=row.get('reason_counts')
        if (row.get('unresolved_count')!=0 or row.get('unresolved')!=[] or not isinstance(counts,dict) or not counts or
                any(key not in allowed or isinstance(n,bool) or not isinstance(n,int) or n<=0 for key,n in counts.items()) or
                sum(counts.values())!=total):
            raise ArithmeticError('空patch仍含未知／模型斜率限制／有效点，不能当有限面外或不配对')
    return dict(empty_patch_retained=True,all_original_patch_fields_retained=True,
        all_queries_completed=True,unknown_count=0,
        reason='全部中心和角点明确有限面外或法向不配对',pressure_zero_certified=False,
        material_domain_certified=False,geometry_only=True)


def build_geometry_templates(loaded,prepared_surface,output,*,query_factory=None,progress=None):
    """逐相逐档fresh查询并流式存模板；不把全周期mask同时常驻。"""
    from .generated_tca import solve_generated_tca
    output=Path(output);output.mkdir(parents=True,exist_ok=False);case=loaded['case'];base=loaded['base_case']
    torque=base['load']['ltca_torque_sign']*base['load']['pinion_torque_n_m']
    ratio=base['manufacturing']['pinion']['request']['tooth_count']/base['manufacturing']['wheel']['request']['tooth_count']
    entries=[]
    for i in loaded['selected_phase_indices']:
        phase=case['phases'][i];template=dict(frames=deepcopy(phase['frames']),driving_angle_rad=phase['driving_angle_rad'],
            reference_follower_angle_rad=phase['reference_follower_angle_rad'],pinion_torque_nm=torque,pinion_to_wheel_ratio=ratio,patches=deepcopy(phase['patches']))
        geometry=dict(tca_results={k:solve_generated_tca(prepared_surface.tca,phase['driving_angle_rad'],tooth_offset=k) for k in loaded['tooth_offsets']}) if query_factory is None else {}
        query=(FiniteCADQueries(prepared_surface,geometry,template,SimpleNamespace(check=lambda:None)) if query_factory is None
            else query_factory(phase))
        for level in case['pressure_grids']:
            current=deepcopy(template);diagnostics=[]
            for patch in current['patches']:
                grid,window=repartition_grid(patch['grid'],*level['spacing_m']);patch['grid']=grid
                mask,evidence=evaluate_grid(query.query,patch['tooth_offset'],grid)
                if evidence['completed'] is not True:
                    write_json(output/'unresolved'/level['name']/f'phase_{i:04d}_k{patch["tooth_offset"]}.json',
                        dict(phase_index=i,level=level['name'],tooth_offset=patch['tooth_offset'],evidence=evidence,unknown_is_not_OUT=True))
                    raise ArithmeticError('fresh CAD查询仍有未决点，不能当OUT或生成完整模板')
                empty=None if mask.any() else empty_geometry_patch_evidence(evidence,grid)
                patch['requested_mask']=mask
                diagnostics.append(dict(tooth_offset=patch['tooth_offset'],window=window,empty_geometry_patch=empty,**evidence))
            if not any(p['requested_mask'].any() for p in current['patches']):
                write_json(output/'rejected_all_empty'/level['name']/f'phase_{i:04d}.json',
                    dict(phase_index=i,level=level['name'],patches=diagnostics,
                        original_tooth_offsets=[p['tooth_offset'] for p in current['patches']],pressure_zero_certified=False))
                raise ArithmeticError('所有原物理齿patch合计请求为空，拒绝创建接触模板')
            relative=Path(level['name'])/f'phase_{i:04d}.json';file=output/relative
            write_json(file,dict(schema='bevel_native_FE_geometry_template_v1',template=current))
            write_json(file.with_name(file.stem+'.geometry_checks.json'),dict(phase_index=i,level=level['name'],patches=diagnostics,
                invalid_is_not_zero_pressure=True,all_queries_fresh=True,K_created=False))
            entries.append(dict(name=level['name'],phase_index=i,driving_angle_rad=phase['driving_angle_rad'],file=relative.as_posix(),sha256=file_sha256(file)))
            if progress:progress(dict(stage='fresh_geometry_template_completed',phase_index=i,level=level['name']))
        del query;gc.collect()
    manifest=dict(schema='bevel_native_FE_geometry_templates_v1',entries=entries,selected_phase_indices=loaded['selected_phase_indices'],
        case_sha256=loaded['case_sha256'],old_gap_pressure_W_C_used=False,K_created=False,full_cycle_certified=False)
    write_json(output/'manifest.json',manifest);return manifest


class FiniteCADQueries:
    """原批量NURBS/裁切分类；Newton未解时以OCP有限面射线补证或重给初值。"""
    def __init__(self,prepared,geometry,template,gate,*,batch_size=256,tolerance_m=2e-11):
        import numpy as np
        from scipy.spatial import ConvexHull
        from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
        from OCP.BRep import BRep_Builder
        from OCP.TopoDS import TopoDS_Compound
        from bevel_gear.generated_tca import generated_tca_frames
        self.np=np;self.prepared=prepared;self.gate=gate;self.batch_size=batch_size;self.tolerance=tolerance_m;self.contexts={};self.engines={}
        self.maximum_residual_m=0.;self.fallback_count=0;self.query_count=0
        phi=template['reference_follower_angle_rad'];alpha=template['driving_angle_rad']
        self.direction=1 if template['pinion_torque_nm']>0 else -1
        for member in MEMBERS:
            compound=TopoDS_Compound();builder=BRep_Builder();builder.MakeCompound(compound)
            for face in prepared.trimmed[member].faces:builder.Add(compound,face)
            engine=IntCurvesFace_ShapeIntersector();engine.Load(compound,1e-7);self.engines[member]=engine
        for patch in template['patches']:
            k=patch['tooth_offset'];frames=generated_tca_frames(prepared.tca,alpha,phi,tooth_offset=k)
            if k==0:
                for member in MEMBERS:
                    if not np.allclose(frames[member+'_rotation_rows'],template['frames'][member+'_rotation_rows'],rtol=0,atol=2e-13):raise ValueError('当前CAD整体R0与显式几何显式参考角不一致')
            tca=geometry['tca_results'][str(k)] if str(k) in geometry['tca_results'] else geometry['tca_results'][k]
            seeds=({m:tca[m+'_uv'] for m in MEMBERS} if tca['converged'] else dict(pinion=tca['candidate']['parameters'][:2],wheel=tca['candidate']['parameters'][2:4]))
            plane=patch['plane'];origin=np.asarray(plane['origin_world_m']);basis=np.asarray(plane['tangent_rows']);normal=np.asarray(plane['normal_world'])
            from bevel_gear.generated_halfspace import _initial_uv
            maps={m:_initial_uv(prepared.trimmed[m].surface,np.asarray(frames[m+'_rotation_rows']),origin,basis,seeds[m]) for m in MEMBERS}
            hulls={};radii={}
            for member in MEMBERS:
                surface=prepared.trimmed[member].surface;values=[]
                for i in range(1,surface.NbUPoles()+1):
                    for j in range(1,surface.NbVPoles()+1):
                        if surface.Weight(i,j)<=0:raise ValueError('凸包排除证明要求原NURBS权重全正')
                        values.append(surface.Pole(i,j).Coord())
                values=np.asarray(values)*.001;rotation=np.asarray(frames[member+'_rotation_rows'])
                projected=(values@rotation.T-origin)@basis.T;hulls[member]=ConvexHull(projected).equations
                radii[member]=float(np.linalg.norm(values,axis=1).max())
            self.contexts[k]=dict(frames=frames,origin=origin,basis=basis,normal=normal,seeds=seeds,maps=maps,hulls=hulls,radii=radii)

    def _finite_face_seed(self,context,member,xy):
        """OCP求实际有限面交点；有界段由正权控制凸包给出完整长度上界。"""
        from OCP.gp import gp_Pnt,gp_Dir,gp_Lin
        np=self.np;rotation=np.asarray(context['frames'][member+'_rotation_rows'])
        point=rotation.T@(context['origin']+np.asarray(xy)@context['basis']);direction=rotation.T@context['normal']
        extent=(np.linalg.norm(point)+context['radii'][member])*1000+1.
        engine=self.engines[member];engine.Perform(gp_Lin(gp_Pnt(*tuple(point*1000)),gp_Dir(*tuple(direction))),-extent,extent)
        if not engine.IsDone():return None,member+':OCP_finite_face_intersection_unresolved'
        selected=[]
        for i in range(1,engine.NbPnt()+1):
            value=float(engine.WParameter(i))*.001
            if not any(abs(value-old[0])<=self.tolerance for old in selected):selected.append((value,i))
        if not selected:return None,member+':OCP_proven_outside_trimmed_face'
        if len(selected)!=1:return None,member+':OCP_multiple_surface_heights'
        trimmed=self.prepared.trimmed[member]
        if not all(trimmed.uv_compatible):return None,member+':OCP_parameter_mapping_unverified'
        i=selected[0][1];b=trimmed.surface.Bounds()
        return [(engine.UParameter(i)-b[0])/(b[1]-b[0]),(engine.VParameter(i)-b[2])/(b[3]-b[2])],None

    def _fallback(self,context,xy):
        from bevel_gear.generated_halfspace import _sample_pair
        seeds={}
        for member in MEMBERS:
            value,reason=self._finite_face_seed(context,member,xy)
            if reason:return None,reason
            seeds[member]=value
        return _sample_pair(self.prepared,context['frames'],context['origin'],context['basis'],context['normal'],
            seeds,xy,self.tolerance,.3,.3)

    def query(self,k,xy):
        """返回每点域/模型状态；未知不当成材料外，不保存CAD间隙作为FE输入。"""
        from bevel_gear.nurbs_batch import sample_pair_batch
        np=self.np;context=self.contexts[k];xy=np.asarray(xy,dtype=float);n=len(xy)
        status=np.full(n,'proven_outside_positive_NURBS_control_hull',dtype=object);valid=np.zeros(n,bool);torque=np.zeros(n,bool)
        candidate=np.ones(n,bool)
        for eq in context['hulls'].values():candidate&=(xy@eq[:,:2].T+eq[:,2]).max(axis=1)<=self.tolerance
        ids=np.flatnonzero(candidate)
        for start in range(0,len(ids),self.batch_size):
            self.gate.check();selected=ids[start:start+self.batch_size];targets=xy[selected]
            samples,reasons=sample_pair_batch(self.prepared,context['frames'],context['origin'],context['basis'],context['normal'],
                context['seeds'],targets,self.tolerance,.3,.3,context['maps'])
            for local,index in enumerate(selected):
                sample,reason=samples[local],reasons[local]
                if sample is None and not any(s in str(reason) for s in ('outside_trimmed_face','halfspace_slope_limit','opposing_normal_limit')):
                    self.fallback_count+=1;sample,reason=self._fallback(context,targets[local])
                if sample is None:status[index]=str(reason);continue
                valid[index]=True;torque[index]=all(self.direction*sample['lever_arms_m'][m]>1e-10 for m in MEMBERS);status[index]='valid'
                self.maximum_residual_m=max(self.maximum_residual_m,*(float(sample[v]['residual_m']) for v in ('first','second')))
            del samples,reasons;gc.collect()
        self.query_count+=n
        known_rejections=('outside_trimmed_face','control_hull','halfspace_slope_limit','opposing_normal_limit')
        unresolved=np.array([s!='valid' and not any(part in s for part in known_rejections) for s in status],dtype=bool)
        return dict(valid=valid,torque=torque,status=status,unresolved=unresolved)


def evaluate_grid(query,k,grid):
    """全部新中心与独立四角重查；返回mask及完整未知/拒绝位置证据。"""
    import numpy as np
    nx,ny=grid['x_count'],grid['y_count'];dx,dy=grid['spacing_x_m'],grid['spacing_y_m'];ox,oy=grid['origin_xy_m']
    def locations(corners):
        x=ox+(np.arange(nx+int(corners))-.5*corners)*dx;y=oy+(np.arange(ny+int(corners))-.5*corners)*dy
        xx,yy=np.meshgrid(x,y);return np.column_stack((xx.ravel(),yy.ravel()))
    corners=query(k,locations(True));centers=query(k,locations(False));cv=corners['valid'].reshape(ny+1,nx+1)
    complete=cv[:-1,:-1]&cv[1:,:-1]&cv[:-1,1:]&cv[1:,1:]
    mask=complete&centers['valid'].reshape(ny,nx)&centers['torque'].reshape(ny,nx)
    def evidence(value,shape):
        status=value['status'];keys,counts=np.unique(status,return_counts=True);unknown=np.flatnonzero(value['unresolved'])
        return dict(reason_counts={str(k):int(v) for k,v in zip(keys,counts)},unresolved_count=len(unknown),
            unresolved=[dict(row=int(i//shape[1]),column=int(i%shape[1]),reason=str(status[i])) for i in unknown])
    corner_e=evidence(corners,(ny+1,nx+1));center_e=evidence(centers,(ny,nx))
    return mask,dict(total_cells=nx*ny,requested_cells=int(mask.sum()),all_center_and_corner_positions_requeried=True,
        torque_direction_rejected_centers=int((centers['valid']&~centers['torque']).sum()),
        incomplete_CAD_corner_cells=int((~complete).sum()),corners=corner_e,centers=center_e,
        completed=corner_e['unresolved_count']==0 and center_e['unresolved_count']==0)

