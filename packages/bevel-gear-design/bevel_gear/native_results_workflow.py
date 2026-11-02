"""正式原生FE结果的分阶段零K后处理；不建立接触刚度或数值因子。

response只读原mesh与物理约束。geometry/domain才恢复同workflow模型档案，
maps/cycle复用本次实际域证据。未知区域与已核零压始终分别保存。
"""
from copy import deepcopy
import gc
import math
from pathlib import Path

from .native_results_archive import (sha,read,relative,fingerprint,load_native_cycle_source,
    load_native_phase_from_cycle,audit_native_phase_response)
from .native_results_pressure import (_write,_bound_domain,export_native_pressure_source,
    audit_native_phase_domain,build_native_cycle_maximum)


MEMBERS=('pinion','wheel')


def load_native_workflow_meshes(source):
    """逐体读取同SHA原JSON，严格检查整数索引后才规范为int64。"""
    import numpy as np
    from .elasticity import TetraMesh
    workflow=source['inputs'].get('workflow_inputs');root=source.get('workflow_root')
    if not isinstance(workflow,dict) or root is None or set(workflow.get('mesh_files',{}))!=set(MEMBERS):raise ValueError('须有同workflow两体原mesh引用')
    meshes={};paths={}
    for member in MEMBERS:
        path=relative(root,workflow['mesh_files'][member]);digest=sha(path)
        if digest!=source['mesh_sha256'][member]:raise ValueError('原mesh字节不属于本态')
        try:
            import orjson
            raw=orjson.loads(path.read_bytes())
        except ImportError:
            raw,_=read(path,digest)
        if set(raw)!={'nodes_m','tetrahedra','boundary_triangles','boundary_tags','metadata'}:raise ValueError('未知TetraMesh档案布局')
        nodes=np.asarray(raw['nodes_m'],dtype=np.float64)
        if nodes.ndim!=2 or nodes.shape[1:]!=(3,) or not len(nodes) or not np.isfinite(nodes).all():raise ValueError('原节点须为非空有限N×3 SI坐标')
        raw['nodes_m']=nodes
        for name,width in (('tetrahedra',4),('boundary_triangles',3),('boundary_tags',None)):
            value=np.asarray(raw[name])
            if value.dtype.kind not in 'iu' or value.dtype.kind=='b' or value.size and (value.min()<0 or value.max()>2**63-1):raise ValueError('原mesh索引/tag须为非负有符号64位范围整数，不能截断float')
            if (width is None and value.ndim!=1) or (width is not None and (value.ndim!=2 or value.shape[1:]!=(width,))):raise ValueError('原mesh数组布局非法：'+name)
            raw[name]=value.astype(np.int64,copy=False)
        if len(raw['boundary_tags'])!=len(raw['boundary_triangles']) or raw['tetrahedra'].size and raw['tetrahedra'].max()>=len(nodes) or raw['boundary_triangles'].size and raw['boundary_triangles'].max()>=len(nodes):raise ValueError('原mesh节点索引或边界tag范围非法')
        meshes[member]=TetraMesh(**raw);paths[member]=path
        if sha(path)!=digest:raise ValueError('读取期间原mesh改变')
        del raw,nodes;gc.collect()
    return meshes,paths


def physical_fixed_nodes(source,meshes,*,boundary_conditions=None):
    """由原孔径/原孔壁容差选约束；不从归档free_dofs反推边界。

    显式plane仅供workflow标为synthetic_FE_fixture_no_CAD的数值块；真实
    齿轮禁止此覆盖。原auditor随后逐位比较真实保存的free DOF。
    """
    import numpy as np
    from .native_fe_operator import safe
    from .distributed_pressure_domain import numeric_array_identity
    workflow=source['inputs'].get('workflow_inputs',{});fixture=workflow.get('cad_binding_state')=='synthetic_FE_fixture_no_CAD'
    fixed={};records={}
    if boundary_conditions is not None and (not fixture or not isinstance(boundary_conditions,dict) or set(boundary_conditions)!=set(MEMBERS)):raise ValueError('显式plane约束只接受明确synthetic无CAD两体fixture，真实齿轮不能覆盖原孔壁')
    for member in MEMBERS:
        nodes=meshes[member].nodes_m
        if boundary_conditions is not None:
            row=boundary_conditions[member]
            if not isinstance(row,dict) or set(row)!={'kind','axis','coordinate_m','tolerance_m'} or row['kind']!='plane':raise ValueError('fixture物理边界须明确plane/axis/coordinate_m/tolerance_m')
            axis=row['axis'];coordinate=row['coordinate_m'];tolerance=row['tolerance_m']
            if isinstance(axis,bool) or not isinstance(axis,int) or axis not in (0,1,2) or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in (coordinate,tolerance)) or not 0<tolerance<=1e-9:raise ValueError('fixture支承平面参数非法')
            selected=np.flatnonzero(np.abs(nodes[:,axis]-coordinate)<=tolerance);record=dict(row,scope='明确无CAD数值块的物理平面，不是齿轮孔壁')
        else:
            from .generated_solid import GeneratedBlankSettings
            raw=workflow.get('base_case',{}).get('manufacturing',{}).get(member,{}).get('blank')
            actual=workflow.get('actual_manufacturing',{}).get(member,{}).get('blank')
            if raw is None:raise ValueError('缺原制造blank；fixture须显式给物理plane，不得反推free DOF')
            try:settings=GeneratedBlankSettings(**raw)
            except TypeError as error:raise ValueError('缺原物理孔径；fixture须显式给物理plane') from error
            if actual is not None and safe(settings)!=safe(GeneratedBlankSettings(**actual)):raise ValueError('原blank与actual规范孔壁设置不同')
            radius=settings.bore_diameter_m*.5;tolerance=max(settings.boolean_tolerance_m*2,radius*1e-8)
            selected=np.flatnonzero(np.abs(np.linalg.norm(nodes[:,:2],axis=1)-radius)<=tolerance)
            record=dict(kind='bore_wall',axis='local_Z',radius_m=radius,tolerance_m=tolerance,blank_settings=safe(settings))
        if len(selected)<3 or len(selected)==len(nodes):raise ValueError('原物理支承没有有效的非整体固定节点集')
        fixed[member]=selected;records[member]=dict(record,selected_node_count=len(selected),selected_node_identity=numeric_array_identity(selected,integer=True),derived_from_saved_free_dofs=False)
    return fixed,dict(schema='native_FE_physical_boundary_selection_v1',members=records,explicit_fixture_plane=boundary_conditions is not None)


def restore_native_phase_response(source,proof):
    """从本态原u/F恢复已独审响应；不加载typed operator或保存的Ku数组。"""
    import numpy as np
    from .native_results_archive import load_native_phase_source
    from .distributed_pressure_domain import distributed_body_identity,distributed_response_fingerprint
    fresh=load_native_phase_source(source['phase_root'],expected_state_sha256=source['state_sha256'],workflow_root=source['workflow_root'])
    if (proof.get('passed') is not True or proof.get('independent_Ku_recheck') is not True or proof.get('source_state_sha256')!=fresh['state_sha256'] or
            proof.get('current_operator_manifest_sha256')!=fresh['identity']['current_operator_manifest_sha256'] or proof.get('operator_fingerprint')!=fresh['identity']['operator_fingerprint']):raise ValueError('独审proof不属于本态原state/operator')
    result=deepcopy(fresh['result']);result['body_fields']={}
    for member,row in fresh['state']['response_receipt']['bodies'].items():
        path=relative(fresh['phase_root'],row['response_file'])
        if sha(path)!=row['response_sha256'] or proof['archived_files_identity'][member]['response_sha256']!=row['response_sha256']:raise ValueError('材料域u/F不是本态原独审响应')
        with np.load(path,allow_pickle=False) as archive:
            u=archive['displacements_m'];force=archive['nodal_forces_n']
        shape=proof['body_checks'][member]['mesh_nodes_identity']['shape']
        if list(u.shape)!=shape or force.shape!=u.shape or not np.isfinite(u).all() or not np.isfinite(force).all():raise ValueError('原u/F布局或有限值非法')
        result['body_fields'][member]=dict(displacements_m=u,nodal_forces_n=force,**{k:result['body_response_checks'][member][k] for k in ('strain_energy_j','external_work_j','maximum_free_residual_n','relative_free_residual')})
    result['distributed_body_identity']=distributed_body_identity(result)
    if distributed_response_fingerprint(result)!=proof['full_response_fingerprint']:raise ValueError('恢复原u/F/φ后的完整响应指纹不同')
    return dict(source=fresh,result=result,response_evidence=proof)


def _workflow_identity(source):
    return dict(source_summary_sha256=source['identity']['source_summary_sha256'],global_workflow_input_fingerprint=source['identity']['source_input_fingerprint'],
        states=[dict(phase_index=e['source_phase_index'],state_sha256=e['native_state_sha256']) for e in source['_entries']])


def _stage_record(directory,name,source):
    if (Path(directory)/'resource_failure.json').exists():raise ValueError('该结果目录已有资源失败记录，不能继续认证；请保留证据并使用新目录/分阶段预算')
    path=Path(directory)/('comparison_fe' if name=='comparison' else name)/'manifest.json';record,_=read(path)
    if record.get('schema')!='bevel_native_FE_postprocessing_stage_v1' or record.get('stage')!=name or record.get('execution_completed') is not True or record.get('workflow_identity')!=_workflow_identity(source):raise ValueError('后处理阶段未完成或不属于同workflow原状态：'+name)
    return record


def _phase_source(source,entry):
    return load_native_phase_from_cycle(source,next(i for i,e in enumerate(source['_entries']) if e['source_phase_index']==entry['source_phase_index']))


def audit_native_responses(source,output,*,boundary_conditions=None,progress=None):
    """逐态原Ku/Wf/能量，释放每态operator；不加载CAD或模型档案。"""
    folder=Path(output)/'response'
    if folder.exists():raise FileExistsError('不覆盖原response阶段证据')
    folder.mkdir(parents=True)
    manifest=dict(schema='bevel_native_FE_postprocessing_stage_v1',stage='response',workflow_identity=_workflow_identity(source),execution_completed=False,phases=[])
    try:
        first=load_native_phase_from_cycle(source,0);meshes,paths=load_native_workflow_meshes(first);fixed,boundaries=physical_fixed_nodes(first,meshes,boundary_conditions=boundary_conditions)
        manifest['physical_boundaries']=boundaries
        for entry in source['_entries']:
            phase=_phase_source(source,entry);audited=audit_native_phase_response(phase,meshes=meshes,mesh_files=paths,fixed_nodes=fixed);proof=audited['response_evidence']
            proof['physical_boundary_selection']=boundaries;proof['stiffness_assembled_by_workflow']=False
            record=_write(folder/f'phase_{phase["phase_index"]:04d}.json',proof)
            manifest['phases'].append(dict(phase_index=phase['phase_index'],source_state_sha256=phase['state_sha256'],proof=record,passed=proof['passed']))
            _write(folder/'manifest.json',manifest)
            if progress:progress(dict(stage='response',phase_index=phase['phase_index'],passed=proof['passed']))
            del audited,proof,phase;gc.collect()
        manifest['execution_completed']=True
    finally:_write(folder/'manifest.json',manifest)
    return manifest


def _response_proof(output,source,phase_index):
    manifest=_stage_record(output,'response',source);row=next((r for r in manifest['phases'] if r['phase_index']==phase_index),None)
    if row is None or row.get('passed') is not True:raise ValueError('缺本态已完成的独立响应')
    proof,_=read(relative(Path(output)/'response',row['proof']['file']),row['proof']['sha256'])
    if proof.get('source_state_sha256')!=row['source_state_sha256'] or proof.get('independent_Ku_recheck') is not True or proof.get('passed') is not True:raise ValueError('本态独审proof内容不完整')
    return proof


def audit_native_domains(source,output,*,base_depth=2,material_depth=3,maximum_regions=150000,progress=None):
    """各态重新建立实际材料域；只共享本mesh/CAD的静态locator。"""
    from .native_results_target import load_native_domain_context
    folder=Path(output)/'domain'
    if folder.exists():raise FileExistsError('不覆盖原domain阶段证据')
    _stage_record(output,'response',source)
    context=load_native_domain_context(source,output);first=load_native_phase_from_cycle(source,0);meshes,paths=load_native_workflow_meshes(first)
    folder.mkdir(parents=True);cache={};manifest=dict(schema='bevel_native_FE_postprocessing_stage_v1',stage='domain',workflow_identity=_workflow_identity(source),execution_completed=False,phases=[])
    try:
        for entry in source['_entries']:
            phase=_phase_source(source,entry);audited=restore_native_phase_response(phase,_response_proof(output,source,phase['phase_index']))
            record=audit_native_phase_domain(audited,meshes=meshes,mesh_files=paths,context=context,output_directory=folder/f'phase_{phase["phase_index"]:04d}',static_geometry_cache=cache,base_depth=base_depth,material_depth=material_depth,maximum_regions=maximum_regions,progress=progress)
            geometry=_stage_record(output,'geometry',source);geometry_folder=Path(output)/'geometry'
            bound=dict(record['manifest']['source_files']);bound[str(geometry_folder/'manifest.json')]=sha(geometry_folder/'manifest.json')
            bound[str(geometry_folder/geometry['target_file'])]=geometry['target_sha256'];bound[str(geometry_folder/geometry['whole_target_file'])]=geometry['whole_target_sha256']
            from . import native_results_target
            bound[str(Path(native_results_target.__file__).resolve())]=sha(native_results_target.__file__)
            record['manifest'].update(source_files=bound,source_files_unchanged=all(sha(p)==h for p,h in bound.items()))
            _write(record['file'],record['manifest']);record['sha256']=sha(record['file'])
            manifest['phases'].append(dict(phase_index=phase['phase_index'],source_state_sha256=phase['state_sha256'],domain=dict(file=str(Path(record['file']).relative_to(folder)),sha256=record['sha256']),
                full_material_domain_certified=record['manifest']['whole_working_flank_zero_pressure_certified']))
            _write(folder/'manifest.json',manifest);del audited,phase,record;gc.collect()
        manifest['execution_completed']=True
    finally:_write(folder/'manifest.json',manifest)
    return manifest


def _domain_record(output,source,phase_index):
    stage=_stage_record(output,'domain',source);row=next((r for r in stage['phases'] if r['phase_index']==phase_index),None)
    if row is None:raise ValueError('缺本态实际材料域')
    return dict(file=relative(Path(output)/'domain',row['domain']['file']),sha256=row['domain']['sha256'])


def native_cycle_window_coverage(source,output):
    """从各原态完整17齿证据聚合实际周期外邻齿证明，不手填passed。"""
    from .contact_cycle import validate_mesh_cycle_angles,cycle_state_sequence_fingerprint,validate_cycle_window_coverage
    from .native_results_identity import source_manufacturing_contract
    rows=[];requested=None;teeth=None
    for index in range(len(source['_entries'])):
        phase=load_native_phase_from_cycle(source,index);proof=_response_proof(output,source,phase['phase_index'])
        if proof.get('source_state_sha256')!=phase['state_sha256']:raise ValueError('周期外齿proof不是同一原状态')
        file,domain=_bound_domain(phase,proof,_domain_record(output,source,phase['phase_index']))
        _,manufacturing,_=source_manufacturing_contract(phase);count=manufacturing['pinion']['request']['tooth_count'];offsets=domain.get('nominal_requested_offsets',domain['requested_offsets']);represented=domain['requested_offsets']
        if teeth is not None and (teeth!=count or requested!=offsets):raise ValueError('周期齿数或请求窗口改变')
        teeth=count;requested=offsets
        if domain.get('whole_working_flank_zero_pressure_certified') is not True or domain.get('all_external_physical_teeth_verified') is not True:raise ValueError('周期需要每态完整实际材料域，未知不得补成蓝色')
        outer,_=read(relative(file.parent,domain['external_evidence']['file']),domain['external_evidence']['sha256'])
        if (outer.get('passed') is not True or outer.get('original_state_unchanged') is not True or outer.get('actual_geometry_checked') is not True or outer.get('response_W_Ku_bound') is not True or outer.get('full_response_fingerprint')!=proof['full_response_fingerprint']):raise ValueError('周期外齿证据未绑定本态真实u/φ')
        physical={(-k)%count for k in represented};external=outer['rows']
        if len(external)!=count-len(physical) or len({r['actual_pinion_tooth'] for r in external})!=len(external) or physical|{r['actual_pinion_tooth'] for r in external}!=set(range(count)):raise ValueError('外齿证据没有覆盖全部其余实际小轮齿')
        for item in external:
            if item.get('passed') is not True or item['actual_pinion_tooth']!=(-item['tooth_offset'])%count:raise ValueError('某个无原field物理齿尚未独立排除')
        outside=[min(requested)-1,max(requested)+1];pairs=[]
        for offset in outside:
            match=[r for r in external if r['tooth_offset']==offset]
            if len(match)!=1:raise ValueError('缺少实际窗口相邻外齿，不以其他物理齿替换')
            item=match[0]
            if item.get('passed') is not True or item['actual_pinion_tooth']!=(-offset)%count:raise ValueError('外邻齿实际物理索引或排除未通过')
            pairs.append({k:item[k] for k in ('tooth_offset','passed','minimum_remaining_gap_m','evidence_kind')})
        rows.append(dict(driving_angle_rad=phase['driving_angle_rad'],source_state_sha256=phase['state_sha256'],outer_pairs=pairs,source_domain_manifest_sha256=sha(file)))
    rows.sort(key=lambda r:r['driving_angle_rad']);angles=[r['driving_angle_rad'] for r in rows];states=[r['source_state_sha256'] for r in rows];validate_mesh_cycle_angles(angles,teeth)
    for index,row in enumerate(rows):row['phase_index']=index
    result=dict(schema='bevel_mesh_cycle_window_coverage_v1',method='full_gear_fe',requested_tooth_offsets=requested,outer_tooth_offsets=[min(requested)-1,max(requested)+1],
        state_sequence_sha256=cycle_state_sequence_fingerprint('full_gear_fe',angles,states,requested),records=rows,
        checks=dict(all_requested_states_accepted=True,outer_pair_search_covered=True,remaining_gap_nonnegative=True,source_states_unchanged=True),passed=True,
        interpretation='每态实际完整材料域及全部其余物理齿已核；此处只聚合原两邻齿余隙，不替代空间/时间收敛')
    return validate_cycle_window_coverage(result,method='full_gear_fe',angles_rad=angles,state_sha256=states,tooth_offsets=requested)


def export_native_instantaneous_maps(source,output,*,levels_pa=None,level_count=33,progress=None):
    """逐态逐承载齿导出原守恒P1，保留完整齿面、unknown和认证zero。"""
    import numpy as np
    import shapely
    from .refined_fe_plotting import _background
    folder=Path(output)/'maps'
    if folder.exists():raise FileExistsError('不覆盖原maps阶段证据')
    _stage_record(output,'domain',source)
    peak=max(r['maximum_pressure_pa'] for r in source['results'])
    if levels_pa is None:
        if isinstance(level_count,bool) or not isinstance(level_count,int) or not 2<=level_count<=257:raise ValueError('压力色阶数须为2至257的整数')
        levels_pa=np.linspace(0.,math.ceil(peak/1e6)*1e6,level_count)
    levels=np.asarray(levels_pa,dtype=float)
    if levels.ndim!=1 or len(levels)<2 or levels[0]!=0 or not np.isfinite(levels).all() or not np.all(np.diff(levels)>0):raise ValueError('绝对Pa色阶须从零严格递增')
    if levels[-1]<peak:raise ValueError('绝对压力色标未覆盖本次原峰值，不允许裁掉高压后作图')
    folder.mkdir(parents=True);manifest=dict(schema='bevel_native_FE_postprocessing_stage_v1',stage='maps',workflow_identity=_workflow_identity(source),execution_completed=False,levels_pa=levels.tolist(),phases=[])
    try:
        for entry in source['_entries']:
            phase=_phase_source(source,entry);proof=_response_proof(output,source,phase['phase_index']);audited=restore_native_phase_response(phase,proof);domain_record=_domain_record(output,source,phase['phase_index']);file,domain=_bound_domain(phase,proof,domain_record)
            background,_=read(relative(file.parent,domain['complete_source_working_flank_background']['file']),domain['complete_source_working_flank_background']['sha256'])
            _,polygons,_=_background(background,domain['pinion_request']);outline=shapely.union_all(shapely.polygons(np.asarray(polygons)))
            records=[]
            for field in audited['result']['fields']:
                if not np.any(np.asarray(field['forces_n'])>0):continue
                offset=field['tooth_offset'];target=folder/f'phase_{phase["phase_index"]:04d}_k_{offset:+d}.json.gz'
                data=export_native_pressure_source(audited,tooth_offset=offset,outline_geometry=outline,levels_pa=levels,domain_manifest=domain_record,output_file=target)
                records.append(dict(tooth_offset=offset,file=target.name,sha256=sha(target),original_P0_peak_pressure_pa=data['source_peak_pressure_pa'],unknown_area_m2=data['unknown_geometry'].area,certified_zero_area_m2=data['certified_zero_geometry'].area))
                del data
            if not records:raise ValueError('本态没有实际非零承载齿，不能生成虚假压力图')
            manifest['phases'].append(dict(phase_index=phase['phase_index'],source_state_sha256=phase['state_sha256'],maps=records));_write(folder/'manifest.json',manifest)
            if progress:progress(dict(stage='maps',phase_index=phase['phase_index'],nonzero_tooth_maps=len(records)))
            del audited,phase,proof;gc.collect()
        manifest['execution_completed']=True
    finally:_write(folder/'manifest.json',manifest)
    return manifest


def export_native_cycle_maximum(source,output):
    """完整闭周期原P1最大包络；缺相、缺非零齿或未知材料域均拒绝。"""
    from .contact_composite import load_triangle_pressure_map,export_triangle_pressure_map
    from .native_results_identity import source_manufacturing_contract
    folder=Path(output)/'cycle'
    if folder.exists():raise FileExistsError('不覆盖原cycle阶段证据')
    maps=_stage_record(output,'maps',source);coverage=native_cycle_window_coverage(source,output);instants=[]
    for phase in maps['phases']:
        for row in phase['maps']:
            path=relative(Path(output)/'maps',row['file'])
            if sha(path)!=row['sha256']:raise ValueError('周期源瞬时图原字节改变')
            instants.append(load_triangle_pressure_map(path,recompute=True))
    _,manufacturing,_=source_manufacturing_contract(load_native_phase_from_cycle(source,0))
    data=build_native_cycle_maximum(instants,pinion_tooth_count=manufacturing['pinion']['request']['tooth_count'],window_coverage=coverage)
    data['source_identity']['cycle_window_coverage_sha256']=fingerprint(coverage)
    folder.mkdir(parents=True);record=_write(folder/'window_coverage.json',coverage);target=folder/'full_gear_fe_sampled_cycle_maximum.json.gz';export_triangle_pressure_map(data,target)
    manifest=dict(schema='bevel_native_FE_postprocessing_stage_v1',stage='cycle',workflow_identity=_workflow_identity(source),execution_completed=True,window_coverage=record,
        pressure_map=dict(file=target.name,sha256=sha(target)),temporal_convergence_certified=False,spatial_convergence_certified=False,source_array_values_modified=False)
    _write(folder/'manifest.json',manifest);return manifest


def export_native_comparison_sources(source,output,*,phase_index=None,tooth_offset=0):
    """导出既有四图CLI所需FE两源；只选真实瞬时图，不制造半解析或周期数据。"""
    import shutil
    from .contact_composite import load_triangle_pressure_map
    folder=Path(output)/'comparison_fe'
    if folder.exists():raise FileExistsError('不覆盖已导出的四图FE源')
    maps=_stage_record(output,'maps',source);cycle=_stage_record(output,'cycle',source)
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('瞬时物理齿偏移须为整数')
    if phase_index is None:
        candidates=sorted((abs(e['driving_angle_rad']),e['source_phase_index']) for e in source['_entries'])
        if len(candidates)>1 and candidates[0][0]==candidates[1][0]:raise ValueError('最近零相位不唯一，请明确instantaneous phase index')
        phase_index=candidates[0][1]
    if isinstance(phase_index,bool) or not isinstance(phase_index,int):raise ValueError('瞬时相位索引须为整数')
    row=next((r for r in maps['phases'] if r['phase_index']==phase_index),None)
    selected=None if row is None else next((r for r in row['maps'] if r['tooth_offset']==tooth_offset),None)
    if selected is None:raise ValueError('所选瞬时相位/物理齿没有实际承载图源，请明确选择真实非零齿')
    original=[('instantaneous',relative(Path(output)/'maps',selected['file']),selected['sha256']),
        ('sampled_cycle_maximum',relative(Path(output)/'cycle',cycle['pressure_map']['file']),cycle['pressure_map']['sha256'])]
    sources=[]
    for kind,file,digest in original:
        if sha(file)!=digest:raise ValueError('四图FE源原字节改变')
        data=load_triangle_pressure_map(file,recompute=True)
        if data['method']!='full_gear_fe' or data['kind']!=kind:raise ValueError('四图FE源方法或图种不符')
        sources.append((kind,file,digest))
    folder.mkdir(parents=True);records=[]
    for kind,file,digest in sources:
        target=folder/f'full_gear_fe_{kind}_triangle_pressure.json.gz';shutil.copyfile(file,target)
        if sha(target)!=digest:raise IOError('四图FE源复制不完整')
        records.append(dict(kind=kind,file=target.name,sha256=digest))
    manifest=dict(schema='bevel_native_FE_postprocessing_stage_v1',stage='comparison',workflow_identity=_workflow_identity(source),execution_completed=True,
        instantaneous_phase_index=phase_index,instantaneous_tooth_offset=tooth_offset,files=records,halfspace_sources_generated=False,
        consumer='examples/replot_contact_composite.py --data <含真实两种方法四源的目录>',source_pressure_values_modified=False)
    _write(folder/'manifest.json',manifest);return manifest


def audit_native_fe_results(summary_path,output,*,stage='all',phase_indices=None,include_cycle=False,boundary_conditions=None,
        source_deflection_m=1e-6,target_deflection_m=1e-6,base_depth=2,material_depth=3,maximum_regions=150000,levels_pa=None,
        instantaneous_phase_index=None,instantaneous_tooth_offset=0,progress=None):
    """正式summary的参数化业务入口；阶段可单独执行，避免operator/CAD同驻。"""
    from .native_results_target import prepare_native_domain_geometry
    stages=('response','geometry','domain','maps','cycle','comparison')
    if stage not in ('all',*stages) or not isinstance(include_cycle,bool):raise ValueError('后处理阶段或周期开关非法')
    if include_cycle and stage not in ('all','cycle'):raise ValueError('include_cycle仅用于all或cycle阶段')
    source=load_native_cycle_source(summary_path,accepted_phase_indices=phase_indices,require_complete_cycle=stage in ('cycle','comparison') or include_cycle)
    output=Path(output).resolve()
    if (output/'resource_failure.json').exists():raise ValueError('结果目录已有资源失败记录，请保留该目录并使用新目录分阶段执行')
    output.mkdir(parents=True,exist_ok=True)
    requested=list(stages if include_cycle else stages[:4]) if stage=='all' else [stage];completed={}
    for name in requested:
        if name=='response':completed[name]=audit_native_responses(source,output,boundary_conditions=boundary_conditions,progress=progress)
        elif name=='geometry':completed[name]=prepare_native_domain_geometry(source,output,source_deflection_m=source_deflection_m,target_deflection_m=target_deflection_m,progress=progress)
        elif name=='domain':completed[name]=audit_native_domains(source,output,base_depth=base_depth,material_depth=material_depth,maximum_regions=maximum_regions,progress=progress)
        elif name=='maps':completed[name]=export_native_instantaneous_maps(source,output,levels_pa=levels_pa,progress=progress)
        elif name=='cycle':completed[name]=export_native_cycle_maximum(source,output)
        else:completed[name]=export_native_comparison_sources(source,output,phase_index=instantaneous_phase_index,tooth_offset=instantaneous_tooth_offset)
        gc.collect()
    return dict(schema='bevel_native_FE_results_postprocessing_v1',stages=completed,workflow_identity=_workflow_identity(source),stiffness_assembled=False,factorization_performed=False,
        real_cycle_maximum_generated='cycle' in completed,temporal_convergence_certified=False,spatial_convergence_certified=False)
