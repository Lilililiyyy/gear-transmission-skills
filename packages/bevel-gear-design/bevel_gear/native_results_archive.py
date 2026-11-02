"""原生FE正式结果档案读取；phase、API与workflow身份分层核对。

读取正式phase档案及distributed工作流，不改变原计算数组。
phase/API/global指纹各自核对；响应物理验证须另调matrix-free审核。
"""
from collections.abc import Sequence
from copy import deepcopy
import gzip
import hashlib
import json
import math
from numbers import Integral
from pathlib import Path

MEMBERS=('pinion','wheel')


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf8')).hexdigest()


def _digest(value,name):
    if not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError(name+'须为明确的SHA256')
    return value


def relative(root,name):
    root=Path(root).resolve();part=Path(name);target=(root/part).resolve()
    if part.is_absolute() or not target.is_relative_to(root):raise ValueError('档案数据路径不能越出明确根目录')
    return target


def read(path,expected=None):
    """同一批原字节同时核SHA和解码，避免状态检查/读取竞态。"""
    path=Path(path);raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if expected is not None and digest!=_digest(expected,str(path)):raise ValueError('原字节SHA不符：'+str(path))
    if path.suffix=='.gz':raw=gzip.decompress(raw)
    return json.loads(raw.decode('utf-8-sig')),digest


def _materials(inputs):
    direct=inputs.get('materials');workflow=inputs.get('workflow_inputs',{}).get('base_case',{}).get('materials')
    def elastic(value):
        if not isinstance(value,dict) or set(value)!=set(MEMBERS):raise ValueError('正式相位缺少明确的两体材料')
        normalized={}
        for member,row in value.items():
            young=row.get('young_modulus_pa');nu=row.get('poisson_ratio')
            if (isinstance(young,bool) or not isinstance(young,(int,float)) or not math.isfinite(young) or young<=0 or
                    isinstance(nu,bool) or not isinstance(nu,(int,float)) or not math.isfinite(nu) or not -1<nu<.5):raise ValueError('原相位材料E/nu非法')
            normalized[member]=dict(young_modulus_pa=young,poisson_ratio=nu)
        return normalized
    # 强度/密度原文保留在inputs；线弹性核只接受明确ElasticMaterial字段。
    if direct is not None and workflow is not None and elastic(direct)!=elastic(workflow):raise ValueError('phase材料与workflow展开材料不同')
    return elastic(direct if direct is not None else workflow)


def load_native_phase_source(phase_root,*,expected_state_sha256,workflow_root=None):
    """先核已接受完成receipt，再读取明确原phase；不信任原残差作物理证明。"""
    root=Path(phase_root).resolve();expected=_digest(expected_state_sha256,'expected_state_sha256');files={}
    def document(name,digest=None):
        path=relative(root,name);value,actual=read(path,digest);files[path]=actual;return value
    def bind(name,digest):
        path=relative(root,name);actual=sha(path)
        if actual!=_digest(digest,name):raise ValueError('原数组SHA不符：'+str(path))
        files[path]=actual;return path
    final=document('final_state_receipt.json')
    if final.get('schema')!='native_FE_final_state_receipt_v1' or final.get('accepted') is not True or final.get('source_identity_unchanged') is not True or final.get('source_state_sha256')!=expected:
        raise ValueError('原phase未完成、被拒绝或不属于明确原state SHA')
    if final.get('state_file')!='state.json.gz':raise ValueError('未知正式phase状态文件布局')
    state=document(final['state_file'],expected);inputs=document('inputs.json');receipt=document(final['response_receipt_file'],final['response_receipt_sha256'])
    if state.get('schema')!='bevel_native_FE_phase_state_v1' or inputs.get('schema')!='bevel_native_FE_phase_inputs_v1':raise ValueError('未知正式原生FE相位协议，禁止回退旧NURBS路线')
    if state.get('accepted') is not True or state.get('source_identity_unchanged') is not True:raise ValueError('原相位未通过接受与来源稳定门控')
    before=state.get('source_identity_before')
    if not isinstance(before,dict) or not before or before!=state.get('source_identity_after') or fingerprint(before)!=final.get('source_identity_sha256'):raise ValueError('原相位源码前后身份不一致')
    if 'input_fingerprint' in inputs:raise ValueError('phase inputs不得嵌自身指纹')
    phase_fp=fingerprint(inputs)
    if state.get('input_fingerprint')!=phase_fp or final.get('input_fingerprint')!=phase_fp:raise ValueError('相位输入指纹不同')
    index=state.get('phase_index');angle=state.get('driving_angle_rad')
    if isinstance(index,bool) or not isinstance(index,int) or index<0 or isinstance(angle,bool) or not isinstance(angle,(int,float)) or not math.isfinite(angle):raise ValueError('原相位索引或rad角度非法')
    if inputs.get('phase_index')!=index or inputs.get('driving_angle_rad')!=angle or receipt.get('phase_index')!=index or receipt.get('input_fingerprint')!=phase_fp or receipt.get('passed') is not True:raise ValueError('原响应的相位或phase指纹不同')
    for key in ('native_contact_api_input_fingerprint','global_workflow_input_fingerprint'):
        value=state.get(key)
        if key=='native_contact_api_input_fingerprint' or value is not None:_digest(value,key)
        if value!=inputs.get(key):raise ValueError('phase/API/global指纹混用：'+key)
    if 'workflow_inputs' in inputs and fingerprint(inputs['workflow_inputs'])!=state.get('global_workflow_input_fingerprint'):
        raise ValueError('phase所存workflow正文指纹与global身份不同')
    result=state['result']
    if (result.get('geometry_source')!='FE_boundary' or result.get('geometry',{}).get('geometry_source')!='FE_boundary' or
            result.get('native_contact_accepted') is not True or result.get('converged') is not True or
            result.get('geometry_linearization_checked') is not True or result.get('equilibrium',{}).get('certificate',{}).get('passed') is not True):
        raise ValueError('原FE_boundary core缺少QP/几何/接受门控')
    if result.get('contact_representation')!='distributed_P0_surface_cells' or result.get('halfspace_compliance_used') is not False:raise ValueError('正式native core须是纯FE分布P0接触')
    from .native_results_geometry import validate_original_native_geometry_gate
    geometry_receipt=validate_original_native_geometry_gate(result,declared_controls=inputs.get('native_contact_solver',{}))
    if result.get('source_input_fingerprint')!=phase_fp or result.get('native_contact_api_input_fingerprint')!=state['native_contact_api_input_fingerprint'] or result['geometry'].get('driving_angle_rad')!=angle:raise ValueError('core与正式phase输入身份不同')
    meshes=result.get('actual_native_mesh_sha256')
    if not isinstance(meshes,dict) or set(meshes)!=set(MEMBERS) or inputs.get('actual_native_mesh_sha256')!=meshes:raise ValueError('两体原mesh身份不同')
    for m,digest in meshes.items():_digest(digest,m+' mesh')
    operator={name:result.get(name) for name in ('current_operator_manifest_file','current_operator_manifest_sha256','operator_fingerprint')}
    if inputs.get('final_operator_binding')!=operator or any(final.get(k)!=v for k,v in operator.items()):raise ValueError('phase最终operator绑定不同')
    manifest=document(operator['current_operator_manifest_file'],operator['current_operator_manifest_sha256'])
    if manifest.get('schema')!='native_FE_typed_operator_manifest_v2' or manifest.get('geometry_source')!='FE_boundary' or manifest.get('operator_fingerprint')!=operator['operator_fingerprint'] or manifest.get('actual_mesh_sha256')!=meshes:raise ValueError('原最终operator不是同mesh的typed FE协议')
    op_root=relative(root,operator['current_operator_manifest_file']).parent
    def operator_file(name,digest):
        target=relative(op_root,name)
        return bind(target.relative_to(root).as_posix(),digest)
    operator_file(manifest['geometry']['file'],manifest['geometry']['sha256'])
    for row in manifest['patches']:
        if row is not None:operator_file(row['arrays_file'],row['arrays_sha256']);operator_file(row['metadata']['file'],row['metadata']['sha256'])
    if set(manifest.get('matrices',{}))!=set(MEMBERS) or set(receipt.get('bodies',{}))!=set(MEMBERS):raise ValueError('原operator或响应缺少两体')
    for member in MEMBERS:
        matrix=manifest['matrices'][member];operator_file(matrix['file'],matrix['sha256']);row=receipt['bodies'][member]
        bind(row['response_file'],row['response_sha256']);bind(row['load_matrix_file'],row['load_matrix_sha256'])
    if state.get('response_receipt')!=receipt or result.get('original_response_receipt')!=receipt:raise ValueError('内嵌与外置原响应receipt不同')
    return dict(method='full_gear_fe',source_kind='full_gear_fe',geometry_source='FE_boundary',state=state,result=result,inputs=inputs,
        phase_root=root,workflow_root=None if workflow_root is None else Path(workflow_root).resolve(),state_file=relative(root,final['state_file']),
        state_sha256=expected,phase_index=index,driving_angle_rad=float(angle),materials=_materials(inputs),mesh_sha256=deepcopy(meshes),
        final_receipt=final,operator_manifest=manifest,operator_manifest_file=relative(root,operator['current_operator_manifest_file']),
        identity=dict(phase_input_fingerprint=phase_fp,native_contact_api_input_fingerprint=state['native_contact_api_input_fingerprint'],
            global_workflow_input_fingerprint=state.get('global_workflow_input_fingerprint'),source_state_sha256=expected,
            current_operator_manifest_sha256=operator['current_operator_manifest_sha256'],operator_fingerprint=operator['operator_fingerprint']),
        original_geometry_receipt=geometry_receipt,response_physics_independently_verified=False,_file_sha256=files,_required_files=sorted(files))


class ReadOnlyNativeResults(Sequence):
    """逐次验证原全局链接和phase，不保存可被改写后复用的压力数组。"""
    def __init__(self,entries):self.entries=deepcopy(entries)
    def __len__(self):return len(self.entries)
    def __getitem__(self,index):
        if isinstance(index,slice):return [self[i] for i in range(*index.indices(len(self)))]
        if isinstance(index,bool) or not isinstance(index,Integral):raise TypeError('结果索引须为整数')
        entry=self.entries[index]
        for file,digest in entry['external_file_sha256'].items():
            if sha(file)!=digest:raise ValueError('原workflow链接在延迟读取前改变')
        return load_native_phase_source(entry['phase_root'],expected_state_sha256=entry['native_state_sha256'],workflow_root=entry['workflow_root'])['result']


def load_native_phase_from_cycle(source,index):
    """高层自动取得正式phase原SHA；不让使用者混用global state SHA。"""
    if isinstance(index,bool) or not isinstance(index,Integral):raise TypeError('相位位置须为整数')
    entry=source['_entries'][index]
    for file,digest in entry['external_file_sha256'].items():
        if sha(file)!=digest:raise ValueError('原workflow链接已改变')
    return load_native_phase_source(entry['phase_root'],expected_state_sha256=entry['native_state_sha256'],workflow_root=entry['workflow_root'])


def _phase_selection(inputs):
    """只读规范化已知两种正式writer字段；不能发明effective后重算原fp。"""
    effective=inputs.get('effective',{});analysis=inputs.get('distributed_fe_case',{}).get('analysis',{})
    primary=('phases' in effective or 'selected' in effective)
    legacy=('driver_angles_rad' in analysis or 'selected_phase_indices' in inputs)
    if primary and not {'phases','selected'}<=set(effective):raise ValueError('effective相位字段不完整')
    if legacy and ('driver_angles_rad' not in analysis or 'selected_phase_indices' not in inputs):raise ValueError('正式distributed相位字段不完整')
    if not primary and not legacy:raise ValueError('未知原生FE workflow相位字段布局')
    if primary and legacy and (effective['phases']!=analysis['driver_angles_rad'] or effective['selected']!=inputs['selected_phase_indices']):raise ValueError('两种正式相位字段声明不一致')
    return (effective['phases'],effective['selected']) if primary else (analysis['driver_angles_rad'],inputs['selected_phase_indices'])


def load_native_cycle_source(summary_path,method=None,*,require_complete_cycle=False,accepted_phase_indices=None):
    """正式distributed summary→既有cycle源同形字典；绝不走旧重建K分支。"""
    if method not in (None,'full_gear_fe') or not isinstance(require_complete_cycle,bool):raise ValueError('原生FE源method或完整周期开关错误')
    path=Path(summary_path).resolve();summary,summary_sha=read(path);root=path.parent
    if summary.get('schema')!='bevel_distributed_fe_workflow_v1':raise ValueError('未知原生FE工作流汇总schema')
    if summary.get('geometry_source')!='FE_boundary':raise ValueError('distributed源须明确geometry_source=FE_boundary')
    before=summary.get('source_identity_before')
    if not isinstance(before,dict) or not before or before!=summary.get('source_identity_after') or summary.get('source_identity_unchanged') is not True:raise ValueError('原workflow计算源前后身份不稳')
    inputs_path=relative(root,summary['input_file']);inputs,input_sha=read(inputs_path)
    if fingerprint(inputs)!=summary.get('input_fingerprint'):raise ValueError('workflow输入指纹不同')
    phases,input_selected=_phase_selection(inputs);selected=summary.get('selected_phase_indices')
    if not isinstance(phases,list) or not phases or any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in phases):raise ValueError('原phase角须为非空有限rad列表')
    if not isinstance(selected,list) or not selected or any(isinstance(i,bool) or not isinstance(i,int) or i<0 or i>=len(phases) for i in selected) or len(set(selected))!=len(selected):raise ValueError('原selected相位列表非法')
    if input_selected!=selected:raise ValueError('workflow输入与summary的原selected不一致')
    chosen=selected if accepted_phase_indices is None else accepted_phase_indices
    if not chosen or len(set(chosen))!=len(chosen) or any(isinstance(i,bool) or not isinstance(i,int) or i not in selected for i in chosen):raise ValueError('显式接受子集须保持原索引')
    if accepted_phase_indices is None and summary.get('analysis_completed') is not True:raise ValueError('未完成workflow只能明确选择已接受子集')
    if require_complete_cycle and (accepted_phase_indices is not None or chosen!=list(range(len(phases)))):raise ValueError('完整周期不能由子集或缺失原相位冒充')
    rows=summary.get('states',[])
    if not isinstance(rows,list) or any(isinstance(row.get('phase_index'),bool) or not isinstance(row.get('phase_index'),int) for row in rows):raise ValueError('原summary状态索引非法')
    lookup={row['phase_index']:row for row in rows}
    if len(lookup)!=len(rows) or any(i not in lookup or lookup[i].get('accepted') is not True for i in chosen):raise ValueError('源summary包含缺失/重复/未接受的所选相位')
    files={path:summary_sha,inputs_path:input_sha};entries=[];materials=None;meshes=None
    for i in chosen:
        row=lookup[i];global_path=relative(root,row['state_file']);global_digest=_digest(row.get('state_sha256') or row.get('source_state_sha256'),'global state SHA');global_state,global_sha=read(global_path,global_digest)
        angle=phases[i]
        if (global_state.get('schema')!='bevel_distributed_fe_state_v1' or global_state.get('phase_index')!=i or global_state.get('accepted') is not True or
                global_state.get('source_identity_unchanged') is not True or global_state.get('source_identity_before')!=before or global_state.get('source_identity_after')!=before or
                global_state.get('input_fingerprint')!=summary['input_fingerprint'] or global_state.get('driving_angle_rad')!=angle):raise ValueError('全局原state索引/输入/来源/接受身份不符')
        archive=row.get('native_phase_archive')
        if not isinstance(archive,dict) or global_state.get('native_phase_archive')!=archive or archive.get('schema')!='native_FE_phase_archive_reference_v1' or archive.get('path_base')!='workflow_root':raise ValueError('global与phase的正式链接缺失或不一致')
        phase_root=relative(root,archive['run_root']);phase=load_native_phase_source(phase_root,expected_state_sha256=archive['state_sha256'],workflow_root=root)
        if phase['inputs'].get('workflow_inputs') is not None and phase['inputs']['workflow_inputs']!=inputs:raise ValueError('phase所存workflow正文与原global inputs不同')
        _materials(dict(materials=phase['materials'],workflow_inputs=dict(base_case=inputs['base_case'])))
        for key,sha_key in (('state_file','state_sha256'),('final_state_receipt_file','final_state_receipt_sha256'),('response_receipt_file','response_receipt_sha256'),('current_operator_manifest_file','current_operator_manifest_sha256')):
            target=relative(root,archive[key]);digest=phase['_file_sha256'].get(target)
            if digest!=archive.get(sha_key):raise ValueError('global正式phase外存引用SHA不符：'+key)
        if phase['phase_index']!=i or phase['driving_angle_rad']!=angle or phase['identity']['global_workflow_input_fingerprint']!=summary['input_fingerprint'] or archive['input_fingerprint']!=phase['identity']['phase_input_fingerprint']:raise ValueError('全局/phase/API指纹或索引混用')
        # global与phase封套的路径/三种fp布局不同；其余物理core须完全相同。
        packaging={'source_input_fingerprint','current_operator_manifest_file','native_contact_api_input_fingerprint','body_fields','original_response_receipt'}
        physical=lambda result:{k:v for k,v in result.items() if k not in packaging}
        if fingerprint(physical(global_state.get('result',{})))!=fingerprint(physical(phase['result'])):raise ValueError('global与正式phase的原压力/TE/几何物理core不同')
        if materials is not None and (materials!=phase['materials'] or meshes!=phase['mesh_sha256']):raise ValueError('各phase实际mesh或材料不一致')
        materials=phase['materials'];meshes=phase['mesh_sha256'];files.update(phase['_file_sha256']);files[global_path]=global_sha
        entries.append(dict(phase_root=phase_root,native_state_sha256=phase['state_sha256'],workflow_root=root,source_phase_index=i,angle_rad=float(angle),
            state_path=global_path,state_sha256=global_sha,summary_path=path,summary_sha256=summary_sha,
            input_path=inputs_path,input_sha256=input_sha,input_fingerprint=summary['input_fingerprint'],
            external_file_sha256={path:summary_sha,inputs_path:input_sha,global_path:global_sha}))
    mesh_files=summary.get('mesh_files')
    if not isinstance(mesh_files,dict) or set(mesh_files)!=set(MEMBERS):raise ValueError('正式workflow须提供两体原mesh相对路径')
    mesh_paths={m:relative(root,name) for m,name in mesh_files.items()}
    for m,file in mesh_paths.items():
        digest=sha(file)
        if digest!=meshes[m] or summary.get('mesh_sha256',{}).get(m)!=digest:raise ValueError('workflow实际mesh原SHA不符')
        files[file]=digest
    cycle=None;manufacturing=inputs['base_case']['manufacturing'];angles=[float(phases[i]) for i in chosen]
    if require_complete_cycle:
        from .contact_cycle import validate_mesh_cycle_angles
        cycle=validate_mesh_cycle_angles(angles,manufacturing['pinion']['request']['tooth_count'])
    physical_inputs=deepcopy(inputs);physical_inputs.get('effective',{}).pop('phases',None);physical_inputs.get('effective',{}).pop('selected',None)
    physical_inputs.pop('selected_phase_indices',None)
    for name in ('distributed_fe_case','refined_fe_case'):
        physical_inputs.get(name,{}).get('analysis',{}).pop('driver_angles_rad',None)
    physics=dict(method='full_gear_fe',geometry_source='FE_boundary',inputs=physical_inputs,materials=materials,mesh_sha256=meshes)
    return dict(summary=summary,inputs=inputs,manufacturing=manufacturing,sampling=inputs['base_case'].get('cad_sampling',{}),
        results=ReadOnlyNativeResults(entries),cycle=cycle,identity=dict(method='full_gear_fe',geometry_source='FE_boundary',source_summary_sha256=summary_sha,
            input_sha256=input_sha,source_input_fingerprint=summary['input_fingerprint'],source_analysis_completed=summary.get('analysis_completed'),
            accepted_phase_indices=None if accepted_phase_indices is None else list(chosen),states=[dict(phase_index=e['source_phase_index'],sha256=e['state_sha256'],native_phase_state_sha256=e['native_state_sha256']) for e in entries]),
        phase_angles_rad=angles,method='full_gear_fe',source_kind='full_gear_fe',geometry_source='FE_boundary',mesh_sha256=meshes,mesh_files=mesh_paths,materials=materials,
        physics_contract=physics,_entries=entries,_source_path=path,_required_files=sorted(files),_file_sha256=files,
        response_physics_independently_verified=False,cycle_convergence_certified=False)


def audit_native_phase_response(source,*,meshes,mesh_files,fixed_nodes):
    """从本态实际u/F/W和原tet独立重算Ku/能量；不读取缓存残差作为证明。"""
    from .native_fe_operator import load_native_operator
    from .native_pressure_domain import audit_native_response_archive
    for m in MEMBERS:
        if sha(mesh_files[m])!=source['mesh_sha256'][m]:raise ValueError('明确实际mesh文件不属于本态')
    fresh=load_native_phase_source(source['phase_root'],expected_state_sha256=source['state_sha256'],workflow_root=source['workflow_root'])
    system=load_native_operator(fresh['operator_manifest_file'].parent)
    result,proof=audit_native_response_archive(fresh['state_file'],run_root=fresh['phase_root'],native_system=system,
        meshes=meshes,materials=fresh['materials'],fixed_nodes=fixed_nodes,mesh_source_sha256=fresh['mesh_sha256'],recheck_body=True)
    if proof.get('passed') is not True or proof.get('independent_Ku_recheck') is not True:raise ArithmeticError('独立原Ku/Wf/能量未通过')
    # 后续材料域只需实际u/F与独立proof；不让最终typed大operator跨阶段常驻。
    return dict(source=fresh,result=result,response_evidence=proof)


def audit_native_phase_domain(*args,**kwargs):
    from .native_results_pressure import audit_native_phase_domain as function
    return function(*args,**kwargs)


def export_native_pressure_source(*args,**kwargs):
    from .native_results_pressure import export_native_pressure_source as function
    return function(*args,**kwargs)
