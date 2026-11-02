"""原生FE原始响应存档：保存原W/u/F，不由彩色云图反推计算结果。"""
from __future__ import annotations
from pathlib import Path
import hashlib


def _canonical_inputs(inputs):
    """规范化实际逐相位输入；不把自身指纹放入被散列的文档。"""
    from dataclasses import asdict,is_dataclass
    import json
    if not isinstance(inputs,dict) or any(name in inputs for name in
            ('input_fingerprint','phase_input_fingerprint','phase_fingerprint')):
        raise ValueError('逐相位inputs须为字典且不能包含自身输入指纹')
    def default(value):
        if is_dataclass(value):return asdict(value)
        if hasattr(value,'tolist'):return value.tolist()
        raise TypeError('逐相位输入包含不能保存为JSON的对象')
    # allow_nan=False拒绝非有限物理输入，不能在保存时静默替换成null。
    raw=json.dumps(inputs,ensure_ascii=False,sort_keys=True,separators=(',',':'),
        allow_nan=False,default=default)
    return json.loads(raw),raw


def native_phase_input_fingerprint(inputs):
    """与正式只读auditor一致的逐相位规范JSON SHA256。"""
    _,raw=_canonical_inputs(inputs)
    return hashlib.sha256(raw.encode('utf8')).hexdigest()


def file_sha256(path):
    """流式绑定原字节，避免大网格/响应复制到内存。"""
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block:=stream.read(4*1024*1024):digest.update(block)
    return digest.hexdigest()


def save_native_fe_response(prepared,system,result,directory,*,phase_index,input_fingerprint):
    """原K独立复算Ku−F、Wf−F及能量，再原样保存两体NPZ/W。

    directory为单相位结果根；输出引用均相对于它。NPZ不含pickle，
    保留原整轮节点、自由度及(k,row,column)列序。固定自由度反力不被
    误当自由度残差。调用时须保留原prepared.systems及body_fields。
    passed仅指响应方程与功，不认证接触几何、网格或全齿面零压。
    """
    import numpy as np
    from scipy import sparse
    if isinstance(phase_index,bool) or not isinstance(phase_index,int) or phase_index<0:
        raise ValueError('原响应相位索引须为非负整数')
    if not isinstance(input_fingerprint,str) or len(input_fingerprint)!=64 or any(c not in '0123456789abcdef' for c in input_fingerprint):
        raise ValueError('原响应需要明确的SHA256输入指纹')
    folder=Path(directory).resolve();target=folder/'responses'/f'phase_{phase_index:04d}'
    if target.exists():raise FileExistsError('响应目录已有原数据，不能覆盖')
    forces=np.asarray(result['equilibrium']['forces_n'],dtype=float)
    records=system['records'];keys=np.asarray([(r['tooth_offset'],r['row'],r['column']) for r in records],dtype=np.int64)
    if forces.shape!=(len(records),) or not np.isfinite(forces).all() or np.any(forces<0):
        raise ValueError('原W列与有限非负平面力布局不同')
    bodies={};total_energy=0.;buffers=[]
    for member in ('pinion','wheel'):
        body=prepared.systems[member];data=result['body_fields'][member]
        u=np.asarray(data['displacements_m'],dtype=float).reshape(-1)
        load=np.asarray(data['nodal_forces_n'],dtype=float).reshape(-1)
        w=system['matrices'][member].tocsc();free=np.asarray(body.free_dofs,dtype=np.int64)
        if (u.shape!=load.shape or len(u)!=body.stiffness.shape[0] or w.shape!=(len(u),len(forces)) or
                not np.isfinite(u).all() or not np.isfinite(load).all() or not len(free) or
                not np.isfinite(w.data).all() or not np.isfinite(body.stiffness.data).all()):
            raise ValueError('原FE位移、节点力、W或自由度布局错误')
        mapped=np.asarray(w@forces).ravel();ku=np.asarray(body.stiffness@u).ravel();residual=ku-load
        relative=float(np.linalg.norm(residual[free])/max(np.linalg.norm(load[free]),1.))
        force_error=float(np.linalg.norm(mapped-load)/max(np.linalg.norm(load),1.))
        energy=.5*float(u@ku);work=float(u@load)
        energy_error=abs(2*energy-work)/max(abs(work),1e-30)
        if (not np.isfinite([relative,force_error,energy,work,energy_error]).all() or energy<0 or
                relative>1e-8 or force_error>1e-12 or energy_error>1e-8):
            raise ArithmeticError('原FE的Ku−F、Wf−F或功/能量未通过')
        bodies[member]=dict(maximum_free_Ku_minus_F_n=float(np.max(np.abs(residual[free]))),
            relative_free_residual=relative,original_Wf_relative_error=force_error,
            strain_energy_j=energy,external_work_j=work,energy_work_relative_error=energy_error)
        buffers.append((member,u,load,residual,free,w));total_energy+=energy
    condensed=float(result['equilibrium']['strain_energy_j'])
    error=abs(total_energy-condensed)/max(abs(condensed),1e-30)
    if not np.isfinite(condensed) or error>1e-8:
        raise ArithmeticError('两体原FE应变能与原凝聚接触能不一致')
    target.mkdir(parents=True)
    for member,u,load,residual,free,w in buffers:
        response=target/(member+'_response.npz');temporary=response.with_name(response.name+'.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream,displacements_m=u.reshape(-1,3),nodal_forces_n=load.reshape(-1,3),
                ku_minus_f_n=residual.reshape(-1,3),free_dofs=free,force_plane_n=forces,contact_cell_keys=keys)
        temporary.replace(response)
        matrix=target/(member+'_load_matrix.npz');temporary=matrix.with_name(matrix.name+'.tmp')
        with temporary.open('wb') as stream:sparse.save_npz(stream,w,compressed=True)
        temporary.replace(matrix)
        bodies[member].update(response_file=response.relative_to(folder).as_posix(),response_sha256=file_sha256(response),
            load_matrix_file=matrix.relative_to(folder).as_posix(),load_matrix_sha256=file_sha256(matrix),load_matrix_shape=list(w.shape))
    return dict(schema='bevel_distributed_FE_response_receipt_v1',passed=True,phase_index=phase_index,
        input_fingerprint=input_fingerprint,body_fields_saved_before_removal=True,pickle_used=False,bodies=bodies,
        total_FE_strain_energy_j=total_energy,condensed_strain_energy_j=condensed,condensed_energy_relative_error=error,
        units=dict(displacement='m',force='N',stiffness='N/m',energy='J',load_matrix='N/N'))


def save_native_fe_phase_state(directory,*,inputs,result,response_receipt,phase_index,
        driving_angle_rad,source_identity_before,source_identity_after,
        native_contact_api_input_fingerprint,global_workflow_input_fingerprint=None,
        include_body_fields=False,operator_manifest_file='operator/manifest.json'):
    """写出已有accepted相位的薄档案封套，不改变压力、u/F或数学接受。

    调用方先保存最终operator和save_native_fe_response返回的原响应，
    后者必须以本inputs规范指纹和directory为根。此函数只核对文件与
    输入身份，写inputs/response_receipt/state/final_receipt；物理独立
    审核仍由audit_native_response_archive执行。状态/manifest/响应路径
    全部相对本相位根，不使用全局workflow路径或outer0代替最终operator。
    """
    import gzip
    import json
    import math
    from .refined_fe_postprocessing import _relative_file
    from .native_fe_operator import safe
    folder=Path(directory).resolve()
    normalized,canonical=_canonical_inputs(inputs)
    fingerprint=hashlib.sha256(canonical.encode('utf8')).hexdigest()
    def checked_sha(value,name):
        if (not isinstance(value,str) or len(value)!=64 or
                any(c not in '0123456789abcdef' for c in value)):
            raise ValueError(name+'须为明确SHA256')
    checked_sha(native_contact_api_input_fingerprint,'native_contact_api_input_fingerprint')
    if global_workflow_input_fingerprint is not None:
        checked_sha(global_workflow_input_fingerprint,'global_workflow_input_fingerprint')
    if isinstance(phase_index,bool) or not isinstance(phase_index,int) or phase_index<0:
        raise ValueError('逐相位档案索引须为非负整数')
    if (isinstance(driving_angle_rad,bool) or not isinstance(driving_angle_rad,(int,float)) or
            not math.isfinite(driving_angle_rad) or not isinstance(include_body_fields,bool)):
        raise ValueError('主动角须为有限rad，内嵌位移开关须为布尔值')
    if (isinstance(normalized.get('phase_index'),bool) or not isinstance(normalized.get('phase_index'),int) or
            isinstance(normalized.get('driving_angle_rad'),bool) or
            normalized.get('phase_index')!=phase_index or normalized.get('driving_angle_rad')!=driving_angle_rad or
            normalized.get('native_contact_api_input_fingerprint')!=native_contact_api_input_fingerprint or
            normalized.get('global_workflow_input_fingerprint')!=global_workflow_input_fingerprint):
        raise ValueError('逐相位输入的索引/角度/API身份/全局身份与显式档案参数不一致')
    if (not isinstance(source_identity_before,dict) or not source_identity_before or
            source_identity_before!=source_identity_after):
        raise ValueError('正式逐相位档案须先通过计算源码前后身份稳定检查')
    if (result.get('native_contact_accepted') is not True or result.get('converged') is not True or result.get('geometry_linearization_checked') is not True or
            result.get('equilibrium',{}).get('certificate',{}).get('passed') is not True or
            result.get('geometry_source')!='FE_boundary'):
        raise ValueError('只能保存已通过原native几何及QP接受的物理core')
    if (response_receipt.get('passed') is not True or response_receipt.get('phase_index')!=phase_index or
            response_receipt.get('input_fingerprint')!=fingerprint):
        raise ValueError('原W/u/F响应receipt未接受或与本相位输入/索引不同')
    if set(response_receipt.get('bodies',{}))!={'pinion','wheel'}:
        raise ValueError('原响应receipt必须包含两体原档案')
    mesh_sha=result.get('actual_native_mesh_sha256')
    if not isinstance(mesh_sha,dict) or set(mesh_sha)!={'pinion','wheel'}:
        raise ValueError('正式逐相位档案须绑定两体实际mesh文件SHA')
    for member,digest in mesh_sha.items():checked_sha(digest,member+'实际mesh SHA')
    if normalized.get('actual_native_mesh_sha256')!=mesh_sha:
        raise ValueError('逐相位实际输入与物理core的mesh字节SHA不同')
    manifest=_relative_file(folder,operator_manifest_file)
    manifest_sha=file_sha256(manifest)
    operator_fp=result.get('operator_fingerprint')
    checked_sha(operator_fp,'operator_fingerprint')
    if result.get('current_operator_manifest_sha256')!=manifest_sha:
        raise ValueError('物理core没有绑定当前最终operator的原manifest SHA')
    document=json.loads(manifest.read_text(encoding='utf8'))
    if document.get('operator_fingerprint')!=operator_fp or document.get('actual_mesh_sha256')!=mesh_sha:
        raise ValueError('当前最终operator的指纹或实际mesh SHA与物理core不同')
    binding=dict(current_operator_manifest_file=operator_manifest_file,
        current_operator_manifest_sha256=manifest_sha,operator_fingerprint=operator_fp)
    if normalized.get('final_operator_binding')!=binding:
        raise ValueError('逐相位输入未绑定同一最终operator路径/SHA/指纹')
    for member,row in response_receipt['bodies'].items():
        for file_key,sha_key in (('response_file','response_sha256'),('load_matrix_file','load_matrix_sha256')):
            file=_relative_file(folder,row[file_key])
            if file_sha256(file)!=row[sha_key]:
                raise ValueError('原响应或W文件SHA不符：'+member)
    names=('inputs.json','response_receipt.json','state.json.gz','final_state_receipt.json')
    protected=[folder/name for name in names]+[folder/(name+'.tmp') for name in names]
    if any(path.exists() or path.is_symlink() for path in protected):
        raise FileExistsError('逐相位正式或未完成档案已存在，不能覆盖')
    phase_result=dict(result)
    phase_result.update(binding,source_input_fingerprint=fingerprint,
        native_contact_api_input_fingerprint=native_contact_api_input_fingerprint,
        original_response_receipt=response_receipt)
    if not include_body_fields:phase_result.pop('body_fields',None)
    state=dict(schema='bevel_native_FE_phase_state_v1',phase_index=phase_index,
        driving_angle_rad=float(driving_angle_rad),accepted=True,input_fingerprint=fingerprint,
        source_identity_before=source_identity_before,source_identity_after=source_identity_after,
        source_identity_unchanged=True,response_receipt=response_receipt,result=phase_result,
        native_contact_api_input_fingerprint=native_contact_api_input_fingerprint,
        global_workflow_input_fingerprint=global_workflow_input_fingerprint,
        release_ready=False,mesh_or_cycle_convergence_certified=False)
    folder.mkdir(parents=True,exist_ok=True)
    def write(name,value,*,packed=False):
        # 最终receipt最后写入；部分写出失败不能形成可审核的完整档案。
        target=folder/name;temporary=target.with_name(target.name+'.tmp')
        raw=(json.dumps(safe(value),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode('utf8')
        with temporary.open('xb') as stream:
            if packed:
                with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as zipped:zipped.write(raw)
            else:stream.write(raw)
        temporary.replace(target)
    write('inputs.json',normalized);write('response_receipt.json',response_receipt)
    write('state.json.gz',state,packed=True)
    final=dict(schema='native_FE_final_state_receipt_v1',state_file='state.json.gz',
        source_state_sha256=file_sha256(folder/'state.json.gz'),input_fingerprint=fingerprint,
        accepted=True,source_identity_unchanged=True,
        source_identity_sha256=hashlib.sha256(json.dumps(source_identity_before,ensure_ascii=False,
            sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf8')).hexdigest(),
        response_receipt_file='response_receipt.json',response_receipt_sha256=file_sha256(folder/'response_receipt.json'),
        **binding)
    write('final_state_receipt.json',final)
    return dict(schema='native_FE_phase_archive_reference_v1',input_fingerprint=fingerprint,
        state_file='state.json.gz',state_sha256=final['source_state_sha256'],
        final_state_receipt_file='final_state_receipt.json',final_state_receipt_sha256=file_sha256(folder/'final_state_receipt.json'),
        response_receipt_file='response_receipt.json',response_receipt_sha256=final['response_receipt_sha256'],
        **binding)
