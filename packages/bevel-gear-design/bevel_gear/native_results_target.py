"""同workflow模型档案的完整CAD目标及保守裁剪；只三角化，不重做Boolean。"""
import gc
from pathlib import Path

from .native_results_archive import sha,read,relative,load_native_phase_from_cycle
from .native_results_pressure import _write,_pose_check
from .native_results_identity import source_manufacturing_contract


def _models(source):
    """只恢复同目录已绑定的原model/grid/BRep，禁止回退重新造实体。"""
    workflow=source['inputs'].get('workflow_inputs');root=source.get('workflow_root')
    source_manufacturing_contract(source)
    if not isinstance(workflow.get('generated_model_archives'),dict) or root is None:raise ValueError('完整材料域需要同workflow的generated_model_archives，不能重Boolean补档')
    try:from .generated_model_archive import load_native_model_archives
    except ImportError as error:raise ImportError('完整材料域缺少generated_model_archive恢复接口，请使用包含该模块的版本') from error
    return load_native_model_archives(workflow,root)


def _context(source,*,source_deflection_m,target_deflection_m,triangles=None,target_identity=None,box=None):
    from .generated_halfspace import prepare_generated_halfspace
    from .distributed_pressure_domain import DistributedPressureDomainAudit
    models=_models(source);surface=prepare_generated_halfspace(models['pinion'].request,models['wheel'].request,cad_models=models)
    workflow,_,cad=source_manufacturing_contract(source)
    return DistributedPressureDomainAudit(surface,source_deflection_m=source_deflection_m,target_deflection_m=target_deflection_m,
        opposing_triangles_m=triangles,opposing_complete_box_m=box,opposing_identity=target_identity,cad_identity=cad)


def crop_complete_triangles(triangles,low,high,*,chunk_size=100000):
    """逐块对全部原三角AABB测试；边界相交者全部保留，不采用近邻抽样。"""
    import numpy as np
    values=np.asarray(triangles);low=np.asarray(low,float);high=np.asarray(high,float)
    if values.ndim!=3 or values.shape[1:]!=(3,3) or not len(values) or low.shape!=(3,) or high.shape!=(3,) or not np.isfinite(low).all() or not np.isfinite(high).all() or not np.all(high>low):raise ValueError('完整三角或裁剪盒布局非法')
    if isinstance(chunk_size,bool) or not isinstance(chunk_size,int) or chunk_size<1:raise ValueError('AABB块大小须为正整数')
    selected=[]
    for first in range(0,len(values),chunk_size):
        block=values[first:first+chunk_size]
        if not np.isfinite(block).all():raise ValueError('完整目标含非有限三角形')
        keep=np.all(block.max(axis=1)>=low,axis=1)&np.all(block.min(axis=1)<=high,axis=1)
        selected.append(np.flatnonzero(keep)+first)
    indices=np.concatenate(selected)
    if not len(indices):raise ValueError('完整CAD目标裁剪为空，不能补造面片')
    return values[indices],dict(all_original_cell_aabbs_tested=True,original_triangle_count=len(values),retained_triangle_count=len(indices))


def prepare_native_domain_geometry(source,output,*,source_deflection_m=1e-6,target_deflection_m=1e-6,angular_deflection_rad=.03,progress=None):
    """单独几何阶段先保存目标；不与大mesh/operator同驻，不加载u/F数组。"""
    import numpy as np
    from .native_results_workflow import _workflow_identity,_response_proof
    from .generated_tca import generated_tca_frames
    from .cad_import import ImportedCadModel,triangulate_imported_model
    for value in (source_deflection_m,target_deflection_m,angular_deflection_rad):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or value<=0:raise ValueError('CAD三角化参数须为正有限SI值')
    folder=Path(output)/'geometry'
    if folder.exists():raise FileExistsError('不覆盖原geometry阶段证据')
    first=load_native_phase_from_cycle(source,0);context=_context(first,source_deflection_m=source_deflection_m,target_deflection_m=target_deflection_m)
    folder.mkdir(parents=True);manifest=dict(schema='bevel_native_FE_postprocessing_stage_v1',stage='geometry',workflow_identity=_workflow_identity(source),execution_completed=False,
        source_deflection_m=source_deflection_m,target_deflection_m=target_deflection_m,angular_deflection_rad=angular_deflection_rad,Boolean_modeling_performed=False,stiffness_assembled=False)
    try:
        low=np.full(3,np.inf);high=-low;source_points=context.triangles.reshape(-1,3);margin=context.tolerance+source_deflection_m+target_deflection_m
        for index in range(len(source['_entries'])):
            phase=load_native_phase_from_cycle(source,index);proof=_response_proof(output,source,phase['phase_index']);_pose_check(phase['result'],context)
            if proof['source_state_sha256']!=phase['state_sha256']:raise ValueError('目标裁剪没有绑定本态真实响应')
            pose=proof['actual_loaded_pose'];whole=generated_tca_frames(context.prepared.tca,pose['driving_angle_rad'],pose['loaded_follower_angle_rad'])
            wheel=np.asarray(whole['wheel_rotation_rows'])
            for offset in phase['result']['requested_tooth_offsets']:
                frames=generated_tca_frames(context.prepared.tca,pose['driving_angle_rad'],pose['loaded_follower_angle_rad'],tooth_offset=offset)
                target_local=(source_points@np.asarray(frames['pinion_rotation_rows']).T)@wheel
                low=np.minimum(low,target_local.min(axis=0));high=np.maximum(high,target_local.max(axis=0))
            margin=max(margin,sum(proof['body_checks'][m]['maximum_displacement_norm_m'] for m in ('pinion','wheel'))+context.tolerance+source_deflection_m+target_deflection_m)
            del phase,proof;gc.collect()
        low-=margin;high+=margin
        mesh=triangulate_imported_model(ImportedCadModel(context.models['wheel'].shape,context.models['wheel'].diagnostics),linear_deflection_m=target_deflection_m,angular_deflection_rad=angular_deflection_rad)
        if mesh['metadata'].get('closed_solid') is not True or mesh['metadata'].get('solid_count')!=1:raise ValueError('目标不是完整单闭CAD实体')
        vertices=np.asarray(mesh.pop('vertices_m'),dtype=float);indices=np.asarray(mesh.pop('triangles'),dtype=np.int64);del mesh;gc.collect()
        whole_triangles=vertices[indices];del vertices,indices;gc.collect()
        whole_file=folder/'whole_wheel_triangles.npy';np.save(whole_file,whole_triangles,allow_pickle=False)
        target,coverage=crop_complete_triangles(whole_triangles,low,high);del whole_triangles;gc.collect()
        target_file=folder/'target_triangles.npz';np.savez_compressed(target_file,triangles_m=target);del target;gc.collect()
        workflow,_,cad=source_manufacturing_contract(first);files={}
        for row in workflow['generated_model_archives'].values():
            for name in ('descriptor_file','gear_brep_file','blank_brep_file'):
                path=relative(first['workflow_root'],row[name]);files[str(path)]=sha(path)
        import sys
        for name in ('native_results_target','generated_model_archive','cad_import','generated_halfspace','distributed_pressure_domain','surface_workflow'):
            module=sys.modules.get(__package__+'.'+name)
            if module is not None:files[str(Path(module.__file__).resolve())]=sha(module.__file__)
        manifest.update(target_file=target_file.name,target_sha256=sha(target_file),whole_target_file=whole_file.name,whole_target_sha256=sha(whole_file),
            target_identity=dict(target_cad_sha256=cad['wheel'],original_whole_mesh_sha256=sha(whole_file),all_original_cell_aabbs_tested=True,whole_mesh_identity_definition='原完整轮体CAD三角坐标NPY文件字节SHA256'),
            complete_box_m=[low.tolist(),high.tolist()],coverage=coverage,background_fingerprint=context.background_fingerprint,
            source_files=files,source_files_unchanged=all(sha(p)==h for p,h in files.items()),
            strict_continuous_CAD_proof=False,geometry_cache_contains_prior_material_decisions=False)
        if progress:progress(dict(stage='geometry',retained_triangles=coverage['retained_triangle_count']))
        manifest['execution_completed']=True
    finally:_write(folder/'manifest.json',manifest)
    return manifest


def load_native_domain_context(source,output):
    """恢复本次纯CAD目标；每次domain仍重新计算实际u/φ材料分类。"""
    import numpy as np
    from .native_results_workflow import _stage_record
    folder=Path(output)/'geometry';record=_stage_record(output,'geometry',source)
    for path,digest in record['source_files'].items():
        if sha(path)!=digest:raise ValueError('原模型档案在材料域审核前改变')
    whole=relative(folder,record['whole_target_file']);target=relative(folder,record['target_file'])
    if sha(whole)!=record['whole_target_sha256'] or sha(target)!=record['target_sha256'] or record['coverage']['all_original_cell_aabbs_tested'] is not True:raise ValueError('本次完整目标/裁剪覆盖身份不符')
    first=load_native_phase_from_cycle(source,0)
    with np.load(target,allow_pickle=False) as archive:
        context=_context(first,source_deflection_m=record['source_deflection_m'],target_deflection_m=record['target_deflection_m'],triangles=archive['triangles_m'],target_identity=record['target_identity'],box=record['complete_box_m'])
    if context.background_fingerprint!=record['background_fingerprint']:raise ValueError('恢复的实际完整背景不同于同模型原几何阶段')
    return context
