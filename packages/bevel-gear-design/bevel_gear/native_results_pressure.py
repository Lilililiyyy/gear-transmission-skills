"""原生FE结果的零K材料域审核与原P0压力图源。

接受态、独立响应、实际u/φ材料域、显示场分层绑定。无域时默认拒绝；
明确诊断模式仅保留未知背景，不能推断零压。周期使用原P1超水平集核。
"""
from dataclasses import asdict,is_dataclass
import gc
import gzip
import json
import math
from pathlib import Path

from .native_results_archive import sha,read,relative,fingerprint,load_native_phase_source


def _write(path,value):
    """流式、固定gzip时间戳；仅更新本次新输出的manifest，不回写原档案。"""
    import numpy as np
    import io
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    def convert(item):
        if isinstance(item,np.ndarray):return item.tolist()
        if isinstance(item,np.generic):return item.item()
        if is_dataclass(item):return asdict(item)
        raise TypeError('未支持的证据值：'+str(type(item)))
    with temporary.open('wb') as raw:
        if path.suffix=='.gz':
            with gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0) as packed:
                with io.TextIOWrapper(packed,encoding='utf8',newline='\n') as stream:
                    json.dump(value,stream,ensure_ascii=False,allow_nan=False,separators=(',',':'),default=convert);stream.write('\n')
        else:
            stream=io.TextIOWrapper(raw,encoding='utf8',newline='\n')
            json.dump(value,stream,ensure_ascii=False,allow_nan=False,separators=(',',':'),default=convert);stream.write('\n');stream.flush();stream.detach()
    temporary.replace(path);return dict(file=path.name,sha256=sha(path))


def _audited(value):
    from .distributed_pressure_domain import distributed_response_fingerprint
    from .native_fe_operator import safe
    source=value['source'];result=value['result'];proof=value['response_evidence']
    fresh=load_native_phase_source(source['phase_root'],expected_state_sha256=source['state_sha256'],workflow_root=source['workflow_root'])
    if proof.get('passed') is not True or proof.get('independent_Ku_recheck') is not True or proof.get('source_state_sha256')!=fresh['state_sha256']:
        raise ValueError('须同一次真实matrix-free Ku/Wf/能量通过的响应')
    if distributed_response_fingerprint(result)!=proof.get('full_response_fingerprint'):raise ValueError('实际u/F/φ/P0与响应指纹不同')
    archived={p['tooth_offset']:p.get('source_material_reference_frames') for p in fresh['result']['geometry']['patches']}
    current={p['tooth_offset']:p.get('source_material_reference_frames') for p in result['geometry']['patches']}
    if safe(archived)!=safe(current):raise ValueError('实际用于显示的材料Rk不同于原state SHA所绑定的Rk')
    return fresh,result,proof


def _pose_check(result,context):
    """限定当前生成齿轮CAD域；不能把材料参考Rk当整体机械姿态。"""
    import numpy as np
    actual=result.get('geometry_check',{}).get('loaded_whole_body_frames');expected=context._loaded_frames(result)
    if not isinstance(actual,dict):raise ValueError('缺少本态有限φ实际整体姿态证据')
    maximum=0.
    for member in ('pinion','wheel'):
        for suffix,shape in (('_rotation_rows',(3,3)),('_axis_world',(3,))):
            a=np.asarray(actual[member+suffix],float);b=np.asarray(expected[member+suffix],float)
            if a.shape!=shape or b.shape!=shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('材料域机械姿态布局非法')
            maximum=max(maximum,float(np.abs(a-b).max()))
    if maximum>2e-12:raise ValueError('原实际R0/φ与CAD域的生成姿态不同')
    for name in ('translation_m','pinion_translation_m','wheel_translation_m','pinion_translation_world_m','wheel_translation_world_m','pinion_axis_origin_world_m','wheel_axis_origin_world_m'):
        if name in actual and (np.asarray(actual[name]).shape!=(3,) or not np.isfinite(actual[name]).all() or np.abs(actual[name]).max()>1e-15):raise ValueError('本生成CAD域尚未支持非零平移或偏置轴')
    return dict(passed=True,maximum_rotation_or_axis_component_difference=maximum,independent_finite_phi_recheck_performed=False)


def audit_native_phase_domain(audited,*,meshes,context,output_directory,static_geometry_cache=None,
        mesh_files=None,cad_files=None,base_depth=2,material_depth=3,maximum_regions=150000,progress=None):
    """完整请求齿和其余实际小轮齿；每态重建DeformedFEMaterialAudit。

    context须由同制造/CAD的prepare_generated_halfspace和完整目标三角域创建；
    这是纯CAD准备，不能传prepare_generated_ltca或旧接受证书替代。
    """
    from .distributed_pressure_domain import (DistributedPressureDomainAudit,
        distributed_pressure_domain_background,distributed_loaded_pose)
    from .native_results_external import native_field_partition,audit_native_unrepresented_teeth
    from .deformed_fe_material import DeformedFEMaterialAudit
    from .native_fe_operator import safe
    from .generated_tca import generated_tca_frames
    import numpy as np
    if not isinstance(context,DistributedPressureDomainAudit):raise TypeError('须明确实际CAD的DistributedPressureDomainAudit上下文')
    if context.target is None:raise ValueError('实际变形材料域需要完整目标三角网及覆盖身份，不能以最近面片替代')
    if not 0<=base_depth<=material_depth<=6 or maximum_regions<10000:raise ValueError('域审核深度/区域预算非法')
    source,result,proof=_audited(audited);folder=Path(output_directory).resolve()
    if folder.exists():raise FileExistsError('不覆盖已有材料域证据')
    from .native_results_identity import verify_native_model_identity,validate_material_projection_frames
    model_identity=verify_native_model_identity(source,meshes=meshes,context=context,mesh_files=mesh_files,cad_files=cad_files)
    partition=native_field_partition(result)
    request=context.models['pinion'].request;teeth=request.tooth_count;offsets=partition['represented_field_offsets']
    if not offsets or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets) or len({(-k)%teeth for k in offsets})!=len(offsets):raise ValueError('请求齿为空或重复实际物理齿')
    physical={(-k)%teeth for k in offsets};external=[k for k in range(-(teeth//2),teeth-teeth//2) if (-k)%teeth not in physical]
    for patch in result['geometry']['patches']:
        reference=patch.get('source_material_reference_frames')
        expected=generated_tca_frames(context.prepared.tca,result['geometry']['driving_angle_rad'],result['geometry']['reference_follower_angle_rad'],tooth_offset=patch['tooth_offset'])
        validate_material_projection_frames(reference,expected)
    cache={} if static_geometry_cache is None else static_geometry_cache
    allowed={'original','closures','tetrahedral_boundary','target_bridge','source_trace'}
    if not isinstance(cache,dict) or any(not isinstance(v,dict) or set(v)-allowed for v in cache.values()):raise ValueError('只允许纯静态CAD/原FE几何缓存，不能携带其他态u或零压证书')
    pose=_pose_check(result,context);folder.mkdir(parents=True)
    proof_record=_write(folder/'response_audit.json',proof);background=_write(folder/'actual_flank_mesh.json.gz',context.flank_mesh)
    sources={str(p):d for p,d in source['_file_sha256'].items()};sources[str(folder/proof_record['file'])]=proof_record['sha256']
    sources.update(model_identity['source_files'])
    identity_module=Path(verify_native_model_identity.__code__.co_filename).resolve();sources[str(identity_module)]=sha(identity_module)
    import sys
    sources[str(Path(__file__).resolve())]=sha(__file__)
    for name in ('bevel_gear.deformed_fe_material','bevel_gear.distributed_pressure_domain','bevel_gear.triangle_surface_distance','bevel_gear.generated_tca'):
        module=sys.modules.get(name)
        if module is not None:sources[str(Path(module.__file__).resolve())]=sha(module.__file__)
    manifest=dict(schema='native_FE_accepted_material_domain_audit_v1',execution_completed=False,
        accepted_source_state_sha256=source['state_sha256'],response_audit_sha256=proof_record['sha256'],
        current_operator_manifest_sha256=proof['current_operator_manifest_sha256'],actual_loaded_pose=distributed_loaded_pose(result),
        pinion_request=safe(request),source_pinion_teeth=teeth,requested_offsets=offsets,external_offsets=external,
        complete_source_working_flank_background=background,generated_domain_mechanical_pose_check=pose,
        material_source_identity=model_identity,
        nominal_requested_offsets=partition['nominal_requested_offsets'],unrepresented_nominal_offsets=partition['unrepresented_nominal_offsets'],
        zero_pressure_field_offsets=partition['zero_pressure_field_offsets'],pressure_fields_synthesized=False,
        certificates=[],whole_working_flank_zero_pressure_certified=False,pressure_mesh_convergence_verified=False,
        phase_convergence_verified=False,stiffness_assembled=False,factorization_performed=False,
        actual_u_phi_recomputed=True,prior_candidate_material_decisions_reused=False,strict_continuous_CAD_proof=False)
    try:
        external_module=Path(audit_native_unrepresented_teeth.__code__.co_filename).resolve();sources[str(external_module)]=sha(external_module)
        outer=audit_native_unrepresented_teeth(context.prepared,result,response_evidence=proof,external_offsets=external)
        manifest.update(external_evidence=_write(folder/'external_all_other_pinion_teeth.json.gz',outer),all_external_physical_teeth_verified=outer['passed'])
        material=DeformedFEMaterialAudit(meshes,result,response_evidence=proof,target_cad_context=context,static_geometry_cache=cache)
        manifest.update(actual_deformed_material_binding=material.binding,static_geometry_cache_identity=material.static_geometry_key)
        for offset in offsets:
            material.source_distance_cache.clear();material.clipped_distance_cache.clear()
            parent=context.audit(result,tooth_offset=offset,response_evidence=proof,maximum_depth=base_depth,maximum_regions=maximum_regions,progress=progress)
            parent_record=_write(folder/f'geometry_domain_k_{offset:+d}.json.gz',parent);certificate=parent
            if parent['unresolved_region_count']:
                certificate=material.refine_certificate(parent,maximum_depth=material_depth,maximum_regions=maximum_regions,parent_sha256=parent_record['sha256'],progress=progress)
            record=_write(folder/f'actual_material_domain_k_{offset:+d}.json.gz',certificate)
            gate=distributed_pressure_domain_background(result,context.background,tooth_offset=offset,certificate=certificate,pinion_request=request)
            manifest['certificates'].append(dict(tooth_offset=offset,certificate=record,parent=parent_record,passed=certificate['passed'],
                unresolved_region_count=certificate['unresolved_region_count'],full_flank_zero_certified=gate['full_flank_zero_certified'],domain_state=gate['domain_state']))
            _write(folder/'manifest.json',manifest);del parent,certificate,gate;gc.collect()
        manifest['execution_completed']=True
        manifest['whole_working_flank_zero_pressure_certified']=bool(outer['passed'] and len(offsets)+len(external)==teeth and all(r['full_flank_zero_certified'] for r in manifest['certificates']))
    finally:
        manifest.update(source_files=sources,source_files_unchanged=all(sha(p)==d for p,d in sources.items()),actual_u_phi_recomputed=True)
        _write(folder/'manifest.json',manifest)
    return dict(file=folder/'manifest.json',sha256=sha(folder/'manifest.json'),manifest=manifest)


def _bound_domain(source,proof,record):
    """新域协议必须具有CAD/blank原文件同源闭环，旧不完整域不能升级认证。"""
    from .native_results_identity import source_manufacturing_contract
    file=Path(record['file']).resolve();domain,_=read(file,record['sha256'])
    if (domain.get('execution_completed') is not True or domain.get('source_files_unchanged') is not True or domain.get('actual_u_phi_recomputed') is not True or
            domain.get('prior_candidate_material_decisions_reused') is not False or domain.get('accepted_source_state_sha256')!=source['state_sha256'] or
            domain.get('actual_deformed_material_binding',{}).get('full_response_fingerprint')!=proof['full_response_fingerprint'] or
            domain.get('current_operator_manifest_sha256')!=proof['current_operator_manifest_sha256']):raise ValueError('域证据不属于本态实际u/φ/operator')
    _,manufacturing,cad=source_manufacturing_contract(source);identity=domain.get('material_source_identity')
    if not isinstance(identity,dict) or identity.get('passed') is not True:raise ValueError('域证据缺少原CAD/blank同源门')
    for member in ('pinion','wheel'):
        row=identity.get('members',{}).get(member,{})
        if row.get('manufacturing')!=manufacturing[member] or row.get('cad_sha256')!=cad[member] or row.get('mesh_sha256')!=source['mesh_sha256'][member]:raise ValueError('域的原CAD/blank/mesh与source不同')
    for name,digest in domain['source_files'].items():
        if sha(name)!=digest:raise ValueError('域证据绑定的原文件改变')
    return file,domain


def _project(field,request,frames,cartesian_material_chart):
    if cartesian_material_chart is None:
        from .pressure_field import project_pressure_field
        return project_pressure_field(field,request,frames)
    # 小TET或一般曲面数值控制可显式声明真实材料平面；不得暗改齿轮默认投影。
    import numpy as np
    chart=cartesian_material_chart;origin=np.asarray(chart['origin_local_m'],float);basis=np.asarray(chart['tangent_rows_local'],float);rotation=np.asarray(frames['pinion_rotation_rows'],float)
    if origin.shape!=(3,) or basis.shape!=(2,3) or rotation.shape!=(3,3) or not np.isfinite(rotation).all() or not np.allclose(rotation.T@rotation,np.eye(3),rtol=0,atol=1e-12) or np.linalg.det(rotation)<=0 or not np.isfinite(origin).all() or not np.isfinite(basis).all() or not np.allclose(basis@basis.T,np.eye(2),rtol=0,atol=1e-12):raise ValueError('显式材料平面和源Rk须为有限正交SI坐标')
    local=(np.asarray(field['points_world_m'])-np.asarray(frames.get('translation_m',[0.,0.,0.])))@rotation
    return dict(field,points_width_height_m=(local-origin)@basis.T,projection_definition='显式实际材料平面；数值接口控制，不称锥齿轮子午投影')


def export_native_pressure_source(audited,*,tooth_offset,outline_geometry,levels_pa,domain_manifest=None,
        pinion_request=None,allow_unknown_background=False,cartesian_material_chart=None,output_file=None):
    """原P0→守恒中心P1→材料Rk→四图可用三角压力源，不改原峰值/压力。

    正式蓝底必须有同实际u/φ完整域；无域仅可明确要求诊断未知背景。
    outline须与域的完整实际背景同一，不能手工给好看轮廓替代。
    """
    import numpy as np
    import shapely
    from .pressure_field import reconstruct_distributed_fe_pressure_field
    from .contact_composite import build_triangle_pressure_map,OriginalP1Contours,export_triangle_pressure_map
    from .distributed_pressure_domain import distributed_pressure_domain_background
    source,result,proof=_audited(audited);outline=outline_geometry
    if not isinstance(allow_unknown_background,bool):raise TypeError('未知背景开关须显式布尔值')
    if domain_manifest is None and not allow_unknown_background:raise ValueError('缺少本态实际材料域，不能外推认证零压背景')
    if cartesian_material_chart is not None and domain_manifest is not None:raise ValueError('Cartesian材料平面仅允许无CAD的unknown诊断，不能产生认证蓝底')
    if cartesian_material_chart is not None and pinion_request is not None:raise ValueError('Cartesian诊断不能混用齿轮pinion_request覆盖')
    bound_domain=None
    if domain_manifest is not None:
        bound_domain=_bound_domain(source,proof,domain_manifest)
        from .native_results_identity import bound_projection_request
        pinion_request=bound_projection_request(source,pinion_request,domain_request=bound_domain[1]['pinion_request'])
    elif cartesian_material_chart is None:
        from .native_results_identity import bound_projection_request
        pinion_request=bound_projection_request(source,pinion_request)
    raw=reconstruct_distributed_fe_pressure_field(result,tooth_offset=tooth_offset)
    original=result['fields'][result['selected_tooth_offsets'].index(tooth_offset)];grid_nodes=raw['grid_node_indices']
    if not np.array_equal(raw['pressure_pa'],np.asarray(original['pressure_pa'])[grid_nodes[:,0],grid_nodes[:,1]]):raise ArithmeticError('中心P1源值不等于原P0中心值')
    patch=next(p for p in result['geometry']['patches'] if p['tooth_offset']==tooth_offset)
    if 'source_material_reference_frames' not in patch:raise ValueError('必须使用该物理齿源材料Rk，不能用机械R0替换')
    field=_project(raw,pinion_request,patch['source_material_reference_frames'],cartesian_material_chart)
    field.update(phase_index=source['phase_index'],tooth_offset=tooth_offset,driver_angle_rad=source['driving_angle_rad'],source_state_sha256=source['state_sha256'])
    source_identity=dict(source_state_sha256=source['state_sha256'],current_operator_manifest_sha256=proof['current_operator_manifest_sha256'],
        operator_fingerprint=proof['operator_fingerprint'],actual_native_mesh_sha256=source['mesh_sha256'],full_response_fingerprint=proof['full_response_fingerprint'],
        geometry_source='FE_boundary',source_kind='full_gear_fe',instantaneous_phase_rad=source['driving_angle_rad'],
        global_workflow_input_fingerprint=source['identity']['global_workflow_input_fingerprint'],materials_fingerprint=fingerprint(source['materials']),
        requested_tooth_offsets=result['requested_tooth_offsets'],nonzero_tooth_offsets=[f['tooth_offset'] for f in result['fields'] if np.any(np.asarray(f['forces_n'])>0)],
        signed_pinion_torque_nm=result['geometry']['pinion_torque_nm'],original_pressure_values_modified=False,
        original_P0_contact_area_plane_m2=float(np.count_nonzero(original['contact_mask'])*original['cell_area_m2']),original_all_pairs_P0_contact_area_plane_m2=result['contact_area_plane_m2'],original_P0_peak_pressure_pa=raw['diagnostics']['source_peak_pressure_pa'],
        original_P0_plane_force_n=float(np.asarray(original['forces_n']).sum()),original_P0_contact_area_actual_first_m2=math.fsum(c['area_actual_first_m2'] for c in result['cells'] if c['tooth_offset']==tooth_offset),
        unknown_background_diagnostic=domain_manifest is None,pressure_mesh_convergence_verified=False,temporal_convergence_certified=False)
    import sys
    source_identity['postprocessing_source_sha256']={Path(__file__).name:sha(__file__)}
    for name in ('bevel_gear.pressure_field','bevel_gear.contact_composite'):
        module=sys.modules[name];source_identity['postprocessing_source_sha256'][name]=sha(module.__file__)
    unknown=outline;zero=shapely.GeometryCollection()
    if domain_manifest is not None:
        file,domain=bound_domain;root=file.parent
        from .refined_fe_plotting import _background
        background_doc,_=read(relative(root,domain['complete_source_working_flank_background']['file']),domain['complete_source_working_flank_background']['sha256'])
        _,background,_=_background(background_doc,domain['pinion_request'])
        actual_outline=shapely.union_all(shapely.polygons(np.asarray(background)))
        if not actual_outline.equals(outline):raise ValueError('显示轮廓不是同一完整真实材料背景')
        row=next(r for r in domain['certificates'] if r['tooth_offset']==tooth_offset)
        certificate,_=read(relative(root,row['certificate']['file']),row['certificate']['sha256'])
        gate=distributed_pressure_domain_background(result,background,tooth_offset=tooth_offset,certificate=certificate,pinion_request=domain['pinion_request'])
        if gate['full_flank_zero_certified'] or gate['domain_state']=='partitioned_pressure_domain':
            polygons=gate.get('_unresolved_polygons_m',[]);unknown=shapely.union_all(shapely.polygons(np.asarray(polygons))).intersection(outline) if len(polygons) else shapely.GeometryCollection()
            support=OriginalP1Contours(field).positive_support.intersection(outline);zero=outline.difference(support.union(unknown))
        source_identity.update(material_domain_manifest_sha256=domain_manifest['sha256'],global_all_pinion_working_flanks_certified=domain['whole_working_flank_zero_pressure_certified'],full_domain_certified=gate['full_flank_zero_certified'],pinion_tooth_count=domain['source_pinion_teeth'])
    data=build_triangle_pressure_map([field],method='full_gear_fe',kind='instantaneous',outline_geometry=outline,levels_pa=levels_pa,
        unknown_geometry=unknown,certified_zero_geometry=zero,source_identity=source_identity)
    if data['source_peak_pressure_pa']!=raw['diagnostics']['source_peak_pressure_pa']:raise ArithmeticError('P1源改变了原P0峰值')
    if output_file is not None:export_triangle_pressure_map(data,output_file)
    return data


def build_native_cycle_maximum(instants,*,pinion_tooth_count,window_coverage):
    """全部真实相位/非零齿映回同材料面，复用原P1超水平集最大包络。"""
    import numpy as np
    import shapely
    from .contact_cycle import validate_mesh_cycle_angles,validate_cycle_window_coverage
    from .contact_composite import build_triangle_pressure_map
    if not isinstance(instants,(list,tuple)) or not instants:raise ValueError('需要本态独立瞬时源列表')
    groups={};first=instants[0]
    for item in instants:
        if item.get('method')!='full_gear_fe' or item.get('kind')!='instantaneous':raise ValueError('需要真实FE瞬时图源')
        identity=item['source_identity']
        reference=first['source_identity']
        if identity.get('pinion_tooth_count')!=pinion_tooth_count or any(identity.get(k)!=reference.get(k) for k in ('global_workflow_input_fingerprint','materials_fingerprint','actual_native_mesh_sha256','signed_pinion_torque_nm')):raise ValueError('周期不能混入不同齿数/实际mesh/材料/载荷或pressure格业务')
        if identity.get('unknown_background_diagnostic') is not False or identity.get('global_all_pinion_working_flanks_certified') is not True:raise ValueError('周期必须有各态实际完整材料域，不接受未知背景诊断图冒充周期')
        if not first['outline_geometry'].equals(item['outline_geometry']) or not np.array_equal(first['levels_pa'],item['levels_pa']):raise ValueError('周期各源须同实际材料轮廓与绝对色标')
        angle=identity['instantaneous_phase_rad'];group=groups.setdefault(angle,dict(state_sha256=identity['source_state_sha256'],fields=set(),expected=set(identity['nonzero_tooth_offsets'])))
        if group['state_sha256']!=identity['source_state_sha256']:raise ValueError('同一真实相位不能混用压力状态')
        for field in item['source_fields']:
            if field['driver_angle_rad']!=angle or field['source_state_sha256']!=identity['source_state_sha256'] or field['tooth_offset'] in group['fields']:raise ValueError('周期存在重复或伪造相位/物理齿源')
            group['fields'].add(field['tooth_offset'])
    if any(not g['expected'] or g['fields']!=g['expected'] for g in groups.values()):raise ValueError('周期不能遗漏任一实际非零承载齿')
    angles=sorted(groups);validate_mesh_cycle_angles(angles,pinion_tooth_count)
    requested=first['source_identity']['requested_tooth_offsets']
    if any(item['source_identity']['requested_tooth_offsets']!=requested for item in instants):raise ValueError('周期实际请求窗口不一致')
    verified=validate_cycle_window_coverage(window_coverage,method='full_gear_fe',angles_rad=angles,state_sha256=[groups[a]['state_sha256'] for a in angles],tooth_offsets=requested)
    if verified.get('passed') is not True:raise ValueError('缺少各实际相位全部外齿覆盖')
    unknown=shapely.union_all([item['unknown_geometry'].union(item['unclassified_geometry']) for item in instants]);zero=first['certified_zero_geometry']
    for item in instants[1:]:zero=zero.intersection(item['certified_zero_geometry'])
    return build_triangle_pressure_map([field for item in instants for field in item['source_fields']],method='full_gear_fe',kind='sampled_cycle_maximum',
        outline_geometry=first['outline_geometry'],levels_pa=first['levels_pa'],unknown_geometry=unknown,certified_zero_geometry=zero,
        source_identity=dict(geometry_source='FE_boundary',source_states=[item['source_identity'] for item in instants],cycle_window_coverage_verified=True,
            original_numerical_arrays_modified=False,temporal_convergence_certified=False))
