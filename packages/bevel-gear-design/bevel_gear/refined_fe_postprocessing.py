"""已接受节点FE状态的独立绘图接线；导入时不加载数值或CAD依赖。

原状态/物理场只读。第一次以同请求实际BRep生成显示网格；后续可仅
读取绑定制造请求与SHA256的网格缓存复绘，无需再次构造CAD或运行FE。
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import shutil


def _read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _relative_file(root,reference):
    if not isinstance(reference,str) or not reference or reference.startswith(('/','\\')) or re.match(r'^[A-Za-z]:',reference):
        raise ValueError('结果引用须为本结果目录内的相对文件路径')
    root=Path(root).resolve();path=(root/reference).resolve()
    if not path.is_relative_to(root):raise ValueError('结果引用越出本结果目录')
    return path


def _validate_models(models,manufacturing):
    """只核对同请求/同齿坯对象，不用任意CAD替代原制造输入。"""
    from .generated_surface import _request
    from .generated_solid import GeneratedBlankSettings
    if not isinstance(models,dict) or set(models)!={'pinion','wheel'}:raise ValueError('cad_models须显式给出pinion/wheel两轮对象')
    for member,model in models.items():
        if model.request!=_request(manufacturing[member]['request']) or model.settings!=GeneratedBlankSettings(**manufacturing[member]['blank']):
            raise ValueError(member+' CAD对象与实际制造请求/齿坯不一致')


def _build_models(manufacturing,sampling):
    from .generated_solid import build_generated_gear
    return {m:build_generated_gear(item['request'],item['blank'],**sampling) for m,item in manufacturing.items()}


def _triangulate_models(models,*,flank_linear_deflection_m=2e-5):
    """同一实际BRep输出整小轮及slot0真实裁切左右工作面；不使用理论矩形。"""
    from .cad_import import ImportedCadModel,triangulate_imported_model
    from .generated_halfspace import _faces
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    model=models['pinion'];compound=TopoDS_Compound();builder=BRep_Builder();builder.MakeCompound(compound)
    for side in ('left_flank','right_flank'):
        for face in _faces(model.shape,model.slot_face_groups[0][side]):builder.Add(compound,face)
    if isinstance(flank_linear_deflection_m,bool) or not isinstance(flank_linear_deflection_m,(int,float)) or not math.isfinite(flank_linear_deflection_m) or flank_linear_deflection_m<=0:
        raise ValueError('工作面显示三角化弦差须为正有限m')
    kwargs=dict(linear_deflection_m=2e-5,angular_deflection_rad=.08)
    result = dict(pinion=triangulate_imported_model(ImportedCadModel(model.shape,{}),**kwargs),
        pinion_flank=triangulate_imported_model(ImportedCadModel(compound,{}),
            linear_deflection_m=flank_linear_deflection_m,angular_deflection_rad=.08))
    patch = model.grid['patches']['left_flank']
    result['pinion_flank']['metadata']['source_surface_sampling'] = dict(
        width_count=patch['width_count'],profile_count=patch['profile_count'])
    return result


def _validate_mesh(mesh,*,flank=False):
    if not isinstance(mesh,dict) or mesh.get('units')!='m' or not mesh.get('vertices_m') or not mesh.get('triangles'):
        raise ValueError('实际CAD缓存须含明确m单位、真实顶点和三角面')
    vertices=mesh['vertices_m'];triangles=mesh['triangles']
    for point in vertices:
        if len(point)!=3 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in point):raise ValueError('实际CAD缓存顶点须为有限三维m坐标')
    for triangle in triangles:
        if len(triangle)!=3 or any(isinstance(i,bool) or not isinstance(i,int) or not 0<=i<len(vertices) for i in triangle):raise ValueError('实际CAD缓存三角面索引越界')
    if flank and len(mesh.get('face_indices',[]))!=len(triangles):raise ValueError('实际裁切齿面缓存缺少BRep面身份')
    if not flank and mesh.get('metadata',{}).get('closed_solid') is not True:raise ValueError('整轮显示缓存须来自实际闭实体CAD')


def _load_or_build_meshes(output,default_output,manufacturing,sampling,cad_models,*,flank_linear_deflection_m=None):
    from .refined_fe_workflow import _fingerprint,_write
    identity=_fingerprint(dict(manufacturing=manufacturing,cad_sampling=sampling))
    directory=output/'cad';manifest_name='cad_mesh_manifest.json'
    if cad_models is not None:_validate_models(cad_models,manufacturing)
    for source in dict.fromkeys((directory,default_output/'cad')):
        manifest_path=source/manifest_name
        if not manifest_path.is_file():continue
        manifest=_read(manifest_path)
        if manifest.get('schema')!='bevel_refined_fe_cad_cache_v1' or manifest.get('geometry_fingerprint')!=identity:
            raise ValueError('已有CAD显示缓存与本次制造请求/齿坯不一致，拒绝混用')
        meshes={}
        if set(manifest.get('files',{}))!={'pinion','pinion_flank'}:raise ValueError('已有CAD显示缓存不完整')
        for member,item in manifest['files'].items():
            path=_relative_file(source,item['file'])
            if _sha(path)!=item['sha256']:raise ValueError('实际CAD显示缓存SHA256变化，拒绝静默替换')
            meshes[member]=_read(path);_validate_mesh(meshes[member],flank=member=='pinion_flank')
        if flank_linear_deflection_m is not None:
            actual=meshes['pinion_flank'].get('metadata',{}).get('linear_deflection_m')
            if not isinstance(actual,(int,float)) or actual>flank_linear_deflection_m*(1+1e-12):
                # 不把粗显示网格误称更细的几何证据；由同一实际CAD重新三角化。
                continue
        directory.mkdir(parents=True,exist_ok=True)
        if source.resolve()!=directory.resolve():
            for item in manifest['files'].values():
                target=directory/item['file'];target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(_relative_file(source,item['file']),target)
            _write(manifest,directory/manifest_name)
        return meshes,manifest,'bound_mesh_cache'
    models=cad_models if cad_models is not None else _build_models(manufacturing,sampling)
    _validate_models(models,manufacturing)
    meshes=_triangulate_models(models,**({} if flank_linear_deflection_m is None else dict(flank_linear_deflection_m=flank_linear_deflection_m)))
    manifest=dict(schema='bevel_refined_fe_cad_cache_v1',geometry_fingerprint=identity,units='m',
        binding='实际同request/blank对象的BRep三角化；缓存以几何输入指纹及文件SHA256绑定',files={})
    directory.mkdir(parents=True,exist_ok=True)
    for member,mesh in meshes.items():
        _validate_mesh(mesh,flank=member=='pinion_flank');path=directory/(member+'_mesh.json');_write(mesh,path)
        manifest['files'][member]=dict(file=path.name,sha256=_sha(path))
    _write(manifest,directory/manifest_name)
    return meshes,manifest,'explicit_same_request_models' if cad_models is not None else 'built_from_effective_inputs'


def _accepted_states(summary_path):
    """校验状态/输入绑定，只读已接受结果；失败相位不补零或绘制占位图。"""
    from .refined_fe_workflow import _fingerprint,_accepted
    summary=_read(summary_path);root=summary_path.parent
    if summary.get('schema')!='bevel_refined_fe_workflow_v1':raise ValueError('需要实际refined FE workflow汇总')
    if summary.get('source_identity_unchanged') is not True:raise ValueError('原业务源码身份未通过，不绘制已接受物理结果')
    input_path=_relative_file(root,summary['input_file']);inputs=_read(input_path)
    effective=inputs['effective'];base=inputs['base_case']
    fingerprint=_fingerprint(dict(base_case=base,**{k:v for k,v in effective.items() if k!='selected'}))
    if inputs.get('input_fingerprint')!=fingerprint or summary.get('input_fingerprint')!=fingerprint:raise ValueError('汇总与实际有效输入指纹不一致')
    accepted=[];excluded=[];used=set()
    for row in summary.get('states',[]):
        if row.get('accepted') is not True:
            excluded.append(dict(phase_index=row.get('phase_index'),status=row.get('status'),reason='原状态未接受，未绘图'));continue
        path=_relative_file(root,row['state_file']);state=_read(path);index=row['phase_index']
        if isinstance(index,bool) or not isinstance(index,int) or index<0 or index in used or index>=len(effective['phases']):raise ValueError('已接受相位索引重复或越界')
        used.add(index)
        result=state.get('result',{})
        if state.get('schema')!='bevel_refined_fe_state_v1' or state.get('accepted') is not True or state.get('source_identity_unchanged') is not True or state.get('input_fingerprint')!=fingerprint or state.get('phase_index')!=index:
            raise ValueError('相位状态身份/接受标志与汇总不一致')
        if not _accepted(result,True):raise ValueError('相位结果未通过原节点FE接受条件，不绘制候选压力')
        angle=effective['phases'][index]
        if not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and math.isclose(v,angle,abs_tol=1e-12,rel_tol=0.)
                for v in (row['driving_angle_rad'],state['driving_angle_rad'],result.get('driving_angle_rad'))):
            raise ValueError('实际相位与有效输入/状态不一致')
        if not math.isclose(result.get('pinion_torque_nm',float('nan')),effective['pinion_torque_nm'],rel_tol=1e-12,abs_tol=1e-10):raise ValueError('实际压力载荷与有效输入不一致')
        accepted.append((index,result,dict(file=row['state_file'],sha256=_sha(path))))
    if not accepted:raise ValueError('没有已接受节点FE状态，未生成接触图占位结果')
    manufacturing={m:base['manufacturing'][m] for m in ('pinion','wheel')}
    return summary,inputs,accepted,excluded,manufacturing,input_path


def plot_refined_fe_workflow(summary_path,output=None,cad_models=None,render3d=False,*,domain_certificates=None,cad_cache=None,tooth_offset=None):
    """从原汇总复绘中英连续接触云图；三维压力叠加已退出公开业务入口。"""
    if not isinstance(render3d,bool):raise ValueError('render3d须为布尔值')
    if render3d:raise ValueError('三维压力叠加已退出公开业务入口；请使用二维连续云图，普通CAD实体显示仍可单独使用')
    summary_path=Path(summary_path).resolve()
    summary,inputs,states,excluded,manufacturing,input_path=_accepted_states(summary_path)
    offset=inputs['refined_fe_case'].get('plot_tooth_offset',0) if tooth_offset is None else tooth_offset
    if isinstance(offset,bool) or not isinstance(offset,int):raise ValueError('绘图物理齿偏移须为整数')
    if domain_certificates is not None:
        if not isinstance(domain_certificates,dict) or any(str(index) not in domain_certificates for index,_,_ in states):
            raise ValueError('完整齿面覆盖证书须按已接受相位的十进制索引字符串逐一提供')
    default_output=summary_path.parent/'plots';output=Path(output).resolve() if output is not None else default_output
    meshes,cache,cad_mode=_load_or_build_meshes(output,Path(cad_cache).resolve() if cad_cache is not None else default_output,
        manufacturing,inputs['effective']['cad_sampling'],cad_models)
    from .refined_fe_plotting import plot_refined_fe_contact_pattern
    from .refined_fe_workflow import _write
    figures=[]
    for index,result,state in states:
        stem=f'phase_{index:04d}_contact_pattern'
        for lang in ('zh','en'):
            # 双语图共用SI源CSV，物理相位标识必须与显示语言无关，避免互相覆盖改变SHA。
            record=plot_refined_fe_contact_pattern([result],output,pinion_request=manufacturing['pinion']['request'],
                flank_mesh=meshes['pinion_flank'],labels=[f'phase_{index:04d}'],lang=lang,stem=stem,
                domain_certificates=None if domain_certificates is None else [domain_certificates[str(index)]],tooth_offset=offset)
            figures.append(dict(phase_index=index,view='complete_flank_and_patch',**record))
    directory=Path(__file__).resolve().parent
    manifest=dict(schema='bevel_refined_fe_plot_manifest_v1',input_fingerprint=summary['input_fingerprint'],
        source_summary_sha256=_sha(summary_path),effective_inputs_sha256=_sha(input_path),
        source_analysis_completed=summary.get('analysis_completed'),accepted_state_files=[state for _,_,state in states],
        excluded_states=excluded,cad_mode=cad_mode,cad_cache=dict(file='cad/cad_mesh_manifest.json',**cache),
        figures=figures,render3d=render3d,pressure_interpolation='C0_P1_scalar_display',
        original_pressure_values_modified=False,physics_states_modified=False,
        source_identity={'bevel_gear/'+name: _sha(directory/name) for name in
            ('refined_fe_postprocessing.py','refined_fe_plotting.py','pressure_field.py','pressure_background.py','halfspace_plotting.py','workflow_plotting.py')})
    _write(manifest,output/'refined_fe_plot_manifest.json')
    return manifest


def compare_contact_workflows(fe_summary_path,surface_summary_path,output,*,phase_index=0,cad_models=None,
        fe_domain_certificate=None,surface_domain_certificate=None,tooth_offset=0,cad_cache=None):
    """比较同制造几何、材料、相位和载荷的已接受FE/半解析连续云图。

    只读已保存物理场；缓存按真实制造输入和SHA绑定。局部半空间闭合量
    不作为整轮TE，显示支撑面积也不替代求解器原接触面积。
    """
    if isinstance(phase_index,bool) or not isinstance(phase_index,int) or phase_index<0:
        raise ValueError('phase_index须为非负整数')
    fe_path=Path(fe_summary_path).resolve();hs_path=Path(surface_summary_path).resolve()
    fe_summary,fe_inputs,fe_states,_,manufacturing,fe_input_path=_accepted_states(fe_path)
    selected=[item for item in fe_states if item[0]==phase_index]
    if len(selected)!=1:raise ValueError('FE结果没有指定的已接受相位')
    _,fe_result,fe_state_identity=selected[0]
    hs_summary=_read(hs_path)
    if hs_summary.get('schema')!='bevel_surface_workflow_v1' or hs_summary.get('analysis_completed') is not True or hs_summary.get('source_identity_unchanged') is not True:
        raise ValueError('半解析业务未接受，不能绘制方法比较云图')
    hs_input_path=_relative_file(hs_path.parent,hs_summary['input_file']);hs_inputs=_read(hs_input_path)
    hs_base=hs_inputs['base_case'];fe_base=fe_inputs['base_case']
    if hs_base.get('manufacturing')!=manufacturing:
        raise ValueError('两方法的双轮制造请求或齿坯不一致，不能混作同几何结果')
    if any(hs_base['materials'][m]!=fe_base['materials'][m] for m in ('pinion','wheel')):
        raise ValueError('两方法材料不一致，不能混作同工况结果')
    phase=fe_result['driving_angle_rad']
    matching=[row for row in hs_summary.get('states',[]) if row.get('accepted') is True and
        isinstance(row.get('driving_angle_rad'),(int,float)) and math.isclose(row['driving_angle_rad'],phase,rel_tol=0.,abs_tol=1e-12)]
    if len(matching)!=1:raise ValueError('半解析结果须恰有一个相同输入角的已接受相位')
    hs_state_path=_relative_file(hs_path.parent,matching[0]['state_file']);hs_state=_read(hs_state_path)
    if hs_state.get('accepted') is not True or hs_state.get('source_identity_unchanged') is not True or hs_state.get('input_fingerprint')!=hs_summary['input_fingerprint']:
        raise ValueError('半解析状态身份或接受标志与业务汇总不一致')
    from .surface_workflow import _accepted as surface_accepted
    hs_result=hs_state['result']
    if not surface_accepted(hs_result,True,hs_inputs['maximum_pinion_torque_relative_residual']):
        raise ValueError('半解析状态未通过原物理结果接受条件')
    if not math.isclose(hs_inputs['pinion_signed_torque_nm'],fe_result['pinion_torque_nm'],rel_tol=1e-12,abs_tol=1e-10):
        raise ValueError('两方法有符号扭矩不一致')
    destination=Path(output).resolve()
    meshes,cache,cad_mode=_load_or_build_meshes(destination,Path(cad_cache).resolve() if cad_cache is not None else fe_path.parent/'plots',
        manufacturing,fe_inputs['effective']['cad_sampling'],cad_models)
    from .refined_fe_plotting import plot_refined_fe_method_comparison
    from .refined_fe_workflow import _fingerprint,_write
    records=[plot_refined_fe_method_comparison(fe_result,hs_result,destination,
        pinion_request=manufacturing['pinion']['request'],flank_mesh=meshes['pinion_flank'],lang=lang,
        tooth_offset=tooth_offset,fe_domain_certificate=fe_domain_certificate,
        surface_domain_certificate=surface_domain_certificate)
        for lang in ('zh','en')]
    manifest=dict(schema='bevel_contact_method_comparison_v1',phase_index=phase_index,
        driving_angle_rad=phase,pinion_torque_nm=fe_result['pinion_torque_nm'],
        geometry_fingerprint=_fingerprint(manufacturing),cad_mode=cad_mode,
        source_state_sha256=dict(full_gear_fe=fe_state_identity['sha256'],local_halfspace=_sha(hs_state_path)),
        source_input_sha256=dict(full_gear_fe=_sha(fe_input_path),local_halfspace=_sha(hs_input_path)),
        pressure_interpolation='C0_P1_scalar_display',original_pressure_values_modified=False,
        comparison_scope='同物理输入的两种柔度模型展示，不等于压力精度对标',figures=records,
        source_identity={'bevel_gear/'+name: _sha(Path(__file__).resolve().parent/name) for name in
            ('refined_fe_postprocessing.py','refined_fe_plotting.py','pressure_field.py','pressure_background.py','halfspace_plotting.py','workflow_plotting.py')})
    _write(manifest,destination/'contact_comparison_manifest.json')
    return manifest
