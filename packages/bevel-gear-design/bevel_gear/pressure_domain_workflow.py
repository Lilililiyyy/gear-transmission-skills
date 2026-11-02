"""由已接受业务状态重建几何，审核完整压力域并保存可绑定复绘的证书。"""
from __future__ import annotations

import json
import math
from pathlib import Path


def _surface_state(path,phase_index):
    """只读原输入及接受状态；失败相位不以零压力替代。"""
    from .refined_fe_postprocessing import _read,_relative_file
    from .surface_workflow import _accepted,_fingerprint
    summary=_read(path);root=path.parent
    if summary.get('schema')!='bevel_surface_workflow_v1' or summary.get('source_identity_unchanged') is not True:
        raise ValueError('需要源身份稳定的半解析业务汇总')
    inputs=_read(_relative_file(root,summary['input_file']))
    fingerprint=_fingerprint(dict(base_case=inputs['base_case'],analysis=inputs['effective_analysis'],
        driver_angles_rad=inputs['surface_case']['analysis']['driver_angles_rad'],solver=inputs['effective_solver'],
        pinion_signed_torque_nm=inputs['pinion_signed_torque_nm'],materials=inputs['materials']))
    if inputs.get('effective_input_fingerprint')!=fingerprint or summary.get('input_fingerprint')!=fingerprint:
        raise ValueError('半解析实际输入指纹已变化，不能套用原压力覆盖证书')
    rows=[r for r in summary['states'] if r['phase_index']==phase_index]
    if len(rows)!=1 or rows[0].get('accepted') is not True:
        raise ValueError('指定半解析相位没有已接受结果')
    state_path=_relative_file(root,rows[0]['state_file']);state=_read(state_path)
    if state.get('accepted') is not True or state.get('source_identity_unchanged') is not True or state.get('input_fingerprint')!=summary['input_fingerprint']:
        raise ValueError('半解析状态身份与汇总不一致')
    result=state['result']
    if not _accepted(result,True,inputs['maximum_pinion_torque_relative_residual']):
        raise ValueError('半解析原结果不满足接受条件')
    angle=rows[0]['driving_angle_rad'];torque=inputs['pinion_signed_torque_nm']
    geometry=result['geometry']
    if not math.isclose(geometry['driving_angle_rad'],angle,abs_tol=1e-12,rel_tol=0.) or not math.isclose(geometry['pinion_torque_nm'],torque,abs_tol=1e-10,rel_tol=1e-12):
        raise ValueError('半解析状态的相位或有符号载荷不一致')
    return summary,inputs,result,state_path


def audit_pressure_workflow(summary_path,output,*,method,phase_index=0,tooth_offset=0,
        cad_models=None,maximum_depth=8,maximum_regions=300000,cad_cache=None,flank_linear_deflection_m=1e-6,progress=None):
    """对一个实际相位和物理齿审核；不重算接触或改写原状态。

    输出certificate.json供对应复绘入口读取。若仍有未决区域，passed为
    false，绘图不可据此把外域铺为零压。CAD/FE可选引擎在执行时才导入。
    """
    from .refined_fe_postprocessing import (_accepted_states,_read,_relative_file,_sha,
        _build_models,_validate_models,_load_or_build_meshes)
    from .refined_fe_workflow import _write
    if method not in ('full_gear_fe','local_halfspace'):raise ValueError('方法须为full_gear_fe或local_halfspace')
    if progress is not None and not callable(progress):raise ValueError('进度回调须可调用')
    if isinstance(flank_linear_deflection_m,bool) or not isinstance(flank_linear_deflection_m,(int,float)) or not math.isfinite(flank_linear_deflection_m) or flank_linear_deflection_m<=0:
        raise ValueError('审核工作面弦差须为正有限m')
    for value,name in ((phase_index,'相位索引'),(tooth_offset,'物理齿偏移')):
        if isinstance(value,bool) or not isinstance(value,int) or name=='相位索引' and value<0:
            raise ValueError(name+'须为有效整数')
    path=Path(summary_path).resolve();destination=Path(output).resolve()
    if method=='full_gear_fe':
        summary,inputs,states,_,manufacturing,_=_accepted_states(path)
        selected=[s for s in states if s[0]==phase_index]
        if len(selected)!=1:raise ValueError('指定FE相位没有唯一已接受结果')
        _,result,identity=selected[0];state_path=_relative_file(path.parent,identity['file'])
        base=inputs['base_case'];sampling=inputs['effective']['cad_sampling']
    else:
        summary,inputs,result,state_path=_surface_state(path,phase_index)
        base=inputs['base_case'];manufacturing=base['manufacturing'];sampling=base.get('cad_sampling',{})
    before=_sha(state_path)
    models=_build_models(manufacturing,sampling) if cad_models is None else cad_models
    _validate_models(models,manufacturing)
    meshes,cache,cad_mode=_load_or_build_meshes(destination,Path(cad_cache).resolve() if cad_cache is not None else destination,
        manufacturing,sampling,models,flank_linear_deflection_m=flank_linear_deflection_m)
    if method=='local_halfspace':
        from .generated_halfspace import prepare_generated_halfspace
        from .pressure_domain_audit import audit_generated_halfspace_pressure_domain
        prepared=prepare_generated_halfspace(manufacturing['pinion']['request'],manufacturing['wheel']['request'],
            pinion_blank_settings=manufacturing['pinion']['blank'],wheel_blank_settings=manufacturing['wheel']['blank'],
            cad_models=models,**sampling)
        certificate=audit_generated_halfspace_pressure_domain(prepared,result,meshes['pinion_flank'],
            manufacturing['pinion']['request'],tooth_offset=tooth_offset,maximum_depth=maximum_depth,
            maximum_regions=maximum_regions,progress=progress)
    else:
        from .elasticity import TetraMesh
        from .generated_ltca import prepare_generated_ltca
        from .pressure_domain_audit import audit_generated_fe_pressure_domain
        if set(summary.get('mesh_files',{}))!={'pinion','wheel'}:
            raise ValueError('FE域审核需要保存的双轮完整SI网格')
        volume_meshes={m:TetraMesh(**_read(_relative_file(path.parent,reference))) for m,reference in summary['mesh_files'].items()}
        materials=inputs['effective']['materials']
        prepared=prepare_generated_ltca(manufacturing['pinion']['request'],manufacturing['wheel']['request'],
            pinion_blank_settings=manufacturing['pinion']['blank'],wheel_blank_settings=manufacturing['wheel']['blank'],
            pinion_material=materials['pinion'],wheel_material=materials['wheel'],directory=destination/'elastic_preparation',
            cad_models=models,meshes=volume_meshes,**sampling)
        certificate=audit_generated_fe_pressure_domain(prepared,result,meshes['pinion_flank'],
            manufacturing['pinion']['request'],tooth_offset=tooth_offset,maximum_depth=maximum_depth,
            maximum_regions=maximum_regions,progress=progress)
    after=_sha(state_path)
    if before!=after:raise ArithmeticError('覆盖审核改变了原物理状态文件')
    _write(certificate,destination/'certificate.json')
    checks=certificate['checks']
    quality_key='included_candidates_kkt_passed' if method=='full_gear_fe' else 'loaded_boundary_checks_passed'
    partial_available=bool(not certificate['passed'] and certificate['unresolved_region_count']>0
        and checks.get('input_geometry_bound') is True and checks.get(quality_key) is True
        and checks.get('audit_inputs_and_source_unchanged') is True and certificate.get('regions',{}).get('regions'))
    manifest=dict(schema='bevel_pressure_domain_workflow_v1',method=method,phase_index=phase_index,tooth_offset=tooth_offset,
        passed=certificate['passed'],certificate_file='certificate.json',source_summary_sha256=_sha(path),
        source_state_sha256=before,original_state_unchanged=before==after,cad_mode=cad_mode,
        cad_cache=cache,unresolved_region_count=certificate['unresolved_region_count'],
        partial_domain_evidence_available=partial_available,
        scope='单相位、单齿的实际CAD三角离散与显式弦差审核；不等于连续CAD严格证明、压力精度或物理验证')
    _write(manifest,destination/'pressure_domain_manifest.json')
    return manifest
