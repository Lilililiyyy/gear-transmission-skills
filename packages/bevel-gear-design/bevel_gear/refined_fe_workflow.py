"""显式制造案例→面上加密带→整轮FE→节点—面接触的独立可重放业务。

本模块只组织输入、计算与输出；不调用半空间柔度，不平滑旧面片压力。
业务编排是本项目新写实现；参与源码以前后SHA256绑定。
导入只需标准库，NumPy/SciPy/OCP/Gmsh在实际实体/有限元阶段按需加载。
节点压力p=f/Σ(A_triangle/3)，真实中位对偶四边形用于显示，力矩积分点
仍是原压力节点；压力网格与物理验证状态不会由迭代成功自动推导。
"""
from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import time


REFINED_FE_CASE_SCHEMA='bevel_refined_fe_contact_case_v1'


def _real(value,name,positive=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or (positive and value<=0):
        raise ValueError(name+('须为正有限SI实数' if positive else '须为有限SI实数'))
    return float(value)


def _integer(value,name,minimum=1):
    if isinstance(value,bool) or not isinstance(value,int) or value<minimum:raise ValueError(name+'须为有效正整数')
    return value


def _json_safe(value):
    """数组转列表、未知数写null；自有引用只允许相对路径，不导出本机路径。"""
    if isinstance(value,dict):return {str(k):_json_safe(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [_json_safe(v) for v in value]
    if isinstance(value,str):
        if re.search(r'(?<![A-Za-z0-9_])[A-Za-z]:[\\/]',value) or value.startswith(('/','\\')):
            raise ValueError('输出含本机绝对路径；须提供相对引用或独立SHA身份')
        return value
    if value is None or isinstance(value,(bool,int)):return value
    if isinstance(value,float):return value if math.isfinite(value) else None
    if callable(getattr(value,'tolist',None)):return _json_safe(value.tolist())
    if callable(getattr(value,'item',None)):return _json_safe(value.item())
    raise TypeError('FE输出含不可序列化对象：'+type(value).__name__)


def _write(value,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(_json_safe(value),ensure_ascii=False,allow_nan=False,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    temporary.replace(path)
    return path


def _fingerprint(value):
    return hashlib.sha256(json.dumps(_json_safe(value),sort_keys=True,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode('utf-8')).hexdigest()


def _source_identity():
    """以公开相对名绑定真实参与计算/输出的源码，不记录本机源码目录。"""
    names=('refined_fe_workflow','generated_nodal_ltca','nodal_contact','fe_contact_solver',
        'generated_ltca','generated_tca','generated_solid','generated_surface','generated_cutter',
        'elasticity','contact_mesh','contact_equilibrium','cad','cad_import','workflow')
    directory=Path(__file__).resolve().parent
    return {'bevel_gear/'+name+'.py':hashlib.sha256((directory/(name+'.py')).read_bytes()).hexdigest() for name in names}


def load_refined_fe_case(path):
    """读取相对base_case并展开制造请求；本入口仅支持实际弧齿制造面。"""
    from .workflow import load_bevel_case
    path=Path(path).resolve();case=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(case,dict) or case.get('schema')!=REFINED_FE_CASE_SCHEMA:
        raise ValueError('refined FE案例schema须为'+REFINED_FE_CASE_SCHEMA)
    relative=case.get('base_case')
    if not isinstance(relative,str) or not relative or relative.startswith(('/','\\')) or Path(relative).is_absolute() or re.match(r'^[A-Za-z]:',relative):
        raise ValueError('refined FE案例须以相对base_case引用完整制造/材料案例')
    base_path=path.parent/relative;base=load_bevel_case(base_path)
    if base.get('kind')!='spiral':raise ValueError('本节点FE业务入口当前仅接入实际spiral弧齿制造面')
    for member in ('pinion','wheel'):
        item=base.get('manufacturing',{}).get(member,{})
        if not isinstance(item.get('request'),dict) or not isinstance(item.get('blank'),dict):
            raise ValueError('基础案例须显式提供两轮制造request与blank')
    return dict(refined_fe_case=deepcopy(case),base_case=base,
        refined_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        base_input_sha256=hashlib.sha256(base_path.read_bytes()).hexdigest())


def _configuration(loaded,phase_indices=None,mesh_size_m=None,pinion_torque_nm=None):
    """实体创建前检查SI网格、真实面加密、载荷/材料与节点接触参数。"""
    case,base=loaded['refined_fe_case'],loaded['base_case']
    analysis=deepcopy(case.get('analysis',{}))
    if not isinstance(analysis,dict) or set(analysis)-{'driver_angles_rad','tooth_offsets','pinion_torque_nm'}:raise ValueError('不支持的refined FE分析参数')
    phases=analysis.get('driver_angles_rad',[0.])
    if not isinstance(phases,list) or not phases:raise ValueError('主动相位须为非空rad列表')
    phases=[_real(v,'主动相位rad') for v in phases]
    if len(set(phases))!=len(phases):raise ValueError('主动相位不能重复')
    selected=list(range(len(phases))) if phase_indices is None else list(phase_indices)
    if not selected or len(set(selected))!=len(selected) or any(isinstance(i,bool) or not isinstance(i,int) or not 0<=i<len(phases) for i in selected):
        raise ValueError('phase_indices须为存在且不重复的零起始索引')
    offsets=analysis.get('tooth_offsets',[-1,0,1])
    if not isinstance(offsets,(list,tuple)) or not offsets or 0 not in offsets or len(set(offsets))!=len(offsets) or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets):
        raise ValueError('齿对偏移须为包含0的不重复整数列表')
    nominal=_real(base['load']['pinion_torque_n_m'],'基础扭矩幅值Nm',True)
    sign=base['load']['ltca_torque_sign']
    if isinstance(sign,bool) or sign not in (-1,1):raise ValueError('基础载荷符号须为-1或1')
    torque=_real(analysis.get('pinion_torque_nm',nominal*sign) if pinion_torque_nm is None else pinion_torque_nm,'有符号小轮扭矩Nm')
    if torque==0:raise ValueError('节点承载接触要求非零有符号Nm')
    fe=dict(element_size_m=.004,contact_refined_size_m=.000075,refinement_radius_m=.0012,maximum_nodes=300000)
    specified=case.get('finite_element',{})
    if not isinstance(specified,dict) or set(specified)-set(fe):raise ValueError('不支持的有限元网格参数')
    fe.update(specified)
    if mesh_size_m is not None:fe['contact_refined_size_m']=mesh_size_m
    for name in ('element_size_m','contact_refined_size_m','refinement_radius_m'):fe[name]=_real(fe[name],name,True)
    _integer(fe['maximum_nodes'],'节点数上限',4)
    if fe['contact_refined_size_m']>=fe['element_size_m']:raise ValueError('局部尺寸须小于主体尺寸')
    if fe['refinement_radius_m']<=fe['contact_refined_size_m']:raise ValueError('加密影响半径须大于局部尺寸')
    refinement=dict(half_widths_m=[.006,.0018],track_spacing_m=None,sample_spacing_m=.0005)
    specified=case.get('refinement',{})
    if not isinstance(specified,dict) or set(specified)-set(refinement):raise ValueError('不支持的实际面加密带参数')
    refinement.update(specified)
    widths=refinement['half_widths_m']
    if not isinstance(widths,(tuple,list)) or len(widths)!=2:raise ValueError('加密半宽须为两个正m')
    refinement['half_widths_m']=[_real(v,'实际面加密半宽m',True) for v in widths]
    refinement['track_spacing_m']=_real(fe['contact_refined_size_m'] if refinement['track_spacing_m'] is None else refinement['track_spacing_m'],'加密轨迹间距m',True)
    refinement['sample_spacing_m']=_real(refinement['sample_spacing_m'],'加密曲线采样间距m',True)
    contact=dict(maximum_distance_m=.00015,maximum_normal_angle_rad=.3,matching_tolerance_m=8e-7,
        mesh_phase_tolerance_m=1e-10,maximum_mesh_phase_iterations=15,maximum_candidates=30000,
        maximum_loaded_angle_rad=.005,maximum_pinion_torque_relative_residual=.02,maximum_working_set=2048)
    specified=case.get('contact',{})
    if not isinstance(specified,dict) or set(specified)-set(contact):raise ValueError('不支持的节点接触参数')
    contact.update(specified)
    for name,value in contact.items():
        if name in ('maximum_mesh_phase_iterations','maximum_candidates','maximum_working_set'):_integer(value,name)
        else:contact[name]=_real(value,name,True)
    if contact['maximum_normal_angle_rad']>=math.pi/2:raise ValueError('接触法向容差须小于pi/2')
    if contact['mesh_phase_tolerance_m']>1e-9:raise ValueError('节点网格初相位容差不得宽于1e-9m')
    materials={}
    for member in ('pinion','wheel'):
        values=base['materials'][member]
        e=_real(values['young_modulus_pa'],member+'弹性模量Pa',True);nu=_real(values['poisson_ratio'],member+'泊松比')
        if not -1<nu<.49:raise ValueError('低阶三维弹性泊松比须位于(-1,.49)')
        materials[member]=dict(young_modulus_pa=e,poisson_ratio=nu)
    sampling=deepcopy(base.get('cad_sampling',{}))
    if set(sampling)-{'width_count','profile_count'}:raise ValueError('不支持的CAD采样参数')
    for key,value in sampling.items():_integer(value,key,3)
    return dict(phases=phases,selected=selected,tooth_offsets=list(offsets),pinion_torque_nm=torque,
        finite_element=fe,refinement=refinement,contact=contact,materials=materials,cad_sampling=sampling)


def _failure(error):
    message=str(error)
    if re.search(r'(?<![A-Za-z0-9_])[A-Za-z]:[\\/]',message) or message.startswith(('/','\\')):
        message='运行阶段发生本机文件/引擎错误，具体路径不导出；请查看异常类型与阶段'
    return dict(type=type(error).__name__,message=message)


def _accepted(result,source_unchanged):
    """接受严格节点FE平衡与真实对偶格；不将该标志解释为峰压收敛。"""
    if not result.get('converged') or not source_unchanged or result.get('contact_representation')!='nodal_dual_cells':return False
    equilibrium=result.get('equilibrium',{})
    if not equilibrium.get('converged') or not equilibrium.get('certificate',{}).get('passed'):return False
    cells=result.get('contacts') or []
    if not cells:return False
    for cell in cells:
        p,a,f=(cell.get(name) for name in ('mean_pressure_pa','area_m2','force_n'))
        if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in (p,a,f)) or p<0 or a<=0 or f<0:return False
        if not math.isclose(p*a,f,rel_tol=1e-10,abs_tol=1e-9):return False
        polygon=cell.get('first_vertices_world_m');point=cell.get('application_point_world_m')
        if not isinstance(polygon,(list,tuple)) or len(polygon)!=4 or not isinstance(point,(list,tuple)) or len(point)!=3:return False
        for xyz in [*polygon,point]:
            if not isinstance(xyz,(list,tuple)) or len(xyz)!=3 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in xyz):return False
    total=result.get('total_normal_force_n')
    if isinstance(total,bool) or not isinstance(total,(int,float)) or not math.isfinite(total):return False
    if not math.isclose(sum(c['force_n'] for c in cells),total,rel_tol=1e-10,abs_tol=1e-8):return False
    energy=result.get('fe_to_condensed_energy_relative_error')
    return isinstance(energy,(int,float)) and math.isfinite(energy) and abs(energy)<=1e-7


def _export_contacts(result,path,phase_index):
    """写真实对偶四边形和原压力节点，不以四边形中心重新定义力矩作用点。"""
    names=('phase_index','tooth_offset','first_pressure_node_index','force_n','mean_pressure_pa','area_m2',
        'remaining_gap_m','polygon_world_m_json','application_point_world_m_json')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=names,lineterminator='\n');writer.writeheader()
        for cell in result['contacts']:
            row={key:cell.get(key) for key in names[:7]};row['remaining_gap_m']=cell.get('remaining_gap_m')
            row.update(phase_index=phase_index,polygon_world_m_json=json.dumps(cell['first_vertices_world_m']),
                application_point_world_m_json=json.dumps(cell['application_point_world_m']))
            writer.writerow(row)


class _StateObserverStateChanged(RuntimeError):
    """只读核查观察器改变原状态时，单独拒绝且保留原字节。"""


def _observe_saved_state(observer,prepared,index,path):
    """在原状态落盘后复用同一FE准备对象核查，禁止观察器改写原结果。"""
    path=Path(path);original=path.read_bytes();digest=hashlib.sha256(original).hexdigest()
    error=None
    try:
        observer(prepared,index,path)
    except Exception as caught:
        error=caught
    current=path.read_bytes() if path.is_file() else None
    if current!=original:
        # 恢复原求解字节，避免核查插件的意外写入污染已接受物理数组。
        temporary=path.with_name(path.name+'.restore.tmp');temporary.write_bytes(original);temporary.replace(path)
        raise _StateObserverStateChanged('只读相位核查改变了原状态；已还原原字节并停止业务') from error
    if error is not None:
        raise error
    return dict(status='completed',source_state_sha256=digest,source_state_unchanged=True)


def run_refined_fe_workflow(case_path,output_dir,*,phase_indices=None,mesh_size_m=None,
        pinion_torque_nm=None,cad_models=None,meshes=None,export_meshes=True,include_fields=False,progress=None,
        state_observer=None):
    """公开冷启动构建真实整轮；显式CAD/网格复用用于调用方已绑定的回放。

    各相位完整JSON先落盘后更新小汇总。网格文件保持SI原数组，诊断明确
    复用状态；无独立CAD来源绑定的复用网格不冒称拓扑身份已认证。可选
    state_observer(prepared,index,state_path)只读核查已落盘相位，以复用
    同一稀疏分解；普通核查失败单列，不冒充力学失败或审核通过，原状态
    改写及源码变化立即停止。图形后处理不在此插值或平滑压力。
    """
    if state_observer is not None and not callable(state_observer):raise TypeError('state_observer须为可调用只读核查函数或None')
    for value,name in ((export_meshes,'export_meshes'),(include_fields,'include_fields')):
        if not isinstance(value,bool):raise ValueError(name+'须为布尔值')
    for value,name in ((cad_models,'cad_models'),(meshes,'meshes')):
        if value is not None and (not isinstance(value,dict) or set(value)!={'pinion','wheel'}):raise ValueError(name+'须显式给出pinion/wheel两项')
    loaded=load_refined_fe_case(case_path);effective=_configuration(loaded,phase_indices,mesh_size_m,pinion_torque_nm)
    fingerprint=_fingerprint(dict(base_case=loaded['base_case'],**{k:v for k,v in effective.items() if k!='selected'}))
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    before=_source_identity()
    summary=dict(schema='bevel_refined_fe_workflow_v1',analysis_completed=False,kind='spiral',
        name=loaded['refined_fe_case'].get('name','实际弧齿整轮高分辨率节点接触'),input_file='refined_fe_inputs.json',
        implementation_identity='本项目新写业务编排；真实制造/整轮FE/节点接触源码以前后SHA256绑定',
        input_fingerprint=fingerprint,source_identity_before=before,source_identity_after=None,source_identity_unchanged=None,
        cad_mode='explicit_same_request_objects' if cad_models is not None else 'built_from_explicit_request_and_blank',
        mesh_mode='explicit_caller_meshes' if meshes is not None else 'generated_by_gmsh_from_actual_solids',
        states=[],selected_phase_indices=effective['selected'],mesh_files={},prepared_geometry_file=None,
        pressure_definition='节点集总压力p=f/Σ(A_triangle/3)；显示为真实FE三角形内中位对偶四边形',
        moment_application='原压力节点积分点；不是对偶四边形几何中心',
        pressure_peak_validated=False,pressure_mesh_converged=None,physical_validation=False,
        state_observer_status='requested' if state_observer is not None else 'not_requested',
        plots_status='independent_postprocessing',contact_model='整轮FE柔度；不调用半空间柔度')
    _write(dict(**loaded,effective=effective,input_fingerprint=fingerprint),output/'refined_fe_inputs.json')
    _write(summary,output/'refined_fe_summary.json')
    def emit(stage,**values):
        if progress is not None:progress(dict(stage=stage,**values))
    manufacturing=loaded['base_case']['manufacturing'];requests={m:manufacturing[m]['request'] for m in ('pinion','wheel')}
    try:
        from .generated_solid import build_generated_gear
        from .generated_tca import prepare_generated_tca
        from .generated_nodal_ltca import generated_refinement_tracks
        from .generated_ltca import prepare_generated_ltca
        emit('prepare_actual_solids',completed=False)
        models=({m:build_generated_gear(requests[m],manufacturing[m]['blank'],**effective['cad_sampling']) for m in requests}
            if cad_models is None else cad_models)
        tca=prepare_generated_tca(requests['pinion'],requests['wheel'],cad_models=models,**effective['cad_sampling'])
        emit('prepare_actual_solids',completed=True)
        tracks={m:[] for m in requests};track_diagnostics=[]
        for angle in effective['phases']:
            control=generated_refinement_tracks(tca,driving_angle_rad=angle,tooth_offsets=effective['tooth_offsets'],**effective['refinement'])
            for member in requests:tracks[member].extend(control['tracks'][member])
            track_diagnostics.append(dict(driving_angle_rad=angle,diagnostics=control['diagnostics']))
        _write(dict(tracks=tracks,diagnostics=track_diagnostics,role='只控制实际FE网格，不作为接触域或载荷分布'),output/'refinement_tracks.json')
        emit('finite_element_preparation',completed=False)
        prepared=prepare_generated_ltca(requests['pinion'],requests['wheel'],
            pinion_blank_settings=manufacturing['pinion']['blank'],wheel_blank_settings=manufacturing['wheel']['blank'],
            pinion_material=effective['materials']['pinion'],wheel_material=effective['materials']['wheel'],
            directory=output/'finite_element',cad_models=models,meshes=meshes,**effective['cad_sampling'],
            refinement_driving_angle_rad=effective['phases'][0],refinement_polylines_m=tracks,**effective['finite_element'])
        _write(prepared.diagnostics,output/'prepared_geometry.json');summary['prepared_geometry_file']='prepared_geometry.json'
        if export_meshes:
            for member in requests:
                mesh=prepared.systems[member].mesh;relative=f'meshes/{member}_mesh.json'
                _write(dict(nodes_m=mesh.nodes_m,tetrahedra=mesh.tetrahedra,boundary_triangles=mesh.boundary_triangles,
                    boundary_tags=mesh.boundary_tags,metadata=mesh.metadata),output/relative)
                summary['mesh_files'][member]=relative
        emit('finite_element_preparation',completed=True)
    except (ArithmeticError,RuntimeError,ValueError,TypeError,ImportError,OSError) as error:
        after=_source_identity();summary.update(status='preparation_failed',failure=_failure(error),
            source_identity_after=after,source_identity_unchanged=after==before)
        _write(summary,output/'refined_fe_summary.json');emit('finite_element_preparation',completed=False,status='preparation_failed')
        return summary
    _write(summary,output/'refined_fe_summary.json')
    from .generated_nodal_ltca import solve_generated_nodal_fe_ltca
    for index in effective['selected']:
        tick=time.perf_counter();angle=effective['phases'][index]
        emit('nodal_contact',phase_index=index,driving_angle_rad=angle,completed=False)
        state=dict(schema='bevel_refined_fe_state_v1',phase_index=index,driving_angle_rad=angle,
            input_fingerprint=fingerprint,source_identity_before=before,accepted=False)
        try:
            result=solve_generated_nodal_fe_ltca(prepared,angle,effective['pinion_torque_nm'],
                tooth_offsets=effective['tooth_offsets'],include_fields=include_fields,progress=progress,**effective['contact'])
            result=dict(result);result['pressure_peak_validated']=False
        except (ArithmeticError,RuntimeError,ValueError,TypeError,ImportError,OSError) as error:
            result=dict(converged=False,status='nodal_computation_failed',contacts=[],pressure_peak_validated=False)
            state['failure']=_failure(error)
        after=_source_identity();unchanged=before==after
        state.update(result=result,source_identity_after=after,source_identity_unchanged=unchanged,
            seconds=time.perf_counter()-tick,accepted=_accepted(result,unchanged))
        relative=f'states/phase_{index:04d}.json';_write(state,output/relative)
        row=dict(phase_index=index,driving_angle_rad=angle,state_file=relative,accepted=state['accepted'],
            converged=bool(result.get('converged')),status=result.get('status'),seconds=state['seconds'],
            candidate_count=result.get('candidate_count'),active_nodal_count=result.get('active_nodal_count') if state['accepted'] else None,
            active_dual_cell_count=result.get('active_contact_count') if state['accepted'] else None,
            maximum_nodal_pressure_pa=result.get('maximum_cell_mean_pressure_pa') if state['accepted'] else None,
            total_normal_force_n=result.get('total_normal_force_n') if state['accepted'] else None,
            elastic_te_rad=result.get('elastic_te_rad') if state['accepted'] else None,
            incremental_mesh_stiffness_nm_per_rad=result.get('incremental_mesh_stiffness_nm_per_rad') if state['accepted'] else None,
            pinion_torque_balance_relative_residual=result.get('pinion_torque_balance_relative_residual'),pressure_peak_validated=False)
        if state['accepted']:
            csv_relative=f'contacts/phase_{index:04d}.csv';_export_contacts(result,output/csv_relative,index);row['contact_csv']=csv_relative
        summary['states'].append(row);summary.update(source_identity_after=after,source_identity_unchanged=unchanged)
        _write(summary,output/'refined_fe_summary.json');emit('nodal_contact',**row,completed=True)
        if not unchanged:raise RuntimeError('实际计算期间参与源码改变；状态已保存，本次业务不验收为可复现完成')
        if state_observer is not None:
            observed_tick=time.perf_counter()
            try:
                row['state_observer']=_observe_saved_state(state_observer,prepared,index,output/relative)
            except _StateObserverStateChanged as error:
                row['state_observer']=dict(status='state_modified_and_restored',failure=_failure(error))
                summary.update(status='state_observer_state_modified',state_observer_status='failed')
                _write(summary,output/'refined_fe_summary.json')
                raise
            except Exception as error:
                row['state_observer']=dict(status='failed',failure=_failure(error),source_state_unchanged=True)
            row['state_observer']['seconds']=time.perf_counter()-observed_tick
            observer_after=_source_identity()
            summary.update(source_identity_after=observer_after,source_identity_unchanged=observer_after==before)
            if observer_after!=before:
                summary.update(status='state_observer_source_changed',state_observer_status='failed')
                _write(summary,output/'refined_fe_summary.json')
                raise RuntimeError('只读相位核查期间计算源码改变；原状态已保存，业务停止')
            summary['state_observer_status']=('contains_failures' if any(
                s.get('state_observer',{}).get('status')!='completed' for s in summary['states']) else 'completed')
            _write(summary,output/'refined_fe_summary.json')
            emit('state_observer',phase_index=index,status=row['state_observer']['status'],completed=True)
    after=_source_identity();summary.update(source_identity_after=after,source_identity_unchanged=before==after)
    summary['analysis_completed']=bool(summary['states']) and all(s['accepted'] for s in summary['states']) and summary['source_identity_unchanged']
    summary['status']='completed' if summary['analysis_completed'] else 'contains_unaccepted_states'
    _write(summary,output/'refined_fe_summary.json')
    return summary
