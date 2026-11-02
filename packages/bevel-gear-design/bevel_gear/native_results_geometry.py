"""只核本态原有限φ门控的完整协议；不冒充独立有限φ重求。"""
import math


def _number(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('原几何门控缺少有限值：'+name)
    return float(value)


def validate_original_native_geometry_gate(result,*,declared_controls=None):
    finite=result.get('geometry_check')
    if (not isinstance(finite,dict) or finite.get('completed') is not True or finite.get('passed') is not True or
            finite.get('source_NURBS_gap_reused') is not False or finite.get('loaded_FE_phi_geometry_check',{}).get('finite_phi_geometry_rebuilt') is not True):raise ValueError('正式phase缺少完整原有限φ几何门控')
    mean=_number(finite['tolerances']['mean_error_m'],'mean_error_m');penetration=_number(finite['tolerances']['penetration_m'],'penetration_m')
    # 只接受本态已绑定输入的显式点间隙容差；缺省仍0.2μm，不从结果反推声明。
    controls={} if declared_controls is None else declared_controls
    if not isinstance(controls,dict):raise ValueError('已绑定的native_contact_solver须为字典')
    declared=_number(controls.get('maximum_piece_penetration_m',2e-7),'declared_piece_penetration_m')
    if not 0<declared<=1e-6:raise ValueError('已声明点间隙容差须为正且不超过1um')
    if not 0<mean<=1e-10 or not 0<penetration<=1e-6:raise ValueError('原有限φ门限超出1e-10m/显式1um')
    if penetration!=declared:raise ValueError('原reported点间隙容差与已绑定声明不同')
    original=finite['initial_piece_vertex_check'];loaded=finite['loaded_FE_phi_geometry_check']
    for branch,name in ((original,'initial'),(loaded,'loaded')):
        if (branch.get('passed') is not True or branch.get('partial_cells')!=0 or
                isinstance(branch.get('checked_vertices'),bool) or not isinstance(branch.get('checked_vertices'),int) or branch['checked_vertices']<=0 or
                _number(branch['minimum_remaining_gap_m'],name+'_minimum') < -penetration):raise ValueError('原有限φ顶点或完整域门控未通过')
    for value,name in ((original['maximum_mean_vs_QP_error_m'],'QP_mean'),(original['maximum_original_q_mean_vs_W_mean_error_m'],'q_W_mean'),(loaded['maximum_cell_mean_error_m'],'loaded_mean')):
        if not 0<=_number(value,name)<=mean:raise ValueError('原有限φ平均间隙门控未通过')
    expected={(field['tooth_offset'],j,i) for field in result['fields'] for j,row in enumerate(field['valid_mask']) for i,valid in enumerate(row) if valid}
    rows=finite['cell_checks'];keys=[(row['tooth_offset'],row['row'],row['column']) for row in rows]
    if not expected or len(keys)!=len(set(keys)) or set(keys)!=expected:raise ValueError('原有限φ逐格证明不完整')
    for row in rows:
        if (row.get('loaded_complete') is not True or
                _number(row['old_minimum_piece_vertex_gap_m'],'old_minimum') < -penetration or
                _number(row['loaded_minimum_piece_vertex_gap_m'],'loaded_minimum') < -penetration or
                abs(_number(row['old_mean_remaining_m'],'old_mean')-_number(row['reference_remaining_gap_m'],'QP_mean'))>mean or
                abs(_number(row['loaded_vs_old_predicted_mean_error_m'],'loaded_error'))>mean):raise ValueError('原有限φ逐格门控未通过')
    step=_number(result['linearization_step_theta_rad'],'末步QP增量')
    if _number(finite['delta_theta_rad'],'geometry_delta')!=step or _number(result['equilibrium']['theta_rad'],'QP_theta')!=step:raise ValueError('原有限φ增量不是本态末步QP增量')
    reference=_number(result['reference_follower_angle_rad'],'reference_phi');angle=_number(result['loaded_follower_angle_rad'],'loaded_phi')
    ratio=_number(result['geometry']['pinion_to_wheel_ratio'],'ratio');driver=_number(result['geometry']['driving_angle_rad'],'driver')
    if ratio<=0 or not math.isclose(reference+step,angle,rel_tol=1e-9,abs_tol=1e-12) or not math.isclose(angle+ratio*driver,_number(result['total_loaded_te_rad'],'total_TE'),rel_tol=1e-9,abs_tol=1e-12):raise ValueError('原受载角/末步增量/总TE簿记不一致')
    return dict(original_geometry_gate_verified=True,independent_finite_phi_recheck_performed=False,checked_cells=len(rows),mean_tolerance_m=mean,penetration_tolerance_m=penetration)
