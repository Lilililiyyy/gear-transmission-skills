"""真实FE柔度列驱动的非负接触QP，不组装全部候选的稠密柔度矩阵。

调用方提供C=WᵀKff⁻¹W的列；本模块只组织受限变量求解与全域KKT复核，
不生成几何、不计算Hertz/BEM柔度、不增加经验对角刚度。未知f为节点合力N，
g为m，ell为m/rad，T为N·m，theta为rad；w=g+C f−ell theta。
局部QP解出后，对全部候选复核w，逐批加入负约化梯度对应变量。只有
全候选平衡/互补通过才返回物理场，候选解与失败状态分别保存。
"""
from __future__ import annotations

import math
from numbers import Integral


def _numpy():
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError('高分辨率FE接触列生成需要可选NumPy') from error
    return np


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(name+'须为正整数')
    return value


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(name+'须为正有限实数')
    return float(value)


def solve_fe_contact_equilibrium(gaps_m, *, lever_arms_m, torque_nm, compliance_columns,
        initial_count=16, add_count=48, maximum_working_set=2048,
        maximum_outer_iterations=100, inner_max_iterations=10000,
        allow_initial_penetration=False, force_tolerance_n=1e-8,
        gap_tolerance_m=1e-10, balance_tolerance_nm=1e-8,
        relative_tolerance=1e-9, symmetry_tolerance=1e-9,
        minimum_normalized_pivot=1e-12, maximum_column_cache_bytes=1073741824,
        initial_indices=None):
    """通过FE柔度列生成求有符号目标扭矩下的全候选接触平衡。

    compliance_columns(indices)须返回有限N×len(indices)数组，严格保持输入
    索引顺序，单位m/N。调用方须保证它来自同一受约束线弹性系统的对称
    半正定柔度；本函数仅检查访问过的互等项与当前受限块的正定性，不能
    将局部块检查称作全部N×N矩阵的独立正定认证。所有已访问列缓存复用。
    负初始gap默认拒绝，显式允许时保留原值而不夹零；零扭矩仅在无初始
    穿入时返回无载解，零扭矩有穿入分支明确拒绝。失败字段为空，最后
    候选留在candidate中，不把局部QP收敛冒充全部候选已通过。
    initial_indices仅提供初始受限变量索引，可复用前一相位活动集；
    不复用旧力或柔度，最终仍对全部当前候选检查非穿透与互补。
    """
    from .contact_equilibrium import solve_contact_equilibrium
    np = _numpy()
    gaps = np.asarray(gaps_m, dtype=float)
    levers = np.asarray(lever_arms_m, dtype=float)
    if gaps.ndim != 1 or not len(gaps) or levers.shape != gaps.shape or not np.isfinite(gaps).all() or not np.isfinite(levers).all():
        raise ValueError('间隙与力臂须为相同非空长度的有限SI向量')
    if isinstance(torque_nm, bool) or not isinstance(torque_nm, (int, float)) or not math.isfinite(torque_nm):
        raise ValueError('目标扭矩须为有限N·m')
    target = float(torque_nm)
    if not callable(compliance_columns):raise TypeError('compliance_columns须为真实FE柔度列回调')
    if not isinstance(allow_initial_penetration, bool):raise ValueError('初始穿入开关须为布尔值')
    if gaps.min() < 0 and not allow_initial_penetration:
        raise ValueError('未加载gap为负；须显式允许，不能静默夹零')
    for value, name in ((initial_count, '初始变量数'), (add_count, '每批新增变量数'),
            (maximum_working_set, '受限变量上限'), (maximum_outer_iterations, '外迭代上限'),
            (inner_max_iterations, '局部QP迭代上限'), (maximum_column_cache_bytes, '柔度列缓存字节上限')):
        _integer(value, name)
    for value, name in ((force_tolerance_n, '接触力容差N'), (gap_tolerance_m, '间隙容差m'),
            (balance_tolerance_nm, '力矩容差N·m'), (relative_tolerance, '相对容差'),
            (symmetry_tolerance, '互等容差'), (minimum_normalized_pivot, '归一化主元门限')):
        _positive(value, name)
    n = len(gaps)
    warm=[]
    if initial_indices is not None:
        if not isinstance(initial_indices,(list,tuple)):
            raise ValueError('初始活动索引须为列表或元组')
        if (any(isinstance(i,bool) or not isinstance(i,Integral) or not 0<=i<n for i in initial_indices)
                or len(set(initial_indices))!=len(initial_indices) or len(initial_indices)>maximum_working_set):
            raise ValueError('初始活动索引越界、重复或超过受限变量上限')
        warm=sorted(map(int,initial_indices))
    history = []
    cache = {}
    queried = []
    linear_solves = 0
    largest = 0
    full_matrix_assembled = False
    candidate = None
    diagnostics = dict(model='真实FE柔度列生成；没有Hertz/BEM或经验刚度',
        candidate_count=n, units=dict(force='N', gap='m', compliance='m/N', lever='m', torque='N·m', angle='rad'),
        full_matrix_assembled=False, maximum_working_set=maximum_working_set,
        maximum_column_cache_bytes=maximum_column_cache_bytes, allow_initial_penetration=allow_initial_penetration,
        minimum_original_gap_m=float(gaps.min()),
        supplied_initial_indices=warm,
        convexity_contract='调用方提供同一真实FE的对称半正定C；只访问部分列，不声称全矩阵独立正定检查')

    def finish(converged, status, failure=None):
        diagnostics.update(queried_column_count=len(cache), queried_column_indices=queried,
            column_cache_bytes=sum(column.nbytes for column in cache.values()),
            largest_restricted_size=largest, full_matrix_assembled=full_matrix_assembled, outer_history=history)
        result = dict(converged=converged, status=status, failure=failure, diagnostics=diagnostics,
            target_load=target, load_unit='N·m', iterations=len(history), linear_solve_count=linear_solves,
            forces_n=[], remaining_gaps_m=[], elastic_deflections_m=[], theta_rad=None,
            normal_approach_m=None, strain_energy_j=None, active_indices=[],
            resultant_load=None, balance_residual=None, relative_balance_residual=None,
            complementarity_residual_nm=None, contact_stationarity_residual_m=None,
            minimum_force_n=None, minimum_remaining_gap_m=None, certificate=None, candidate=None)
        if converged:
            result.update(candidate)
        else:
            result['candidate'] = candidate
        return result

    if target == 0:
        if gaps.min() < 0:
            raise NotImplementedError('零扭矩且有初始穿入的自平衡接触分支尚未接入')
        candidate = dict(forces_n=[0.]*n, remaining_gaps_m=gaps.tolist(), elastic_deflections_m=[0.]*n,
            theta_rad=0., normal_approach_m=None, strain_energy_j=0., active_indices=[],
            resultant_load=0., balance_residual=0., relative_balance_residual=0.,
            complementarity_residual_nm=0., contact_stationarity_residual_m=0.,
            minimum_force_n=0., minimum_remaining_gap_m=float(gaps.min()),
            certificate=dict(passed=True, energy_balance_residual_j=0., zero_load=True))
        return finish(True, 'converged_zero_load')
    feasible = np.flatnonzero(levers*target > 0)
    if not len(feasible):raise ValueError('非负接触力不能形成所给有符号目标扭矩')
    if any(levers[i]*target<=0 for i in warm):
        raise ValueError('初始活动索引的非负力不能形成目标扭矩方向')
    score = gaps[feasible]/abs(levers[feasible])
    order = np.lexsort((feasible, score))
    working = sorted(feasible[order[:min(initial_count,maximum_working_set,len(feasible))]].tolist())
    if warm:
        extra=[i for i in working if i not in warm]
        working=sorted(warm+extra[:maximum_working_set-len(warm)])
    diagnostics['initial_working_indices']=working.copy()

    def get_columns(indices):
        fresh = [int(i) for i in indices if int(i) not in cache]
        if fresh:
            if (len(cache)+len(fresh))*n*8 > maximum_column_cache_bytes:
                raise MemoryError('真实FE柔度列缓存超过显式字节上限；未组装全部矩阵或丢弃物理候选')
            values = np.asarray(compliance_columns(fresh), dtype=float)
            if values.shape != (n,len(fresh)) or not np.isfinite(values).all():
                raise ValueError('真实FE柔度列回调的形状/有限数检查失败')
            block = values[np.asarray(fresh), :]
            scale = max(float(np.abs(block).max()), 1e-300)
            if np.any(np.diag(block) <= 0) or np.max(abs(block-block.T)) > symmetry_tolerance*scale:
                raise ValueError('新访问的真实FE柔度列缺少正自影响或不满足Maxwell互等')
            if cache:
                old = list(cache)
                left = values[np.asarray(old), :]
                right = np.asarray([cache[i][fresh] for i in old])
                scale = max(scale, float(np.abs(left).max()), float(np.abs(right).max()), 1e-300)
                if np.max(abs(left-right)) > symmetry_tolerance*scale:
                    raise ValueError('跨批次真实FE柔度列不满足Maxwell互等；可能混入不同系统/索引')
            for j, index in enumerate(fresh):
                cache[index] = values[:,j].copy()
                queried.append(index)
        return np.column_stack([cache[int(i)] for i in indices])

    for outer in range(1, maximum_outer_iterations+1):
        largest = max(largest, len(working))
        try:
            columns = get_columns(working)
            full_matrix_assembled = full_matrix_assembled or len(working)==n
            restricted = columns[np.asarray(working), :]
            local = solve_contact_equilibrium(gaps[working].tolist(), restricted.tolist(),
                lever_arms_m=levers[working].tolist(), torque_nm=target,
                allow_initial_penetration=allow_initial_penetration, max_iterations=inner_max_iterations,
                force_tolerance_n=force_tolerance_n, gap_tolerance_m=gap_tolerance_m,
                balance_tolerance_nm=balance_tolerance_nm, relative_tolerance=relative_tolerance,
                symmetry_tolerance=symmetry_tolerance, minimum_normalized_pivot=minimum_normalized_pivot)
        except (ValueError, ArithmeticError, MemoryError) as error:
            return finish(False, 'restricted_system_failed', str(error))
        linear_solves += local.linear_solve_count
        force = np.zeros(n); force[np.asarray(working)] = local.forces_n
        elastic = columns @ np.asarray(local.forces_n)
        theta = local.theta_rad
        if theta is None:
            return finish(False,'restricted_angle_missing','局部QP未给出可复核的转角')
        remaining = gaps+elastic-levers*theta
        force_limit = local.diagnostics['force_acceptance_n']
        gap_limit = local.diagnostics['gap_acceptance_m']
        balance_limit = local.diagnostics['balance_acceptance']
        comp_limit = local.diagnostics['complementarity_acceptance_nm']
        active = np.flatnonzero(force > force_limit)
        resultant = float(levers @ force)
        balance = resultant-target
        comp = float(np.max(abs(force*remaining)))
        stationarity = float(np.max(abs(remaining[active]))) if len(active) else 0.
        strain = .5*float(force @ elastic)
        energy_error = float(gaps @ force+2*strain-theta*target)
        energy_limit = max(comp_limit*max(1,len(active)),
            128*np.finfo(float).eps*max(abs(float(gaps @ force)),2*abs(strain),abs(theta*target),1e-300))
        certificate = dict(minimum_force_n=float(force.min()), minimum_remaining_gap_m=float(remaining.min()),
            active_stationarity_error_m=stationarity, complementarity_residual_nm=comp,
            balance_residual_nm=balance, relative_balance_residual=abs(balance)/max(abs(target),balance_limit),
            strain_energy_j=strain, energy_balance_residual_j=energy_error,
            force_acceptance_n=force_limit, gap_acceptance_m=gap_limit,
            balance_acceptance_nm=balance_limit, complementarity_acceptance_nm=comp_limit,
            energy_acceptance_j=energy_limit, checked_candidate_count=n)
        certified = (local.converged and force.min() >= -force_limit and remaining.min() >= -gap_limit
            and stationarity <= gap_limit and abs(balance) <= balance_limit and comp <= comp_limit
            and abs(energy_error) <= energy_limit)
        certificate['passed'] = bool(certified)
        candidate = dict(forces_n=force.tolist(),remaining_gaps_m=remaining.tolist(),
            elastic_deflections_m=elastic.tolist(),theta_rad=float(theta),normal_approach_m=None,
            strain_energy_j=strain,active_indices=active.tolist(),resultant_load=resultant,
            balance_residual=balance,relative_balance_residual=certificate['relative_balance_residual'],
            complementarity_residual_nm=comp,contact_stationarity_residual_m=stationarity,
            minimum_force_n=float(force.min()),minimum_remaining_gap_m=float(remaining.min()),certificate=certificate)
        in_working = np.zeros(n,dtype=bool);in_working[np.asarray(working)] = True
        violated = np.flatnonzero((remaining < -gap_limit)&~in_working)
        history.append(dict(iteration=outer,restricted_count=len(working),active_count=len(active),
            queried_column_count=len(cache),violated_outside_working_count=len(violated),
            minimum_global_remaining_gap_m=float(remaining.min()),restricted_converged=local.converged,
            restricted_iterations=local.iterations))
        if certified:return finish(True,'converged')
        if not local.converged:return finish(False,'restricted_qp_not_converged',local.failure)
        if not len(violated):return finish(False,'global_certificate_failed','没有遗漏负约化梯度，但全候选残差证书未通过')
        preserved = active.tolist()
        if not preserved:
            preserved=[working[int(np.argmax(local.forces_n))]]
        capacity = maximum_working_set-len(preserved)
        if capacity <= 0:
            return finish(False,'working_set_limit','承载活动变量达到显式上限，仍有全域未闭合接触约束')
        selected = violated[np.lexsort((violated,remaining[violated]))[:min(add_count,capacity,len(violated))]]
        working=sorted(set(preserved)|set(map(int,selected)))
    return finish(False,'outer_iteration_limit','FE接触列生成达到显式外迭代上限，全候选KKT尚未通过')
