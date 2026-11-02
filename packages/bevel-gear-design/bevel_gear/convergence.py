"""按相同物理输入评估直齿FE网格序列；数值平衡与离散收敛分开。"""
from __future__ import annotations

import json
import math


def assess_straight_ltca_convergence(levels, *, te_relative_tolerance=0.05,
        stiffness_relative_tolerance=0.05, pressure_relative_tolerance=0.1,
        te_absolute_tolerance_rad=1e-8, minimum_levels=3):
    """比较每对相邻加密级；不把最后两级接近解释为严格误差上界。

    levels每项含contact_size_m、case_controls、result。case_controls须
    为除局部网格尺寸外完全相同的几何/材料/载荷/接触/边界配置；调用
    方负责将真实配置随结果保存，不能用一个相同标签代替输入对应。
    当前检查限连续折线加密模式，主体尺寸、CAD精度、接触姿态均固定。
    TE为弹性从动滞后rad，刚度N·m/rad，压力为面片平均值Pa。
    """
    if isinstance(minimum_levels, bool) or not isinstance(minimum_levels, int) or minimum_levels < 3:
        raise ValueError('网格检查至少需要三级，minimum_levels须为不小于3的整数')
    for value in (te_relative_tolerance, stiffness_relative_tolerance, pressure_relative_tolerance,
                  te_absolute_tolerance_rad):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value <= 0:
            raise ValueError('各收敛门槛须为正有限数')
    samples=list(levels)
    if len(samples)<minimum_levels:
        raise ValueError('完成的网格级数不足，不能评估收敛')
    fields=(('elastic_follower_lag_rad',te_relative_tolerance,te_absolute_tolerance_rad),
            ('incremental_mesh_stiffness_nm_per_rad',stiffness_relative_tolerance,0.),
            ('maximum_cell_mean_pressure_pa',pressure_relative_tolerance,0.))
    common=None; last_h=math.inf
    rows=[]
    for sample in samples:
        h=sample['contact_size_m']
        if isinstance(h,bool) or not isinstance(h,(int,float)) or not math.isfinite(h) or not 0<h<last_h:
            raise ValueError('加密尺寸须按正有限值严格递减，不自动重排或合并重复级')
        if not isinstance(sample['case_controls'],dict) or not sample['case_controls']:
            raise ValueError('须保存实际case_controls，不能只提供问题标签')
        if sample['case_controls'].get('refinement_geometry')!='continuous_polyline':
            raise ValueError('当前收敛检查只接受明确的continuous_polyline控制场')
        controls=json.dumps(sample['case_controls'],sort_keys=True,allow_nan=False)
        if common is None:common=controls
        elif controls!=common:raise ValueError('物理/CAD/候选或支承配置不同，不能混为同一网格序列')
        result=sample['result']
        if not result.get('converged'):
            raise ValueError('有网格级尚未通过离散平衡，不能评估其收敛')
        preparation=result['preparation']
        if not preparation['full_gears']:
            raise ValueError('当前检查限整轮，不能把自由侧面的扇区混入')
        if preparation['contact_refined_size_m']!=h:
            raise ValueError('给定尺寸与实际结果的加密元数据不一致')
        row=dict(contact_size_m=h,driving_angle_rad=result['driving_angle_rad'],
                 pinion_torque_nm=result['pinion_torque_nm'])
        for field,_,_ in fields:
            value=result[field]
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
                raise ValueError('所评估TE、刚度、面片平均压力须为正有限数')
            row[field]=value
        if rows and (row['driving_angle_rad'],row['pinion_torque_nm'])!=(rows[0]['driving_angle_rad'],rows[0]['pinion_torque_nm']):
            raise ValueError('所测姿态或实际扭矩不同')
        rows.append(row);last_h=h
    comparisons=[]
    for coarse,fine in zip(rows,rows[1:]):
        changes={}
        for field,relative,absolute in fields:
            difference=fine[field]-coarse[field]
            allowed=max(absolute,relative*abs(fine[field]))
            changes[field]=dict(difference=difference,relative_change=difference/abs(fine[field]),
                absolute_threshold=allowed,passed=abs(difference)<=allowed)
        comparisons.append(dict(coarse_size_m=coarse['contact_size_m'],fine_size_m=fine['contact_size_m'],
            metrics=changes,passed=all(v['passed'] for v in changes.values())))
    return dict(checked_levels=len(rows),rows=rows,comparisons=comparisons,
        passed=all(c['passed'] for c in comparisons),
        pressure_peak_validated=False,physical_validation=False,
        interpretation='所测连续加密序列的相邻变化检查；非Richardson误差估计、非物理验证',
        limitations=['当前CAD和主体网格保持固定，不能据局部序列推断全部离散误差',
            '面片平均压力最大值不是Hertz峰值；通过门槛不自动使pressure_peak_validated变为真',
            '接触面积/活跃集可变，未估计渐近收敛阶或完整周期精度'])
