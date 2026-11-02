"""真实间隙曲率主轴用于压力网格布置；不以二次函数代替原几何。"""
from __future__ import annotations

import math


def rotate_plane_basis(basis,angle_rad):
    """绕参考外法向旋转正交投影坐标；两行仍为单位切向基。"""
    import numpy as np
    if isinstance(angle_rad,bool) or not isinstance(angle_rad,(int,float)) or not math.isfinite(angle_rad):
        raise ValueError('投影网格旋转角须为有限rad')
    value=np.asarray(basis,dtype=float)
    if value.shape!=(2,3) or not np.isfinite(value).all() or not np.allclose(value@value.T,np.eye(2),rtol=0,atol=1e-12):
        raise ValueError('投影基须为正交单位2×3行向量')
    c=math.cos(angle_rad);s=math.sin(angle_rad)
    return np.array([[c,s],[-s,c]])@value


def principal_gap_direction(sample_gap,*,step_m=1e-4,maximum_angle_difference_rad=.01,
        maximum_relative_curvature_difference=.05):
    """用两档中心差分求真实间隙Hessian，长接触方向取小正曲率主轴。

    sample_gap(x,y)直接查询实际双齿面的有符号法向间隙m。差分仅选择
    网格方向；后续每一个中心、角点、Gauss点仍查询真实曲面。有限面
    stencil失败、Hessian非正定、方向不稳定均明确拒绝自动选择。
    """
    import numpy as np
    for value,name in ((step_m,'差分步长m'),(maximum_angle_difference_rad,'方向差容限rad'),
            (maximum_relative_curvature_difference,'曲率相对差容限')):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
            raise ValueError(name+'须为正有限数')
    samples={}
    def gap(x,y):
        key=(x,y)
        if key not in samples:
            value=sample_gap(x,y)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ArithmeticError('真实间隙差分点未求解，不能生成主轴')
            samples[key]=float(value)
        return samples[key]
    studies=[]
    for h in (step_m,step_m*.5):
        center=gap(0.,0.)
        hxx=(gap(h,0.)-2*center+gap(-h,0.))/h**2
        hyy=(gap(0.,h)-2*center+gap(0.,-h))/h**2
        hxy=(gap(h,h)-gap(h,-h)-gap(-h,h)+gap(-h,-h))/(4*h**2)
        matrix=np.array([[hxx,hxy],[hxy,hyy]]);values,vectors=np.linalg.eigh(matrix)
        if values[0]<=0 or not np.isfinite(values).all():raise ArithmeticError('真实间隙局部Hessian不是正定，自动主轴布网不适用')
        if (values[1]-values[0])/values[1]<.05:raise ArithmeticError('真实间隙主曲率接近，主轴方向不唯一')
        vector=vectors[:,0]
        angle=(math.atan2(float(vector[1]),float(vector[0]))+math.pi/2)%math.pi-math.pi/2
        studies.append(dict(step_m=h,hessian_1_m=matrix.tolist(),principal_curvatures_1_m=values.tolist(),angle_rad=angle))
    first,last=studies;delta=abs((last['angle_rad']-first['angle_rad']+math.pi/2)%math.pi-math.pi/2)
    curvature_error=float(np.max(np.abs(np.asarray(first['principal_curvatures_1_m'])-last['principal_curvatures_1_m'])/np.asarray(last['principal_curvatures_1_m'])))
    if delta>maximum_angle_difference_rad or curvature_error>maximum_relative_curvature_difference:
        raise ArithmeticError('真实间隙主轴的差分加密检查未通过')
    return dict(angle_rad=last['angle_rad'],studies=studies,angle_refinement_error_rad=delta,
        relative_curvature_refinement_error=curvature_error,actual_gap_queries=len(samples),
        purpose='仅用于真实曲面的压力网格主轴；不替代间隙、不拟合压力、不改变材料几何')


def generated_patch_gap_direction(prepared,patch,*,step_m=1e-4):
    """从实际有限NURBS面查询主轴；输入patch必须为明确的参考姿态。"""
    import numpy as np
    from .generated_halfspace import _sample_pair,_initial_uv
    frames=patch['frames'];origin=np.asarray(patch['plane_origin_world_m']);basis=np.asarray(patch['plane_tangent_rows'])
    normal=np.asarray(patch['plane_outward_normal_world']);seeds=patch['seed_uv']
    maps={m:_initial_uv(prepared.trimmed[m].surface,np.asarray(frames[m+'_rotation_rows']),origin,basis,seeds[m]) for m in ('pinion','wheel')}
    def sample(x,y):
        value,reason=_sample_pair(prepared,frames,origin,basis,normal,seeds,(x,y),1e-13,.3,.3,maps)
        if value is None:raise ArithmeticError('实际有限齿面主轴差分查询失败：'+str(reason))
        return value['gap_m']
    result=principal_gap_direction(sample,step_m=step_m)
    result.update(tooth_offset=patch['tooth_offset'],reference_frames=frames)
    return result


def build_principal_aligned_patches(prepared,driving_angle_rad,pinion_torque_nm,*,step_m=1e-4,
        fallback_to_surface_basis=False,**parameters):
    """每个新相位先查询真实间隙主轴，再按固定旋转方向建立原压力格。

    无法求主轴时默认报错。显式开启fallback仅保留原曲面投影基并逐齿
    记录原因；不会删掉该齿、压力格或未知区域，也不构成接触排除证据。
    固定方向进入build_parameters，受载外迭代因此复用同一网格约定。
    """
    from .generated_halfspace import build_generated_halfspace_patches
    if not isinstance(fallback_to_surface_basis,bool):raise ValueError('主轴失败后保留曲面基的开关须为布尔值')
    if 'plane_basis_angles_rad' in parameters:raise ValueError('自动主轴和显式投影角不能同时提供')
    # 方向探针不生成压力自由度。即使相切点贴近有限面边缘，也先保留
    # 真实中心供方向查询；正式网格仍遵守调用者的完整四角域检查。
    probe=dict(parameters,grid_shape=(3,3),spacing_m=(step_m,step_m),cover_complete_working_flank=False,
        check_cell_corners=False)
    initial=build_generated_halfspace_patches(prepared,driving_angle_rad,pinion_torque_nm,**probe)
    directions=[];angles={}
    for patch in initial['patches']:
        try:
            direction=generated_patch_gap_direction(prepared,patch,step_m=step_m)
            directions.append(dict(accepted=True,**direction));angles[patch['tooth_offset']]=direction['angle_rad']
        except ArithmeticError as error:
            if not fallback_to_surface_basis:raise
            angles[patch['tooth_offset']]=0.
            directions.append(dict(tooth_offset=patch['tooth_offset'],accepted=False,angle_rad=0.,
                fallback='显式保留原曲面投影基；该齿全部原候选仍保留',reason=str(error)))
    full=dict(parameters,tca_results=initial['tca_results'],plane_basis_angles_rad=angles)
    result=build_generated_halfspace_patches(prepared,driving_angle_rad,pinion_torque_nm,**full)
    result['grid_alignment']=dict(method='actual_gap_hessian_principal_axes',directions=directions,
        step_m=step_m,fallback_to_surface_basis=fallback_to_surface_basis,
        pressure_or_geometry_replaced=False)
    return result
