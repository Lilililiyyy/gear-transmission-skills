"""理想球面渐开线直齿副的零间隙相切TCA，新写解析几何实现。

联立两齿面位置重合与外法向反向，不用正间隙距离带代替相切。
仅共同锥顶、无装配误差的理想直齿；不等同于octoid齿形或实际
加工工艺齿面。角度rad、点坐标m，弹性接触载荷不属于本模块。
"""
from __future__ import annotations

import math

from .geometry import GeometryResult
from .tooth_surface import build_spherical_tooth_surface


def _numeric():
    try:
        import numpy as np
        from scipy.optimize import least_squares
    except ImportError as error:
        raise ImportError("解析直齿相切TCA需要可选NumPy/SciPy依赖") from error
    return np, least_squares


def straight_pair_frames(geometry, driving_angle_rad, follower_angle_rad, *, tooth_offset=0, flank="right"):
    """返回单齿局部到共同世界的保向矩阵；两轮同名侧构成工作面。

    小轮轴为+Z，大轮轴在XZ平面；大轮切向基为-Y，故不能将左右侧
    标签直接套入另一套装配约定。初相位按节锥节点使所选工作侧相切。
    tooth_offset在主动为+k齿、从动为−k齿，增量角是整轮轴角。
    """
    np, _ = _numeric()
    if not isinstance(geometry, GeometryResult):
        raise TypeError("直齿副入口要求GeometryResult")
    p, w = [build_spherical_tooth_surface(geometry,m) for m in ("pinion","wheel")]
    if flank not in ("left","right"):
        raise ValueError("工作侧须为left或right")
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):
        raise ValueError("配对齿偏移须为整数")
    for angle in (driving_angle_rad,follower_angle_rad):
        if isinstance(angle,bool) or not isinstance(angle,(int,float)) or not math.isfinite(angle):
            raise ValueError("两轮增量角须为有限弧度")
    sign=1 if flank=="right" else -1
    def rz(a):
        c,s=math.cos(a),math.sin(a)
        return np.array(((c,-s,0.),(s,c,0.),(0.,0.,1.)))
    sigma=geometry.inputs.shaft_angle_rad
    c,s=math.cos(sigma),math.sin(sigma)
    frame=np.array(((-c,0.,s),(0.,-1.,0.),(s,0.,c)))
    a=rz(-sign*p.profile.pitch_angular_tooth_thickness_rad/2 + driving_angle_rad + tooth_offset*2*math.pi/p.teeth)
    b=frame@rz(-sign*w.profile.pitch_angular_tooth_thickness_rad/2 + follower_angle_rad - tooth_offset*2*math.pi/w.teeth)
    return dict(pinion_rotation_rows=a.tolist(),wheel_rotation_rows=b.tolist(),
        pinion_axis_world=(0.,0.,1.),wheel_axis_world=(s,0.,c),origin_m=(0.,0.,0.),
        initial_pinion_angle_rad=-sign*p.profile.pitch_angular_tooth_thickness_rad/2,
        initial_wheel_angle_rad=-sign*w.profile.pitch_angular_tooth_thickness_rad/2,
        flank=flank,tooth_offset=tooth_offset,
        convention="列向量world=R@local；共同锥顶，无平移或安装失准")


def solve_straight_tca(geometry, driving_angle_rad, *, tooth_offset=0, flank="right",
                       contact_line_samples=9, cone_distance_bounds_m=None,
                       residual_tolerance=1e-9, maximum_evaluations=180):
    """在有效齿廓域求解位置/法向联合残差，返回相切或明确的未啮合状态。

    以两极角和从动增量为未知量。位置先除球面半径，与单位法向共同
    构造无量纲残差。沿母线的有限接触段另检查两轮边界，不把无界
    锥面相切自动当作整个实际齿宽都在接触。
    """
    np, least_squares = _numeric()
    if not isinstance(geometry,GeometryResult):
        raise TypeError("TCA需要已计算的GeometryResult")
    p,w=[build_spherical_tooth_surface(geometry,m) for m in ("pinion","wheel")]
    straight_pair_frames(geometry,driving_angle_rad,0.,tooth_offset=tooth_offset,flank=flank)
    if isinstance(contact_line_samples,bool) or not isinstance(contact_line_samples,int) or not 2<=contact_line_samples<=501:
        raise ValueError("接触线采样数须为2至501的整数")
    if isinstance(maximum_evaluations,bool) or not isinstance(maximum_evaluations,int) or maximum_evaluations<1:
        raise ValueError("maximum_evaluations须为正整数")
    if not isinstance(residual_tolerance,(int,float)) or isinstance(residual_tolerance,bool) or not math.isfinite(residual_tolerance) or residual_tolerance<=0:
        raise ValueError("相切无量纲残差门槛须为正有限数")
    radius=geometry.mean_cone_distance_m
    low1,high1=p.profile_bounds(radius)
    low2,high2=w.profile_bounds(radius)
    # 基锥处导数奇异，仅留极小的开放域余量；不修改齿面本身。
    low1=max(low1,p.profile.base_cone_angle_rad+1e-9)
    low2=max(low2,w.profile.base_cone_angle_rad+1e-9)
    nominal=-driving_angle_rad*p.teeth/w.teeth
    half_pitch=math.pi/w.teeth
    def evaluate(x, details=False):
        delta1,delta2,follower=x
        frames=straight_pair_frames(geometry,driving_angle_rad,float(follower),tooth_offset=tooth_offset,flank=flank)
        a,b=np.array(frames['pinion_rotation_rows']),np.array(frames['wheel_rotation_rows'])
        pp=a@p.profile.evaluate(1.,delta1,flank)
        pw=b@w.profile.evaluate(1.,delta2,flank)
        n1=a@p.profile.normal(delta1,flank)
        n2=b@w.profile.normal(delta2,flank)
        if details:return pp,pw,n1,n2,frames
        return np.r_[pp-pw,n1+n2]
    initial=[min(high1,max(low1,p.profile.pitch_cone_angle_rad)),
             min(high2,max(low2,w.profile.pitch_cone_angle_rad)),nominal]
    fit=least_squares(evaluate,initial,bounds=([low1,low2,nominal-half_pitch],[high1,high2,nominal+half_pitch]),
        xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=maximum_evaluations)
    pp,pw,n1,n2,frames=evaluate(fit.x,True)
    position_residual=float(np.linalg.norm(pp-pw)*radius)
    normal_residual=float(np.linalg.norm(n1+n2))
    angular_residual=float(np.linalg.norm(fit.fun))
    converged=bool(fit.success and angular_residual<=residual_tolerance)
    if cone_distance_bounds_m is None:
        lower,upper=geometry.inner_cone_distance_m,geometry.outer_cone_distance_m
    else:
        if len(cone_distance_bounds_m)!=2:
            raise ValueError("球面锥距边界须含两个SI长度")
        lower,upper=cone_distance_bounds_m
        if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in (lower,upper)) or not 0<lower<upper:
            raise ValueError("球面锥距须满足0<下界<上界，单位m")
    line=[]
    for r in np.linspace(lower,upper,contact_line_samples):
        within=True
        for surface,delta in ((p,fit.x[0]),(w,fit.x[1])):
            # 默认球面端界与CAD轴向端界不同；显式CAD边界由调用者给定。
            try:
                if cone_distance_bounds_m is None:
                    lo,hi=surface.profile_bounds(float(r))
                else:
                    root=surface._boundary_polar_angle(float(r),surface.mean_dedendum_m,surface.dedendum_taper_angle_rad,-1)
                    tip=surface._boundary_polar_angle(float(r),surface.mean_addendum_m,surface.addendum_taper_angle_rad,1)
                    lo,hi=max(root,surface.profile.base_cone_angle_rad),tip
                within=within and lo-1e-10<=delta<=hi+1e-10
            except ValueError:within=False
        if converged and within:
            line.append(dict(cone_distance_m=float(r),point_world_m=(.5*(pp+pw)*r).tolist(),
                normal_world=n1.tolist(),position_residual_m=float(np.linalg.norm(pp-pw)*r)))
    if converged and len(line)<2:
        converged=False
        status='outside_finite_tooth'
    else:status='tangent' if converged else ('evaluation_limit' if not fit.success else 'outside_active_profile')
    normal=.5*(n1-n2)
    normal=normal/np.linalg.norm(normal)
    point=.5*(pp+pw)*radius
    lever1=float(normal@np.cross(np.array(frames['pinion_axis_world']),point))
    lever2=float(-normal@np.cross(np.array(frames['wheel_axis_world']),point))
    return dict(converged=converged,status=status,driving_angle_rad=driving_angle_rad,
        follower_angle_rad=float(fit.x[2]),nominal_follower_angle_rad=nominal,
        unloaded_transmission_error_rad=float(fit.x[2]-nominal) if converged else None,
        pinion_polar_angle_rad=float(fit.x[0]),wheel_polar_angle_rad=float(fit.x[1]),
        position_residual_m=position_residual,opposing_normal_residual=normal_residual,
        dimensionless_residual=angular_residual,function_evaluations=int(fit.nfev),
        mean_contact_point_world_m=point.tolist() if converged else None,
        mean_contact_normal_world=normal.tolist() if converged else None,
        pinion_normal_lever_m=lever1,wheel_normal_lever_m=lever2,
        contact_line=line if converged else [],frames=frames,
        model='理想球面渐开线直齿；位置重合与法向相切；有限边界另检验',
        limitations=['共同锥顶、无安装失准；不代表octoid或实际切削齿形',
            '未求弹性载荷、压力或材料变形；残差不是制造精度',
            '接触线是解析锥面相切母线的有效采样，非承载接触斑'])
