"""实际刀盘展成工作面的局部相切TCA；独立求位置与法向，不用正距离代替。

输入为显式双轮制造参数。只支持共同锥顶、零安装平移、主动左/从动右
工作面；不是全齿面侵入排除、全周期接触域或承载接触分析。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .generated_surface import _request,build_generated_surface_grid
from .generated_cutter import evaluate_generated_point


@dataclass
class PreparedGeneratedTca:
    """显式双轮机床输入及独立NURBS工作面，CAD内部长度mm。"""
    requests: dict
    shells: dict
    initial_angles_rad: dict
    shaft_angle_rad: float
    cad_models: dict | None
    diagnostics: dict


def _numeric():
    try:
        import numpy as np
        from scipy.optimize import least_squares
    except ImportError as error:
        raise ImportError('展成曲面相切需要可选NumPy/SciPy依赖') from error
    return np,least_squares


def _initial(request,side,opposite):
    machine=request.machine
    station=machine.reference_cone_distance_m-request.face_width_m/2
    diameter=request.tooth_count*request.module_m*station/machine.reference_cone_distance_m
    point=evaluate_generated_point(machine,diameter,station)['flanks'][side]['point_m']
    return math.atan2(point[1],point[0])+(math.pi if opposite else 0.)


def prepare_generated_tca(pinion_request,wheel_request,*,width_count=25,profile_count=25,cad_models=None):
    """建立实际工作面；可复用同制造输入模型的shell，避免重建。

    未自动设计双轮机床设置。参考锥距和节锥角和只验证共同锥顶条件，
    不据此保证输入工艺构成共轭齿面。模型几何和宏观标准评级分别使用。
    """
    from .cad import build_generated_shell
    requests=dict(pinion=_request(pinion_request),wheel=_request(wheel_request))
    first,second=requests.values()
    if not math.isclose(first.machine.reference_cone_distance_m,second.machine.reference_cone_distance_m,rel_tol=1e-9,abs_tol=1e-12):
        raise ValueError('本相切入口要求双轮共享参考锥距，不支持偏置或安装平移')
    shaft=first.machine.pitch_cone_angle_rad+second.machine.pitch_cone_angle_rad
    if not 0<shaft<math.pi:
        raise ValueError('两节锥角和须在0至pi之间')
    shells={}
    if cad_models is not None and (not isinstance(cad_models,dict) or set(cad_models)!=set(requests)):
        raise ValueError('复用CAD须为pinion/wheel两模型字典')
    for member,request in requests.items():
        if cad_models is None:
            shells[member]=build_generated_shell(build_generated_surface_grid(request,width_count=width_count,profile_count=profile_count))
        else:
            model=cad_models[member]
            if getattr(model,'request',None)!=request:
                raise ValueError('复用CAD的实际制造输入与TCA输入不同')
            shells[member]=model.shell
    initial=dict(pinion=_initial(first,'left',False),wheel=_initial(second,'right',True))
    return PreparedGeneratedTca(requests,shells,initial,shaft,cad_models,
        dict(units='m/rad',cad_internal_length_unit='mm',working_patches=['left_flank','right_flank'],
             installation='共同锥顶，制造局部框架显式换轴，无安装平移',cad_reused=cad_models is not None))


def generated_tca_frames(prepared,driving_angle_rad,follower_angle_rad,*,tooth_offset=0):
    """制造局部坐标至公共坐标旋转矩阵，含初装与同一齿对的齿距索引。"""
    np,_=_numeric()
    if not isinstance(prepared,PreparedGeneratedTca):raise TypeError('需要PreparedGeneratedTca')
    if isinstance(tooth_offset,bool) or not isinstance(tooth_offset,int):raise ValueError('齿对偏移须为整数')
    for value in (driving_angle_rad,follower_angle_rad):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError('主从动增量角须为有限弧度数')
    def rz(angle):
        c,s=math.cos(angle),math.sin(angle)
        return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
    c,s=math.cos(-prepared.shaft_angle_rad),math.sin(-prepared.shaft_angle_rad)
    rx=np.array([[1.,0.,0.],[0.,c,-s],[0.,s,c]])
    q=np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,-1.]])
    first,second=prepared.requests.values()
    a=rx@rz(driving_angle_rad+prepared.initial_angles_rad['pinion']+tooth_offset*2*math.pi/first.tooth_count)@q
    b=rz(follower_angle_rad+prepared.initial_angles_rad['wheel']-tooth_offset*2*math.pi/second.tooth_count)@q
    return dict(pinion_rotation_rows=a.tolist(),wheel_rotation_rows=b.tolist(),
        pinion_axis_world=(rx@np.array([0.,0.,1.])).tolist(),wheel_axis_world=[0.,0.,1.],
        translation_m=[0.,0.,0.])


def _surface_value(surface,u,v):
    np,_=_numeric()
    from OCP.gp import gp_Pnt,gp_Vec
    u0,u1,v0,v1=surface.Bounds()
    point=gp_Pnt();a=gp_Vec();b=gp_Vec()
    surface.D1(u0+u*(u1-u0),v0+v*(v1-v0),point,a,b)
    normal=np.cross(a.Coord(),b.Coord())
    magnitude=float(np.linalg.norm(normal))
    if magnitude<=1e-25 or not math.isfinite(magnitude):raise ArithmeticError('工作面局部参数导数退化')
    return np.array(point.Coord())*.001,normal/magnitude


def solve_generated_tca(prepared,driving_angle_rad,*,tooth_offset=0,position_tolerance_m=1e-9,
        normal_tolerance=1e-7,maximum_evaluations=200,seeds=(.25,.5,.75)):
    """五未知量u1/v1/u2/v2/从动角，联立三位置与法向平行条件。

    六残差中法向叉积有两个独立约束；最后另检验外法向反向与参数域
    内点。没有强制从动齿数比、零TE或把小残差当作全域无干涉证明。
    未复用闭实体时法向仅按工作面参数约定，结果显式未做材料内外检查。
    """
    np,least_squares=_numeric()
    if not isinstance(prepared,PreparedGeneratedTca):raise TypeError('需要PreparedGeneratedTca')
    for value in (position_tolerance_m,normal_tolerance):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
            raise ValueError('位置与法向容差须为正有限数')
    if isinstance(maximum_evaluations,bool) or not isinstance(maximum_evaluations,int) or maximum_evaluations<1:
        raise ValueError('每个初值的迭代评价上限须为正整数')
    seeds=list(seeds)
    if not seeds or any(isinstance(s,bool) or not isinstance(s,(int,float)) or not 0<s<1 for s in seeds):
        raise ValueError('宽向初值须为0和1之间的非空有限数列')
    first,second=prepared.requests.values()
    # 先通过统一框架校验主动角/索引，再建立有界从动相位窗。
    generated_tca_frames(prepared,driving_angle_rad,0.,tooth_offset=tooth_offset)
    nominal=-driving_angle_rad*first.tooth_count/second.tooth_count
    radius=first.machine.reference_cone_distance_m
    half_window=.25*2*math.pi/second.tooth_count
    surfaces=[prepared.shells[m].surfaces[p] for m,p in (('pinion','left_flank'),('wheel','right_flank'))]
    def evaluate(x):
        frames=generated_tca_frames(prepared,driving_angle_rad,float(x[4]),tooth_offset=tooth_offset)
        a,b=(np.array(frames[m+'_rotation_rows']) for m in ('pinion','wheel'))
        p,n=_surface_value(surfaces[0],*x[:2]);w,m=_surface_value(surfaces[1],*x[2:4])
        return a@p,b@w,a@n,b@m
    residual_evaluations=0
    def residual(x):
        nonlocal residual_evaluations
        residual_evaluations+=1
        p,w,n,m=evaluate(x)
        return np.r_[(p-w)/radius,np.cross(n,m)]
    trials=[]
    lower=[1e-9]*4+[nominal-half_window];upper=[1-1e-9]*4+[nominal+half_window]
    for u in seeds:
        trial=least_squares(residual,[u,.5,u,.5,nominal],bounds=(lower,upper),
            xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=maximum_evaluations)
        trials.append(trial)
    best=min(trials,key=lambda t:float(t.cost))
    p,w,n,m=evaluate(best.x)
    gap=float(np.linalg.norm(p-w));cross=float(np.linalg.norm(np.cross(n,m)))
    dot=float(n@m);boundary=bool(np.any(np.minimum(best.x[:4],1-best.x[:4])<1e-6))
    normal_checks={}
    if prepared.cad_models is not None:
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.TopAbs import TopAbs_IN,TopAbs_OUT
        from OCP.gp import gp_Pnt
        normals=[]
        for member,surface,uv in (('pinion',surfaces[0],best.x[:2]),('wheel',surfaces[1],best.x[2:4])):
            local,normal=_surface_value(surface,*uv)
            model=prepared.cad_models[member]
            states=[]
            for sign in (1,-1):
                point=(local+sign*5e-7*normal)*1000
                states.append(BRepClass3d_SolidClassifier(model.shape,gp_Pnt(*point),1e-7).State())
            if states==[TopAbs_OUT,TopAbs_IN]:factor=1
            elif states==[TopAbs_IN,TopAbs_OUT]:factor=-1
            else:factor=None
            normal_checks[member]=dict(classified=factor is not None,parameter_normal_to_outward=factor)
            normals.append(factor)
        if all(v is not None for v in normals):dot*=normals[0]*normals[1]
    oriented=(not normal_checks or all(v['classified'] for v in normal_checks.values())) and dot<0
    converged=gap<=position_tolerance_m and cross<=normal_tolerance and oriented and not boundary
    status='tangent' if converged else ('boundary_contact_unresolved' if boundary else 'tangency_not_satisfied')
    frames=generated_tca_frames(prepared,driving_angle_rad,float(best.x[4]),tooth_offset=tooth_offset)
    candidate=dict(parameters=best.x.tolist(),follower_increment_angle_rad=float(best.x[4]),
        unloaded_te_rad=float(best.x[4]-nominal),pinion_contact_point_world_m=p.tolist(),
        wheel_contact_point_world_m=w.tolist(),frames=frames)
    return dict(converged=converged,status=status,driving_angle_rad=driving_angle_rad,tooth_offset=tooth_offset,
        follower_increment_angle_rad=float(best.x[4]) if converged else None,nominal_follower_increment_angle_rad=nominal,
        unloaded_te_rad=float(best.x[4]-nominal) if converged else None,position_residual_m=gap,normal_cross_residual=cross,
        opposing_normal_dot=dot,pinion_uv=best.x[:2].tolist() if converged else None,wheel_uv=best.x[2:4].tolist() if converged else None,
        pinion_contact_point_world_m=p.tolist() if converged else None,wheel_contact_point_world_m=w.tolist() if converged else None,
        frames=frames if converged else None,candidate=candidate,
        normal_material_checks=normal_checks,solid_interference_checked=False,
        optimizer_success=bool(best.success),total_residual_evaluations=residual_evaluations,
        optimizer_reported_nfev=sum(int(t.nfev) for t in trials),
        jacobian_singular_values=np.linalg.svd(best.jac,compute_uv=False).tolist(),
        definition='实际NURBS展成面的位置/法向局部相切；从动相位为未知量',
        limitations=['只检查所选工作面与齿对；没有穷尽其他局部相切或全轮侵入',
            'NURBS插值与机床输入的适用性须另验；不是已收敛制造误差保证',
            '未加载，不含接触压力、柔度、材料或齿根应力'])
