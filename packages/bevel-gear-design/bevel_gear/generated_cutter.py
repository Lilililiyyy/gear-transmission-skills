"""双面刀盘展成包络：直刃与刀尖圆弧共同参与逐点高度搜索。

机床设置显式输入，公开接口使用m/rad；内层保留mm公式及1e-7 mm停止条件。
本模块给出指定机床模型的两侧包络，不自动完成机床设置设计、成对接触或LTCA。
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class CutterMachineSettings:
    """单轮机床几何；设置量不是轴偏置或宏观齿轮参数的同义字段。"""
    reference_cone_distance_m: float
    dedendum_m: float
    vertical_setting_m: float
    axial_setting_m: float
    horizontal_setting_m: float
    machine_root_angle_rad: float
    pitch_cone_angle_rad: float
    roll_ratio: float
    cradle_horizontal_m: float
    cradle_vertical_m: float
    direction: int
    cutter_diameter_m: float
    blade_point_width_m: float
    inner_pressure_angle_rad: float
    outer_pressure_angle_rad: float
    inner_tip_radius_m: float
    outer_tip_radius_m: float

    def __post_init__(self):
        for key,value in self.__dict__.items():
            if not isinstance(value,(float,int)) or isinstance(value,bool) or not math.isfinite(value):
                raise ValueError(f'{key}须为有限实数')
        if self.direction not in (-1,1):raise ValueError('direction须为−1或1')
        for key in ('reference_cone_distance_m','roll_ratio','cutter_diameter_m','cradle_horizontal_m'):
            if getattr(self,key)<=0:raise ValueError(f'{key}须为正')
        for key in ('dedendum_m','blade_point_width_m','inner_tip_radius_m','outer_tip_radius_m'):
            if getattr(self,key)<0:raise ValueError(f'{key}不得为负')
        for key in ('machine_root_angle_rad','pitch_cone_angle_rad','inner_pressure_angle_rad','outer_pressure_angle_rad'):
            if not 0<getattr(self,key)<math.pi/2:raise ValueError(f'{key}须在(0,π/2)内')
        if self.cradle_vertical_m<0:raise ValueError('摇台竖直量须为非负幅值')

    def _mm(self):
        return {key:(value*1000 if key.endswith('_m') else value) for key,value in self.__dict__.items()}


def _acos(value):
    """只容纳浮点舍入越界，不把实际无交点改造成交点。"""
    if not math.isfinite(value) or abs(value)>1+2e-12:
        raise ValueError('刀具与所求截面没有有效实数交点，请检查机床设置和直径范围')
    return math.acos(max(-1.,min(1.,value)))


def _coordinates_mm(xcenter,zcenter,diameter,station,pitch,setting,height):
    """工件截面与刀盘局部坐标；保留分象限旋转及原离散常数。"""
    radius=diameter*.5
    length=math.hypot(xcenter,zcenter)
    quarter=1.5707963
    angle=(quarter if zcenter>=0 else -quarter) if xcenter==0 else math.atan(zcenter/xcenter)
    angle+=setting
    x=z=0.
    if zcenter<=0 and xcenter>=0:
        x=math.cos(angle)*length;z=math.sin(angle)*length
    if zcenter<0 and xcenter<0:
        if angle>=quarter:
            x=math.cos(3.1415926-angle)*length;z=-math.sin(3.1415926-angle)*length
        else:
            x=-math.cos(angle)*length;z=-math.sin(angle)*length
    if zcenter>=0 and xcenter<0:
        x=-math.cos(angle)*length;z=-math.sin(angle)*length
    if zcenter>0 and xcenter>=0:
        if angle>=quarter:
            angle=3.1415926-angle;x=-math.cos(angle)*length
        else:x=math.cos(angle)*length
        z=math.sin(angle)*length
    x+=(radius-(radius/math.cos(pitch)-math.tan(pitch)*station)/math.cos(pitch))/math.tan(pitch)
    return (height/math.cos(setting)+math.tan(setting)*x-z,
            (radius-math.tan(setting)*x+z)*math.cos(setting),x)


def _tip_correction_mm(angle,radius,height):
    """刀尖圆弧与直刃相切高度及径向回退量，不是经验圆角插值。"""
    limit=radius*(1-math.sin(angle))
    if height>=limit:return 0.,limit
    if radius<=0:raise ValueError('负高度处不能用零刀尖圆角生成齿根')
    effective=max(0.,height)
    arc=_acos((radius-effective)/radius)
    return radius*(math.cos(angle)-math.sin(arc))-math.tan(angle)*(limit-effective),limit


def _angles_at_height_mm(settings,xcenter,zcenter,diameter,station,height):
    s=settings
    setting=s['machine_root_angle_rad'];pitch=s['pitch_cone_angle_rad']
    q,_,axial=_coordinates_mm(xcenter,zcenter,diameter,station,pitch,setting,height)
    radicand=diameter*diameter*.25-q*q
    if q>diameter*.5:cross=0.
    elif radicand<0:raise ValueError('机床截面离开有效圆柱截交范围')
    else:cross=math.sqrt(radicand)
    denominator=axial/math.cos(setting)+math.tan(setting)*height
    if denominator==0:raise ValueError('刀盘局部截面极角分母为零')
    psi=math.atan(cross/denominator)
    roll_radius=denominator if abs(psi)<1e-16 else cross/math.sin(psi)
    correction_in,_=_tip_correction_mm(s['inner_pressure_angle_rad'],s['inner_tip_radius_m'],height)
    correction_out,_=_tip_correction_mm(s['outer_pressure_angle_rad'],s['outer_tip_radius_m'],height)
    inner=s['cutter_diameter_m']-2*(math.tan(s['inner_pressure_angle_rad'])*height+s['blade_point_width_m']*.5)+2*correction_in
    outer=s['cutter_diameter_m']+2*(math.tan(s['outer_pressure_angle_rad'])*height+s['blade_point_width_m']*.5)-2*correction_out
    if min(inner,outer)<=0:raise ValueError('刀尖截面出现非正有效刀盘直径')
    center_radius=math.hypot(s['cradle_horizontal_m'],s['cradle_vertical_m'])
    theta=math.atan(s['cradle_vertical_m']/s['cradle_horizontal_m'])
    offset=s['vertical_setting_m'];direction=s['direction'];ratio=s['roll_ratio']
    radii=[];phases=[]
    for first in (True,False):
        included=(math.pi/2-psi if offset<0 else psi+math.pi/2) if first else (psi+math.pi/2 if offset<0 else math.pi/2-psi)
        radius=math.sqrt(offset*offset+roll_radius*roll_radius-2*offset*roll_radius*math.cos(included))
        phase=(psi if first else -psi) if abs(offset)<=1e-15 else math.pi/2-_acos((radius*radius+offset*offset-roll_radius*roll_radius)/(2*radius*offset))
        radii.append(radius);phases.append(phase)
    base=math.asin(2*cross/diameter)
    boundaries=[]
    for side,radius,phase in zip((1,-1),radii,phases):
        angles=[]
        for cutter in (inner,outer):
            triangle=_acos((center_radius*center_radius+radius*radius-cutter*cutter*.25)/(2*radius*center_radius))
            angles.append(side*base+((theta-triangle)*direction-phase)*ratio)
        boundaries.append(angles)
    inward,outward=boundaries
    if direction>0:
        left=min(math.pi/2-inward[0],math.pi/2-outward[0])
        right=min(math.pi/2+inward[1],math.pi/2+outward[1])
    else:
        left=min(math.pi/2+inward[0],math.pi/2+outward[0])
        right=min(math.pi/2-inward[1],math.pi/2-outward[1])
    return left,right


def _envelope_minimum(evaluate,upper,side,iteration_limit,height_tolerance):
    """四点评价的包络搜索；保留0.382/0.618及最先最小值的分支规则。"""
    a,b=0.,upper
    c=a+.382*(b-a);d=a+.618*(b-a)
    fa,fc,fd,fb=(evaluate(h)[side] for h in (a,c,d,b))
    iterations=1
    while b-a>height_tolerance and iterations<iteration_limit:
        index=min(range(4),key=lambda i:(fa,fc,fd,fb)[i])
        if index<2:
            b,fb=d,fd;d,fd=c,fc;c=a+.382*(b-a);fc=evaluate(c)[side]
        else:
            a,fa=c,fc;c,fc=d,fd;d=a+.618*(b-a);fd=evaluate(d)[side]
        iterations+=1
    if b-a>height_tolerance:raise ArithmeticError('刀具包络高度搜索达到迭代上限')
    height=(a+b)*.5
    return height,evaluate(height)[side],iterations,b-a


def evaluate_generated_point(settings:CutterMachineSettings|dict,diameter_m:float,station_m:float,
                             *,iteration_limit=1000,height_tolerance_m=1e-10,
                             height_interval_policy='source_compatible')->dict:
    """计算两侧展成包络点及刀具区域。局部+Z从锥顶朝外，+X为径向。

    station_m沿模型的宽向锥距定义；两侧高度分别求解，不能复用一侧结果。
    region='cutter_tip_fillet'来自实际圆弧刀刃的包络高度判据。
    source_compatible保留边界处非正高度区间的直接中点行为，并标记未进行有效区间搜索；
    positive_only拒绝这些输入。不把源边界外推视作经过物理验证的刀尖包络。
    """
    if isinstance(settings,dict):settings=CutterMachineSettings(**settings)
    if not isinstance(settings,CutterMachineSettings):raise TypeError('settings须为CutterMachineSettings或具名字典')
    for key,value in (('diameter_m',diameter_m),('station_m',station_m),('height_tolerance_m',height_tolerance_m)):
        if isinstance(value,bool) or not math.isfinite(value) or value<=0:raise ValueError(f'{key}须为正有限实数')
    if isinstance(iteration_limit,bool) or not isinstance(iteration_limit,int) or iteration_limit<1:raise ValueError('iteration_limit须为正整数')
    if height_interval_policy not in ('source_compatible','positive_only'):raise ValueError('未知高度区间策略')
    s=settings._mm();diameter=diameter_m*1000;station=station_m*1000
    root=s['machine_root_angle_rad'];pitch=s['pitch_cone_angle_rad']
    radial=math.sin(pitch)*s['reference_cone_distance_m']-math.cos(pitch)*s['dedendum_m']
    axial=math.cos(pitch)*s['reference_cone_distance_m']+math.sin(pitch)*s['dedendum_m']
    xc=math.cos(root)*s['axial_setting_m']
    zc=-s['horizontal_setting_m']-math.sin(root)*s['axial_setting_m']-(radial/math.tan(root)-axial)*math.sin(root)
    _,upper,_=_coordinates_mm(xc,zc,diameter,station,pitch,root,0.)
    if not math.isfinite(upper):raise ValueError('机床包络高度边界不是有限值')
    if upper<=0 and height_interval_policy=='positive_only':raise ValueError('所求直径没有正的刀具高度搜索区间')
    evaluate=lambda h:_angles_at_height_mm(s,xc,zc,diameter,station,h)
    outputs={}
    for index,(flank,kind) in enumerate((('left','inner'),('right','outer'))):
        height,angle,iterations,width=_envelope_minimum(evaluate,upper,index,iteration_limit,height_tolerance_m*1000)
        limit=s[kind+'_tip_radius_m']*(1-math.sin(s[kind+'_pressure_angle_rad']))
        azimuth=(math.pi/2-angle if s['direction']>0 else -math.pi/2+angle) if index==0 else (-math.pi/2+angle if s['direction']>0 else math.pi/2-angle)
        radius=diameter_m*.5
        outputs[flank]=dict(point_m=(radius*math.cos(azimuth),radius*math.sin(azimuth),station_m/math.cos(pitch)-radius*math.tan(pitch)),
                            azimuth_rad=azimuth,envelope_angle_rad=angle,selected_height_m=height/1000,
                            cutter_tip_height_limit_m=limit/1000,
                            region='straight_blade' if height+1e-10>=limit else 'cutter_tip_fillet',
                            iterations=iterations,final_height_interval_m=abs(width)/1000,
                            search_state='nonpositive_source_interval' if upper<=0 else 'converged')
    return dict(diameter_m=diameter_m,station_m=station_m,maximum_cutter_height_m=upper/1000,
                flanks=outputs,geometry_model='two_sided_disk_cutter_generated_envelope',
                height_interval_policy=height_interval_policy,
                diagnostics=['非正高度区间按兼容策略取中点，属于源边界外推；未执行正区间包络搜索'] if upper<=0 else [],
                limitations=['机床参数是显式输入，不保证对应可制造配对设计','只计算几何包络，不求成对接触、载荷或压力'])
