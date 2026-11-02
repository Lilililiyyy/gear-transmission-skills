"""实际展成工作面的显式法向齿顶修缘；位移进入NURBS与实体，不修改压力。"""
import math


def tip_relief_value(v, amplitude_m, start):
    """C2五次渐增材料去除；长度m，起点与齿顶的一、二阶导数均为零。"""
    if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in (v,amplitude_m,start)):
        raise ValueError('法向修缘参数须为有限实数')
    if not 0 <= v <= 1 or amplitude_m < 0 or not 0 <= start < 1:
        raise ValueError('齿廓参数须在[0,1]、修缘量非负、起点在[0,1)')
    t = max(0.,(v-start)/(1-start))
    return amplitude_m*t*t*t*(10-15*t+6*t*t)


def _cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def _unit(vector):
    length = math.sqrt(sum(x*x for x in vector))
    if length <= 1e-30:
        raise ArithmeticError('修缘工作面法向退化')
    return tuple(x/length for x in vector)


def apply_normal_tip_relief(grid,request):
    """沿原CAD材料外法向反方向去除材料，保留根端接边和完整五面布局。"""
    if request.tip_relief_m == 0:
        return grid
    if 'normal_tip_relief' in grid['metadata']:
        raise ValueError('不能对已修缘网格再次施加同一法向位移')
    from .cad import build_generated_shell
    from OCP.gp import gp_Pnt,gp_Vec
    baseline = build_generated_shell(grid)
    patches = grid['patches']
    for name,other in (('left_flank','right_flank'),('right_flank','left_flank')):
        surface = baseline.surfaces[name]; bounds = surface.Bounds()
        width = patches[name]['width_count']; profile = patches[name]['profile_count']
        def normal(u,v):
            p=gp_Pnt();du=gp_Vec();dv=gp_Vec()
            surface.D1(bounds[0]+u*(bounds[1]-bounds[0]),bounds[2]+v*(bounds[3]-bounds[2]),p,du,dv)
            return _unit(_cross(du.Coord(),dv.Coord()))
        middle = (profile//2)*width+width//2
        u,v = patches[name]['uv'][middle]
        toward_slot = tuple(b-a for a,b in zip(patches[name]['points_m'][middle],patches[other]['points_m'][middle]))
        sense = 1 if sum(a*b for a,b in zip(normal(u,v),toward_slot)) > 0 else -1
        # 根端零位移，刀尖过渡及槽底保持原几何；实际法向由同一基准NURBS提供。
        points=[]
        for (u,v),point in zip(patches[name]['uv'],patches[name]['points_m']):
            delta=tip_relief_value(v,request.tip_relief_m,request.tip_relief_start)
            points.append(tuple(p-delta*sense*n for p,n in zip(point,normal(u,v))))
        patches[name]['points_m']=points
    grid['metadata']['normal_tip_relief']=dict(amplitude_m=request.tip_relief_m,
        start_profile_fraction=request.tip_relief_start,law='C2五次渐增法向材料去除',
        original_generated_formula_modified=False,actual_surface_geometry_modified=True,
        normal_definition='原NURBS的D1，按槽空腔方向识别材料外法向')
    return grid


def independent_relieved_point(request,sections,name,u,v):
    """由原刀盘点的中央差分独立评价修缘点，用于离插值节点的几何核查。"""
    from .generated_surface import _patch_point
    def point(x,y,patch=name):
        return _patch_point(sections.section(sections.station(x)),patch,y)[0]
    base=point(u,v);delta=tip_relief_value(v,request.tip_relief_m,request.tip_relief_start)
    if delta == 0:
        return base
    h=1e-3
    if not h < u < 1-h or not h < v < 1-h:
        raise ValueError('独立差分精度点须避开参数边界')
    du=tuple((a-b)/(2*h) for a,b in zip(point(u+h,v),point(u-h,v)))
    dv=tuple((a-b)/(2*h) for a,b in zip(point(u,v+h),point(u,v-h)))
    normal=_unit(_cross(du,dv))
    other='right_flank' if name == 'left_flank' else 'left_flank'
    toward_slot=tuple(b-a for a,b in zip(base,point(u,v,other)))
    sense=1 if sum(a*b for a,b in zip(normal,toward_slot)) > 0 else -1
    return tuple(p-delta*sense*n for p,n in zip(base,normal))
