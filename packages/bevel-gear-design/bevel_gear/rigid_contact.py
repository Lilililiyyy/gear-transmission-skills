"""可选BRep刚性间隙查询，为后续TCA提供几何基础。

距离由Open CASCADE计算，实体侵入另做布尔交集体积检查。开放曲面零距离
只能报告未决；本模块不求从动相位、传动误差、弹性变形或LTCA压力。
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class RigidPose:
    """绕全局轴线旋转后平移；轴心和平移单位m，角度rad。"""
    axis_origin_m: tuple=(0.,0.,0.)
    axis_direction: tuple=(0.,0.,1.)
    angle_rad: float=0.
    translation_m: tuple=(0.,0.,0.)

    def __post_init__(self):
        for key in ('axis_origin_m','axis_direction','translation_m'):
            vector=getattr(self,key)
            if len(vector)!=3 or any(isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) for v in vector):
                raise ValueError(f'{key}须含三个有限数')
        if sum(v*v for v in self.axis_direction)==0:raise ValueError('旋转轴方向不得为零向量')
        if isinstance(self.angle_rad,bool) or not math.isfinite(self.angle_rad):raise ValueError('转角须为有限弧度值')


@dataclass(frozen=True)
class RigidTransform:
    """列向量作用的正交旋转：world=rotation_rows@local+translation_m。"""
    rotation_rows:tuple=((1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    translation_m:tuple=(0.,0.,0.)

    def __post_init__(self):
        r=self.rotation_rows
        if len(r)!=3 or any(len(row)!=3 or any(isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) for v in row) for row in r):
            raise ValueError('rotation_rows须为3×3有限旋转矩阵，按行存储')
        RigidPose(translation_m=self.translation_m)
        for i in range(3):
            for j in range(3):
                if abs(sum(r[i][k]*r[j][k] for k in range(3))-(1 if i==j else 0))>1e-10:
                    raise ValueError('旋转矩阵须正交，不能包含缩放或剪切')
        determinant=sum(r[0][i]*(r[1][(i+1)%3]*r[2][(i+2)%3]-r[1][(i+2)%3]*r[2][(i+1)%3]) for i in range(3))
        if abs(determinant-1)>1e-10:raise ValueError('旋转矩阵须保向，不能把镜像当刚体旋转')


def _transform(shape,pose,unit):
    from OCP.gp import gp_Trsf,gp_Pnt,gp_Dir,gp_Ax1,gp_Vec
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    if isinstance(pose,dict):pose=RigidTransform(**pose) if 'rotation_rows' in pose else RigidPose(**pose)
    if not isinstance(pose,(RigidPose,RigidTransform)):raise TypeError('位姿须为RigidPose、RigidTransform或具名字典')
    transform=gp_Trsf()
    if isinstance(pose,RigidTransform):
        values=[value for row,shift in zip(pose.rotation_rows,pose.translation_m) for value in (*row,shift/unit)]
        transform.SetValues(*values)
    else:
        transform.SetRotation(gp_Ax1(gp_Pnt(*(v/unit for v in pose.axis_origin_m)),gp_Dir(*pose.axis_direction)),pose.angle_rad)
        transform.SetTranslationPart(gp_Vec(*(transform.TranslationPart().Coord(i+1)+pose.translation_m[i]/unit for i in range(3))))
    operation=BRepBuilderAPI_Transform(shape,transform,True)
    if not operation.IsDone():raise ArithmeticError('BRep刚体变换失败')
    return operation.Shape()


def query_clearance(first,second,*,kernel_length_unit_m,first_pose=None,second_pose=None,
                    contact_tolerance_m=1e-7,volume_tolerance_m3=1e-18,
                    verify_solid_interference=True)->dict:
    """查询两形状的距离和候选最近点；内核长度单位须显式提供。

    本项目cad生成的shape内部为mm，应传kernel_length_unit_m=0.001。
    不能仅从非负最短距离判断实体没有相互穿透，也不能据此声称已完成TCA。
    """
    for key,value in (('kernel_length_unit_m',kernel_length_unit_m),('contact_tolerance_m',contact_tolerance_m),('volume_tolerance_m3',volume_tolerance_m3)):
        if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or value<=0:raise ValueError(f'{key}须为正有限数')
    if not isinstance(verify_solid_interference,bool):raise TypeError('verify_solid_interference须为布尔值')
    try:
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.TopExp import TopExp
        from OCP.TopAbs import TopAbs_SOLID
        from OCP.TopTools import TopTools_IndexedMapOfShape
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
    except ImportError as error:raise ImportError('刚性BRep查询需要可选OCP依赖') from error
    for shape in (first,second):
        if shape.IsNull() or not BRepCheck_Analyzer(shape).IsValid():raise ValueError('刚性查询需要非空且拓扑有效的形状')
    a=_transform(first,first_pose or RigidPose(),kernel_length_unit_m)
    b=_transform(second,second_pose or RigidPose(),kernel_length_unit_m)
    distance=BRepExtrema_DistShapeShape(a,b)
    distance.Perform()
    if not distance.IsDone():raise ArithmeticError('BRep最短距离查询未完成')
    gap=distance.Value()*kernel_length_unit_m
    solutions=[]
    for index in range(1,distance.NbSolution()+1):
        p,q=distance.PointOnShape1(index),distance.PointOnShape2(index)
        solutions.append(dict(first_point_m=tuple(p.Coord(i)*kernel_length_unit_m for i in (1,2,3)),
                              second_point_m=tuple(q.Coord(i)*kernel_length_unit_m for i in (1,2,3))))
    solid_counts=[]
    for shape in (a,b):
        mapped=TopTools_IndexedMapOfShape();TopExp.MapShapes_s(shape,TopAbs_SOLID,mapped);solid_counts.append(mapped.Extent())
    volume=None;penetrating=None;boolean_state='not_requested'
    if verify_solid_interference:
        if all(solid_counts):
            common=BRepAlgoAPI_Common(a,b);common.Build()
            if common.IsDone():
                properties=GProp_GProps();BRepGProp.VolumeProperties_s(common.Shape(),properties)
                volume=abs(properties.Mass())*kernel_length_unit_m**3
                penetrating=volume>volume_tolerance_m3;boolean_state='computed'
            else:boolean_state='failed'
        else:boolean_state='open_geometry_volume_undefined'
    if penetrating is True:state='interference'
    elif gap>contact_tolerance_m:state='separated'
    elif penetrating is False:state='touching_within_tolerance'
    else:state='zero_distance_unresolved'
    return dict(minimum_distance_m=gap,candidate_pairs=solutions,candidate_pair_count=len(solutions),
                state=state,intersection_volume_m3=volume,is_penetrating=penetrating,
                solid_counts=solid_counts,boolean_state=boolean_state,
                contact_tolerance_m=contact_tolerance_m,kernel_length_unit_m=kernel_length_unit_m,
                interpretation='刚性给定位姿的BRep几何查询，非接触相位求解或LTCA',
                limitations=['未求载荷、接触压力、材料变形或传动误差','开放曲面零距离不能判定实体接触或侵入'])
