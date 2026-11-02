"""实际展成有限工作面的半空间接触几何，公开长度m、角度rad。

在共同齿轮位姿的局部切平面建立规则压力格，每格重新对两张真实
NURBS求射线交点并检查裁切BRep域。几何不由旧FE三角片或二次拟合
替代；局部半空间不包含轮体、轴承或跨齿结构柔度。
"""
from __future__ import annotations

from dataclasses import dataclass
import math


def _numeric():
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError('实际半空间接触几何需要NumPy与OCP') from error
    return np


def _real(value,name,positive=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or (positive and value<=0):
        raise ValueError(name+'须为'+('正' if positive else '')+'有限实数')
    return float(value)


def _surface_sample(surface,uv):
    """归一化UV的点和两个导数，所有长度从CAD毫米显式转为m。"""
    from OCP.gp import gp_Pnt,gp_Vec
    np=_numeric()
    bounds=surface.Bounds();u,v=uv
    p,du,dv=gp_Pnt(),gp_Vec(),gp_Vec()
    surface.D1(bounds[0]+u*(bounds[1]-bounds[0]),bounds[2]+v*(bounds[3]-bounds[2]),p,du,dv)
    return (np.array(p.Coord())*.001,
            np.array(du.Coord())*((bounds[1]-bounds[0])*.001),
            np.array(dv.Coord())*((bounds[3]-bounds[2])*.001))


def _normal(du,dv):
    np=_numeric();value=np.cross(du,dv);length=float(np.linalg.norm(value))
    if length<1e-22:raise ArithmeticError('实际曲面参数导数退化')
    return value/length


def _faces(shape,indices):
    from OCP.TopTools import TopTools_IndexedMapOfShape
    from OCP.TopExp import TopExp
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopoDS import TopoDS
    mapping=TopTools_IndexedMapOfShape();TopExp.MapShapes_s(shape,TopAbs_FACE,mapping)
    result=[]
    for index in indices:
        if not isinstance(index,int) or isinstance(index,bool) or not 0<=index<mapping.Extent():
            raise ValueError('裁切BRep面索引无效')
        result.append(TopoDS.Face_s(mapping.FindKey(index+1)))
    if not result:raise ValueError('没有实际工作侧的裁切面')
    return result


class _TrimmedSurface:
    """持有真实裁切面；仅在已证明参数一致时用二维分类器加速。"""
    def __init__(self,surface,faces,shape):
        from OCP.BRep import BRep_Tool
        from OCP.BRepClass import BRepClass_FaceClassifier
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.BRepTools import BRepTools
        from OCP.TopAbs import TopAbs_IN,TopAbs_ON,TopAbs_OUT
        from OCP.gp import gp_Pnt,gp_Pnt2d
        np=_numeric();self.surface=surface;self.faces=faces;self.classifier=BRepClass_FaceClassifier()
        self._Point=gp_Pnt;self._Point2=gp_Pnt2d;self._inside=(TopAbs_IN,TopAbs_ON)
        self.uv_compatible=[]
        bounds=surface.Bounds()
        for face in faces:
            actual=BRep_Tool.Surface_s(face);error=0.
            for uv in ((.2,.3),(.5,.5),(.8,.7)):
                source,*_=_surface_sample(surface,uv)
                point=actual.Value(bounds[0]+uv[0]*(bounds[1]-bounds[0]),bounds[2]+uv[1]*(bounds[3]-bounds[2]))
                error=max(error,float(np.linalg.norm(source-np.array(point.Coord())*.001)))
            self.uv_compatible.append(error<1e-10)
        sense=None
        for face in faces:
            fb=BRepTools.UVBounds_s(face)
            for a,b in ((.5,.5),(.25,.5),(.75,.5),(.5,.25),(.5,.75)):
                real=(fb[0]+a*(fb[1]-fb[0]),fb[2]+b*(fb[3]-fb[2]))
                uv=((real[0]-bounds[0])/(bounds[1]-bounds[0]),(real[1]-bounds[2])/(bounds[3]-bounds[2]))
                point,du,dv=_surface_sample(surface,uv)
                self.classifier.Perform(face,gp_Pnt(*tuple(point*1000)),1e-7)
                if self.classifier.State() not in self._inside:continue
                normal=_normal(du,dv)
                states=[BRepClass3d_SolidClassifier(shape,gp_Pnt(*tuple((point+sign*5e-7*normal)*1000)),1e-7).State() for sign in (1,-1)]
                if states==[TopAbs_OUT,TopAbs_IN]:sense=1
                elif states==[TopAbs_IN,TopAbs_OUT]:sense=-1
                if sense is not None:break
            if sense is not None:break
        if sense is None:raise ArithmeticError('实际裁切工作面的材料内外法向未能分类')
        self.outward_sense=sense

    def inside(self,uv,point):
        bounds=self.surface.Bounds()
        for face,compatible in zip(self.faces,self.uv_compatible):
            if compatible:
                location=self._Point2(bounds[0]+uv[0]*(bounds[1]-bounds[0]),bounds[2]+uv[1]*(bounds[3]-bounds[2]))
            else:location=self._Point(*tuple(point*1000))
            self.classifier.Perform(face,location,1e-7)
            if self.classifier.State() in self._inside:return True
        return False


@dataclass
class PreparedGeneratedHalfspace:
    """双轮实际裁切工作面与共享初装，不持有FE或厂商组件。"""
    tca: object
    trimmed: dict
    diagnostics: dict


def prepare_generated_halfspace(pinion_request,wheel_request,*,cad_models=None,
        pinion_blank_settings=None,wheel_blank_settings=None,width_count=25,profile_count=25):
    """从显式双轮制造请求准备实际有限工作面，可复用相同请求的完整CAD。

    未给cad_models时必须给两轮齿坯设置并实际创建裁切整轮，不能把开放
    齿槽壳的整个UV矩形自动认作有效材料工作面。
    """
    from .generated_surface import _request
    from .generated_tca import prepare_generated_tca
    requests=dict(pinion=_request(pinion_request),wheel=_request(wheel_request))
    if cad_models is None:
        from .generated_solid import build_generated_gear
        if pinion_blank_settings is None or wheel_blank_settings is None:
            raise ValueError('须给相同请求的实际cad_models，或显式两轮齿坯设置')
        cad_models={member:build_generated_gear(requests[member],settings,width_count=width_count,profile_count=profile_count)
                    for member,settings in (('pinion',pinion_blank_settings),('wheel',wheel_blank_settings))}
    if not isinstance(cad_models,dict) or set(cad_models)!=set(requests):raise ValueError('cad_models须含pinion和wheel')
    for member,model in cad_models.items():
        if model.request!=requests[member] or model.diagnostics.get('solids')!=1 or not model.diagnostics.get('valid'):
            raise ValueError('实际CAD请求/闭实体与本半空间几何输入不一致')
    tca=prepare_generated_tca(*requests.values(),cad_models=cad_models)
    trimmed={}
    for member,patch in (('pinion','left_flank'),('wheel','right_flank')):
        model=cad_models[member]
        faces=_faces(model.shape,model.slot_face_groups[0][patch])
        trimmed[member]=_TrimmedSurface(model.shell.surfaces[patch],faces,model.shape)
    return PreparedGeneratedHalfspace(tca,trimmed,dict(units='m/rad',cad_internal_unit='mm',
        trimmed_material_domain=True,material_normal_sense={m:t.outward_sense for m,t in trimmed.items()},
        source='实际展成NURBS及同请求闭实体的裁切BRep工作面',
        model='局部半空间表面接触几何；不生成轮体/跨齿结构柔度'))


def _ray(surface,rotation,origin,basis,target_xy,seed_uv,tolerance_m,maximum_iterations=15):
    """有界Newton求固定切平面坐标的实际NURBS射线交点。"""
    np=_numeric();uv=np.clip(np.asarray(seed_uv,dtype=float),1e-12,1-1e-12)
    last=math.inf
    for iteration in range(maximum_iterations):
        local,du,dv=_surface_sample(surface,uv);world=rotation@local
        projected=basis@(world-origin);residual=projected-np.asarray(target_xy)
        error=float(np.linalg.norm(residual))
        if error<=tolerance_m:
            return dict(converged=True,uv=uv,local=local,world=world,du_world=rotation@du,dv_world=rotation@dv,
                        residual_m=error,iterations=iteration+1)
        jac=np.column_stack((basis@(rotation@du),basis@(rotation@dv)))
        if not np.isfinite(jac).all() or abs(float(np.linalg.det(jac)))<1e-16:
            return dict(converged=False,status='ray_graph_degenerate',residual_m=error)
        step=np.linalg.solve(jac,residual)
        fraction=min(1.,.2/max(float(np.max(np.abs(step))),1e-20))
        trial=np.clip(uv-fraction*step,0.,1.)
        if np.linalg.norm(trial-uv)<1e-13 or (error>=last and iteration>10):
            return dict(converged=False,status='outside_or_unresolved_parameter_domain',residual_m=error)
        uv=trial;last=error
    return dict(converged=False,status='ray_iteration_limit',residual_m=last)


def _initial_uv(surface,rotation,origin,basis,seed_uv):
    np=_numeric();point,du,dv=_surface_sample(surface,seed_uv)
    coordinate=basis@(rotation@point-origin)
    jac=np.column_stack((basis@(rotation@du),basis@(rotation@dv)))
    if abs(float(np.linalg.det(jac)))<1e-16:raise ArithmeticError('切平面图的参数导数退化')
    return coordinate,np.linalg.inv(jac)


def _sample_pair(prepared,frames,origin,basis,normal,seed_uv,target_xy,tolerance_m,maximum_slope,maximum_normal_angle_rad,seed_maps=None):
    np=_numeric();values={};failure=None
    for member in ('pinion','wheel'):
        trimmed=prepared.trimmed[member];rotation=np.asarray(frames[member+'_rotation_rows'])
        coordinate,inverse=(seed_maps[member] if seed_maps is not None else
                            _initial_uv(trimmed.surface,rotation,origin,basis,seed_uv[member]))
        guess=np.asarray(seed_uv[member])+inverse@(np.asarray(target_xy)-coordinate)
        value=_ray(trimmed.surface,rotation,origin,basis,target_xy,guess,tolerance_m)
        if not value['converged']:return None,member+':'+value['status']
        if not trimmed.inside(value['uv'],value['local']):return None,member+':outside_trimmed_face'
        outward=trimmed.outward_sense*_normal(value['du_world'],value['dv_world'])
        cosine=float(outward@normal)
        slope=float(np.linalg.norm(basis@outward)/abs(cosine)) if abs(cosine)>1e-12 else math.inf
        if slope>maximum_slope:return None,member+':halfspace_slope_limit'
        value.update(outward=outward,slope=slope);values[member]=value
    first,second=values['pinion'],values['wheel']
    cosine=float(-first['outward']@second['outward'])
    if cosine<math.cos(maximum_normal_angle_rad):return None,'opposing_normal_limit'
    gap=float((second['world']-first['world'])@normal)
    arms={};moment_arms={}
    for member,sign in (('pinion',1),('wheel',-1)):
        value=values[member];velocity=np.cross(frames[member+'_axis_world'],value['world'])
        # 固定投影格的高度导数须除以N·n；不能把斜齿面射线重投影漏掉。
        arms[member]=sign*float(value['outward']@velocity/(value['outward']@normal))
        moment_arms[member]=sign*float(normal@velocity)
    return dict(gap_m=gap,first=first,second=second,lever_arms_m=arms,normal_moment_arms_m=moment_arms,
                maximum_slope=max(first['slope'],second['slope']),normal_alignment=cosine),None


def _rejection_kind(reason):
    """区分材料裁切、半空间门限与未解查询，避免混称真实齿端。"""
    if 'halfspace_slope_limit' in reason or 'opposing_normal_limit' in reason:return 'model_validity'
    if 'outside_trimmed_face' in reason:return 'trimmed_domain'
    if 'torque_direction' in reason:return 'torque_direction'
    # 有界Newton未解并不能证明射线已落在材料边缘以外。
    return 'unresolved_query'


def build_generated_halfspace_patches(prepared,driving_angle_rad,pinion_torque_nm,*,
        tooth_offsets=(-1,0,1),grid_shape=(161,81),spacing_m=(1e-4,2.5e-5),
        position_tolerance_m=2e-11,maximum_surface_slope=.3,maximum_normal_angle_rad=.3,
        check_cell_corners=True,tca_results=None,follower_correction_rad=0.,cover_complete_working_flank=False,
        ray_backend='scalar',plane_basis_angles_rad=None):
    """以共同首接触相位求真实规则间隙格，返回每齿独立patch和全部质量状态。

    grid_shape为(nx,ny)，数组形状(ny,nx)；spacing_m为(dx,dy)。每格
    压力作用于参考切平面的dx*dy投影面积，原三维点/法向另存。严格域
    策略要求中心及四角均在双轮真实裁切面内，部分裁切格不冒充整面积。
    """
    from .generated_tca import solve_generated_tca,generated_tca_frames
    from .halfspace_contact import RectangularContactGrid
    np=_numeric()
    if not isinstance(prepared,PreparedGeneratedHalfspace):raise TypeError('需要PreparedGeneratedHalfspace')
    angle=_real(driving_angle_rad,'小轮角');torque=_real(pinion_torque_nm,'小轮有符号扭矩')
    if torque==0:raise ValueError('本接触入口要求非零有符号扭矩')
    if len(grid_shape)!=2 or any(isinstance(n,bool) or not isinstance(n,int) or n<3 for n in grid_shape):raise ValueError('grid_shape须为两个不小于3的整数')
    if len(spacing_m)!=2:raise ValueError('spacing_m须为两个正间距')
    dx,dy=(_real(v,'压力格距',True) for v in spacing_m)
    requested_nx,requested_ny=grid_shape
    if requested_nx*requested_ny>500_000:raise ValueError('单接触patch超过显式五十万格上限')
    tolerance=_real(position_tolerance_m,'射线残差容差',True)
    slope_limit=_real(maximum_surface_slope,'半空间表面斜率门限',True)
    normal_limit=_real(maximum_normal_angle_rad,'两面法向夹角门限',True)
    offsets=tuple(tooth_offsets)
    if not offsets or len(set(offsets))!=len(offsets) or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets):raise ValueError('齿对偏移须为不重复整数')
    basis_angles=None
    if plane_basis_angles_rad is not None:
        if not isinstance(plane_basis_angles_rad,dict):raise ValueError('投影网格方向须按每个原齿对索引提供rad字典')
        basis_angles={}
        for key,value in plane_basis_angles_rad.items():
            if isinstance(key,bool) or not isinstance(key,(int,str)):raise ValueError('网格方向齿对索引须为整数或其JSON字符串')
            try:index=int(key)
            except ValueError as error:raise ValueError('网格方向齿对索引不可解析') from error
            if str(index)!=str(key) or index in basis_angles:raise ValueError('网格方向齿对索引须为规范且不重复的整数')
            basis_angles[index]=_real(value,'实际投影网格旋转角rad')
        if set(basis_angles)!=set(offsets):raise ValueError('每个所选齿对都必须有明确的投影网格方向')
    direction=1 if torque>0 else -1
    if ray_backend not in ('scalar','batch'):
        raise ValueError('实际NURBS射线后端须为scalar或batch')
    if ray_backend == 'batch':
        from .nurbs_batch import sample_pair_batch
    if tca_results is None:tca_results={k:solve_generated_tca(prepared.tca,angle,tooth_offset=k) for k in offsets}
    else:tca_results={k:tca_results[k] if k in tca_results else tca_results[str(k)] for k in offsets}
    valid_reference=[value for value in tca_results.values() if value['converged']]
    if not valid_reference:raise ArithmeticError('所选齿对没有有效实际未加载相切，不能以候选点冒充参考')
    reference=min(valid_reference,key=lambda value:(direction*value['follower_increment_angle_rad'],abs(value['tooth_offset'])))
    follower=reference['follower_increment_angle_rad']+_real(follower_correction_rad,'共同相位修正')
    patches=[]
    for offset in offsets:
        value=tca_results[offset]
        if value['converged']:
            seeds=dict(pinion=value['pinion_uv'],wheel=value['wheel_uv']);seed_source='accepted_tca_parameters'
        else:
            # 只作为迭代初值；任何候选点都不直接进入物理网格。
            parameters=value['candidate']['parameters']
            seeds=dict(pinion=parameters[:2],wheel=parameters[2:4]);seed_source='unaccepted_tca_parameters_only_as_newton_seed'
        frames=generated_tca_frames(prepared.tca,angle,follower,tooth_offset=offset)
        rotation=np.asarray(frames['pinion_rotation_rows'])
        point,du,dv=_surface_sample(prepared.trimmed['pinion'].surface,seeds['pinion'])
        origin=rotation@point
        normal=prepared.trimmed['pinion'].outward_sense*_normal(rotation@du,rotation@dv)
        first_tangent=rotation@du;first_tangent-=normal*float(normal@first_tangent);first_tangent/=np.linalg.norm(first_tangent)
        second_tangent=np.cross(normal,first_tangent)
        basis=np.array([first_tangent,second_tangent])
        if basis_angles is not None:
            from .contact_plane import rotate_plane_basis
            basis=rotate_plane_basis(basis,basis_angles[offset])
        seed_maps={member:_initial_uv(prepared.trimmed[member].surface,np.asarray(frames[member+'_rotation_rows']),origin,basis,seeds[member])
                   for member in ('pinion','wheel')}
        if cover_complete_working_flank:
            surface=prepared.trimmed['pinion'].surface;poles=[]
            for u in range(1,surface.NbUPoles()+1):
                for v in range(1,surface.NbVPoles()+1):
                    if surface.Weight(u,v)<=0:raise ValueError('正权NURBS凸包覆盖不接受非正权')
                    poles.append(basis@(rotation@(np.asarray(surface.Pole(u,v).Coord())*.001)-origin))
            poles=np.asarray(poles)
            low=np.floor(poles.min(axis=0)/[dx,dy])*[dx,dy]-[dx,dy]
            high=np.ceil(poles.max(axis=0)/[dx,dy])*[dx,dy]+[dx,dy]
            nx,ny=map(int,np.ceil((high-low)/[dx,dy]))
            center_origin=low+np.array([dx,dy])/2
        else:
            nx,ny=requested_nx,requested_ny;center_origin=np.array([-(nx-1)*dx/2,-(ny-1)*dy/2])
        if nx*ny>500_000:raise ValueError('完整工作面按本格距超过显式五十万格上限')
        grid=RectangularContactGrid(nx,ny,dx,dy,origin_xy_m=tuple(center_origin))
        shape=(ny,nx);gaps=np.full(shape,np.nan);valid=np.zeros(shape,dtype=bool)
        first_points=np.full(shape+(3,),np.nan);second_points=first_points.copy()
        first_uv=np.full(shape+(2,),np.nan);second_uv=first_uv.copy()
        wheel_arms=np.full(shape,np.nan);pinion_arms=wheel_arms.copy()
        normal_wheel=wheel_arms.copy();normal_pinion=wheel_arms.copy();slopes=wheel_arms.copy();angles=wheel_arms.copy()
        first_jacobian=wheel_arms.copy();second_jacobian=wheel_arms.copy()
        first_normals=first_points.copy();second_normals=first_points.copy()
        failures={};maximum_error=0.
        rejection_masks={kind:np.zeros(shape,dtype=bool) for kind in ('model_validity','trimmed_domain','unresolved_query','torque_direction')}
        corner_valid=None
        corner_first=np.full((ny+1,nx+1,3),np.nan);corner_second=corner_first.copy();corner_failures={}
        if check_cell_corners:
            corner_valid=np.zeros((ny+1,nx+1),dtype=bool)
            corner_rejections={kind:np.zeros(corner_valid.shape,dtype=bool) for kind in rejection_masks}
            corner_samples=corner_reasons=None
            if ray_backend == 'batch':
                gx,gy=np.meshgrid(center_origin[0]+(np.arange(nx+1)-.5)*dx,
                    center_origin[1]+(np.arange(ny+1)-.5)*dy)
                corner_samples,corner_reasons=sample_pair_batch(prepared,frames,origin,basis,normal,seeds,
                    np.column_stack((gx.ravel(),gy.ravel())),tolerance,slope_limit,normal_limit,seed_maps)
            for j in range(ny+1):
                for i in range(nx+1):
                    if corner_samples is None:
                        sample,failure=_sample_pair(prepared,frames,origin,basis,normal,seeds,
                            (center_origin[0]+(i-.5)*dx,center_origin[1]+(j-.5)*dy),tolerance,slope_limit,normal_limit,seed_maps)
                    else:
                        index=j*(nx+1)+i;sample=corner_samples[index];failure=corner_reasons[index]
                    if sample is not None:
                        corner_valid[j,i]=True;corner_first[j,i]=sample['first']['world'];corner_second[j,i]=sample['second']['world']
                    else:
                        corner_failures[failure]=corner_failures.get(failure,0)+1
                        corner_rejections[_rejection_kind(failure)][j,i]=True
            del corner_samples,corner_reasons
        centre_samples=centre_reasons=centre_lookup=None
        if ray_backend == 'batch':
            candidates=(corner_valid[:-1,:-1]&corner_valid[1:,:-1]&corner_valid[:-1,1:]&corner_valid[1:,1:]
                if check_cell_corners else np.ones((ny,nx),dtype=bool))
            rows,columns=np.nonzero(candidates)
            xy=np.column_stack((center_origin[0]+columns*dx,center_origin[1]+rows*dy))
            centre_samples,centre_reasons=sample_pair_batch(prepared,frames,origin,basis,normal,seeds,
                xy,tolerance,slope_limit,normal_limit,seed_maps)
            centre_lookup=np.full((ny,nx),-1,dtype=int);centre_lookup[rows,columns]=np.arange(len(rows))
        for j in range(ny):
            for i in range(nx):
                if check_cell_corners and not corner_valid[j:j+2,i:i+2].all():
                    for kind,mask in corner_rejections.items():rejection_masks[kind][j,i]=bool(mask[j:j+2,i:i+2].any())
                    failures['partial_or_outside_cell']=failures.get('partial_or_outside_cell',0)+1;continue
                if centre_samples is None:
                    sample,failure=_sample_pair(prepared,frames,origin,basis,normal,seeds,
                        (center_origin[0]+i*dx,center_origin[1]+j*dy),tolerance,slope_limit,normal_limit,seed_maps)
                else:
                    index=centre_lookup[j,i];sample=centre_samples[index];failure=centre_reasons[index]
                if sample is None:
                    rejection_masks[_rejection_kind(failure)][j,i]=True
                    failures[failure]=failures.get(failure,0)+1;continue
                if direction*sample['lever_arms_m']['wheel']<=1e-10 or direction*sample['lever_arms_m']['pinion']<=1e-10:
                    rejection_masks['torque_direction'][j,i]=True
                    failures['torque_direction_not_supported']=failures.get('torque_direction_not_supported',0)+1;continue
                valid[j,i]=True;gaps[j,i]=sample['gap_m']
                for prefix,member in (('first','pinion'),('second','wheel')):
                    item=sample[prefix]
                    (first_points if prefix=='first' else second_points)[j,i]=item['world']
                    (first_uv if prefix=='first' else second_uv)[j,i]=item['uv']
                    maximum_error=max(maximum_error,item['residual_m'])
                wheel_arms[j,i]=sample['lever_arms_m']['wheel'];pinion_arms[j,i]=sample['lever_arms_m']['pinion']
                normal_wheel[j,i]=sample['normal_moment_arms_m']['wheel'];normal_pinion[j,i]=sample['normal_moment_arms_m']['pinion']
                slopes[j,i]=sample['maximum_slope']
                angles[j,i]=math.acos(max(-1.,min(1.,sample['normal_alignment'])))
                first_normals[j,i]=sample['first']['outward'];second_normals[j,i]=sample['second']['outward']
                first_jacobian[j,i]=1/abs(float(sample['first']['outward']@normal))
                second_jacobian[j,i]=1/abs(float(sample['second']['outward']@normal))
        patches.append(dict(tooth_offset=offset,grid=grid,gaps_m=gaps,valid_mask=valid,
            model_validity_rejected_mask=rejection_masks['model_validity'],trimmed_domain_rejected_mask=rejection_masks['trimmed_domain'],
            unresolved_query_rejected_mask=rejection_masks['unresolved_query'],torque_direction_rejected_mask=rejection_masks['torque_direction'],
            wheel_lever_arms_m=wheel_arms,pinion_lever_arms_m=pinion_arms,
            wheel_normal_moment_arms_m=normal_wheel,pinion_normal_moment_arms_m=normal_pinion,
            first_points_world_m=first_points,second_points_world_m=second_points,
            first_corner_points_world_m=corner_first,second_corner_points_world_m=corner_second,corner_valid_mask=corner_valid,
            surface_slopes=slopes,opposing_normal_angles_rad=angles,
            first_outward_normals_world=first_normals,second_outward_normals_world=second_normals,
            first_surface_area_jacobian=first_jacobian,second_surface_area_jacobian=second_jacobian,
            first_uv=first_uv,second_uv=second_uv,frames=frames,
            plane_origin_world_m=origin.tolist(),plane_tangent_rows=basis.tolist(),plane_outward_normal_world=normal.tolist(),
            diagnostics=dict(initial_seed_source=seed_source,tca_seed_accepted=bool(value['converged']),
                valid_cells=int(valid.sum()),total_cells=nx*ny,rejected_cells=failures,
                rejected_corners=corner_failures,complete_first_flank_projected_convex_hull_covered=bool(cover_complete_working_flank),
                rejected_by_category={kind:int(mask.sum()) for kind,mask in rejection_masks.items()},
                minimum_signed_gap_m=float(gaps[valid].min()) if valid.any() else None,
                maximum_ray_position_residual_m=maximum_error,maximum_surface_slope=float(slopes[valid].max()) if valid.any() else None,
                cell_area_m2=dx*dy,area_definition='参考切平面矩形投影面积；不是任意裁切三角形的原面积',
                lever_definition='负的固定投影格间隙对参数轴转角导数；保留局部斜率分母',
                normal_moment_levers_separately_stored=True,partial_cells_excluded=bool(check_cell_corners)),seed_uv=seeds))
    if not any(p['valid_mask'].any() for p in patches):raise ArithmeticError('没有双面实际裁切域内的有效压力格')
    result=dict(patches=patches,driving_angle_rad=angle,pinion_torque_nm=torque,
        pinion_to_wheel_ratio=prepared.tca.requests['pinion'].tooth_count/prepared.tca.requests['wheel'].tooth_count,
        reference_tca=reference,tca_results=tca_results,reference_follower_angle_rad=follower,
        surface_reference_phase_correction_rad=follower_correction_rad,
        build_parameters=dict(tooth_offsets=offsets,grid_shape=grid_shape,spacing_m=spacing_m,
            position_tolerance_m=position_tolerance_m,maximum_surface_slope=maximum_surface_slope,
            maximum_normal_angle_rad=maximum_normal_angle_rad,check_cell_corners=check_cell_corners,
            cover_complete_working_flank=cover_complete_working_flank,ray_backend=ray_backend),
        diagnostics=dict(units='m/rad/N/Pa',common_wheel_phase=True,
            geometry='实际NURBS逐格射线求交及实际裁切BRep域；没有拟合二次间隙',
            compliance_model='每齿局部弹性半空间；不含轮体/跨齿结构柔度',
            closure_definition='该局部半空间模型的从动参数轴闭合量，不是整轮FE全局TE'))
    if basis_angles is not None:result['build_parameters']['plane_basis_angles_rad']=basis_angles
    return result


def check_generated_halfspace_linearization(prepared,geometry,result,*,maximum_error_m=1e-10,
        derivative_epsilon_rad=1e-6):
    """以加载转角重新求真实射线交点，并独立复核功共轭力臂。

    全部活跃格检查g_true与g0−ell2*theta；最多12个分布活跃点另作
    对从动角的中心差分。未通过时只返回诊断，不把几何误差当弹性位移。
    """
    from .generated_tca import generated_tca_frames
    np=_numeric();limit=_real(maximum_error_m,'几何线性化容差m',True)
    epsilon=_real(derivative_epsilon_rad,'导数差分角rad',True)
    theta=result.get('theta_rad')
    if not result.get('converged') or theta is None:return dict(passed=False,status='contact_solver_not_converged')
    selected=[p for p in geometry['patches'] if p['valid_mask'].any()]
    maximum=0.;derivative_error=0.;checked=0;fd_count=0;failed=[];fd_failed=[]
    actual_stationarity=0.;actual_minimum=math.inf;actual_complementarity=0.
    for patch,field in zip(selected,result['fields']):
        origin=np.asarray(patch['plane_origin_world_m']);basis=np.asarray(patch['plane_tangent_rows']);normal=np.asarray(patch['plane_outward_normal_world'])
        angle=geometry['driving_angle_rad'];follower=geometry['reference_follower_angle_rad'];offset=patch['tooth_offset']
        frames=generated_tca_frames(prepared.tca,angle,follower+theta,tooth_offset=offset)
        active=np.argwhere(field['contact_mask'])
        fd_positions=set(np.linspace(0,max(0,len(active)-1),min(12,len(active)),dtype=int).tolist())
        grid=patch['grid'];x0,y0=grid.origin_xy_m
        for number,(j,i) in enumerate(active):
            xy=(x0+i*grid.spacing_x_m,y0+j*grid.spacing_y_m)
            seeds=dict(pinion=patch['first_uv'][j,i],wheel=patch['second_uv'][j,i])
            value,reason=_sample_pair(prepared,frames,origin,basis,normal,seeds,xy,1e-13,
                geometry['build_parameters']['maximum_surface_slope'],geometry['build_parameters']['maximum_normal_angle_rad'])
            checked+=1
            if value is None:
                failed.append(dict(tooth_offset=offset,row=int(j),column=int(i),reason=reason));continue
            predicted=patch['gaps_m'][j,i]-patch['wheel_lever_arms_m'][j,i]*theta
            maximum=max(maximum,abs(value['gap_m']-predicted))
            actual_remaining=value['gap_m']+field['elastic_deflection_m'][j,i]
            actual_stationarity=max(actual_stationarity,abs(float(actual_remaining)))
            actual_minimum=min(actual_minimum,float(actual_remaining))
            actual_complementarity+=abs(float(field['forces_n'][j,i]*actual_remaining))
            if number in fd_positions:
                values=[]
                for sign in (-1,1):
                    frame=generated_tca_frames(prepared.tca,angle,follower+sign*epsilon,tooth_offset=offset)
                    sample,_=_sample_pair(prepared,frame,origin,basis,normal,seeds,xy,1e-13,
                        geometry['build_parameters']['maximum_surface_slope'],geometry['build_parameters']['maximum_normal_angle_rad'])
                    values.append(sample)
                if all(sample is not None for sample in values):
                    observed=(values[1]['gap_m']-values[0]['gap_m'])/(2*epsilon)
                    derivative_error=max(derivative_error,abs(observed+patch['wheel_lever_arms_m'][j,i]));fd_count+=1
                else:fd_failed.append(dict(tooth_offset=offset,row=int(j),column=int(i)))
    actual_acceptance=float(result['certificate']['gap_acceptance_m'])+limit
    return dict(passed=bool(not failed and not fd_failed and maximum<=limit and derivative_error<=5e-7 and actual_stationarity<=actual_acceptance),
        maximum_nonlinear_gap_error_m=maximum,maximum_accepted_error_m=limit,
        active_cells_reprojected=checked,failed_loaded_ray_queries=failed,
        actual_loaded_active_stationarity_error_m=actual_stationarity,
        actual_loaded_active_minimum_remaining_gap_m=actual_minimum if checked else None,
        actual_loaded_active_complementarity_residual_nm=actual_complementarity,
        actual_loaded_active_gap_acceptance_m=actual_acceptance,
        inactive_cells_nonlinearly_reprojected=False,
        finite_difference_points=fd_count,maximum_wheel_lever_difference_m=derivative_error,
        failed_finite_difference_queries=fd_failed,
        finite_difference_acceptance_m=5e-7,definition='真实重新投影的g(theta)对比g0−ell2*theta；力臂对固定投影格功共轭')


def solve_generated_halfspace_contact(geometry,*,pinion_young_modulus_pa,pinion_poisson_ratio,
        wheel_young_modulus_pa,wheel_poisson_ratio,maximum_iterations=20000,
        gap_tolerance_m=1e-10,relative_tolerance=1e-7,allow_initial_penetration=False,
        maximum_initial_penetration_m=1e-10,prepared=None,maximum_outer_iterations=4,
        maximum_linearization_error_m=1e-10):
    """对真实间隙patch求半解析压力，可外层重投影消除几何线性化误差。

    负的未加载gap仅能在显式允许且不超过maximum_initial_penetration_m
    时作为数值误差保留；加载外迭代中的刚性穿入另行记录，由弹性开口
    抵消。没有给prepared时几何重投影状态保留未检查，不称非线性验收。
    """
    from .halfspace_contact import prepare_halfspace_operator,solve_halfspace_contact
    np=_numeric();base_geometry=geometry;total_theta=0.;history=[]
    initial_min=min(float(p['gaps_m'][p['valid_mask']].min()) for p in geometry['patches'] if p['valid_mask'].any())
    penetration_limit=_real(maximum_initial_penetration_m,'未加载负间隙数值门限m',True)
    if initial_min < -penetration_limit:raise ArithmeticError('未加载实际工作面有超限负gap，须重求共同首接触相位，不能以允许穿透掩盖')
    if isinstance(maximum_outer_iterations,bool) or not isinstance(maximum_outer_iterations,int) or maximum_outer_iterations<1:raise ValueError('maximum_outer_iterations须为正整数')
    # 两轴输入的名义力矩换算与既有副定义一致；主动轴重构仍独立报告。
    ratio=geometry.get('pinion_to_wheel_ratio')
    if ratio is None or ratio<=0:raise ValueError('几何结果须给正的pinion_to_wheel_ratio')
    target=geometry['pinion_torque_nm']/ratio
    for outer in range(maximum_outer_iterations):
        requests=[];operators=[];arms=[];masks=[];selected=[]
        for patch in geometry['patches']:
            if not patch['valid_mask'].any():continue
            operators.append(prepare_halfspace_operator(patch['grid'],pinion_young_modulus_pa,pinion_poisson_ratio,
                                                       wheel_young_modulus_pa,wheel_poisson_ratio))
            requests.append(patch['gaps_m']);arms.append(patch['wheel_lever_arms_m']);masks.append(patch['valid_mask']);selected.append(patch)
        if not operators:raise ArithmeticError('没有有效半空间压力格')
        result=solve_halfspace_contact(requests,operators,lever_arms_m=arms,torque_nm=target,valid_masks=masks,
            maximum_iterations=maximum_iterations,gap_tolerance_m=gap_tolerance_m,relative_tolerance=relative_tolerance,
            allow_initial_penetration=allow_initial_penetration or outer>0)
        if not result['converged']:break
        step=result['theta_rad']
        check=(check_generated_halfspace_linearization(prepared,geometry,result,maximum_error_m=maximum_linearization_error_m)
               if prepared is not None else dict(passed=None,status='not_checked_no_prepared_geometry'))
        history.append(dict(iteration=outer+1,reference_angle_offset_rad=total_theta,step_theta_rad=step,
                            loaded_rigid_overlap_allowed=outer>0,geometry_check=check))
        if prepared is None or check['passed']:
            total_theta+=step;break
        if outer+1==maximum_outer_iterations:
            candidate=dict(result,geometry_check=check)
            result=dict(converged=False,status='geometry_linearization_limit',candidate=candidate,fields=[],theta_rad=None,
                        maximum_pressure_pa=None,incremental_stiffness_nm_per_rad=None)
            break
        total_theta+=step
        geometry=build_generated_halfspace_patches(prepared,base_geometry['driving_angle_rad'],base_geometry['pinion_torque_nm'],
            **base_geometry['build_parameters'],tca_results=base_geometry['tca_results'],
            follower_correction_rad=base_geometry['surface_reference_phase_correction_rad']+total_theta)
    if result.get('converged'):
        result['linearization_step_theta_rad']=result['theta_rad'];result['theta_rad']=total_theta
        result['loaded_follower_angle_rad']=geometry['reference_follower_angle_rad']+result['linearization_step_theta_rad']
        unloaded={p['tooth_offset']:p for p in base_geometry['patches']}
        for patch in geometry['patches']:
            original=unloaded[patch['tooth_offset']]
            patch['unloaded_gaps_m']=original['gaps_m'].copy();patch['unloaded_valid_mask']=original['valid_mask'].copy()
            patch['affine_total_angle_intercept_gaps_m']=patch['gaps_m']+patch['wheel_lever_arms_m']*(total_theta-result['linearization_step_theta_rad'])
    result.update(model='实际有限展成面＋每齿局部半空间半解析接触；非整轮FE',
        geometry=geometry,selected_tooth_offsets=[p['tooth_offset'] for p in selected],
        unloaded_reference_follower_angle_rad=base_geometry['reference_follower_angle_rad'],
        geometry_outer_history=history,geometry_linearization_checked=prepared is not None,
        minimum_unloaded_gap_m=initial_min,maximum_initial_penetration_m=penetration_limit,
        initial_penetration_explicitly_allowed=allow_initial_penetration,
        wheel_target_torque_nm=target,pressure_definition='参考接触切平面的法向压力Pa',
        closure_definition='局部半空间的从动参数轴闭合量；不含轮体/跨齿结构柔度',
        geometry_angle_definition='geometry.reference_follower_angle_rad是最后线性化参考；加载角=该参考+linearization_step_theta_rad；theta_rad是相对原无载参考的总量',
        physical_validation=False,limitations=['局部小斜率、无摩擦线弹性半空间',
            '不含轮体/轴承及跨齿结构耦合；闭合量不能冒充原整轮FE的全局TE',
            '实际裁切边缘或窗口边界承载须独立检查，不能强制绘成单连通斑点'])
    if result.get('converged'):
        reconstructed=0.;shares={};edge=[];active_slope=0.;active_angle=0.;projected_area=0.;first_area=0.;second_area=0.
        for patch,field in zip(selected,result['fields']):
            forces=np.asarray(field['forces_n']);active=np.asarray(field['contact_mask'])
            reconstructed+=float(np.nansum(forces*patch['pinion_lever_arms_m']))
            shares[str(patch['tooth_offset'])]=float(forces.sum())
            projected_area+=float(active.sum())*patch['grid'].cell_area_m2
            first_area+=float(np.nansum(patch['first_surface_area_jacobian'][active]))*patch['grid'].cell_area_m2
            second_area+=float(np.nansum(patch['second_surface_area_jacobian'][active]))*patch['grid'].cell_area_m2
            if active.any():
                active_slope=max(active_slope,float(np.nanmax(patch['surface_slopes'][active])))
                active_angle=max(active_angle,float(np.nanmax(patch['opposing_normal_angles_rad'][active])))
            outer=bool(active[0,:].any() or active[-1,:].any() or active[:,0].any() or active[:,-1].any())
            padded=np.pad(patch['valid_mask'],1,constant_values=False)
            interior=(padded[1:-1,:-2]&padded[1:-1,2:]&padded[:-2,1:-1]&padded[2:,1:-1])
            boundaries={}
            for name in ('model_validity','trimmed_domain','unresolved_query','torque_direction'):
                rejected=np.pad(patch[name+'_rejected_mask'],1,constant_values=False)
                adjacent=(rejected[1:-1,:-2]|rejected[1:-1,2:]|rejected[:-2,1:-1]|rejected[2:,1:-1])
                boundaries['touches_'+name+'_edge']=bool((active&adjacent).any())
            edge.append(dict(tooth_offset=patch['tooth_offset'],touches_window_edge=outer,
                             touches_trimmed_or_model_validity_edge=bool((active&~interior).any()),**boundaries))
        result.update(pinion_torque_reconstructed_nm=reconstructed,
            pinion_torque_relative_residual=(reconstructed-geometry['pinion_torque_nm'])/abs(geometry['pinion_torque_nm']),
            pair_plane_normal_force_shares_n=shares,pair_normal_force_shares_n=shares,contact_boundary_checks=edge,
            active_region_maximum_surface_slope=active_slope,active_region_maximum_opposing_normal_angle_rad=active_angle,
            contact_projected_area_m2=projected_area,contact_actual_first_surface_area_m2=first_area,
            contact_actual_second_surface_area_m2=second_area,
            local_parallel_plane_approximation='参考平面法向合力f=p*Aplane；实际面法向合力与面积含各自1/abs(N·n)雅可比，未混作旧FE面片面积')
        domain_ok=not any(item['touches_window_edge'] or item['touches_model_validity_edge'] or
                          item['touches_unresolved_query_edge'] or item['touches_torque_direction_edge'] for item in edge)
        result['active_contact_domain_accepted']=domain_ok
        if not domain_ok:
            candidate=dict(result)
            result.update(converged=False,status='active_contact_reaches_unaccepted_domain_boundary',candidate=candidate,
                          fields=[],theta_rad=None,maximum_pressure_pa=None,incremental_stiffness_nm_per_rad=None)
            if prepared is not None and any(item['touches_unresolved_query_edge'] for item in edge):
                from .halfspace_boundary import certify_halfspace_contact_boundaries
                result=certify_halfspace_contact_boundaries(prepared,result)
    return result
