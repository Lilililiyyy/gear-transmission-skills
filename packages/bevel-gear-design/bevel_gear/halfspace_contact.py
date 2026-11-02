"""矩形常压力单元的半解析弹性半空间接触，按Boussinesq方程新写。

两体合成柔度为(1/E*)/π乘源矩形内1/r的解析面积积分。未知量为
单元力f=p*dx*dy；规则网格用零填充FFT实现非周期线性卷积。多patch
各自局部平面、半空间柔度块互不耦合，不含齿弯曲、轮体或系统柔度。
非负力与给定合力/有符号力矩同时求解，不插值已有粗网格压力。

数学定义见Popov等《Normal Contact Without Adhesion》§2.2和Wang等
《FFT-Based Methods for Computational Contact Mechanics》§2/§4：
https://doi.org/10.1007/978-3-662-58709-6_2
https://doi.org/10.3389/fmech.2020.00061
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math


def _numpy():
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError('半解析表面接触需要可选NumPy依赖') from error
    return np


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int,float)) or not math.isfinite(value):
        raise ValueError(name+'须为有限实数')
    if positive and value<=0:raise ValueError(name+'须为正数')
    return float(value)


@dataclass(frozen=True)
class RectangularContactGrid:
    """局部接触平面规则网格；origin为第一个单元中心，长度m。

    数组布局始终为(y_count,x_count)，每格物理参考平面面积为dx*dy。
    该面积不能自动当成有明显坡度的曲面真实面积。
    """
    x_count: int
    y_count: int
    spacing_x_m: float
    spacing_y_m: float
    origin_xy_m: tuple = (0.,0.)

    def __post_init__(self):
        for name in ('x_count','y_count'):
            value=getattr(self,name)
            if isinstance(value,bool) or not isinstance(value,int) or value<1:
                raise ValueError(name+'须为正整数')
        for name in ('spacing_x_m','spacing_y_m'):_number(getattr(self,name),name,positive=True)
        if not math.isfinite(self.spacing_x_m*self.spacing_y_m) or self.spacing_x_m*self.spacing_y_m<=0:
            raise ValueError('单元面积须为有限正值')
        if len(self.origin_xy_m)!=2:raise ValueError('origin_xy_m须为二维坐标')
        object.__setattr__(self,'origin_xy_m',tuple(_number(v,'origin_xy_m') for v in self.origin_xy_m))

    @property
    def shape(self):return (self.y_count,self.x_count)

    @property
    def cell_area_m2(self):return self.spacing_x_m*self.spacing_y_m

    def centres_m(self):
        """返回x、y二维中心坐标；只在调用时加载NumPy。"""
        np=_numpy()
        return np.meshgrid(self.origin_xy_m[0]+np.arange(self.x_count)*self.spacing_x_m,
                           self.origin_xy_m[1]+np.arange(self.y_count)*self.spacing_y_m)

    def as_dict(self):
        return dict(x_count=self.x_count,y_count=self.y_count,spacing_x_m=self.spacing_x_m,
                    spacing_y_m=self.spacing_y_m,origin_xy_m=list(self.origin_xy_m))


def rectangular_potential_integral(dx_m,dy_m,half_width_m,half_height_m):
    """解析积分∫[-a,a]∫[-b,b]1/sqrt((x-dx)^2+(y-dy)^2) dxdy，单位m。

    采用四角原函数 x*asinh(y/|x|)+y*asinh(x/|y|)，坐标为0时取
    连续极限，不以小圆半径替换自作用奇点；输入dx/dy支持广播数组。
    """
    np=_numpy()
    a=_number(half_width_m,'half_width_m',positive=True)
    b=_number(half_height_m,'half_height_m',positive=True)
    rawx,rawy=np.asarray(dx_m),np.asarray(dy_m)
    if rawx.dtype.kind not in 'iuf' or rawy.dtype.kind not in 'iuf':raise ValueError('矩形相对坐标须为实数，不接受布尔/复数/字符串')
    dx,dy=np.broadcast_arrays(rawx.astype(float),rawy.astype(float))
    if not np.isfinite(dx).all() or not np.isfinite(dy).all():raise ValueError('矩形相对坐标须有限')
    def primitive(x,y):
        first=np.zeros_like(x);second=np.zeros_like(x)
        np.divide(y,np.abs(x),out=first,where=x!=0)
        np.divide(x,np.abs(y),out=second,where=y!=0)
        return x*np.arcsinh(first)+y*np.arcsinh(second)
    value=(primitive(dx+a,dy+b)-primitive(dx-a,dy+b)
           -primitive(dx+a,dy-b)+primitive(dx-a,dy-b))
    if np.any(value<=0) or not np.isfinite(value).all():
        raise ArithmeticError('矩形积分失去正值或有限性，未夹零修补')
    return value


@dataclass
class HalfSpaceOperator:
    """两体合成半空间柔度；apply输入N，返回m，不装配全部稠密矩阵。"""
    grid: RectangularContactGrid
    inverse_effective_modulus_per_pa: float
    kernel_m_per_n: object = field(repr=False)
    kernel_fft: object = field(repr=False)
    self_compliance_m_per_n: float
    matvec_count: int = 0

    def apply(self, forces_n):
        np=_numpy();raw=np.asarray(forces_n)
        if raw.dtype.kind not in 'iuf':raise ValueError('半空间作用力须为实数数组')
        values=raw.astype(float)
        if values.shape!=self.grid.shape or not np.isfinite(values).all():
            raise ValueError('半空间作用力须为与grid一致的有限二维数组')
        padded=np.zeros(self.kernel_m_per_n.shape)
        padded[:self.grid.y_count,:self.grid.x_count]=values
        answer=np.fft.irfft2(np.fft.rfft2(padded)*self.kernel_fft,s=padded.shape)
        self.matvec_count+=1
        return answer[:self.grid.y_count,:self.grid.x_count]

    def dense_matrix(self, *, maximum_cells=2048):
        """小系统独立对照用的精确Toeplitz矩阵，不用于大网格生产求解。"""
        np=_numpy();count=self.grid.x_count*self.grid.y_count
        if count>maximum_cells:raise ValueError('稠密对照超过显式单元上限')
        yy,xx=np.indices(self.grid.shape);x=xx.ravel();y=yy.ravel()
        return self.kernel_m_per_n[(y[:,None]-y[None,:])%(2*self.grid.y_count),
                                   (x[:,None]-x[None,:])%(2*self.grid.x_count)]

    def as_dict(self):
        return dict(grid=self.grid.as_dict(),inverse_effective_modulus_per_pa=self.inverse_effective_modulus_per_pa,
            self_compliance_m_per_n=self.self_compliance_m_per_n,matvec_count=self.matvec_count,
            model='两体弹性半空间、矩形常压力源单元中心配点',
            units=dict(force='N',pressure='Pa',gap='m',compliance='m/N'))


def prepare_halfspace_operator(grid,pinion_young_modulus_pa,pinion_poisson_ratio,
        wheel_young_modulus_pa,wheel_poisson_ratio,*,maximum_cells=1048576):
    """由两体E/ν构造半解析影响系数，不使用经验刚度或对角正则。"""
    np=_numpy()
    if isinstance(grid,dict):grid=RectangularContactGrid(**grid)
    if not isinstance(grid,RectangularContactGrid):raise TypeError('grid须为RectangularContactGrid或其字典')
    if isinstance(maximum_cells,bool) or not isinstance(maximum_cells,int) or maximum_cells<1:
        raise ValueError('maximum_cells须为正整数')
    if grid.x_count*grid.y_count>maximum_cells:raise ValueError('规则接触网格超过显式单元上限')
    inverse=0.
    for young,poisson in ((pinion_young_modulus_pa,pinion_poisson_ratio),(wheel_young_modulus_pa,wheel_poisson_ratio)):
        young=_number(young,'弹性模量Pa',positive=True);poisson=_number(poisson,'泊松比')
        if not -1<poisson<.5:raise ValueError('各向同性弹性泊松比须在(-1,.5)')
        inverse+=(1-poisson*poisson)/young
    if not math.isfinite(inverse) or inverse<=0:raise ValueError('材料合成柔度超出有限数值域')
    ix=np.arange(2*grid.x_count);iy=np.arange(2*grid.y_count)
    ix=np.where(ix<grid.x_count,ix,ix-2*grid.x_count)
    iy=np.where(iy<grid.y_count,iy,iy-2*grid.y_count)
    x,y=np.meshgrid(ix*grid.spacing_x_m,iy*grid.spacing_y_m)
    kernel=rectangular_potential_integral(x,y,grid.spacing_x_m/2,grid.spacing_y_m/2)
    kernel*=inverse/(math.pi*grid.cell_area_m2)
    if not np.isfinite(kernel).all() or np.any(kernel<=0):
        raise ArithmeticError('影响系数超出有限正数域，未正则修补')
    # 差值±count不属于真实中心对；置零的仅是线性卷积零填充隔离行/列。
    kernel[grid.y_count,:]=0.;kernel[:,grid.x_count]=0.
    return HalfSpaceOperator(grid,inverse,kernel,np.fft.rfft2(kernel),float(kernel[0,0]))


def _project_weighted_simplex(values, weights):
    """欧氏投影到q>=0、aᵀq=1，权重a严格为正；不改变给定力臂。"""
    np=_numpy()
    ratio=values/weights
    order=np.argsort(ratio)[::-1]
    sorted_ratio=ratio[order]
    limits=(np.cumsum(weights[order]*values[order])-1.)/np.cumsum(weights[order]**2)
    eligible=np.flatnonzero(sorted_ratio>limits)
    if not len(eligible):raise ArithmeticError('非负等式投影未找到可行支集')
    threshold=limits[eligible[-1]]
    projected=np.maximum(values-threshold*weights,0.)
    if abs(float(weights@projected)-1.)>1e-10:
        raise ArithmeticError('非负等式投影载荷残差超限')
    return projected


def solve_halfspace_contact(gaps_m,operators,*,load_n=None,lever_arms_m=None,torque_nm=None,
        valid_masks=None,maximum_iterations=20000,gap_tolerance_m=1e-10,
        relative_tolerance=1e-7,force_tolerance_n=1e-8,allow_initial_penetration=False,
        compute_stiffness=True,check_interval=20):
    """矩阵自由凸接触QP，求非负单元力与合力或有符号力矩平衡。

    g+C*f−ell*theta>=0、f>=0、fᵀ(g+C*f−ell*theta)=0。总法向力
    模式ell=1，闭合量单位m；力矩模式ell单位m，theta单位rad，必须由
    同一几何/轴符号给出。多patch仅共享一个广义闭合量，不含跨patch
    结构柔度。有效域外gap可以未知NaN，必须提供False掩码并保留未知。
    失败压力/力不作为物理输出，只在candidate中保留末次试算。
    """
    np=_numpy()
    single=isinstance(operators,HalfSpaceOperator)
    ops=[operators] if single else list(operators)
    if not ops or any(not isinstance(op,HalfSpaceOperator) for op in ops):raise TypeError('operators须为半空间算子或其非空列表')
    gaps=[gaps_m] if single else list(gaps_m)
    if len(gaps)!=len(ops):raise ValueError('gap与operator的patch数须一致')
    if (load_n is None)==(torque_nm is None):raise ValueError('须且仅须给load_n或torque_nm')
    if load_n is not None:
        target=_number(load_n,'法向合力N',positive=True);direction=1;mode='normal_force';lever_scale=1.
        if lever_arms_m is not None:raise ValueError('法向合力模式不接收力臂')
        levers=[np.ones(op.grid.shape) for op in ops]
    else:
        target=_number(torque_nm,'有符号力矩Nm')
        if target==0:raise ValueError('力矩须非零')
        if lever_arms_m is None:raise ValueError('力矩模式须显式给同符号/同几何力臂')
        direction=1 if target>0 else -1;mode='torque'
        levers=[lever_arms_m] if single else list(lever_arms_m)
        if len(levers)!=len(ops):raise ValueError('力臂patch数不一致')
    masks=[np.ones(op.grid.shape,dtype=bool) for op in ops] if valid_masks is None else ([valid_masks] if single else list(valid_masks))
    if len(masks)!=len(ops):raise ValueError('有效掩码patch数不一致')
    if not isinstance(allow_initial_penetration,bool) or not isinstance(compute_stiffness,bool):raise ValueError('策略开关须为布尔值')
    for value,name in ((maximum_iterations,'maximum_iterations'),(check_interval,'check_interval')):
        if isinstance(value,bool) or not isinstance(value,int) or value<1:raise ValueError(name+'须为正整数')
    for value,name in ((gap_tolerance_m,'gap_tolerance_m'),(relative_tolerance,'relative_tolerance'),(force_tolerance_n,'force_tolerance_n')):
        _number(value,name,positive=True)
    prepared=[]
    for gap,lever,mask,op in zip(gaps,levers,masks,ops):
        rawgap,rawlever=np.asarray(gap),np.asarray(lever)
        if rawgap.dtype.kind not in 'iuf' or rawlever.dtype.kind not in 'iuf':
            raise ValueError('gap和力臂须为实数数组，不接受布尔/复数/字符串')
        gap=rawgap.astype(float);lever=rawlever.astype(float);mask=np.asarray(mask)
        if gap.shape!=op.grid.shape or lever.shape!=op.grid.shape or mask.shape!=op.grid.shape or mask.dtype!=np.bool_:
            raise ValueError('gap/力臂/布尔掩码须与规则网格二维布局一致')
        if not np.isfinite(gap[mask]).all() or not np.isfinite(lever[mask]).all():raise ValueError('有效接触单元的gap和力臂须有限')
        if np.any(direction*lever[mask]<=0):raise ValueError('有效力臂须能承受给定力矩方向；不能取绝对值修补')
        if not allow_initial_penetration and np.any(gap[mask]<0):raise ValueError('负初始gap须显式允许，未夹零')
        prepared.append((gap,lever,mask))
    if not any(mask.any() for _,_,mask in prepared):raise ValueError('没有有效接触单元')
    if mode=='torque':lever_scale=max(float(np.max(np.abs(lever[mask]))) if mask.any() else 0. for _,lever,mask in prepared)
    force_scale=abs(target)/lever_scale
    compliance_scale=max(op.self_compliance_m_per_n for op in ops)
    displacement_scale=compliance_scale*force_scale
    if not math.isfinite(force_scale) or not math.isfinite(displacement_scale) or displacement_scale<=0 or not math.isfinite(force_scale*displacement_scale):
        raise ValueError('载荷/柔度的位移或能量尺度超出有限数值域')
    lengths=[int(mask.sum()) for _,_,mask in prepared];ends=np.cumsum([0,*lengths])
    a=np.concatenate([direction*lever[mask]/lever_scale for _,lever,mask in prepared])
    if np.any(a<1e-12):raise ValueError('有向力臂相对最大值低于1e-12，等式投影病态')
    g=np.concatenate([gap[mask]/displacement_scale for gap,_,mask in prepared])
    if not np.isfinite(g).all():raise ValueError('初始gap相对于载荷位移尺度超出有限数值域')
    def operate(vector):
        answer=[]
        for k,(op,(_,_,mask)) in enumerate(zip(ops,prepared)):
            forces=np.zeros(op.grid.shape);forces[mask]=vector[ends[k]:ends[k+1]]
            answer.append(op.apply(forces)[mask]/compliance_scale)
        return np.concatenate(answer)
    bound=float(np.max(operate(np.ones(len(g)))))
    if not math.isfinite(bound) or bound<=0:raise ArithmeticError('半空间柔度算子谱界不合法')
    q=_project_weighted_simplex(-g/bound,a);y=q.copy();acceleration=1.
    certificate=None;converged=False;iterations=0
    for iteration in range(1,maximum_iterations+1):
        gradient=operate(y)+g
        updated=_project_weighted_simplex(y-gradient/bound,a)
        next_acceleration=(1+math.sqrt(1+4*acceleration*acceleration))/2
        extrapolated=updated+(acceleration-1)/next_acceleration*(updated-q)
        # 梯度方向重启只影响迭代路径，不增加柔度或修改原gap。
        if float((y-updated)@(updated-q))>0:
            next_acceleration=1.;extrapolated=updated.copy()
        q=updated;y=extrapolated;acceleration=next_acceleration;iterations=iteration
        if iteration%check_interval and iteration!=maximum_iterations:continue
        h=operate(q)+g
        multiplier=float(q@h)  # aᵀq=1，对应当前接触力加权的共同闭合量。
        remaining=(h-a*multiplier)*displacement_scale
        forces=q*force_scale
        active=forces>max(force_tolerance_n,float(forces.max())*relative_tolerance*.1)
        active_error=float(np.max(np.abs(remaining[active]))) if active.any() else math.inf
        minimum=float(remaining.min())
        work=float(forces@(h*displacement_scale))
        complementarity=float(np.max(np.abs(forces*remaining)))
        twice_energy=float(forces@((h-g)*displacement_scale))
        minimum_force=float(forces.min())
        # 虚功恒等式：2U+fᵀg=Q*θ；只检查，不藉此修改压力或闭合量。
        energy_residual=twice_energy+float(forces@(g*displacement_scale))-multiplier*force_scale*displacement_scale
        gap_limit=gap_tolerance_m+relative_tolerance*max(float(np.max(np.abs((h-g)*displacement_scale))),abs(multiplier*displacement_scale),1e-20)
        complementarity_limit=gap_limit*max(float(forces.max()),1.)+relative_tolerance*max(abs(work),1e-20)
        balance=float(a@q)-1.
        certificate=dict(minimum_remaining_gap_m=minimum,active_stationarity_error_m=active_error,
            complementarity_residual_nm=complementarity,relative_balance_residual=balance,
            gap_acceptance_m=gap_limit,complementarity_acceptance_nm=complementarity_limit,
            minimum_force_n=minimum_force,strain_energy_j=.5*twice_energy,
            energy_balance_residual_j=energy_residual)
        energy_ok=math.isfinite(twice_energy) and twice_energy>0 and math.isfinite(energy_residual)
        if minimum>=-gap_limit and active_error<=gap_limit and complementarity<=complementarity_limit and abs(balance)<=relative_tolerance and minimum_force>=-force_tolerance_n and energy_ok:
            converged=True;break
    h=operate(q)+g;multiplier=float(q@h);forces=q*force_scale
    deflections=(h-g)*displacement_scale
    remaining=(h-a*multiplier)*displacement_scale
    theta=direction*multiplier*displacement_scale/lever_scale
    stiffness=None;stiffness_status='contact_not_converged' if compute_stiffness and not converged else 'not_requested';stiffness_iterations=0
    active=forces>max(force_tolerance_n,float(forces.max())*relative_tolerance*.1)
    if converged and compute_stiffness:
        try:
            from scipy.sparse.linalg import LinearOperator,cg
        except ImportError:
            stiffness_status='requires_scipy'
        else:
            positions=np.flatnonzero(active)
            def restricted(value):
                whole=np.zeros(len(g));whole[positions]=value
                return operate(whole)[positions]
            matrix=LinearOperator((len(positions),len(positions)),matvec=restricted,dtype=float)
            def callback(_):
                nonlocal stiffness_iterations
                stiffness_iterations+=1
            import inspect
            tolerance={'rtol':1e-9} if 'rtol' in inspect.signature(cg).parameters else {'tol':1e-9}
            solution,info=cg(matrix,a[positions],atol=0.,maxiter=3000,callback=callback,**tolerance)
            if info==0:
                stiffness=lever_scale**2/compliance_scale*float(a[positions]@solution)
                stiffness_status='converged' if stiffness>0 and math.isfinite(stiffness) else 'nonpositive_or_nonfinite'
                if stiffness_status!='converged':stiffness=None
            else:stiffness_status='iteration_limit' if info>0 else 'linear_solver_failure'
    fields=[]
    for k,(op,(gap,lever,mask)) in enumerate(zip(ops,prepared)):
        force_grid=np.zeros(op.grid.shape);force_grid[mask]=forces[ends[k]:ends[k+1]]
        gap_grid=np.full(op.grid.shape,np.nan);gap_grid[mask]=remaining[ends[k]:ends[k+1]]
        displacement=np.zeros(op.grid.shape);displacement[mask]=deflections[ends[k]:ends[k+1]]
        contact=np.zeros(op.grid.shape,dtype=bool);contact[mask]=active[ends[k]:ends[k+1]]
        fields.append(dict(grid=op.grid.as_dict(),forces_n=force_grid,pressure_pa=force_grid/op.grid.cell_area_m2,
            remaining_gap_m=gap_grid,elastic_deflection_m=displacement,contact_mask=contact,valid_mask=mask.copy(),
            active_cell_count=int(contact.sum()),cell_area_m2=op.grid.cell_area_m2))
    resultant=float(sum(np.sum(field['forces_n'][mask]*lever[mask]) for field,(_,lever,mask) in zip(fields,prepared)))
    energy=.5*float(forces@deflections)
    candidate=None if converged else dict(fields=fields,theta=theta,certificate=certificate,
                                          interpretation='末次候选，不能作为已接受接触压力')
    return dict(converged=converged,status='converged' if converged else 'iteration_limit',
        model='半解析Boussinesq矩形常压力、无摩擦弹性半空间、共享广义闭合量',
        mode=mode,fields=fields if converged else [],candidate=candidate,iterations=iterations,
        target_load=target,resultant_load=resultant if converged else None,
        balance_residual=resultant-target if converged else None,load_unit='N' if mode=='normal_force' else 'N·m',
        normal_approach_m=theta if converged and mode=='normal_force' else None,
        theta_rad=theta if converged and mode=='torque' else None,
        strain_energy_j=energy if converged else None,
        total_normal_force_n=float(forces.sum()) if converged else None,
        maximum_pressure_pa=max(float(field['pressure_pa'].max()) for field in fields) if converged else None,
        active_cell_count=int(active.sum()) if converged else None,
        incremental_stiffness_n_per_m=stiffness if mode=='normal_force' else None,
        incremental_stiffness_nm_per_rad=stiffness if mode=='torque' else None,
        stiffness_status=stiffness_status,stiffness_iterations=stiffness_iterations,
        certificate=certificate,operators=[op.as_dict() for op in ops],
        diagnostics=dict(solver='FISTA非负加权等式投影，正谱行和上界；活跃支集CG切线',
            normalized_operator_bound=bound,force_scale_n=force_scale,compliance_scale_m_per_n=compliance_scale,
            displacement_scale_m=displacement_scale,allow_initial_penetration=allow_initial_penetration,
            minimum_initial_gap_m=min(float(gap[mask].min()) if mask.any() else math.inf for gap,_,mask in prepared)),
        limitations=['局部平面、小坡度、各向同性线弹性半空间；不含有限齿厚/齿根/轮体柔度',
            '多个patch的柔度为块对角，不包含跨齿和系统结构耦合',
            '规则格压力为常单元配点离散；迭代通过不代表几何/压力网格已收敛'])
