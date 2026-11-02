"""压力的工程后处理：源面元法向力守恒恢复、随接触轨迹的相位插值。

不求解新的接触平衡，不覆盖原始压力。恢复场与插值包络均须单独标注，
不能用显示峰值替代原求解峰值或据此声称网格收敛。
"""
from __future__ import annotations

import math
import numpy as np


def recover_grid_pressure(pressure_pa, valid_mask, cell_area_m2, spacing_m,
        smoothing_length_m):
    """在原参考面元上作非负高斯恢复，并严格保持该齿原法向合力。

    spacing_m为(dx,dy)，单位m；数组索引为[y,x]。核宽采用物理长度，
    不是图像RGB模糊半径。有效域边缘留一格零压支承，保证中心P1积分。
    这会改变显示压力及其峰值，不改变原解、位移、TE或刚度。
    """
    from scipy.ndimage import gaussian_filter, binary_erosion
    p=np.asarray(pressure_pa,dtype=float); valid=np.asarray(valid_mask,dtype=bool)
    dx,dy=map(float,spacing_m); area=float(cell_area_m2); length=float(smoothing_length_m)
    if p.ndim!=2 or valid.shape!=p.shape or not np.isfinite(p).all() or np.any(p<0):
        raise ValueError('压力须为有限非负二维SI数组，有效域布局须一致')
    if any(not math.isfinite(v) or v<=0 for v in (dx,dy,area,length)):
        raise ValueError('格距、面元面积和恢复长度须为正有限SI数值')
    if np.any(p[~valid]!=0): raise ValueError('原有效域外不得有非零压力')
    force=float(p.sum()*area)
    if force==0:return p.copy(),dict(original_force_n=0.,recovered_force_n=0.,relative_force_error=0.,raw_peak_pa=0.,recovered_peak_pa=0.)
    inside=binary_erosion(valid,structure=np.ones((3,3)),border_value=0)
    q=gaussian_filter(p,sigma=(length/dy,length/dx),mode='constant',cval=0.,truncate=3.)
    q=np.where(inside,np.maximum(q,0.),0.)
    before=float(q.sum()*area)
    if before<=0:raise ValueError('恢复后没有完整支承域，须保留原图或扩展实际输入域')
    q*=force/before
    yy,xx=np.indices(p.shape); raw_center=np.array([(p*xx).sum(),(p*yy).sum()])/p.sum()
    new_center=np.array([(q*xx).sum(),(q*yy).sum()])/q.sum()
    return q,dict(original_force_n=force,recovered_force_n=float(q.sum()*area),
        relative_force_error=float(q.sum()*area/force-1.),raw_peak_pa=float(p.max()),
        recovered_peak_pa=float(q.max()),smoothing_length_m=length,
        sigma_rows_columns=[length/dy,length/dx],force_centroid_shift_m=((new_center-raw_center)*[dx,dy]).tolist(),
        original_pressure_modified=False,recovered_field_is_new_contact_solution=False)


def grid_triangle_field(source, pressure_pa):
    """在真实材料坐标上重建中心P1；参考积分面积仍用原cell_area/2。"""
    p=np.asarray(pressure_pa,float); valid=np.asarray(source['valid_mask'],bool)
    xy=np.asarray(source['points_width_height_m'],float)
    if xy.shape!=p.shape+(2,) or valid.shape!=p.shape:raise ValueError('材料坐标/压力格布局不同')
    ny,nx=p.shape; triangles=[]
    for j in range(ny-1):
        for i in range(nx-1):
            for corners in (((j,i),(j,i+1),(j+1,i+1)),((j,i),(j+1,i+1),(j+1,i))):
                if all(valid[q] for q in corners) and any(p[q]>0 for q in corners):
                    triangles.append([y*nx+x for y,x in corners])
    if not triangles:return None
    full=np.asarray(triangles,int);used=np.unique(full);inverse=np.full(ny*nx,-1,int);inverse[used]=np.arange(len(used))
    t=inverse[full];values=p.ravel()[used];points=xy.reshape(-1,2)[used]
    a=points[t[:,1]]-points[t[:,0]];b=points[t[:,2]]-points[t[:,0]]
    signed=a[:,0]*b[:,1]-a[:,1]*b[:,0]
    t[signed<0]=t[signed<0][:,[0,2,1]]
    integral=float(values[t].mean(axis=1).sum()*source['cell_area_m2']/2)
    original=float(p.sum()*source['cell_area_m2'])
    if not math.isclose(integral,original,rel_tol=2e-10,abs_tol=1e-9):
        raise ValueError('恢复压力贴近输入域边界，中心P1显示积分不守恒')
    weights=p*source['cell_area_m2'];center=(xy*weights[...,None]).sum(axis=(0,1))/weights.sum()
    return dict(points_width_height_m=points,triangles=t,pressure_pa=values,
        driver_angle_rad=float(source['driver_angle_rad']),phase_index=int(source['phase_index']),
        tooth_offset=int(source['tooth_offset']),source_state_sha256=source['source_state_sha256'],
        center_m=center,reference_force_n=original,reconstructed_reference_integral_n=integral)


def raster_grid(outline, shape=(161,641)):
    """规则SI显示格；实际轮廓以几何点包含关系裁切，齿面外始终掩码。"""
    import shapely
    xmin,ymin,xmax,ymax=outline.bounds
    x=np.linspace(xmin,xmax,shape[1]);y=np.linspace(ymin,ymax,shape[0]);xx,yy=np.meshgrid(x,y)
    mask=shapely.intersects_xy(outline,xx,yy)
    return dict(x_m=x,y_m=y,xx_m=xx,yy_m=yy,flank_mask=mask)


def sample_triangle_field(field, grid):
    """仅在输入原三角形内线性采样，不跨原孔洞或缺失三角形外推。"""
    from matplotlib.tri import Triangulation,LinearTriInterpolator
    xy=np.asarray(field['points_width_height_m'],float);tri=np.asarray(field['triangles'],int)
    values=np.asarray(field['pressure_pa'],float)
    mesh=Triangulation(xy[:,0],xy[:,1],tri)
    result=np.zeros_like(grid['xx_m'])
    # 只访问该真实三角域的包围盒，不在整个齿面反复搜索每个小接触斑。
    low=xy.min(axis=0);high=xy.max(axis=0)
    x0=max(0,int(np.searchsorted(grid['x_m'],low[0])-1));x1=min(result.shape[1],int(np.searchsorted(grid['x_m'],high[0]))+1)
    y0=max(0,int(np.searchsorted(grid['y_m'],low[1])-1));y1=min(result.shape[0],int(np.searchsorted(grid['y_m'],high[1]))+1)
    sampled=LinearTriInterpolator(mesh,values)(grid['xx_m'][y0:y1,x0:x1],grid['yy_m'][y0:y1,x0:x1])
    result[y0:y1,x0:x1]=np.maximum(np.asarray(sampled.filled(0.)),0.)
    return result


def raw_raster(data, grid):
    """离散原相位最大值；不做时间插值或压力恢复。"""
    result=np.zeros_like(grid['xx_m'])
    for field in data['source_fields']:np.maximum(result,sample_triangle_field(field,grid),out=result)
    return result


def advected_cycle_maximum(fields, phases_rad, grid, *, subdivisions=16):
    """邻接真实相位间作接触中心平移及压力凸组合，输出工程包络。

    p(t,x)=(1-t)p_i(x-tΔc)+t p_j(x+(1-t)Δc)。只连接同一物理齿、
    在相邻原相位均承载的场；不跨缺失/零载相位连接，也不补造求解状态。
    除中心平移和线性压力插值外，不拟合椭圆、不钳制峰值。
    """
    from scipy.ndimage import shift
    phases=list(map(float,phases_rad))
    if len(phases)<2 or np.any(np.diff(phases)<=0) or subdivisions<1:
        raise ValueError('原相位须严格升序，区间细分数须为正')
    dx=grid['x_m'][1]-grid['x_m'][0];dy=grid['y_m'][1]-grid['y_m'][0]
    result=np.zeros_like(grid['xx_m']);lookup={};source_peaks=[]
    for field in fields:
        angle=field['driver_angle_rad']
        if angle not in phases:raise ValueError('承载场不在原真实相位表内')
        key=(field['tooth_offset'],angle)
        if key in lookup:raise ValueError('同齿同相位场重复')
        raster=sample_triangle_field(field,grid);np.maximum(result,raster,out=result)
        lookup[key]=(raster,np.asarray(field['center_m'],float));source_peaks.append(float(raster.max()))
    segments=0;interpolated=0
    for tooth in sorted({k[0] for k in lookup}):
        for a,b in zip(phases[:-1],phases[1:]):
            if (tooth,a) not in lookup or (tooth,b) not in lookup:continue
            left,c0=lookup[(tooth,a)];right,c1=lookup[(tooth,b)];delta=(c1-c0)/[dx,dy]
            segments+=1
            for i in range(1,subdivisions):
                t=i/subdivisions
                first=shift(left,shift=(t*delta[1],t*delta[0]),order=1,mode='constant',cval=0.,prefilter=False)
                second=shift(right,shift=(-(1-t)*delta[1],-(1-t)*delta[0]),order=1,mode='constant',cval=0.,prefilter=False)
                np.maximum(result,(1-t)*first+t*second,out=result);interpolated+=1
    if source_peaks and result.max()>max(source_peaks)*(1+1e-12):raise ArithmeticError('凸插值产生了超出端点的峰值')
    return result,dict(solved_phase_count=len(phases),display_phase_count=(len(phases)-1)*subdivisions+1,
        connected_same_tooth_segments=segments,interpolated_field_count=interpolated,
        new_FE_solves=0,pressure_geometry_interpolation='接触中心随相位平移，端点标量场线性凸组合',
        continuous_phase_extremum_certified=False,interpolated_equilibrium_verified=False)


def recover_sequence(sequence, grid, *, subdivisions=16, smoothing_length_m=None):
    """恢复原9态而非借用半解析结果；返回瞬时恢复图与相位插值包络。"""
    fields=[];diagnostics=[];instant=np.zeros_like(grid['xx_m']);length=sequence['smoothing_length_m'] if smoothing_length_m is None else smoothing_length_m
    for source in sequence['fields']:
        q,proof=recover_grid_pressure(source['pressure_pa'],source['valid_mask'],source['cell_area_m2'],sequence['grid_spacing_m'],length)
        proof.update(phase_index=source['phase_index'],tooth_offset=source['tooth_offset'])
        diagnostics.append(proof)
        if not np.any(q):continue
        field=grid_triangle_field(source,q);fields.append(field)
        if source['driver_angle_rad']==0. and source['tooth_offset']==0:
            instant=sample_triangle_field(field,grid)
    envelope,temporal=advected_cycle_maximum(fields,sequence['driver_angles_rad'],grid,subdivisions=subdivisions)
    return instant,envelope,dict(spatial_recovery=diagnostics,temporal_recovery=temporal,
        smoothing_length_m=length,original_source_fields_unchanged=True,display_is_postprocessed=True,
        interpretation='源面元合力守恒；相位插值包络仅供工程展示，不是新增FE平衡解或精度验证')
