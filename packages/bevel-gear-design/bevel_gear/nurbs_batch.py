"""对同一OCP有理B样条作批量精确评价；复用原控制点、权重和结点，不拟合齿面。"""
from __future__ import annotations


class NurbsSurfaceBatch:
    """SciPy B样条基函数批量计算原NURBS及归一化UV导数，返回m。"""
    def __init__(self,surface):
        import numpy as np
        from scipy.interpolate import BSpline
        self.np=np;self.bounds=surface.Bounds()
        nu,nv=surface.NbUPoles(),surface.NbVPoles()
        ku=[surface.UKnot(i) for i in range(1,surface.NbUKnots()+1) for _ in range(surface.UMultiplicity(i))]
        kv=[surface.VKnot(i) for i in range(1,surface.NbVKnots()+1) for _ in range(surface.VMultiplicity(i))]
        self.u=BSpline(ku,np.eye(nu),surface.UDegree(),extrapolate=False)
        self.v=BSpline(kv,np.eye(nv),surface.VDegree(),extrapolate=False)
        self.du=self.u.derivative();self.dv=self.v.derivative()
        self.knots_u=np.asarray(ku);self.knots_v=np.asarray(kv)
        self.degree_u=surface.UDegree();self.degree_v=surface.VDegree()
        self.nu=nu;self.nv=nv
        poles=np.asarray([[surface.Pole(i,j).Coord() for j in range(1,nv+1)] for i in range(1,nu+1)])
        weights=np.asarray([[surface.Weight(i,j) for j in range(1,nv+1)] for i in range(1,nu+1)])
        if not np.isfinite(weights).all() or (weights<=0).any():
            raise ValueError('批量NURBS要求有限正权重')
        self.homogeneous=np.concatenate((poles*weights[:,:,None],weights[:,:,None]),axis=2)

    def sample(self,uv):
        """批量SI点及两个归一化参数导数；不改变原表面的定义域。"""
        np=self.np;uv=np.asarray(uv,dtype=float)
        if uv.ndim!=2 or uv.shape[1]!=2 or not np.isfinite(uv).all() or ((uv<0)|(uv>1)).any():
            raise ValueError('批量NURBS输入须为[0,1]内的有限N×2参数')
        b=self.bounds;u=b[0]+uv[:,0]*(b[1]-b[0]);v=b[2]+uv[:,1]*(b[3]-b[2])
        outputs=[]
        for start in range(0,len(uv),4096):
            us=u[start:start+4096];vs=v[start:start+4096]
            # B样条每方向只有degree+1个非零基函数；不对整张控制网作稠密收缩。
            iu=np.minimum(np.searchsorted(self.knots_u,us,side='right')-1,self.nu-1)[:,None]-self.degree_u+np.arange(self.degree_u+1)
            iv=np.minimum(np.searchsorted(self.knots_v,vs,side='right')-1,self.nv-1)[:,None]-self.degree_v+np.arange(self.degree_v+1)
            indices=np.arange(len(us))[:,None]
            bu=self.u(us)[indices,iu];bv=self.v(vs)[indices,iv]
            du=self.du(us)[indices,iu];dv=self.dv(vs)[indices,iv]
            control=self.homogeneous[iu[:,:,None],iv[:,None,:]]
            value=np.einsum('ia,iabk,ib->ik',bu,control,bv,optimize=False)
            a=np.einsum('ia,iabk,ib->ik',du,control,bv,optimize=False)
            c=np.einsum('ia,iabk,ib->ik',bu,control,dv,optimize=False)
            p=value[:,:3]/value[:,3,None]
            du=(a[:,:3]-p*a[:,3,None])/value[:,3,None]*(b[1]-b[0])
            dv=(c[:,:3]-p*c[:,3,None])/value[:,3,None]*(b[3]-b[2])
            outputs.append((p*.001,du*.001,dv*.001))
        if not outputs:
            empty=np.empty((0,3));return empty,empty.copy(),empty.copy()
        return tuple(np.concatenate([part[k] for part in outputs],axis=0) for k in range(3))


def ray_batch(surface,rotation,origin,basis,target_xy,seed_uv,tolerance_m,maximum_iterations=15):
    """同一有界Newton射线方程的批量实现；每个点独立接受、失败和停止。"""
    import numpy as np
    evaluator=surface if isinstance(surface,NurbsSurfaceBatch) else NurbsSurfaceBatch(surface)
    xy=np.asarray(target_xy,dtype=float);uv=np.clip(np.asarray(seed_uv,dtype=float),1e-12,1-1e-12).copy()
    rotation=np.asarray(rotation);origin=np.asarray(origin);basis=np.asarray(basis)
    n=len(xy);pending=np.ones(n,bool);accepted=np.zeros(n,bool);last=np.full(n,np.inf)
    errors=np.full(n,np.inf);iterations=np.zeros(n,int);statuses=np.full(n,'ray_iteration_limit',dtype=object)
    for iteration in range(maximum_iterations):
        indices=np.flatnonzero(pending)
        if not len(indices):break
        local,du,dv=evaluator.sample(uv[indices]);world=local@rotation.T
        residual=(world-origin)@basis.T-xy[indices]
        error=np.linalg.norm(residual,axis=1);errors[indices]=error;iterations[indices]=iteration+1
        done=error<=tolerance_m;good=indices[done]
        accepted[good]=True;pending[good]=False;statuses[good]='converged'
        keep=~done;ids=indices[keep]
        if not len(ids):continue
        a=(du[keep]@rotation.T)@basis.T;b=(dv[keep]@rotation.T)@basis.T
        det=a[:,0]*b[:,1]-b[:,0]*a[:,1]
        singular=~np.isfinite(det)|(np.abs(det)<1e-16)
        bad=ids[singular];pending[bad]=False;statuses[bad]='ray_graph_degenerate'
        ok=~singular;ids=ids[ok];a=a[ok];b=b[ok];det=det[ok];r=residual[keep][ok]
        if not len(ids):continue
        step=np.column_stack(((b[:,1]*r[:,0]-b[:,0]*r[:,1])/det,
            (-a[:,1]*r[:,0]+a[:,0]*r[:,1])/det))
        fraction=np.minimum(1.,.2/np.maximum(np.max(np.abs(step),axis=1),1e-20))
        trial=np.clip(uv[ids]-fraction[:,None]*step,0.,1.)
        stalled=(np.linalg.norm(trial-uv[ids],axis=1)<1e-13)|((errors[ids]>=last[ids])&(iteration>10))
        bad=ids[stalled];pending[bad]=False;statuses[bad]='outside_or_unresolved_parameter_domain'
        moving=ids[~stalled];uv[moving]=trial[~stalled];last[moving]=errors[moving]
    local,du,dv=evaluator.sample(uv)
    return dict(converged=accepted,status=statuses,residual_m=errors,iterations=iterations,uv=uv,
        local=local,world=local@rotation.T,du_world=du@rotation.T,dv_world=dv@rotation.T)


def sample_pair_batch(prepared,frames,origin,basis,normal,seed_uv,target_xy,tolerance_m,
        maximum_slope,maximum_normal_angle_rad,seed_maps):
    """批量求原双面射线，实际BRep裁切仍由同一个OCP分类器逐点判定。"""
    import math
    import numpy as np
    xy=np.asarray(target_xy,dtype=float);normal=np.asarray(normal);basis=np.asarray(basis)
    members={};reasons=[None]*len(xy)
    for member in ('pinion','wheel'):
        trimmed=prepared.trimmed[member];rotation=np.asarray(frames[member+'_rotation_rows'])
        if not hasattr(trimmed,'_batch_evaluator'):
            trimmed._batch_evaluator=NurbsSurfaceBatch(trimmed.surface)
        coordinate,inverse=seed_maps[member]
        guesses=np.asarray(seed_uv[member])+(xy-coordinate)@inverse.T
        values=ray_batch(trimmed._batch_evaluator,rotation,origin,basis,xy,guesses,tolerance_m)
        outward=np.cross(values['du_world'],values['dv_world']);length=np.linalg.norm(outward,axis=1)
        if np.any(length<1e-22):raise ArithmeticError('批量实际曲面参数导数退化')
        outward*=trimmed.outward_sense/length[:,None]
        cosine=outward@normal
        slopes=np.linalg.norm(outward@basis.T,axis=1)/np.maximum(np.abs(cosine),1e-300)
        for i in range(len(xy)):
            if reasons[i] is not None:continue
            if not values['converged'][i]:reasons[i]=member+':'+values['status'][i]
            elif not trimmed.inside(values['uv'][i],values['local'][i]):reasons[i]=member+':outside_trimmed_face'
            elif slopes[i]>maximum_slope:reasons[i]=member+':halfspace_slope_limit'
        values.update(outward=outward,slope=slopes);members[member]=values
    samples=[]
    for i in range(len(xy)):
        if reasons[i] is not None:samples.append(None);continue
        first,second=members['pinion'],members['wheel']
        cosine=float(-first['outward'][i]@second['outward'][i])
        if cosine<math.cos(maximum_normal_angle_rad):
            reasons[i]='opposing_normal_limit';samples.append(None);continue
        arms={};moments={}
        for member,sign in (('pinion',1),('wheel',-1)):
            value=members[member];outward=value['outward'][i]
            velocity=np.cross(frames[member+'_axis_world'],value['world'][i])
            arms[member]=sign*float(outward@velocity/(outward@normal))
            moments[member]=sign*float(normal@velocity)
        def one(value):
            return {key:array[i] for key,array in value.items() if key!='status'}
        samples.append(dict(gap_m=float((second['world'][i]-first['world'][i])@normal),
            first=one(first),second=one(second),lever_arms_m=arms,normal_moment_arms_m=moments,
            maximum_slope=max(first['slope'][i],second['slope'][i]),normal_alignment=cosine))
    return samples,reasons
