"""新写实现：真实曲面P0分布牵引与整轮FE功共轭接触。

压力为参考投影矩形内P0未知量，以Gauss积分重新查询实际有限NURBS。
每单位平面合力按FE形函数分布，C=W.T K^-1 W，不使用半空间柔度、
不平滑旧节点压力。弱非穿透约束为同一压力小片的积分平均，逐Gauss点
剩余间隙另报；离散KKT、Gauss残差、几何重投影与网格收敛分别报告。
"""
from __future__ import annotations

import math


def _positive_real(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
        raise ValueError(name+'须为正有限数')
    return float(value)


def _working_facets(prepared,member,tooth_offset):
    """同一不可变整轮体网格复用边界拓扑；物理齿索引不换成参考齿。"""
    from .contact_mesh import boundary_facets
    from .generated_ltca import select_generated_flank_facets
    cache=prepared._facet_cache;boundary_key=('distributed_boundary',member)
    if boundary_key not in cache:cache[boundary_key]=boundary_facets(prepared.systems[member].mesh)
    model=prepared.models[member];tooth=(-tooth_offset if member=='pinion' else tooth_offset)%model.request.tooth_count
    key=('distributed_working_flank',member,tooth)
    if key not in cache:
        side='left_flank' if member=='pinion' else 'right_flank'
        cache[key]=select_generated_flank_facets(prepared.systems[member].mesh,model.shell.surfaces[side],
            tooth_count=model.request.tooth_count,tooth_index=tooth,boundary=cache[boundary_key],matching_tolerance_m=8e-7)
    return cache[key]


class _SurfaceTraceLocator:
    """把已核实的实际CAD点映至同物理齿的真实FE边界形函数。"""
    def __init__(self, facets, rotation_rows):
        import numpy as np
        import vtk
        from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray
        self.np=np;self.vtk=vtk;self.rotation=np.asarray(rotation_rows)
        self.triangles=np.asarray(facets['triangles'],dtype=int)
        self.nodes=np.asarray(facets['nodes_m']);self.world_nodes=self.nodes@self.rotation.T
        self.world_vertices=self.world_nodes[self.triangles]
        self.normals=np.asarray(facets['normals'])@self.rotation.T
        self.face_edge_lengths_m=np.linalg.norm(self.world_vertices-np.roll(self.world_vertices,1,axis=1),axis=2)
        points=vtk.vtkPoints();points.SetData(numpy_to_vtk(self.world_nodes,deep=True))
        cells=vtk.vtkCellArray();cells.SetCells(len(self.triangles),numpy_to_vtkIdTypeArray(
            np.column_stack((np.full(len(self.triangles),3),self.triangles)).astype(np.int64).ravel(),deep=True))
        poly=vtk.vtkPolyData();poly.SetPoints(points);poly.SetPolys(cells)
        self.poly=poly;self.locator=vtk.vtkStaticCellLocator();self.locator.SetDataSet(poly);self.locator.BuildLocator()

    def locate(self, point, outward_normal, *, tolerance_m):
        from bevel_gear.contact_mesh import _closest_on_triangles
        np=self.np;vtk=self.vtk;closest=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);squared=vtk.mutable(0.)
        self.locator.FindClosestPoint(list(point),closest,cell,sub,squared);index=int(cell)
        if not 0<=index<len(self.triangles):raise ArithmeticError('实际FE面点定位失败')
        mapped,bary,distance2=_closest_on_triangles(np.asarray(point),self.world_vertices[index:index+1],np)
        error=math.sqrt(float(distance2[0]))
        if error>tolerance_m:raise ArithmeticError(f'CAD→FE形函数映射偏差{error:g}m超过{tolerance_m:g}m')
        if float(np.asarray(outward_normal)@self.normals[index])<math.cos(.3):
            raise ArithmeticError('CAD点定位落在不对应的FE法向面')
        return dict(nodes=self.triangles[index],weights=bary[0],mapped_world=mapped[0],distance_m=error,
            face_edge_lengths_m=self.face_edge_lengths_m[index])

    def locate_many(self,points,outward_normals,*,tolerance_m):
        """保留逐点VTK最近面选择，随后批量执行同一有限三角投影公式。

        返回同顺序的映射及失败原因列表。失败点不补零、不吞掉；调用方
        可按原点索引报错。该批量接口不采用近邻插值或ProbeFilter。
        """
        from .contact_mesh import _closest_on_triangles
        np=self.np;vtk=self.vtk;queries=np.asarray(points,dtype=float);normals=np.asarray(outward_normals,dtype=float)
        if queries.ndim!=2 or queries.shape[1]!=3 or normals.shape!=queries.shape or not np.isfinite(queries).all() or not np.isfinite(normals).all():
            raise ValueError('批量实际FE定位须给同形有限N×3点与外法向')
        count=len(queries);result=[None]*count;reasons=[None]*count;ids=np.full(count,-1,dtype=int)
        closest=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);squared=vtk.mutable(0.)
        for i,query in enumerate(queries):
            self.locator.FindClosestPoint(query.tolist(),closest,cell,sub,squared);ids[i]=int(cell)
        valid=(ids>=0)&(ids<len(self.triangles));indices=np.flatnonzero(valid)
        for i in np.flatnonzero(~valid):reasons[int(i)]='实际FE面点定位失败'
        if not len(indices):return result,reasons
        selected=ids[indices]
        mapped,bary,distance2=_closest_on_triangles(queries[indices],self.world_vertices[selected],np)
        errors=np.sqrt(distance2);cosines=np.einsum('ij,ij->i',normals[indices],self.normals[selected])
        for row,i in enumerate(indices):
            index=int(i);facet=int(selected[row]);error=float(errors[row])
            if error>tolerance_m:
                reasons[index]=f'CAD→FE形函数映射偏差{error:g}m超过{tolerance_m:g}m';continue
            if cosines[row]<math.cos(.3):
                reasons[index]='CAD点定位落在不对应的FE法向面';continue
            result[index]=dict(nodes=self.triangles[facet],weights=bary[row],mapped_world=mapped[row],distance_m=error,
                face_edge_lengths_m=self.face_edge_lengths_m[facet])
        return result,reasons


def build_distributed_contact_system(prepared_fe, prepared_surface, geometry, *,
        mapping_tolerance_m=2e-4, maximum_cells=30000,ray_backend='scalar',quadrature_order=4,progress=None):
    """对实际有限曲面积分，以投影高度的真实导数建立分布力列。

    设参考平面法向n、两面外法向N。固定投影坐标的间隙导数为
    dg=(N2·du2)/(N2·n)−(N1·du1)/(N1·n)，故两侧单位平面合力
    的节点牵引分别为−N1/(N1·n)、+N2/(N2·n)。曲面面积Jac单列。
    """
    import numpy as np
    from time import perf_counter
    started=perf_counter();timing=dict(facet_and_locator_seconds=0.,batch_ray_seconds=0.,scalar_ray_seconds=0.,native_FE_mapping_seconds=0.)
    from scipy import sparse
    from bevel_gear.generated_tca import generated_tca_frames
    from bevel_gear.generated_ltca import select_generated_flank_facets
    from bevel_gear.contact_mesh import boundary_facets
    from bevel_gear.generated_halfspace import _sample_pair,_initial_uv
    mapping_tolerance_m=_positive_real(mapping_tolerance_m,'CAD至FE映射容差m')
    if isinstance(maximum_cells,bool) or not isinstance(maximum_cells,int) or maximum_cells<1:raise ValueError('压力小片上限须为正整数')
    if ray_backend not in ('scalar','batch'):raise ValueError('真实曲面射线后端须为scalar或batch')
    if isinstance(quadrature_order,bool) or not isinstance(quadrature_order,int) or not 2<=quadrature_order<=12:
        raise ValueError('每方向Gauss积分阶数须为2至12的整数')
    angle=geometry['driving_angle_rad'];phi=geometry['reference_follower_angle_rad']
    frames=generated_tca_frames(prepared_fe.tca,angle,phi)
    bounds={m:None for m in ('pinion','wheel')}
    rotations={m:np.asarray(frames[m+'_rotation_rows']) for m in bounds}
    locators={};records=[];matrix_rows={m:[] for m in bounds};matrix_cols={m:[] for m in bounds};matrix_values={m:[] for m in bounds}
    max_mapping=0.;failures=[]
    if quadrature_order==2:
        q=1/math.sqrt(3);gauss=[(-q,-q,.25),(q,-q,.25),(q,q,.25),(-q,q,.25)]
    else:
        locations,weights=np.polynomial.legendre.leggauss(quadrature_order)
        gauss=[(float(a),float(b),float(wa*wb/4)) for b,wb in zip(locations,weights) for a,wa in zip(locations,weights)]
    for patch in geometry['patches']:
        k=patch['tooth_offset'];normal=np.asarray(patch['plane_outward_normal_world'])
        origin=np.asarray(patch['plane_origin_world_m']);basis=np.asarray(patch['plane_tangent_rows'])
        for m,side in (('pinion','left_flank'),('wheel','right_flank')):
            phase_tick=perf_counter()
            facets=_working_facets(prepared_fe,m,k)
            locators[m,k]=_SurfaceTraceLocator(facets,rotations[m])
            timing['facet_and_locator_seconds']+=perf_counter()-phase_tick
        grid=patch['grid'];dx=grid.spacing_x_m;dy=grid.spacing_y_m;ox,oy=grid.origin_xy_m
        valid=np.asarray(patch['valid_mask'],dtype=bool)
        seed_maps={m:_initial_uv(prepared_surface.trimmed[m].surface,np.asarray(patch['frames'][m+'_rotation_rows']),origin,basis,patch['seed_uv'][m]) for m in bounds}
        cell_indices=list(zip(*np.nonzero(valid)));batch_samples=batch_reasons=None;sample_index=0;batch_mappings={}
        if ray_backend=='batch' and cell_indices:
            from .nurbs_batch import sample_pair_batch
            queries=np.asarray([(ox+(i+a*.5)*dx,oy+(j+b*.5)*dy) for j,i in cell_indices for a,b,weight in gauss])
            phase_tick=perf_counter()
            batch_samples,batch_reasons=sample_pair_batch(prepared_surface,patch['frames'],origin,basis,normal,
                patch['seed_uv'],queries,2e-11,.3,.3,seed_maps)
            timing['batch_ray_seconds']+=perf_counter()-phase_tick
            if any(sample is None for sample in batch_samples):
                first=next(i for i,sample in enumerate(batch_samples) if sample is None)
                raise ArithmeticError('已核实矩形内部Gauss点查询失败：'+str(batch_reasons[first]))
            phase_tick=perf_counter()
            for m,prefix in (('pinion','first'),('wheel','second')):
                mapped,errors=locators[m,k].locate_many(np.asarray([s[prefix]['world'] for s in batch_samples]),
                    np.asarray([s[prefix]['outward'] for s in batch_samples]),tolerance_m=mapping_tolerance_m)
                if any(error is not None for error in errors):raise ArithmeticError(next(error for error in errors if error is not None))
                batch_mappings[m]=mapped
            timing['native_FE_mapping_seconds']+=perf_counter()-phase_tick
        for j,i in cell_indices:
            if len(records)>=maximum_cells:raise ValueError('分布接触小片数超过显式上限')
            samples=[];gap=0.;levers={m:0. for m in bounds};geometric_levers={m:0. for m in bounds};actual_area={m:0. for m in bounds}
            for a,b,weight in gauss:
                target=(ox+(i+a*.5)*dx,oy+(j+b*.5)*dy)
                if batch_samples is None:
                    phase_tick=perf_counter()
                    sample,failure=_sample_pair(prepared_surface,patch['frames'],origin,basis,normal,patch['seed_uv'],
                        target,2e-11,.3,.3,seed_maps)
                    timing['scalar_ray_seconds']+=perf_counter()-phase_tick
                else:sample,failure=batch_samples[sample_index],batch_reasons[sample_index]
                sample_index+=1
                if sample is None:raise ArithmeticError('已核实矩形内部Gauss点查询失败：'+str(failure))
                if batch_mappings:mappings={m:batch_mappings[m][sample_index-1] for m in bounds}
                else:
                    phase_tick=perf_counter()
                    mappings={m:locators[m,k].locate(sample['first' if m=='pinion' else 'second']['world'],
                        sample['first' if m=='pinion' else 'second']['outward'],tolerance_m=mapping_tolerance_m) for m in bounds}
                    timing['native_FE_mapping_seconds']+=perf_counter()-phase_tick
                coefficients={}
                for m,prefix,sign in (('pinion','first',-1.),('wheel','second',1.)):
                    n=np.asarray(sample[prefix]['outward']);cosine=float(n@normal)
                    vector_world=sign*n/cosine;coefficients[m]=rotations[m].T@vector_world
                    # 力矩以实际FE形函数映射点求取，与节点载荷矩阵严格共轭。
                    mapped=mappings[m];velocity=np.cross(frames[m+'_axis_world'],mapped['mapped_world'])
                    levers[m]-=weight*float(vector_world@velocity)
                    geometric_levers[m]+=weight*sample['lever_arms_m'][m]
                    actual_area[m]+=weight*dx*dy/abs(cosine)
                    max_mapping=max(max_mapping,mapped['distance_m'])
                gap+=weight*sample['gap_m']
                samples.append(dict(gap_m=sample['gap_m'],integration_weight=weight,target_xy_m=target,
                    seed_uv={m:sample[p]['uv'] for m,p in (('pinion','first'),('wheel','second'))},
                    geometric_levers_m=sample['lever_arms_m'],mappings=mappings,coefficients_local=coefficients,
                    levers_m={m:-float((rotations[m]@coefficients[m])@np.cross(frames[m+'_axis_world'],mappings[m]['mapped_world'])) for m in bounds}))
            column=len(records)
            for sample in samples:
                for m in bounds:
                    mapping=sample['mappings'][m];vector=sample['coefficients_local'][m]
                    for node,weight in zip(mapping['nodes'],mapping['weights']):
                        for d in range(3):
                            matrix_rows[m].append(3*int(node)+d);matrix_cols[m].append(column)
                            matrix_values[m].append(sample['integration_weight']*float(weight)*float(vector[d]))
            # pinion的输入反力矩与wheel相同，均取负的接触节点力矩。
            records.append(dict(tooth_offset=k,row=int(j),column=int(i),gap_m=gap,area_plane_m2=dx*dy,
                actual_area_m2=actual_area,levers_m=levers,geometric_levers_m=geometric_levers,quadrature=samples,
                maximum_mapping_error_m=max(q['mappings'][m]['distance_m'] for q in samples for m in bounds),
                first_corner_points_world_m=[patch['first_corner_points_world_m'][j,i].tolist(),patch['first_corner_points_world_m'][j,i+1].tolist(),patch['first_corner_points_world_m'][j+1,i+1].tolist(),patch['first_corner_points_world_m'][j+1,i].tolist()]))
        if progress:progress(dict(stage='distributed_channels',tooth_offset=k,assembled_cells=len(records)))
    if not records:raise ArithmeticError('没有实际有限面压力小片')
    matrices={m:sparse.coo_matrix((matrix_values[m],(matrix_rows[m],matrix_cols[m])),shape=(3*len(prepared_fe.systems[m].mesh.nodes_m),len(records))).tocsc() for m in bounds}
    timing['total_assembly_seconds']=perf_counter()-started
    return dict(records=records,matrices=matrices,frames=frames,geometry=geometry,
        ray_backend=ray_backend,quadrature_order=quadrature_order,
        assembly_timing=timing,
        maximum_cad_fe_mapping_error_m=max_mapping,mapping_tolerance_m=mapping_tolerance_m,
        pressure_basis=f'真实投影小片P0法向压力，{quadrature_order}×{quadrature_order} Gauss分布牵引；非节点点载压力')


def solve_distributed_contact_system(prepared, system, pinion_torque_nm, *, maximum_working_set=2048,
        maximum_active_mapping_error_m=2e-6,gap_tolerance_m=1e-10,relative_tolerance=1e-9,
        allow_loaded_rigid_overlap=False,maximum_initial_penetration_m=1e-10,
        initial_indices=None,include_body_fields=False,progress=None):
    """解分布牵引的FE互补问题，原QP未知为每小片平面合力N。"""
    import numpy as np
    from time import perf_counter
    from bevel_gear.fe_contact_solver import solve_fe_contact_equilibrium
    if isinstance(pinion_torque_nm,bool) or not isinstance(pinion_torque_nm,(int,float)) or not math.isfinite(pinion_torque_nm) or pinion_torque_nm==0:
        raise ValueError('分布受载接触入口要求有限非零主动扭矩N·m')
    maximum_active_mapping_error_m=_positive_real(maximum_active_mapping_error_m,'承载CAD至FE映射容差m')
    maximum_initial_penetration_m=_positive_real(maximum_initial_penetration_m,'无载负间隙数值门限m')
    records=system['records'];matrices=system['matrices'];count=len(records)
    gaps=np.asarray([r['gap_m'] for r in records]);levers=np.asarray([r['levers_m']['wheel'] for r in records])
    if not allow_loaded_rigid_overlap and float(gaps.min()) < -maximum_initial_penetration_m:
        raise ArithmeticError('真实无载有限面存在超限负间隙，不能以分布接触掩盖几何侵入')
    pinion_levers=np.asarray([r['levers_m']['pinion'] for r in records]);free={m:w[prepared.systems[m].free_dofs] for m,w in matrices.items()};cache={}
    timing=dict(compliance_factor_solve_seconds=0.,compliance_projection_seconds=0.,FE_RHS_columns_solved=0)
    def columns(indices):
        indices=list(map(int,indices));missing=[i for i in indices if i not in cache]
        if missing:
            values=np.zeros((count,len(missing)))
            for m,w in free.items():
                for start in range(0,len(missing),24):
                    selected=missing[start:start+24];phase_tick=perf_counter()
                    u=prepared.systems[m].factorization.solve(w[:,selected].toarray())
                    timing['compliance_factor_solve_seconds']+=perf_counter()-phase_tick
                    timing['FE_RHS_columns_solved']+=len(selected);phase_tick=perf_counter()
                    values[:,start:start+len(selected)]+=w.T@u
                    timing['compliance_projection_seconds']+=perf_counter()-phase_tick
            for j,index in enumerate(missing):cache[index]=values[:,j].copy()
        return np.column_stack([cache[i] for i in indices])
    ratio=prepared.models['pinion'].request.tooth_count/prepared.models['wheel'].request.tooth_count
    phase_tick=perf_counter()
    equilibrium=solve_fe_contact_equilibrium(gaps,lever_arms_m=levers,torque_nm=pinion_torque_nm/ratio,
        compliance_columns=columns,maximum_working_set=maximum_working_set,
        allow_initial_penetration=bool(gaps.min()<0),gap_tolerance_m=gap_tolerance_m,relative_tolerance=relative_tolerance,
        initial_indices=initial_indices)
    timing['qp_total_seconds']=perf_counter()-phase_tick
    if not equilibrium['converged']:return dict(converged=False,status=equilibrium['status'],equilibrium=equilibrium)
    forces=np.asarray(equilibrium['forces_n']);active=forces>equilibrium['certificate']['force_acceptance_n'];areas=np.asarray([r['area_plane_m2'] for r in records]);pressures=forces/areas
    active_mapping=max((r['maximum_mapping_error_m'] for r,a in zip(records,active) if a),default=math.inf)
    if active_mapping>maximum_active_mapping_error_m:
        return dict(converged=False,status='active_CAD_FE_mapping_requires_refinement',equilibrium=equilibrium,
            maximum_active_cad_fe_mapping_error_m=active_mapping,maximum_allowed_m=maximum_active_mapping_error_m)
    phase_tick=perf_counter();responses={}
    for m,w in matrices.items():
        body=prepared.systems[m];load=np.asarray(w@forces).ravel();u=np.zeros(len(load))
        u[body.free_dofs]=body.factorization.solve(load[body.free_dofs])
        internal=np.asarray(body.stiffness@u).ravel();residual=internal-load
        relative_residual=float(np.linalg.norm(residual[body.free_dofs])/max(np.linalg.norm(load[body.free_dofs]),1.))
        if not np.isfinite(u).all() or relative_residual>1e-8:
            raise ArithmeticError('分布牵引FE完整响应不满足原刚度方程残差要求')
        responses[m]=dict(displacements_m=u.reshape((-1,3)),nodal_forces_n=load.reshape((-1,3)),
            strain_energy_j=.5*float(u@internal),external_work_j=float(u@load),
            maximum_free_residual_n=float(np.abs(residual[body.free_dofs]).max()),
            relative_free_residual=relative_residual)
    timing['body_response_seconds']=perf_counter()-phase_tick;phase_tick=perf_counter()
    displacements={m:r['displacements_m'] for m,r in responses.items()};point_gaps=[];cells=[]
    for index,record in enumerate(records):
        for q in record['quadrature']:
            elastic=0.
            for m in matrices:
                mapping=q['mappings'][m];u=displacements[m]
                elastic+=float(mapping['weights']@u[mapping['nodes']]@q['coefficients_local'][m])
            point_gaps.append(q['gap_m']+elastic-q['levers_m']['wheel']*equilibrium['theta_rad'])
        if active[index]:cells.append(dict(tooth_offset=record['tooth_offset'],row=record['row'],column=record['column'],
            pressure_pa=float(pressures[index]),force_plane_n=float(forces[index]),area_plane_m2=float(areas[index]),
            area_actual_first_m2=record['actual_area_m2']['pinion'],first_corner_points_world_m=record['first_corner_points_world_m']))
    energy=sum(r['strain_energy_j'] for r in responses.values())
    stiffness=float(levers[np.flatnonzero(active)]@np.linalg.solve(columns(np.flatnonzero(active))[active],levers[active]))
    timing['quadrature_diagnostic_seconds']=perf_counter()-phase_tick;timing['queried_C_columns']=len(cache)
    result=dict(converged=True,status='converged_discrete_patch_average',equilibrium=equilibrium,cells=cells,
        maximum_pressure_pa=float(pressures.max()),contact_area_plane_m2=float(areas[active].sum()),
        contact_area_actual_first_m2=sum(r['actual_area_m2']['pinion'] for r,a in zip(records,active) if a),
        total_force_plane_n=float(forces.sum()),active_cell_count=int(active.sum()),candidate_cell_count=count,
        pinion_torque_relative_residual=float((pinion_levers@forces-pinion_torque_nm)/abs(pinion_torque_nm)),
        elastic_te_rad=equilibrium['theta_rad'],incremental_stiffness_nm_per_rad=stiffness,
        fe_strain_energy_j=energy,condensed_energy_relative_error=energy/equilibrium['strain_energy_j']-1.,
        minimum_gauss_remaining_gap_m=float(min(point_gaps)),maximum_cad_fe_mapping_error_m=system['maximum_cad_fe_mapping_error_m'],
        maximum_active_cad_fe_mapping_error_m=active_mapping,
        pressure_filter_applied=False,halfspace_compliance_used=False,
        pressure_basis=system['pressure_basis'],limitations=['弱约束为小片积分平均，逐Gauss点非穿透须另检并收敛','单状态离散KKT不代表网格收敛、全域和外齿覆盖已验证'])
    result.update(contact_representation='distributed_P0_surface_cells',
        assembly_timing=system.get('assembly_timing'),numerics_timing=timing,
        quadrature_order=system.get('quadrature_order',2),
        active_cell_keys=[(r['tooth_offset'],r['row'],r['column']) for r,a in zip(records,active) if a],
        requested_tooth_offsets=[p['tooth_offset'] for p in system['geometry']['patches']],
        selected_tooth_offsets=[],fields=[],
        body_response_checks={m:{k:v for k,v in r.items() if k not in ('displacements_m','nodal_forces_n')} for m,r in responses.items()},
        numerical_force_acceptance_n=equilibrium['certificate']['force_acceptance_n'],
        geometry_linearization_checked=False,mesh_convergence_verified=False,external_tooth_coverage_verified=False)
    opening=sum(np.asarray(w.T@displacements[m].ravel()).ravel() for m,w in matrices.items())
    cursor=0;boundary_checks=[]
    for patch in system['geometry']['patches']:
        valid=np.asarray(patch['valid_mask'],dtype=bool);count=int(valid.sum())
        if not count:continue
        area=patch['grid'].spacing_x_m*patch['grid'].spacing_y_m;field=np.zeros(valid.shape);field[valid]=forces[cursor:cursor+count]
        elastic=np.zeros(valid.shape);elastic[valid]=opening[cursor:cursor+count]
        mask=np.zeros(valid.shape,dtype=bool);mask[valid]=active[cursor:cursor+count];cursor+=count
        result['selected_tooth_offsets'].append(patch['tooth_offset'])
        result['fields'].append(dict(tooth_offset=patch['tooth_offset'],pressure_pa=field/area,
            forces_n=field,valid_mask=valid.copy(),contact_mask=mask,cell_area_m2=area,elastic_deflection_m=elastic))
        padded=np.pad(valid,1,constant_values=False)
        interior=padded[1:-1,:-2]&padded[1:-1,2:]&padded[:-2,1:-1]&padded[2:,1:-1]
        boundary_checks.append(dict(tooth_offset=patch['tooth_offset'],touches_window_edge=bool(mask[0].any() or mask[-1].any() or mask[:,0].any() or mask[:,-1].any()),
            touches_invalid_candidate_neighbor=bool((mask&~interior).any())))
    if cursor!=len(forces):raise ArithmeticError('原P0力向量与实际有效网格索引不一致')
    result['contact_boundary_checks']=boundary_checks
    for m in matrices:
        edges=np.asarray([q['mappings'][m]['face_edge_lengths_m'] for r,a in zip(records,active) if a for q in r['quadrature']]).ravel()
        result['body_response_checks'][m].update(loaded_FE_face_edge_p95_m=float(np.quantile(edges,.95)),loaded_FE_face_edge_max_m=float(edges.max()))
    if include_body_fields:result['body_fields']=responses
    return result


def check_distributed_loaded_geometry(prepared,prepared_surface,system,result,*,
        maximum_error_m=1e-10,check_all_candidates=True):
    """实际受载姿态重查NURBS及FE位移插值，不把刚体误差当弹性开口。

    以本次小角步更新从动轮姿态，并在重新查询的材料点计算FE位移。
    接触弱约束仍为同一Gauss积分平均；逐Gauss最小值另报。完整有限面及
    外齿排除不是本函数的覆盖范围，须另行验证。
    """
    import numpy as np
    from time import perf_counter
    from .generated_tca import generated_tca_frames
    from .generated_halfspace import _sample_pair,_initial_uv
    maximum_error_m=_positive_real(maximum_error_m,'实际加载几何误差容差m')
    if not result.get('converged') or 'body_fields' not in result:
        return dict(passed=False,status='需要已解出的P0压力及同一次FE响应')
    started=perf_counter()
    geometry=system['geometry'];step=result['equilibrium']['theta_rad']
    alpha=geometry['driving_angle_rad'];phi=geometry['reference_follower_angle_rad']+step
    frames0=generated_tca_frames(prepared.tca,alpha,phi)
    rotations={m:np.asarray(frames0[m+'_rotation_rows']) for m in prepared.systems}
    offsets={r['tooth_offset'] for r in system['records']};locators={}
    for m in prepared.systems:
        for k in offsets:locators[m,k]=_SurfaceTraceLocator(_working_facets(prepared,m,k),rotations[m])
    contexts={}
    for patch in geometry['patches']:
        k=patch['tooth_offset'];frames=generated_tca_frames(prepared.tca,alpha,phi,tooth_offset=k)
        origin=np.asarray(patch['plane_origin_world_m']);basis=np.asarray(patch['plane_tangent_rows'])
        contexts[k]=dict(patch=patch,frames=frames,origin=origin,basis=basis,
            normal=np.asarray(patch['plane_outward_normal_world']),
            seed_maps={m:_initial_uv(prepared_surface.trimmed[m].surface,np.asarray(frames[m+'_rotation_rows']),
                origin,basis,patch['seed_uv'][m]) for m in prepared.systems})
    forces=np.asarray(result['equilibrium']['forces_n']);active=forces>result['numerical_force_acceptance_n']
    locator_seconds=perf_counter()-started;batched={};batched_mappings={};batch_calls=0;batch_points=0;scalar_points=0;scalar_seconds=0.;mapping_seconds=0.;ray_precompute_seconds=0.
    if system.get('ray_backend')=='batch':
        from .nurbs_batch import sample_pair_batch
        for k,context in contexts.items():
            items=[(index,qi,q['target_xy_m']) for index,r in enumerate(system['records']) if r['tooth_offset']==k
                and (check_all_candidates or active[index]) for qi,q in enumerate(r['quadrature'])]
            if not items:continue
            batch_calls+=1;batch_points+=len(items)
            phase_tick=perf_counter()
            values,reasons=sample_pair_batch(prepared_surface,context['frames'],context['origin'],context['basis'],context['normal'],
                context['patch']['seed_uv'],np.asarray([x[2] for x in items]),1e-13,.3,.3,context['seed_maps'])
            ray_precompute_seconds+=perf_counter()-phase_tick
            for item,value,reason in zip(items,values,reasons):batched[item[0],item[1]]=(value,reason)
            valid_positions=[i for i,value in enumerate(values) if value is not None]
            if valid_positions:
                phase_tick=perf_counter()
                for m,prefix in (('pinion','first'),('wheel','second')):
                    mapped,errors=locators[m,k].locate_many(np.asarray([values[i][prefix]['world'] for i in valid_positions]),
                        np.asarray([values[i][prefix]['outward'] for i in valid_positions]),tolerance_m=system['mapping_tolerance_m'])
                    for position,mapping,error in zip(valid_positions,mapped,errors):
                        item=items[position];batched_mappings[item[0],item[1],m]=(mapping,error)
                mapping_seconds+=perf_counter()-phase_tick
    maximum_rigid=0.;actual_averages=[];actual_forces=[];actual_active=[];point_min=math.inf;failures=[];checked=0
    for index,record in enumerate(system['records']):
        if not check_all_candidates and not active[index]:continue
        k=record['tooth_offset'];context=contexts[k];patch=context['patch'];remaining=[];integration_weights=[]
        for qi,q in enumerate(record['quadrature']):
            checked+=1
            if batched:value,reason=batched[index,qi]
            else:
                phase_tick=perf_counter();scalar_points+=1
                value,reason=_sample_pair(prepared_surface,context['frames'],context['origin'],context['basis'],
                    context['normal'],q['seed_uv'],q['target_xy_m'],1e-13,.3,.3)
                scalar_seconds+=perf_counter()-phase_tick
            if value is None:
                failures.append(dict(tooth_offset=k,row=record['row'],column=record['column'],gauss=qi,active=bool(active[index]),reason=reason));continue
            predicted=q['gap_m']-q['geometric_levers_m']['wheel']*step
            maximum_rigid=max(maximum_rigid,abs(value['gap_m']-predicted));elastic=0.
            for m,prefix,sign in (('pinion','first',-1.),('wheel','second',1.)):
                sample=value[prefix]
                try:
                    if batched_mappings:
                        mapping,error=batched_mappings[index,qi,m]
                        if error is not None:raise ArithmeticError(error)
                    else:
                        phase_tick=perf_counter()
                        mapping=locators[m,k].locate(sample['world'],sample['outward'],tolerance_m=system['mapping_tolerance_m'])
                        mapping_seconds+=perf_counter()-phase_tick
                except ArithmeticError as error:
                    failures.append(dict(tooth_offset=k,row=record['row'],column=record['column'],gauss=qi,active=bool(active[index]),reason=str(error)));break
                n=np.asarray(sample['outward']);vector=rotations[m].T@(sign*n/float(n@context['normal']))
                u=result['body_fields'][m]['displacements_m']
                elastic+=float(mapping['weights']@u[mapping['nodes']]@vector)
            else:
                remaining.append(value['gap_m']+elastic)
                integration_weights.append(q.get('integration_weight',.25))
        if len(remaining)==len(record['quadrature']):
            point_min=min(point_min,min(remaining));actual_averages.append(float(np.asarray(integration_weights)@remaining))
            actual_forces.append(float(forces[index]));actual_active.append(bool(active[index]))
    averages=np.asarray(actual_averages);mask=np.asarray(actual_active,dtype=bool);f=np.asarray(actual_forces)
    stationarity=float(np.abs(averages[mask]).max()) if mask.any() else math.inf
    minimum=float(averages.min()) if len(averages) else -math.inf
    acceptance=float(result['equilibrium']['certificate']['gap_acceptance_m'])+maximum_error_m
    return dict(passed=bool(not failures and maximum_rigid<=maximum_error_m and stationarity<=acceptance and minimum>=-acceptance),
        ray_backend_used='batch' if batched else 'scalar',batch_ray_calls=batch_calls,batch_ray_points=batch_points,scalar_ray_points=scalar_points,
        performance=dict(total_seconds=perf_counter()-started,locator_and_context_seconds=locator_seconds,
            batch_ray_seconds=ray_precompute_seconds,scalar_ray_seconds=scalar_seconds,native_FE_mapping_seconds=mapping_seconds),
        maximum_rigid_linearization_error_m=maximum_rigid,maximum_accepted_geometry_error_m=maximum_error_m,
        actual_loaded_active_average_stationarity_error_m=stationarity,
        actual_loaded_minimum_average_gap_m=minimum,actual_loaded_minimum_gauss_gap_m=point_min if math.isfinite(point_min) else None,
        actual_loaded_complementarity_residual_nm=float(np.sum(np.abs(f*averages))),
        average_gap_acceptance_m=acceptance,checked_gauss_points=checked,check_all_candidates=check_all_candidates,
        failed_loaded_queries=failures,
        definition='真实受载NURBS重投影＋同物理齿FE材料点位移；同权Gauss弱间隙与点值分别报告')


def solve_generated_distributed_fe_ltca(prepared,driving_angle_rad,pinion_torque_nm,*,
        prepared_surface=None,initial_geometry=None,tooth_offsets=(-2,-1,0,1,2),
        grid_shape=(161,41),spacing_m=(3e-4,1.5e-4),follower_correction_rad=1e-10,
        maximum_outer_iterations=5,maximum_geometry_error_m=1e-10,
        maximum_gauss_penetration_m=None,require_interior_contact=True,
        mapping_tolerance_m=2e-4,maximum_active_mapping_error_m=2e-6,
        maximum_working_set=2048,ray_backend='scalar',quadrature_order=4,
        initial_total_te_rad=None,initial_active_cell_keys=None,include_body_fields=False,progress=None,state_observer=None,
        maximum_candidate_cells=30000):
    """实际展成整轮的分布牵引FE，外迭代重建几何而复用同一K分解。

    输出原始P0二维压力，不把中心P1显示场当作压力未知量。未提供
    独立网格与外齿证明时对应验证状态始终为False。
    """
    from .generated_halfspace import prepare_generated_halfspace,build_generated_halfspace_patches
    if isinstance(maximum_outer_iterations,bool) or not isinstance(maximum_outer_iterations,int) or maximum_outer_iterations<1:
        raise ValueError('几何外迭代上限须为正整数')
    if isinstance(maximum_candidate_cells,bool) or not isinstance(maximum_candidate_cells,int) or maximum_candidate_cells<1:
        raise ValueError('实际压力小片数上限须为正整数')
    maximum_geometry_error_m=_positive_real(maximum_geometry_error_m,'几何误差容差m')
    if maximum_gauss_penetration_m is not None and (isinstance(maximum_gauss_penetration_m,bool) or not isinstance(maximum_gauss_penetration_m,(int,float)) or not math.isfinite(maximum_gauss_penetration_m) or maximum_gauss_penetration_m<0):
        raise ValueError('Gauss点允许穿入须为有限非负m或None')
    if initial_total_te_rad is not None and (isinstance(initial_total_te_rad,bool) or not isinstance(initial_total_te_rad,(int,float)) or not math.isfinite(initial_total_te_rad)):
        raise ValueError('初始总TE须为有限弧度数；仅作为迭代初值')
    warm_keys=None
    if initial_active_cell_keys is not None:
        warm_keys=[]
        for key in initial_active_cell_keys:
            if not isinstance(key,(list,tuple)) or len(key)!=3 or any(isinstance(x,bool) or not isinstance(x,int) for x in key):
                raise ValueError('初始承载格索引须为不重复的整数(k,row,column)')
            warm_keys.append(tuple(key))
        if len(set(warm_keys))!=len(warm_keys):raise ValueError('初始承载格索引不能重复')
    if prepared_surface is None:
        prepared_surface=prepare_generated_halfspace(prepared.models['pinion'].request,prepared.models['wheel'].request,cad_models=prepared.models)
    if prepared_surface.tca.requests!=prepared.tca.requests:raise ValueError('FE和有限曲面几何制造请求不一致')
    geometry=(build_generated_halfspace_patches(prepared_surface,driving_angle_rad,pinion_torque_nm,
        tooth_offsets=tooth_offsets,grid_shape=grid_shape,spacing_m=spacing_m,
        follower_correction_rad=follower_correction_rad,ray_backend=ray_backend) if initial_geometry is None else initial_geometry)
    if geometry['driving_angle_rad']!=driving_angle_rad or geometry['pinion_torque_nm']!=pinion_torque_nm:
        raise ValueError('复用的实际几何与本次相位或输入扭矩不一致')
    base=geometry;total=0. if initial_total_te_rad is None else float(initial_total_te_rad);history=[];result=None;system=None
    rebuild_parameters=dict(base['build_parameters'])
    # 初始缓存可以来自scalar；之后的重建明确遵守调用者选择的后端。
    rebuild_parameters['ray_backend']=ray_backend
    minimum_unloaded=min(float(p['gaps_m'][p['valid_mask']].min()) for p in base['patches'] if p['valid_mask'].any())
    if minimum_unloaded < -1e-10:raise ArithmeticError('初始无载几何存在超限负间隙；不能用受载初值掩盖')
    if total:
        geometry=build_generated_halfspace_patches(prepared_surface,driving_angle_rad,pinion_torque_nm,
            **rebuild_parameters,tca_results=base['tca_results'],
            follower_correction_rad=base['surface_reference_phase_correction_rad']+total)
    for outer in range(maximum_outer_iterations):
        if progress:progress(dict(stage='distributed_FE_geometry_iteration',iteration=outer+1,reference_angle_offset_rad=total))
        system=build_distributed_contact_system(prepared,prepared_surface,geometry,mapping_tolerance_m=mapping_tolerance_m,
            ray_backend=ray_backend,quadrature_order=quadrature_order,progress=progress,maximum_cells=maximum_candidate_cells)
        keys={ (r['tooth_offset'],r['row'],r['column']):i for i,r in enumerate(system['records']) }
        warm_indices=([keys[key] for key in warm_keys if key in keys] if warm_keys is not None else None)
        result=solve_distributed_contact_system(prepared,system,pinion_torque_nm,
            maximum_working_set=maximum_working_set,maximum_active_mapping_error_m=maximum_active_mapping_error_m,
            allow_loaded_rigid_overlap=bool(total or outer>0),initial_indices=warm_indices,include_body_fields=True,progress=progress)
        if not result.get('converged'):break
        step=result['equilibrium']['theta_rad']
        check=check_distributed_loaded_geometry(prepared,prepared_surface,system,result,maximum_error_m=maximum_geometry_error_m)
        history.append(dict(iteration=outer+1,reference_angle_offset_rad=total,step_theta_rad=step,
            loaded_rigid_overlap_explicitly_allowed=bool(total or outer>0),warm_working_set_size=len(warm_indices) if warm_indices is not None else 0,geometry_check=check))
        warm_keys=result['active_cell_keys']
        total+=step
        if check['passed']:break
        if outer+1==maximum_outer_iterations:
            result=dict(converged=False,status='loaded_geometry_reprojection_not_converged',candidate=result,geometry_check=check)
            break
        geometry=build_generated_halfspace_patches(prepared_surface,driving_angle_rad,pinion_torque_nm,
            **rebuild_parameters,tca_results=base['tca_results'],
            follower_correction_rad=base['surface_reference_phase_correction_rad']+total)
    if result is None:raise ArithmeticError('没有执行实际分布FE接触迭代')
    if result.get('converged'):
        result.update(linearization_step_theta_rad=result['equilibrium']['theta_rad'],elastic_te_rad=total,theta_rad=total,
            unloaded_reference_follower_angle_rad=base['reference_follower_angle_rad'],
            loaded_follower_angle_rad=geometry['reference_follower_angle_rad']+result['equilibrium']['theta_rad'],
            geometry_linearization_checked=True,geometry_check=history[-1]['geometry_check'])
        boundary_ok=not any(x['touches_window_edge'] or x['touches_invalid_candidate_neighbor'] for x in result['contact_boundary_checks'])
        if require_interior_contact and not boundary_ok:
            result=dict(converged=False,status='loaded_pressure_reaches_candidate_boundary',candidate=result)
        elif maximum_gauss_penetration_m is not None and result['minimum_gauss_remaining_gap_m'] < -maximum_gauss_penetration_m:
            result=dict(converged=False,status='Gauss_nonpenetration_requires_refinement',candidate=result)
    result.update(geometry=geometry,geometry_outer_history=history,
        initial_total_te_rad=initial_total_te_rad,minimum_unloaded_reference_gap_m=minimum_unloaded,
        model='实际有限工作面P0分布牵引＋整轮三维FE结构柔度；无半空间柔度',
        mesh_convergence_verified=False,external_tooth_coverage_verified=False)
    if result.get('converged') and state_observer is not None:state_observer(prepared,system,result)
    if not include_body_fields:
        result.pop('body_fields',None)
        if isinstance(result.get('candidate'),dict):result['candidate'].pop('body_fields',None)
    return result
