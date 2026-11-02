"""实际刀盘展成整轮的三维线弹性FE接触；SI接口，CAD内部mm。

工作面来自有限UV的实际NURBS，柔度来自整轮四面体，不使用理想直齿
曲面或经验刚度。槽底根圆连接、孔壁固定与低阶面片压力均有适用限制。
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from pathlib import Path

from .elasticity import (ElasticMaterial, ElasticSystem, TetraMesh, ContactLoadChannel,
    mesh_step_tetra, condense_contact_flexibility, contact_load_matrix, _numeric)
from .generated_tca import (PreparedGeneratedTca, prepare_generated_tca, solve_generated_tca,
    generated_tca_frames)


@dataclass
class PreparedGeneratedLtca:
    """真实整轮CAD、工作面TCA、两轮FE系统；同轮所有齿共享同一弹性体。"""
    models: dict
    tca: PreparedGeneratedTca
    systems: dict
    diagnostics: dict
    _contact_cache: dict = field(default_factory=dict, repr=False)
    _facet_cache: dict = field(default_factory=dict, repr=False)


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int,float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(name + "须为正有限实数")
    return float(value)


def _rz(angle, np):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])


def _check_mesh_metadata(mesh):
    """拒绝与SI节点/单实体相矛盾的复用元数据；无CAD绑定仍明确记unknown。"""
    np, _, _ = _numeric()
    metadata=mesh.metadata
    if not isinstance(metadata,dict):raise ValueError("体网格metadata须为字典")
    if metadata.get('units') not in (None,'m'):
        raise ValueError("体网格单位元数据与SI节点m矛盾")
    if 'volume_entities' in metadata and (isinstance(metadata['volume_entities'],bool) or metadata['volume_entities']!=1):
        raise ValueError("体网格须对应单一整轮实体")
    nodes=np.asarray(mesh.nodes_m,dtype=float)
    if nodes.ndim!=2 or nodes.shape[1]!=3 or not np.isfinite(nodes).all() or not len(nodes):
        raise ValueError("体网格节点须为有限N×3 SI数组")
    actual=np.array([nodes.min(axis=0),nodes.max(axis=0)])
    if 'bounds_m' in metadata:
        stated=np.asarray(metadata['bounds_m'],dtype=float)
        if stated.shape!=(2,3) or not np.isfinite(stated).all() or not np.allclose(stated,actual,rtol=1e-10,atol=1e-12):
            raise ValueError("体网格bounds_m与实际节点边界不一致")
    return dict(units_declared=metadata.get('units'),volume_entities_declared=metadata.get('volume_entities'),
        actual_node_bounds_m=actual.tolist(),cad_binding_state='unknown',
        source_step_sha256=metadata.get('source_step_sha256'),
        binding_limit='元数据单位/边界/体积一致不证明全轮拓扑相同；外部网格须另保留源CAD绑定')


def select_generated_flank_facets(mesh, surface, *, tooth_count, tooth_index=0,
        matching_tolerance_m=2e-7, maximum_surface_normal_angle_rad=.5, boundary=None):
    """三顶点逐个投影至有限UV工作面，旋回周期齿距，保留原FE朝外法向。

    surface为实际展成shell中的Geom_BSplineSurface，单位mm。包围盒取
    正权NURBS控制极点的凸包，避免用稀疏采样漏筛。三顶点全匹配并检查
    面片法向与局部曲面切平面相符；不把端面/根面标签直接当工作面。
    """
    from .contact_mesh import boundary_facets
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.gp import gp_Pnt, gp_Vec
    np, _, _ = _numeric()
    tolerance = _positive(matching_tolerance_m, "工作面匹配容差m")
    angle_limit = _positive(maximum_surface_normal_angle_rad, "面片切平面角容差rad")
    if angle_limit >= math.pi/2:
        raise ValueError("面片切平面角容差须小于pi/2")
    if isinstance(tooth_count,bool) or not isinstance(tooth_count,int) or tooth_count<3:
        raise ValueError("齿数须为不小于3的整数")
    if isinstance(tooth_index,bool) or not isinstance(tooth_index,int) or not 0<=tooth_index<tooth_count:
        raise ValueError("实际齿索引须在0至齿数减1之间")
    if not isinstance(mesh,TetraMesh):
        raise TypeError("mesh须为TetraMesh")
    facets = boundary_facets(mesh) if boundary is None else boundary
    nodes = facets["nodes_m"]
    rotation = _rz(tooth_index*2*math.pi/tooth_count,np)
    reference_nodes = nodes @ rotation
    poles = []
    for u in range(1,surface.NbUPoles()+1):
        for v in range(1,surface.NbVPoles()+1):
            if surface.Weight(u,v)<=0:
                raise ValueError("NURBS控制极点存在非正权，不能使用凸包预筛")
            poles.append(surface.Pole(u,v).Coord())
    poles = np.asarray(poles)*.001
    lo, hi = poles.min(axis=0)-tolerance, poles.max(axis=0)+tolerance
    boundary_nodes = np.unique(facets["triangles"])
    candidates = boundary_nodes[np.all((reference_nodes[boundary_nodes]>=lo)&(reference_nodes[boundary_nodes]<=hi),axis=1)]
    valid = np.zeros(len(nodes),dtype=bool)
    distances = np.full(len(nodes),math.inf)
    normals = np.zeros((len(nodes),3))
    parameters = np.full((len(nodes),2),math.nan)
    bounds = surface.Bounds()
    projector = GeomAPI_ProjectPointOnSurf()
    projector.Init(surface,*bounds,min(tolerance*1000*.1,1e-7))
    failures = 0
    for index in candidates:
        point = reference_nodes[index]*1000
        projector.Perform(gp_Pnt(*point))
        if not projector.NbPoints():
            failures += 1
            continue
        distance = projector.LowerDistance()*.001
        distances[index] = distance
        if distance>tolerance:
            continue
        u,v = projector.LowerDistanceParameters()
        if not (bounds[0]-1e-12<=u<=bounds[1]+1e-12 and bounds[2]-1e-12<=v<=bounds[3]+1e-12):
            continue
        p,du,dv = gp_Pnt(),gp_Vec(),gp_Vec()
        surface.D1(u,v,p,du,dv)
        normal = np.cross(du.Coord(),dv.Coord())
        magnitude = float(np.linalg.norm(normal))
        if magnitude<=1e-25:
            continue
        valid[index] = True
        normals[index] = rotation@(normal/magnitude)
        parameters[index] = (u,v)
    triangles = facets["triangles"]
    mask = np.all(valid[triangles],axis=1)
    preliminary = np.flatnonzero(mask)
    for index in preliminary:
        analytic = normals[triangles[index]].sum(axis=0)
        magnitude = float(np.linalg.norm(analytic))
        if magnitude<=1e-12 or abs(float(analytic@facets['normals'][index]/magnitude))<math.cos(angle_limit):
            mask[index] = False
    result = {key:(value if key in ("nodes_m","metadata") else value[mask]) for key,value in facets.items()}
    result["tooth_indices"] = np.full(int(mask.sum()),tooth_index,dtype=np.int64)
    result["maximum_vertex_matching_distance_m"] = np.max(distances[triangles[mask]],axis=1) if mask.any() else np.empty(0)
    result["vertex_uv"] = parameters[triangles[mask]]
    result["metadata"] = {**facets["metadata"],"tooth_index":tooth_index,"tooth_count":tooth_count,
        "selection":"实际NURBS有限UV逐FE三顶点投影；保留四面体朝外法向和全局节点索引",
        "matching_tolerance_m":tolerance,"surface_uv_bounds":bounds,"projection_candidates":len(candidates),
        "matched_vertices":int(valid.sum()),"projection_failures":failures,
        "selected_facets":int(mask.sum()),"maximum_surface_normal_angle_rad":angle_limit}
    return result


def prepare_generated_ltca(pinion_request,wheel_request,*,pinion_blank_settings,wheel_blank_settings,
        pinion_material,wheel_material,directory,element_size_m=.001,contact_refined_size_m=.0003,
        refinement_radius_m=.0012,refinement_driving_angle_rad=0.,maximum_nodes=120000,
        width_count=25,profile_count=25,cad_models=None,meshes=None,
        refinement_driving_angles_rad=None,refinement_tooth_offsets=(-1,0,1),refinement_polylines_m=None,
        factorization_factory=None,refinement_core_radius_m=None,refinement_regions=None):
    """实际整轮CAD→Gmsh→完全固定孔壁→两体线弹性；可显式复用同模型CAD/网格。

    meshes仅接收调用方明确准备的SI四面体网格；拒绝矛盾单位/边界/实体数，
    与CAD体积差超过2%也拒绝，但这不是网格收敛判据。无源CAD绑定的复用
    网格记unknown，不冒称同一CAD已验证。孔壁固定不含轴/轴承/配合柔度。
    factorization_factory可明确接入完成SPD分解并提供solve的专业数值引擎；
    默认SciPy分解，不改变Kff、材料、单元或边界，响应另检相对残差。
    refinement_regions可按pinion/wheel提供额外局部SI折线场；只在自动
    划网时与原全轨迹场取Min。复用网格不能再应用新的网格控制场。
    """
    from .generated_surface import _request
    from .generated_solid import build_generated_gear,GeneratedBlankSettings
    from .cad_import import ImportedCadModel,export_model_step
    np, _, _ = _numeric()
    regions={} if refinement_regions is None else refinement_regions
    if (not isinstance(regions,dict) or set(regions)-{'pinion','wheel'} or
            any(not isinstance(value,(list,tuple)) for value in regions.values())):
        raise ValueError('附加加密区须按pinion/wheel给出区域列表')
    if meshes is not None and any(regions.values()):raise ValueError('复用体网格不能再应用新的加密区；须先核验已有网格的区域身份')
    if refinement_core_radius_m is not None and contact_refined_size_m is None:
        raise ValueError('近场细网格核心半径须与局部单元尺寸同时提供')
    requests = dict(pinion=_request(pinion_request),wheel=_request(wheel_request))
    settings = dict(pinion=pinion_blank_settings,wheel=wheel_blank_settings)
    materials = dict(pinion=pinion_material,wheel=wheel_material)
    directory = Path(directory);directory.mkdir(parents=True,exist_ok=True)
    for mapping,name in ((cad_models,"CAD模型"),(meshes,"体网格")):
        if mapping is not None and (not isinstance(mapping,dict) or set(mapping)!=set(requests)):
            raise ValueError(name+"须为pinion/wheel两项字典")
    models = {}
    for member in requests:
        s=settings[member]
        if isinstance(s,dict):s=GeneratedBlankSettings(**s)
        if not isinstance(s,GeneratedBlankSettings):raise TypeError("齿坯设置类型无效")
        settings[member]=s
        model=(build_generated_gear(requests[member],s,width_count=width_count,profile_count=profile_count)
               if cad_models is None else cad_models[member])
        if model.request!=requests[member] or model.settings!=s:
            raise ValueError("复用CAD的制造输入或孔壁/齿坯设置与请求不符")
        if model.diagnostics.get('solids')!=1 or not model.diagnostics.get('valid') or model.diagnostics.get('object_scope')!='full_gear':
            raise ValueError("实际展成模型不是单一有效整轮实体")
        models[member]=model
    tca=prepare_generated_tca(*requests.values(),cad_models=models)
    reference=solve_generated_tca(tca,refinement_driving_angle_rad)
    if not reference['converged']:
        raise ArithmeticError("局部加密参考姿态没有有效实际工作面相切解")
    tracks={}
    if refinement_polylines_m is not None:
        if refinement_driving_angles_rad is not None:raise ValueError('显式加密带与相位轨迹不可同时指定')
        if not isinstance(refinement_polylines_m,dict) or set(refinement_polylines_m)!={'pinion','wheel'}:
            raise ValueError('显式加密带须按pinion/wheel给出真实SI折线')
        tracks=refinement_polylines_m
    if refinement_driving_angles_rad is not None:
        angles=list(refinement_driving_angles_rad)
        offsets=list(refinement_tooth_offsets)
        if len(angles)<2 or any(isinstance(a,bool) or not isinstance(a,(int,float)) or not math.isfinite(a) for a in angles):
            raise ValueError("齿距加密轨迹须提供至少两个有限主动角rad")
        if any(b<=a for a,b in zip(angles,angles[1:])):
            raise ValueError("齿距加密主动角须严格递增，不可重复或乱序")
        if not offsets or len(set(offsets))!=len(offsets) or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets):
            raise ValueError("加密齿号须为不重复整数序列")
        samples={k:[solve_generated_tca(tca,a,tooth_offset=k) for a in angles] for k in offsets}
        if not all(any(samples[k][i]['converged'] for k in offsets) for i in range(len(angles))):
            raise ValueError("加密轨迹存在所有候选齿对均无有效相切的姿态；未连接失败候选点")
        # 齿对进入/退出实体裁切边界时，k=0不一定仍有效。只将每一齿对的
        # 连续真实有效段连成独立轨迹，不跨越失败姿态，也不连接齿对之间的空档。
        segments=[]
        for k in offsets:
            current=[]
            for sample in samples[k]:
                if sample['converged']:
                    current.append(sample)
                else:
                    if len(current)>=2:segments.append(current)
                    current=[]
            if len(current)>=2:segments.append(current)
        if not segments:
            raise ValueError("真实有效相切点不足以形成连续加密轨迹")
        for member in requests:
            z=requests[member].tooth_count
            tracks[member]=[]
            for segment in segments:
                local=np.array([np.array(sample['frames'][member+'_rotation_rows']).T @
                                np.array(sample[member+'_contact_point_world_m']) for sample in segment])
                for k in offsets:
                    rotation=_rz(((-k if member=='pinion' else k)%z)*2*math.pi/z,np)
                    tracks[member].append((local@rotation.T).tolist())
    systems={}
    for member in requests:
        model=models[member];material=materials[member]
        if isinstance(material,dict):material=ElasticMaterial(**material)
        path=directory/(member+'_elastic_body.step')
        if meshes is None:
            export_model_step(ImportedCadModel(model.shape,model.diagnostics),path)
            frame=np.array(reference['frames'][member+'_rotation_rows'])
            local=frame.T@np.array(reference[member+'_contact_point_world_m'])
            controls=(dict(refinement_polylines_m=tracks[member]) if tracks else dict(refinement_points_m=[local.tolist()]))
            options={} if contact_refined_size_m is None else dict(**controls,
                refined_size_m=contact_refined_size_m,refinement_radius_m=refinement_radius_m,
                refinement_core_radius_m=refinement_core_radius_m)
            if regions.get(member):options['refinement_regions']=regions[member]
            mesh=mesh_step_tetra(path,element_size_m=element_size_m,maximum_nodes=maximum_nodes,**options)
        else:
            mesh=meshes[member]
            if not isinstance(mesh,TetraMesh):raise TypeError("复用体网格须为TetraMesh")
        binding=_check_mesh_metadata(mesh)
        if meshes is None:
            import hashlib
            binding.update(cad_binding_state='generated_from_this_exported_cad',source_step_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        nodes=np.asarray(mesh.nodes_m)
        bore=model.settings.bore_diameter_m/2
        tolerance=max(model.settings.boolean_tolerance_m*2,bore*1e-8)
        fixed=np.flatnonzero(np.abs(np.linalg.norm(nodes[:,:2],axis=1)-bore)<=tolerance)
        system=ElasticSystem(mesh,material,fixed,maximum_dofs=3*maximum_nodes,
            factorization_factory=factorization_factory)
        volume_error=system.diagnostics['volume_m3']/model.diagnostics['volume_m3']-1
        if abs(volume_error)>.02:
            raise ValueError("体网格与实际CAD体积差超过2%；拒绝未匹配轮体或过粗体积分辨率")
        system.diagnostics.update(bore_radius_m=bore,bore_selection_tolerance_m=tolerance,
            cad_volume_m3=model.diagnostics['volume_m3'],
            relative_mesh_volume_error=volume_error,maximum_mesh_volume_error_relative=.02,mesh_cad_binding=binding)
        systems[member]=system
    diagnostics=dict(full_gears=True,element_size_m=element_size_m,contact_refined_size_m=contact_refined_size_m,
        refinement_radius_m=refinement_radius_m,refinement_core_radius_m=refinement_core_radius_m,reference_tca=reference,
        refinement_driving_angles_rad=refinement_driving_angles_rad,
        refinement_tooth_offsets=list(refinement_tooth_offsets) if tracks else [],
        refinement_control_geometry=('explicit_surface_patch_polylines' if refinement_polylines_m is not None else
            'independent_contact_track_polylines' if tracks else 'single_tangent_point'),
        cad_reused=cad_models is not None,mesh_reused=meshes is not None,
        members={m:s.diagnostics for m,s in systems.items()},
        support_model="两轮实际轴孔孔壁完全固定；不包含轴/轴承/箱体柔度",
        root_model="实际展成复合根过渡；槽底仍为显式根圆连接约定",
        neighbour_coupling="整轮三维弹性，保留同轮多齿之间全部非对角耦合")
    if any(regions.values()):
        diagnostics['additional_refinement_region_counts']={m:len(regions.get(m,())) for m in requests}
        diagnostics['additional_refinement_scope']='附加网格控制；不自动证明单相位或全周期空间收敛'
    return PreparedGeneratedLtca(models,tca,systems,diagnostics)


def solve_generated_fe_ltca(prepared,driving_angle_rad,pinion_torque_nm,*,tooth_offsets=(-1,0,1),
        maximum_distance_m=2e-4,maximum_normal_angle_rad=.3,maximum_channels=300,
        matching_tolerance_m=2e-7,mesh_phase_tolerance_m=1e-9,maximum_mesh_phase_iterations=10,
        maximum_loaded_angle_rad=.005,maximum_pinion_torque_relative_residual=.02,include_fields=True):
    """实际工作面有限元承载接触；T1为相对TCA参数增量轴的有符号扭矩。

    当前工作侧为主动left/从动right，须由实际法向和力臂判定正/负受载。
    KKT在从动轴满足T2=T1/(z1/z2)，主动轴重构扭矩另报残差；它不受第二
    独立力矩约束。未加载TCA相位、FE离散初相位、弹性角TE分别保留。
    """
    from .contact_mesh import boundary_facets,pair_facet_contacts
    from .contact_equilibrium import solve_contact_equilibrium
    np, _, _ = _numeric()
    if not isinstance(prepared,PreparedGeneratedLtca):
        raise TypeError("需要已准备的PreparedGeneratedLtca")
    for value,name in ((driving_angle_rad,"主动角rad"),(pinion_torque_nm,"小轮有符号扭矩Nm")):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError(name+"须为有限实数")
    if pinion_torque_nm==0:
        raise ValueError("本承载入口要求非零有符号扭矩")
    for value,name in ((maximum_distance_m,"候选距离m"),(matching_tolerance_m,"曲面匹配容差m"),
            (mesh_phase_tolerance_m,"网格相位残差m"),(maximum_loaded_angle_rad,"弹性小转角上限rad")):
        _positive(value,name)
    if isinstance(maximum_channels,bool) or not isinstance(maximum_channels,int) or not 1<=maximum_channels<=2048:
        raise ValueError("接触通道上限须为1至2048整数")
    _positive(maximum_pinion_torque_relative_residual, '主动轴重构扭矩相对残差上限')
    if isinstance(maximum_mesh_phase_iterations,bool) or not isinstance(maximum_mesh_phase_iterations,int) or maximum_mesh_phase_iterations<1:
        raise ValueError("网格初相位迭代上限须为正整数")
    offsets=tuple(tooth_offsets)
    if not offsets or 0 not in offsets or len(set(offsets))!=len(offsets) or any(isinstance(k,bool) or not isinstance(k,int) for k in offsets):
        raise ValueError("齿对偏移须为不重复整数序列并包含参考齿对0")
    direction=1 if pinion_torque_nm>0 else -1
    key=(driving_angle_rad,direction,offsets,maximum_distance_m,maximum_normal_angle_rad,maximum_channels,
         matching_tolerance_m,mesh_phase_tolerance_m,maximum_mesh_phase_iterations)
    cache_hit=key in prepared._contact_cache
    if not cache_hit:
        tca_results={k:solve_generated_tca(prepared.tca,driving_angle_rad,tooth_offset=k) for k in offsets}
        candidates=[value for value in tca_results.values() if value['converged']]
        if not candidates:
            raise ArithmeticError("所选齿对没有有效的实际NURBS未加载相切解；不能以域外候选代替实际接触")
        # 载荷方向决定最先闭合的齿对；不能在整个齿距内强制齿对0始终是参考。
        # 所有候选从动角均相对同一整轮轴，选出一个参考后全部齿共用同一相位。
        reference=min(candidates,key=lambda item:(direction*item['follower_increment_angle_rad'],abs(item['tooth_offset'])))
        follower=reference['follower_increment_angle_rad']
        facets={};selection={}
        for member,side in (('pinion','left_flank'),('wheel','right_flank')):
            model=prepared.models[member];mesh=prepared.systems[member].mesh
            boundary_key=('boundary',member)
            if boundary_key not in prepared._facet_cache:
                prepared._facet_cache[boundary_key]=boundary_facets(mesh)
            boundary=prepared._facet_cache[boundary_key]
            for k in offsets:
                # q将局部z反向；参数齿距+k对应小轮实体槽-k、大轮实体槽+k。
                tooth=(-k if member=='pinion' else k)%model.request.tooth_count
                selection_key=(member,tooth,matching_tolerance_m)
                if selection_key not in prepared._facet_cache:
                    prepared._facet_cache[selection_key]=select_generated_flank_facets(mesh,model.shell.surfaces[side],tooth_count=model.request.tooth_count,
                        tooth_index=tooth,matching_tolerance_m=matching_tolerance_m,boundary=boundary)
                f=prepared._facet_cache[selection_key]
                facets[member,k]=f;selection[member,k]=f['metadata']
        sign_diagnostics=[]
        def at_phase(correction):
            frames=generated_tca_frames(prepared.tca,driving_angle_rad,follower+correction)
            a,b=(frames[m+'_rotation_rows'] for m in ('pinion','wheel'))
            records=[];signs=[]
            for k in offsets:
                pairs=pair_facet_contacts(facets['pinion',k],facets['wheel',k],a,b,
                    maximum_distance_m=maximum_distance_m,maximum_normal_angle_rad=maximum_normal_angle_rad)
                for pair in pairs:
                    normal=np.array(pair['normal_world'])
                    lever1=float(normal@np.cross(frames['pinion_axis_world'],pair['first_point_world_m']))
                    lever2=float(-normal@np.cross(frames['wheel_axis_world'],pair['second_point_world_m']))
                    signs.append((lever1,lever2))
                    if direction*lever1<=1e-10 or direction*lever2<=1e-10:
                        continue
                    records.append(dict(pair,tooth_offset=k,pinion_lever_m=lever1,wheel_lever_m=lever2))
            sign_diagnostics.append(dict(follower_correction_rad=correction,all_opposing_candidates=len(signs),
                signed_load_candidates=len(records),pinion_lever_range_m=[min(v[0] for v in signs),max(v[0] for v in signs)] if signs else None,
                wheel_lever_range_m=[min(v[1] for v in signs),max(v[1] for v in signs)] if signs else None))
            if not records:
                raise ArithmeticError("所选实际工作侧没有可承受该有符号扭矩的候选；须检查轴/工作侧与扭矩方向，不能对力臂取绝对值")
            if len(records)>maximum_channels:
                raise ValueError(f"真实候选{len(records)}个超过显式通道上限{maximum_channels}；不静默删点")
            return records
        correction=0.;history=[]
        for _ in range(maximum_mesh_phase_iterations):
            records=at_phase(correction)
            gaps=np.array([r['signed_gap_m'] for r in records])
            lever=np.array([r['wheel_lever_m'] for r in records])
            minimum=float(gaps.min())
            history.append(dict(follower_correction_rad=correction,minimum_signed_gap_m=minimum,candidate_count=len(records)))
            if abs(minimum)<=mesh_phase_tolerance_m:
                break
            # dg/dtheta2=-ell2；方向乘积为正只用于选择首个闭合，不改原signed ell。
            correction+=direction*float(np.min(gaps/(direction*lever)))
            if abs(correction)>.02:
                raise ArithmeticError("FE网格初相位修正超出0.02rad，不能用大相位覆盖齿形失配")
        else:
            raise ArithmeticError("实际FE面片的未加载相位没有达到给定间隙残差")
        if float(gaps.min()) < -mesh_phase_tolerance_m:
            raise ArithmeticError("未加载FE工作面仍有超限负间隙")
        count=len(records);compliance=np.zeros((count,count));channels={};condensation={}
        for member,prefix in (('pinion','first'),('wheel','second')):
            rows=[ContactLoadChannel(tuple(r[prefix+'_node_indices']),tuple(r[prefix+'_weights']),
                                    tuple(r[prefix+'_force_normal_local'])) for r in records]
            channels[member]=rows
            # 同轮全部齿一起凝聚；单轮可半正定，最终两体C仍由KKT严格检验正定。
            condensed=condense_contact_flexibility(prepared.systems[member],rows,require_positive_definite=False)
            compliance+=np.array(condensed['compliance_m_per_n'])
            condensation[member]={k:v for k,v in condensed.items() if k!='compliance_m_per_n'}
        prepared._contact_cache[key]=dict(reference=reference,tca_results=tca_results,follower=follower,
            records=records,gaps=gaps,lever=lever,compliance=compliance,channels=channels,
            condensation=condensation,correction=correction,history=history,
            selection=[dict(member=m,tooth_offset=k,**v) for (m,k),v in selection.items()],
            sign_diagnostics=sign_diagnostics)
    c=prepared._contact_cache[key]
    ratio=prepared.tca.requests['pinion'].tooth_count/prepared.tca.requests['wheel'].tooth_count
    target=pinion_torque_nm/ratio
    equilibrium=solve_contact_equilibrium(c['gaps'].tolist(),c['compliance'].tolist(),
        lever_arms_m=c['lever'].tolist(),torque_nm=target,allow_initial_penetration=bool(c['gaps'].min()<0),
        gap_tolerance_m=max(1e-10,mesh_phase_tolerance_m*2))
    forces=np.array(equilibrium.forces_n)
    mesh_frames=generated_tca_frames(prepared.tca,driving_angle_rad,c['follower']+c['correction'])
    frame_arrays={m:np.array(mesh_frames[m+'_rotation_rows']) for m in ('pinion','wheel')}
    contacts=[]
    for index,record in enumerate(c['records']):
        contacts.append(dict(record,force_n=float(forces[index]),mean_pressure_pa=float(forces[index]/record['area_m2']),
            remaining_gap_m=equilibrium.remaining_gaps_m[index],active=index in equilibrium.active_indices,
            first_vertices_world_m=(prepared.systems['pinion'].nodes[list(record['first_node_indices'])]@frame_arrays['pinion'].T).tolist(),
            second_vertices_world_m=(prepared.systems['wheel'].nodes[list(record['second_node_indices'])]@frame_arrays['wheel'].T).tolist()))
    summaries={};fields={}
    for member in ('pinion','wheel'):
        system=prepared.systems[member]
        load=contact_load_matrix(system,c['channels'][member],matrix_format='csc')@forces
        response=system.solve(load)
        summaries[member]=dict(maximum_element_von_mises_pa=max(response['element_von_mises_pa']),
            maximum_displacement_m=float(np.linalg.norm(response['displacements_m'],axis=1).max()),
            maximum_free_residual_n=response['maximum_free_residual_n'],total_force_balance_n=response['total_force_balance_n'],
            strain_energy_j=response['strain_energy_j'])
        if include_fields:fields[member]=response
    active=equilibrium.active_indices
    stiffness=None
    if active:
        ca=c['compliance'][np.ix_(active,active)];la=c['lever'][active]
        stiffness=float(la@np.linalg.solve(ca,la))
    pinion_reconstructed=float(sum(r['pinion_lever_m']*f for r,f in zip(c['records'],forces)))
    theta=equilibrium.theta_rad
    model_valid=theta is not None and abs(theta)<=maximum_loaded_angle_rad
    pinion_relative_residual=(pinion_reconstructed-pinion_torque_nm)/abs(pinion_torque_nm)
    pinion_balance_valid=abs(pinion_relative_residual)<=maximum_pinion_torque_relative_residual
    accepted=equilibrium.converged and model_valid and pinion_balance_valid
    candidate=None if accepted else dict(theta_rad=theta,incremental_mesh_stiffness_nm_per_rad=stiffness,
        contacts=contacts,stress_summary=summaries,elastic_fields=fields,equilibrium=asdict(equilibrium),
        pinion_torque_reconstructed_nm=pinion_reconstructed,wheel_torque_reconstructed_nm=equilibrium.resultant_load,
        interpretation='未通过的候选，不作为已接受加载物理结果')
    loaded_theta=theta if accepted else None
    status=(equilibrium.status if not equilibrium.converged else 'small_rotation_limit' if not model_valid
            else 'pinion_torque_residual_limit' if not pinion_balance_valid else 'converged')
    return dict(converged=accepted,status=status,contact_representation='facets',
        model='实际展成整轮三维线弹性FE＋有限UV工作面三角面片无摩擦接触',
        driving_angle_rad=driving_angle_rad,pinion_torque_nm=pinion_torque_nm,torque_magnitude_nm=abs(pinion_torque_nm),
        torque_direction=direction,wheel_target_torque_nm=target,pinion_torque_reconstructed_nm=pinion_reconstructed if accepted else None,
        pinion_torque_balance_residual_nm=pinion_reconstructed-pinion_torque_nm if accepted else None,
        pinion_torque_balance_relative_residual=(pinion_reconstructed-pinion_torque_nm)/abs(pinion_torque_nm) if accepted else None,
        wheel_torque_reconstructed_nm=equilibrium.resultant_load if accepted else None,
        torque_constraint='从动真实轴力臂严格KKT；目标按名义齿数比由输入T1换算，主动重构扭矩单独核验',
        unloaded_tca_follower_angle_rad=c['follower'],unloaded_tca_te_rad=c['reference']['unloaded_te_rad'],
        mesh_unloaded_correction_rad=c['correction'],elastic_follower_lag_rad=loaded_theta,
        actual_unloaded_te_rad=c['reference']['unloaded_te_rad'],mesh_correction_rad=c['correction'],elastic_te_rad=loaded_theta,
        mesh_reference_follower_angle_rad=c['follower']+c['correction'],mesh_reference_frames=mesh_frames,
        loaded_follower_angle_rad=None if not accepted else c['follower']+c['correction']+theta,
        total_relative_te_rad=None if not accepted else c['reference']['unloaded_te_rad']+c['correction']+theta,
        elastic_follower_lag_magnitude_rad=None if not accepted else abs(theta),
        total_follower_te_rad=None if not accepted else c['reference']['unloaded_te_rad']+c['correction']+theta,
        incremental_mesh_stiffness_nm_per_rad=stiffness if accepted else None,geometry_flexibility_cache_hit=cache_hit,
        candidate_count=len(contacts),active_contact_count=len(active) if accepted else None,total_normal_force_n=float(forces.sum()) if accepted else None,
        active_area_m2=float(sum(c['records'][k]['area_m2'] for k in active)) if accepted else None,
        maximum_cell_mean_pressure_pa=max((p['mean_pressure_pa'] for p in contacts),default=0.) if accepted else None,
        pair_force_shares={str(k):float(sum(f for r,f in zip(c['records'],forces) if r['tooth_offset']==k)) for k in offsets} if accepted else None,
        contacts=contacts if accepted else [],equilibrium=asdict(equilibrium),stress_summary=summaries if accepted else {},elastic_fields=fields if accepted else {},candidate=candidate,
        tca={str(k):v for k,v in c['tca_results'].items()},selection=c['selection'],
        unloaded_reference_tooth_offset=c['reference'].get('tooth_offset',0),
        quality_checks={'small_rotation_passed':model_valid,'pinion_torque_balance_passed':pinion_balance_valid,
            'maximum_pinion_torque_relative_residual':maximum_pinion_torque_relative_residual,
            'candidate_pinion_torque_relative_residual':pinion_relative_residual},
        mesh_phase_history=c['history'],signed_lever_checks=c['sign_diagnostics'],condensation=c['condensation'],
        preparation=prepared.diagnostics,discretization_converged=None,pressure_peak_validated=False,
        limitations=['面片平均压力不是已收敛Hertz峰值；常应力单元峰值不是标准齿根应力',
            '根底仍含根圆连接约定；固定孔壁未包含轴/轴承/箱体及配合柔度',
            '未加载TCA、网格初相位与弹性TE分别报告，不把数值齿形误差称加工TE',
            '相切与候选仅覆盖显式齿对/角域，未作所有轮齿的实体侵入穷尽检查',
            '目标从动扭矩采用名义齿比；主动扭矩的离散重构残差须独立审查，不宣称两独立轴约束同时成立'])
