"""直齿实体FE承载图与可追溯源数据；仅绘制真实活跃三角片。

失败结果拒绝成功图；面片平均压力不是Hertz峰值。网格数值初相位
与弹性滞后分开，载荷扫描不把二者之和称为加工传动误差。
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path


_ARCSEC_PER_RAD = 180/math.pi*3600


def _get(value,key):
    return value[key] if isinstance(value,dict) else getattr(value,key)


def _number(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError(f"{name}须为有限实数")
    return float(value)


def _converged(result):
    if not isinstance(result,dict) or result.get("converged") is not True or not isinstance(result.get("equilibrium"),dict) or result["equilibrium"].get("converged") is not True:
        raise ValueError("未收敛的FE接触结果不能导出成功图或成功源数据")
    return result


def _write_csv(path,rows,header=None):
    with path.open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else header)
        writer.writeheader()
        writer.writerows(rows)
    return str(path.resolve())


def _contact_data(result,prepared):
    _converged(result)
    nodes=_get(_get(_get(prepared,"systems")["pinion"],"mesh"),"nodes_m")
    pitch=_number(_get(_get(_get(prepared,"geometry"),"pinion"),"pitch_cone_angle_rad"),"节锥角")
    if not 0<pitch<math.pi/2:raise ValueError("小轮节锥角须为锐角")
    contacts=result.get("contacts",[])
    if not contacts:raise ValueError("成功承载图需要实际接触候选面片")
    rows,vertices,polygons=[],[],[]
    force_by_pair={str(k):0. for k in result.get("pair_force_shares",{})}
    for index,contact in enumerate(contacts):
        ids=contact["first_node_indices"]
        if len(ids)!=3 or any(isinstance(i,bool) or not isinstance(i,int) or not 0<=i<len(nodes) for i in ids):
            raise ValueError("第一接触面须为实际FE全局节点三角形")
        points=[tuple(float(x) for x in nodes[i]) for i in ids]
        if any(len(p)!=3 or not all(math.isfinite(v) for v in p) for p in points):
            raise ValueError("FE面片节点坐标无效")
        ab=[points[1][j]-points[0][j] for j in range(3)]
        ac=[points[2][j]-points[0][j] for j in range(3)]
        cross=(ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0])
        geometric_area=.5*math.sqrt(sum(v*v for v in cross))
        area=_number(contact["area_m2"],"面片面积")
        force=_number(contact["force_n"],"面片合力")
        pressure=_number(contact["mean_pressure_pa"],"面片平均压力")
        if area<=0 or not math.isclose(area,geometric_area,rel_tol=1e-6,abs_tol=1e-20):
            raise ValueError("压力面片面积与实际FE节点三角形不一致")
        if not math.isclose(pressure,force/area,rel_tol=1e-8,abs_tol=1e-6):
            raise ValueError("面片平均压力与合力/面积不一致")
        if not isinstance(contact.get("active"),bool):raise ValueError("面片active须为明确布尔状态")
        active=contact["active"]
        if active and force<=0:raise ValueError("活跃压力面片须有正合力")
        offset=contact["tooth_offset"]
        if isinstance(offset,bool) or not isinstance(offset,int):raise ValueError("齿偏移须为整数")
        force_by_pair[str(offset)]=force_by_pair.get(str(offset),0.)+force
        row={"contact_index":index,"active":active,"tooth_offset":offset,"area_m2":area,
             "force_n":force,"mean_pressure_pa":pressure,"mean_pressure_mpa":pressure/1e6}
        for key in ("signed_gap_m","remaining_gap_m","distance_m","normal_angle_rad","pinion_lever_m","wheel_lever_m"):
            row[key]=_number(contact[key],key)
        for prefix in ("first","second"):
            if len(contact[prefix+"_node_indices"])!=3 or len(contact[prefix+"_weights"])!=3:
                raise ValueError("接触节点与权重须为完整三角形的三个值")
            if any(_number(w,"节点权重")<0 for w in contact[prefix+"_weights"]) or not math.isclose(sum(contact[prefix+"_weights"]),1.,abs_tol=1e-10):
                raise ValueError("接触权重须非负且和为1")
            for j,node in enumerate(contact[prefix+"_node_indices"]):row[f"{prefix}_node_{j}"]=int(node)
            for j,weight in enumerate(contact[prefix+"_weights"]):row[f"{prefix}_weight_{j}"]=_number(weight,"节点权重")
            for field in ("point_world_m","force_normal_local"):
                for axis,value in zip("xyz",contact[prefix+"_"+field]):row[prefix+"_"+field+"_"+axis]=_number(value,field)
        for axis,value in zip("xyz",contact["normal_world"]):row["normal_world_"+axis]=_number(value,"世界法向")
        row["second_closest_feature"]=contact.get("second_closest_feature","")
        polygon=[]
        for local,node_index in enumerate(ids):
            point=points[local]
            radius=math.sqrt(sum(v*v for v in point))
            if radius<=0:raise ValueError("FE齿面节点不能位于锥顶原点")
            polar=math.atan2(math.hypot(point[0],point[1]),point[2])
            x,y=radius*1000,radius*(polar-pitch)*1000
            polygon.append((x,y))
            vertices.append(dict(contact_index=index,active=active,tooth_offset=offset,local_vertex=local,
                node_index=node_index,x_m=point[0],y_m=point[1],z_m=point[2],cone_distance_m=radius,
                polar_angle_rad=polar,pitch_cone_angle_rad=pitch,plot_x_mm=x,plot_y_mm=y))
        row["plot_centroid_x_mm"]=sum(p[0] for p in polygon)/3
        row["plot_centroid_y_mm"]=sum(p[1] for p in polygon)/3
        rows.append(row)
        if active:polygons.append(dict(contact_index=index,tooth_offset=offset,vertices_mm=polygon,mean_pressure_mpa=pressure/1e6))
    if not polygons:raise ValueError("没有真实活跃面片，不能生成成功压力图")
    if len(polygons)!=result["active_contact_count"] or len(rows)!=result["candidate_count"]:raise ValueError("面片数与结果摘要不一致")
    total=sum(row["force_n"] for row in rows)
    if not math.isclose(total,_number(result["total_normal_force_n"],"总法向力"),rel_tol=1e-8,abs_tol=1e-8):
        raise ValueError("源数据面片合力和与结果摘要不一致")
    for key,value in result.get("pair_force_shares",{}).items():
        if not math.isclose(force_by_pair[str(key)],float(value),rel_tol=1e-8,abs_tol=1e-8):raise ValueError("每齿对合力与候选面片不一致")
    pair_rows=[dict(tooth_offset=int(key),normal_force_n=value) for key,value in sorted(force_by_pair.items(),key=lambda item:int(item[0]))]
    return rows,vertices,polygons,pair_rows


def export_straight_ltca_source_data(result,prepared,output_dir):
    """只导出真实物理源数据，不加载绘图后端；失败或不一致数据拒绝导出。"""
    rows,vertices,polygons,pair_rows=_contact_data(result,prepared)
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    files=[_write_csv(output/"straight_ltca_contacts.csv",rows),
        _write_csv(output/"straight_ltca_facet_vertices.csv",vertices),
        _write_csv(output/"straight_ltca_pair_forces.csv",pair_rows)]
    summary=[]
    for key in ("driving_angle_rad","pinion_torque_nm","wheel_target_torque_nm","pinion_torque_reconstructed_nm",
        "wheel_torque_reconstructed_nm","analytic_unloaded_te_rad","mesh_unloaded_correction_rad",
        "elastic_follower_lag_rad","total_follower_te_rad","incremental_mesh_stiffness_nm_per_rad",
        "total_normal_force_n","active_contact_count","candidate_count","active_area_m2","maximum_cell_mean_pressure_pa"):
        unit=("N·m/rad" if key.endswith("nm_per_rad") else "rad" if key.endswith("_rad") else
              "N·m" if key.endswith("_nm") else "N" if key.endswith("_n") else
              "m²" if key.endswith("_m2") else "Pa" if key.endswith("_pa") else "1")
        summary.append(dict(field=key,value=result.get(key),unit=unit,source_field="result."+key))
    for key in ("balance_residual","relative_balance_residual","complementarity_residual_nm","contact_stationarity_residual_m","iterations","status"):
        unit="N·m" if key in ("balance_residual","complementarity_residual_nm") else "m" if key.endswith("_m") else "1"
        summary.append(dict(field=key,value=result["equilibrium"].get(key),unit=unit,source_field="result.equilibrium."+key))
    files.append(_write_csv(output/"straight_ltca_summary.csv",summary))
    sources={}
    for key in rows[0]:
        origin="result.contacts[i]."+key
        if key=="contact_index":origin="result.contacts中的0起始序号"
        elif key=="mean_pressure_mpa":origin="result.contacts[i].mean_pressure_pa/1e6"
        elif key.startswith("plot_centroid_"):origin="实际第一三角片三个绘图顶点坐标的算术平均"
        for prefix in ("first","second"):
            if key.startswith(prefix+"_node_"):origin=f"result.contacts[i].{prefix}_node_indices[{key.rsplit('_',1)[1]}]"
            elif key.startswith(prefix+"_weight_"):origin=f"result.contacts[i].{prefix}_weights[{key.rsplit('_',1)[1]}]"
            elif key.startswith(prefix+"_point_world_m_") or key.startswith(prefix+"_force_normal_local_"):
                origin="result.contacts[i]."+key[:-2]+f"[{ 'xyz'.index(key[-1])}]"
        if key.startswith("normal_world_"):origin=f"result.contacts[i].normal_world[{'xyz'.index(key[-1])}]"
        sources[key]=origin
    manifest={"converged":True,"coordinate_definition":{"plot_x_mm":"norm(prepared.systems.pinion.mesh.nodes_m[node])*1000",
        "plot_y_mm":"norm(node)*(atan2(hypot(x,y),z)-prepared.geometry.pinion.pitch_cone_angle_rad)*1000"},
        "contact_field_sources":sources,
        "facet_vertex_sources":{"node_index":"result.contacts[i].first_node_indices[local_vertex]",
            "x_m/y_m/z_m":"prepared.systems.pinion.mesh.nodes_m[node_index]",
            "cone_distance_m":"实际节点向量模长","polar_angle_rad":"atan2(hypot(x_m,y_m),z_m)",
            "pitch_cone_angle_rad":"prepared.geometry.pinion.pitch_cone_angle_rad"},
        "pair_force_source":"同一tooth_offset的result.contacts[i].force_n求和，并核验result.pair_force_shares",
        "units":{"nodes":"m","force":"N","pressure":"Pa，绘图MPa","area":"m²","angles":"rad"},
        "active_polygon_count":len(polygons),"candidate_count":len(rows),"plotting":"仅活跃实际FE三角PolyCollection，无插值/平滑",
        "limitations":result.get("limitations",[]),"phase_definition":"网格数值初相位独立于弹性滞后，不作为加工TE"}
    path=output/"straight_ltca_source_manifest.json"
    path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    return files+[str(path.resolve())]


def plot_straight_ltca(result,prepared,output_dir,*,publication=True):
    """183mm双栏真实活跃面片平均压力与每齿对合力，默认600dpi及SVG/PDF。"""
    rows,vertices,polygons,pair_rows=_contact_data(result,prepared)
    from .plotting import _pyplot,_size,_panel,_save
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    plt=_pyplot()
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    offsets=sorted({item["tooth_offset"] for item in polygons})
    fig=plt.figure(figsize=_size(183,max(85,31*len(offsets)+22)),constrained_layout=True)
    outer=fig.add_gridspec(1,2,width_ratios=[2.2,1.])
    left=outer[0].subgridspec(len(offsets),1)
    right=outer[1].subgridspec(2,1,height_ratios=[3,1.3])
    pressure_axes=[]
    from .halfspace_plotting import pressure_colour_scale
    normalization,pressure_cmap,_=pressure_colour_scale([p['mean_pressure_mpa']*1e6 for p in polygons])
    all_x=[p[0] for item in polygons for p in item["vertices_mm"]]
    all_y=[p[1] for item in polygons for p in item["vertices_mm"]]
    for index,offset in enumerate(offsets):
        ax=fig.add_subplot(left[index]);pressure_axes.append(ax)
        cells=[item for item in polygons if item["tooth_offset"]==offset]
        collection=PolyCollection([item["vertices_mm"] for item in cells],
            array=[item["mean_pressure_mpa"] for item in cells],cmap=pressure_cmap,norm=normalization,
            edgecolors="#70777B",linewidths=.18)
        ax.add_collection(collection)
        ax.update_datalim(list(zip(all_x,all_y)));ax.autoscale_view();ax.margins(.06)
        ax.set_ylabel("齿廓坐标 / mm")
        if index==len(offsets)-1:ax.set_xlabel("锥距 / mm")
        else:ax.tick_params(labelbottom=False)
        heading=f"实际活跃面片 · 齿偏移 {offset:+d}"
        if index==0:_panel(ax,"a",heading)
        else:ax.set_title(heading,loc="left",fontsize=8)
    colorbar=fig.colorbar(collection,ax=pressure_axes,orientation="horizontal",pad=.06,aspect=35,shrink=.85)
    colorbar.set_label("面片平均压力 / MPa（非 Hertz 峰值）")
    ax=fig.add_subplot(right[0])
    positions=list(range(len(pair_rows)))
    forces=[row["normal_force_n"] for row in pair_rows]
    ax.bar(positions,forces,color="#2166AC",width=.6)
    ax.set(xticks=positions,xticklabels=[f"{row['tooth_offset']:+d}" for row in pair_rows],
           xlabel="齿对偏移",ylabel="法向合力 / N")
    ax.margins(y=.18);_panel(ax,"b","齿对承载")
    note=fig.add_subplot(right[1]);note.axis("off")
    eq=result["equilibrium"]
    note.text(0,1,"收敛残差\n"+f"力矩：{eq['balance_residual']:.2e} N·m\n"+
        f"互补：{eq['complementarity_residual_nm']:.2e} N·m\n"+
        f"接触间隙：{eq['contact_stationarity_residual_m']:.2e} m\n"+
        f"活跃/候选：{result['active_contact_count']}/{result['candidate_count']}",va="top",fontsize=7)
    files=_save(fig,output,"straight_fe_ltca",publication=publication,title="真实FE面片的承载接触")
    plt.close(fig)
    return files+export_straight_ltca_source_data(result,prepared,output)


def _sweep_rows(results):
    if len(results)<2:raise ValueError("载荷扫描至少需要两个实际收敛结果")
    for result in results:_converged(result)
    reference=results[0]
    if not isinstance(reference.get("preparation"),dict) or not reference["preparation"]:
        raise ValueError("载荷扫描需要模型准备信息，不能省略工况一致性核对")
    rows=[]
    for index,result in enumerate(results):
        if result.get("preparation")!=reference.get("preparation") or not math.isclose(float(result["driving_angle_rad"]),float(reference["driving_angle_rad"]),abs_tol=1e-12,rel_tol=0):
            raise ValueError("载荷曲线须使用相同准备模型和主动角，不能混合不同离散工况")
        torque=_number(result["pinion_torque_nm"],"小轮扭矩")
        lag=_number(result["elastic_follower_lag_rad"],"弹性从动滞后")
        stiffness=_number(result["incremental_mesh_stiffness_nm_per_rad"],"增量刚度")
        if torque<=0 or stiffness<=0:raise ValueError("载荷扫描要求正扭矩和正增量刚度")
        rows.append(dict(input_index=index,pinion_torque_nm=torque,wheel_target_torque_nm=result["wheel_target_torque_nm"],
            elastic_follower_lag_rad=lag,elastic_follower_lag_arcsec=lag*_ARCSEC_PER_RAD,
            incremental_mesh_stiffness_nm_per_rad=stiffness,incremental_mesh_stiffness_knm_per_rad=stiffness/1000,
            mesh_unloaded_correction_rad=_number(result["mesh_unloaded_correction_rad"],"网格初相位"),
            analytic_unloaded_te_rad=_number(result["analytic_unloaded_te_rad"],"解析无载TE"),
            active_contact_count=result["active_contact_count"],active_area_m2=result["active_area_m2"],
            total_normal_force_n=result["total_normal_force_n"],balance_residual_nm=result["equilibrium"]["balance_residual"],
            complementarity_residual_nm=result["equilibrium"]["complementarity_residual_nm"]))
    rows=sorted(rows,key=lambda row:row["pinion_torque_nm"])
    if len({row["pinion_torque_nm"] for row in rows})!=len(rows):raise ValueError("重复扭矩不自动平均或伪装成重复统计")
    return rows


def export_ltca_load_sweep_source_data(results,output_dir):
    """载荷扫描CSV保留原rad和N·m单位，单独列网格数值初相位。"""
    rows=_sweep_rows(list(results))
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    files=[_write_csv(output/"straight_ltca_load_sweep.csv",rows)]
    manifest={"converged":True,"source":"每行对应实际收敛result，input_index保留原输入顺序",
        "derived_fields":{"elastic_follower_lag_arcsec":"elastic_follower_lag_rad*180/pi*3600",
            "incremental_mesh_stiffness_knm_per_rad":"incremental_mesh_stiffness_nm_per_rad/1000"},
        "phase_definition":"只绘弹性从动滞后；mesh_unloaded_correction_rad单列，不混为加工TE",
        "statistics":"确定性载荷样本，无统计误差线；线段只连接实算点"}
    manifest["field_sources"]={key:("输入列表中的0起始位置" if key=="input_index" else
        manifest["derived_fields"][key] if key in manifest["derived_fields"] else
        "result.equilibrium.balance_residual" if key=="balance_residual_nm" else
        "result.equilibrium.complementarity_residual_nm" if key=="complementarity_residual_nm" else "result."+key) for key in rows[0]}
    path=output/"straight_ltca_load_sweep_manifest.json"
    path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    return files+[str(path.resolve())]


def plot_ltca_load_sweep(results,output_dir,*,publication=True):
    """183mm双栏绘实际扭矩—弹性从动滞后及冻结活跃集的增量刚度。"""
    results=list(results)
    rows=_sweep_rows(results)
    from .plotting import _pyplot,_size,_panel,_save
    plt=_pyplot();output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=_size(183,75),constrained_layout=True)
    torque=[row["pinion_torque_nm"] for row in rows]
    axes[0].plot(torque,[row["elastic_follower_lag_arcsec"] for row in rows],"o-",color="#2166AC")
    axes[1].plot(torque,[row["incremental_mesh_stiffness_knm_per_rad"] for row in rows],"s-",color="#616970")
    for ax in axes:ax.set_xlabel("小轮扭矩 / (N·m)")
    axes[0].set_ylabel("弹性从动滞后 / arcsec")
    axes[1].set_ylabel("增量扭转刚度 / (kN·m/rad)")
    _panel(axes[0],"a","弹性角响应");_panel(axes[1],"b","活跃集增量刚度")
    files=_save(fig,output,"straight_ltca_load_sweep",publication=publication,title="实算载荷扫描的弹性响应")
    plt.close(fig)
    return files+export_ltca_load_sweep_source_data(results,output)


def _phase_cycle_rows(results,pinion_teeth):
    """按真实主动角排序，保留两个齿距端点；不取模或补点。"""
    if isinstance(pinion_teeth,bool) or not isinstance(pinion_teeth,int) or pinion_teeth<=0:
        raise ValueError("pinion_teeth须为正整数")
    if len(results)<2:raise ValueError("离散齿距图至少需要两个实际收敛姿态")
    for result in results:_converged(result)
    reference=results[0]
    preparation=reference.get("preparation")
    if not isinstance(preparation,dict) or not preparation:raise ValueError("离散齿距图需要准备网格信息")
    torque=_number(reference["pinion_torque_nm"],"小轮扭矩")
    wheel_torque=_number(reference["wheel_target_torque_nm"],"从动目标扭矩")
    if torque<=0 or wheel_torque<=0:raise ValueError("离散齿距图要求正扭矩")
    pitch=2*math.pi/pinion_teeth
    rows=[]
    for index,result in enumerate(results):
        if result.get("preparation")!=preparation:
            raise ValueError("离散齿距图不能混合不同准备网格")
        if not math.isclose(_number(result["pinion_torque_nm"],"小轮扭矩"),torque,rel_tol=1e-10,abs_tol=1e-12) or not math.isclose(_number(result["wheel_target_torque_nm"],"从动扭矩"),wheel_torque,rel_tol=1e-10,abs_tol=1e-12):
            raise ValueError("离散齿距图不能混合不同扭矩")
        angle=_number(result["driving_angle_rad"],"主动角")
        elastic=_number(result["elastic_follower_lag_rad"],"弹性从动滞后")
        mesh=_number(result["mesh_unloaded_correction_rad"],"网格数值初相位")
        stiffness=_number(result["incremental_mesh_stiffness_nm_per_rad"],"增量刚度")
        if stiffness<=0:raise ValueError("增量刚度须为正")
        rows.append(dict(input_index=index,driving_angle_rad=angle,driving_angle_deg=math.degrees(angle),
            absolute_pitch_coordinate=angle/pitch,pinion_torque_nm=torque,wheel_target_torque_nm=wheel_torque,
            elastic_follower_lag_rad=elastic,elastic_follower_lag_arcsec=elastic*_ARCSEC_PER_RAD,
            mesh_unloaded_correction_rad=mesh,mesh_unloaded_correction_arcsec=mesh*_ARCSEC_PER_RAD,
            analytic_unloaded_te_rad=_number(result["analytic_unloaded_te_rad"],"解析无载TE"),
            incremental_mesh_stiffness_nm_per_rad=stiffness,incremental_mesh_stiffness_knm_per_rad=stiffness/1000,
            active_contact_count=result["active_contact_count"],candidate_count=result["candidate_count"],
            active_area_m2=result["active_area_m2"],total_normal_force_n=result["total_normal_force_n"],
            balance_residual_nm=result["equilibrium"]["balance_residual"],
            complementarity_residual_nm=result["equilibrium"]["complementarity_residual_nm"]))
    rows=sorted(rows,key=lambda row:row["driving_angle_rad"])
    if len({row["driving_angle_rad"] for row in rows})!=len(rows):raise ValueError("重复主动角不平均、不折叠")
    span=rows[-1]["driving_angle_rad"]-rows[0]["driving_angle_rad"]
    if span>pitch+max(1e-12,pitch*1e-10):raise ValueError("主动角跨度超过一个齿距，须显式分组后再绘图")
    first=rows[0]["driving_angle_rad"]
    for index,row in enumerate(rows):
        row["sample_index"]=index
        row["relative_pitch_phase"]=(row["driving_angle_rad"]-first)/pitch
    return rows,pitch


def export_straight_phase_cycle_source_data(results,output_dir,*,pinion_teeth):
    """离散周期CSV保留绝对主动角、两个端点及两类角度，不混为加工TE。"""
    rows,pitch=_phase_cycle_rows(list(results),pinion_teeth)
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    files=[_write_csv(output/"straight_ltca_phase_cycle.csv",rows)]
    derived={"driving_angle_deg":"degrees(result.driving_angle_rad)",
        "absolute_pitch_coordinate":"result.driving_angle_rad/(2*pi/pinion_teeth)",
        "relative_pitch_phase":"(主动角-本批最小主动角)/(2*pi/pinion_teeth)",
        "elastic_follower_lag_arcsec":"result.elastic_follower_lag_rad*180/pi*3600",
        "mesh_unloaded_correction_arcsec":"result.mesh_unloaded_correction_rad*180/pi*3600",
        "incremental_mesh_stiffness_knm_per_rad":"result.incremental_mesh_stiffness_nm_per_rad/1000"}
    sources={key:(derived[key] if key in derived else "输入列表0起始序号" if key=="input_index" else
        "按真实主动角排序后的0起始序号" if key=="sample_index" else
        "result.equilibrium.balance_residual" if key=="balance_residual_nm" else
        "result.equilibrium.complementarity_residual_nm" if key=="complementarity_residual_nm" else "result."+key) for key in rows[0]}
    endpoint_keys=("elastic_follower_lag_rad","elastic_follower_lag_arcsec","mesh_unloaded_correction_rad",
        "mesh_unloaded_correction_arcsec","incremental_mesh_stiffness_nm_per_rad")
    span=(rows[-1]["driving_angle_rad"]-rows[0]["driving_angle_rad"])/pitch
    manifest={"converged":True,"sample_count":len(rows),"pinion_teeth":pinion_teeth,"pinion_pitch_rad":pitch,
        "reference_driving_angle_rad":rows[0]["driving_angle_rad"],"coverage_in_pitch":span,
        "covers_one_pitch":math.isclose(span,1.,rel_tol=0,abs_tol=1e-10),"periodicity_enforced":False,
        "endpoint_differences_end_minus_start":{key:rows[-1][key]-rows[0][key] for key in endpoint_keys},
        "field_sources":sources,"grid_convergence_verified":False,
        "interpolation":"只连接实际离散点，未插值补点、平滑、闭合或FFT",
        "caption":f"{len(rows)}个实际FE收敛姿态，覆盖{span:.6g}个小轮齿距。弹性从动滞后与网格数值初相位分线显示，后者不作为加工TE。端点差原样保留；采样分辨率有限，尚未验证网格收敛。"}
    path=output/"straight_ltca_phase_cycle_manifest.json"
    path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    return files+[str(path.resolve())]


def plot_straight_phase_cycle(results,output_dir,*,pinion_teeth,publication=True):
    """183mm双栏离散一齿距角响应与刚度；两个端点不强制相等。"""
    results=list(results)
    rows,pitch=_phase_cycle_rows(results,pinion_teeth)
    from .plotting import _pyplot,_size,_panel,_save
    plt=_pyplot();output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=_size(183,82),constrained_layout=True)
    fig.set_layout_engine("constrained",rect=(0,.06,1,.94))
    phase=[row["relative_pitch_phase"] for row in rows]
    axes[0].plot(phase,[row["elastic_follower_lag_arcsec"] for row in rows],"o-",color="#2166AC",label="弹性从动滞后")
    axes[0].plot(phase,[row["mesh_unloaded_correction_arcsec"] for row in rows],"s--",color="#616970",label="网格数值初相位")
    axes[0].set_ylabel("角响应 / arcsec")
    axes[0].legend(loc="best")
    axes[1].plot(phase,[row["incremental_mesh_stiffness_knm_per_rad"] for row in rows],"o-",color="#616970")
    axes[1].set_ylabel("增量扭转刚度 / (kN·m/rad)")
    for ax in axes:
        ax.set_xlabel("主动角距起点 / 小轮齿距")
        if len(phase)<=9:ax.set_xticks(phase)
    _panel(axes[0],"a","离散角响应（非加工TE）")
    _panel(axes[1],"b","离散啮合刚度")
    fig.text(.5,.012,f"{len(rows)}个实算点；保留端点差；未验证网格收敛",ha="center",fontsize=7)
    files=_save(fig,output,"straight_ltca_phase_cycle",publication=publication,title="离散齿距内的弹性与网格角响应")
    plt.close(fig)
    return files+export_straight_phase_cycle_source_data(results,output,pinion_teeth=pinion_teeth)
