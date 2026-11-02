"""真实FE边界上的节点—面接触与集总压力对偶网格，单位统一为SI。"""
from __future__ import annotations
import math


def pair_nodal_contacts(first,second,first_rotation_rows,second_rotation_rows,*,
        maximum_distance_m=2e-4,maximum_normal_angle_rad=.3,closest_backend='auto'):
    """每个独立第一面节点只设一个接触未知量，第二面使用实际三角形重心权重。

    a_i=Σ邻接三角形面积/3，节点压力p_i=f_i/a_i。显示网格是每个
    原三角形的三个中位对偶四边形，不是把旧面片压力插值成光滑图。
    可选VTK只加速真实最近点查询，法向不合格时退回完整精确搜索。
    """
    import numpy as np
    from .contact_mesh import _facet_arrays,_rotation,_closest_on_triangles
    if not math.isfinite(maximum_distance_m) or maximum_distance_m<=0:raise ValueError('候选距离须为正m')
    if not 0<maximum_normal_angle_rad<math.pi/2:raise ValueError('法向容差须为正锐角rad')
    if closest_backend not in ('auto','vtk','exact'):raise ValueError('最近点后端须为auto/vtk/exact')
    r1,r2=_rotation(first_rotation_rows,np),_rotation(second_rotation_rows,np)
    tri1,verts1,cent1,norm1,areas1=_facet_arrays(first,np)
    tri2,verts2,cent2,norm2,areas2=_facet_arrays(second,np)
    if not len(tri1) or not len(tri2):return []
    ids=np.unique(tri1);nodes=np.asarray(first['nodes_m']);second_nodes=np.asarray(second['nodes_m'])
    tributary=np.zeros(len(nodes));normals=np.zeros_like(nodes)
    for vertex in range(3):
        np.add.at(tributary,tri1[:,vertex],areas1/3)
        np.add.at(normals,tri1[:,vertex],norm1*areas1[:,None])
    normals[ids]/=np.linalg.norm(normals[ids],axis=1)[:,None]
    world_nodes=nodes@r1.T;world_vertices=verts2@r2.T;world_normals2=norm2@r2.T
    locator=None
    if closest_backend!='exact':
        try:
            import vtk
            from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray
            points=vtk.vtkPoints();points.SetData(numpy_to_vtk(second_nodes@r2.T,deep=True))
            cells=vtk.vtkCellArray();cells.SetCells(len(tri2),numpy_to_vtkIdTypeArray(
                np.column_stack((np.full(len(tri2),3),tri2)).astype(np.int64).ravel(),deep=True))
            poly=vtk.vtkPolyData();poly.SetPoints(points);poly.SetPolys(cells)
            locator=vtk.vtkStaticCellLocator();locator.SetDataSet(poly);locator.BuildLocator()
        except ImportError:
            if closest_backend=='vtk':raise ImportError('所选VTK最近点后端未安装')
    rows=[];cos_limit=math.cos(maximum_normal_angle_rad)
    for node in ids:
        point=world_nodes[node];normal1=r1@normals[node];chosen=None
        if locator is not None:
            closest=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);distance2=vtk.mutable(0.)
            locator.FindClosestPoint(point.tolist(),closest,cell,sub,distance2)
            index=int(cell)
            if 0<=index<len(tri2) and -world_normals2[index]@normal1>=cos_limit:
                if float(distance2)>maximum_distance_m**2:continue
                nearest,weights,squared=_closest_on_triangles(point,world_vertices[index:index+1],np)
                chosen=index;which=0
        if chosen is None:
            eligible=np.flatnonzero(-world_normals2@normal1>=cos_limit)
            if not len(eligible):continue
            nearest,weights,squared=_closest_on_triangles(point,world_vertices[eligible],np)
            which=int(np.argmin(squared));chosen=int(eligible[which])
        distance=math.sqrt(max(0.,float(squared[which])))
        if distance>maximum_distance_m:continue
        normal=normal1-world_normals2[chosen];normal/=np.linalg.norm(normal)
        bary=weights[which];second_point=nearest[which]
        rows.append(dict(first_node_indices=(int(node),),first_weights=(1.,),
            second_node_indices=tuple(map(int,tri2[chosen])),second_weights=tuple(map(float,bary)),
            first_point_world_m=point.tolist(),second_point_world_m=second_point.tolist(),
            normal_world=normal.tolist(),first_force_normal_local=(-r1.T@normal).tolist(),
            second_force_normal_local=(r2.T@normal).tolist(),area_m2=float(tributary[node]),
            first_pressure_node_index=int(node),second_facet_index=chosen,
            signed_gap_m=float((second_point-point)@normal),distance_m=distance,
            normal_angle_rad=math.acos(float(np.clip(-world_normals2[chosen]@normal1,-1,1))),
            contact_discretization='独立节点—实际FE三角面，节点面积集总',
            geometry_gap_model='实际FE边界节点到对面三角形，不使用半空间柔度'))
    return rows


def nodal_pressure_cells(facets,nodal_values,rotation_rows,*,tooth_offset=0):
    """由已求解节点力输出实际FE面上的中位对偶单元；逐格pA=f严格保持。

    每个四边形面积为原三角形的1/3；载荷的数值积分点仍为该节点，
    不把多边形几何中心误当有限元力矩作用点。压力为集总节点压力。
    """
    import numpy as np
    from .contact_mesh import _facet_arrays,_rotation
    triangles,vertices,centres,normals,areas=_facet_arrays(facets,np)
    rotation=_rotation(rotation_rows,np);rows=[]
    for triangle,points,centre,area in zip(triangles,vertices,centres,areas):
        for j,node in enumerate(triangle):
            value=nodal_values.get(int(node))
            if value is None or not value['active'] or value['force_n']<=0:continue
            pressure=value['force_n']/value['area_m2']
            polygon=np.array([points[j],(points[j]+points[(j+1)%3])/2,centre,(points[j]+points[(j+2)%3])/2])
            cell_area=float(area/3)
            rows.append(dict(tooth_offset=tooth_offset,first_pressure_node_index=int(node),
                first_vertices_world_m=(polygon@rotation.T).tolist(),area_m2=cell_area,
                force_n=float(pressure*cell_area),mean_pressure_pa=float(pressure),active=True,
                application_point_world_m=(points[j]@rotation.T).tolist(),
                pressure_representation='实际FE三角形的节点集总中位对偶四边形',
                integration_point='原节点，非四边形中心',remaining_gap_m=value['remaining_gap_m']))
    return rows
