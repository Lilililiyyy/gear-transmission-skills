"""展成五面网格的标准库导出；CAD内核与图形依赖由调用方显式选择。"""
from pathlib import Path
from .io import write_json


def generated_grid_mesh(grid:dict)->dict:
    """五张参数网格连接为开放三角网格，不推断材料外法向或补造实体。"""
    if grid.get('units')!='m' or not grid.get('patches'):raise ValueError('需要单位m的完整展成网格')
    points=[];triangles=[];patch_ids=[];names=[]
    for name,patch in grid['patches'].items():
        width,profile=patch['width_count'],patch['profile_count']
        if width<2 or profile<2 or len(patch['points_m'])!=width*profile:raise ValueError('展成网格维度与点数不一致')
        offset=len(points);index=len(names);names.append(name);points.extend(patch['points_m'])
        for v in range(profile-1):
            for u in range(width-1):
                a=offset+v*width+u;b=a+1;c=a+width;d=c+1
                triangles.extend(((a,b,d),(a,d,c)));patch_ids.extend((index,index))
    return dict(units='m',vertices_m=points,triangles=triangles,face_indices=patch_ids,patch_names=names,
                metadata=dict(grid['metadata'],normal_scope='参数网格绕序，不保证材料外法向'))


def export_generated_grid(grid:dict,output_dir)->list:
    """保存完整输入/边界/区域JSON与分组OBJ，OBJ使用m和1起始顶点索引。"""
    directory=Path(output_dir);directory.mkdir(parents=True,exist_ok=True)
    files=[write_json(grid,directory/'generated_grid.json')]
    mesh=generated_grid_mesh(grid)
    target=directory/'generated_surfaces.obj'
    with target.open('w',encoding='utf-8',newline='\n') as stream:
        stream.write('# 单位m；单齿槽开放五面；根过渡由刀尖包络产生，槽底为声明的连接圆弧\n')
        for point in mesh['vertices_m']:stream.write('v '+' '.join(format(value,'.17g') for value in point)+'\n')
        previous=None
        for face,indices in zip(mesh['face_indices'],mesh['triangles']):
            if face!=previous:stream.write('g '+mesh['patch_names'][face]+'\n');previous=face
            stream.write('f '+' '.join(str(i+1) for i in indices)+'\n')
    files.append(str(target.resolve()))
    return files
