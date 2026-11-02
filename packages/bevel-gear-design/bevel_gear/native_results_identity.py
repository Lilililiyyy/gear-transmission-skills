"""原生FE材料域同源门：原mesh字节→metadata→规范制造→实际CAD文件。

本模块只校验输入身份，不把这些校验称为CAD几何或材料OUT证明。
"""
from pathlib import Path
import gc
from .native_results_archive import sha,relative,fingerprint


def normalize_manufacturing(value):
    from .generated_surface import _request
    from .generated_solid import GeneratedBlankSettings
    from .native_fe_operator import safe
    if not isinstance(value,dict) or set(value)!= {'request','blank'}:raise ValueError('规范制造须明确request/blank')
    return dict(request=safe(_request(value['request'])),blank=safe(GeneratedBlankSettings(**value['blank'])))


def source_manufacturing_contract(source):
    workflow=source['inputs'].get('workflow_inputs')
    if not isinstance(workflow,dict):raise ValueError('真实CAD域需要原workflow制造/CAD档案')
    if fingerprint(workflow)!=source['identity']['global_workflow_input_fingerprint']:raise ValueError('workflow正文没有绑定原global指纹')
    actual=workflow.get('actual_manufacturing');raw=workflow.get('base_case',{}).get('manufacturing')
    cad=workflow.get('actual_cad_sha256');cad_files=workflow.get('cad_files')
    members={'pinion','wheel'}
    if any(not isinstance(value,dict) or set(value)!=members for value in (actual,raw,cad,cad_files)):raise ValueError('缺规范actual_manufacturing/actual_cad_sha256/安全cad_files，不能认证真实CAD域')
    normalized={}
    for member in ('pinion','wheel'):
        normalized[member]=normalize_manufacturing(actual[member])
        if normalized[member]!=actual[member]:raise ValueError('actual_manufacturing必须为完整规范值')
        if normalize_manufacturing(dict(request=raw[member]['request'],blank=raw[member].get('blank',{})))!=normalized[member]:raise ValueError('原base_case显式制造/blank与规范值不同')
        digest=cad[member]
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):raise ValueError('原CAD须为明确SHA256')
        name=cad_files[member]
        if isinstance(name,dict):
            if name.get('sha256',digest)!=digest:raise ValueError('cad_files记录与actual_cad_sha256不同')
            name=name['file']
        part=Path(name)
        if part.is_absolute() or '..' in part.parts:raise ValueError('原CAD引用须为安全相对路径')
        for key in ('actual_manufacturing','actual_cad_sha256','cad_files'):
            if key in source['inputs'] and source['inputs'][key]!=workflow[key]:raise ValueError('phase重复声明的CAD/制造与workflow不同')
    return workflow,normalized,cad


def verify_native_model_identity(source,*,meshes,context,mesh_files=None,cad_files=None):
    """逐体解析同SHA原mesh的metadata；绝不只相信内存对象附带的声明。"""
    from .native_fe_operator import safe
    import numpy as np
    workflow,manufacturing,cad=source_manufacturing_contract(source);root=source.get('workflow_root')
    if mesh_files is None:
        if root is None or set(workflow.get('mesh_files',{}))!={'pinion','wheel'}:raise ValueError('须明确原mesh_files或带合法workflow根')
        mesh_files={m:relative(root,name) for m,name in workflow['mesh_files'].items()}
    if cad_files is None:
        if root is None:raise ValueError('须明确原cad_files或带合法workflow根')
        cad_files={m:relative(root,row['file'] if isinstance(row,dict) else row) for m,row in workflow['cad_files'].items()}
    if set(meshes)!={'pinion','wheel'} or set(mesh_files)!={'pinion','wheel'} or set(cad_files)!={'pinion','wheel'}:raise ValueError('两体mesh/CAD须完整明确')
    bindings={};records={}
    for member in ('pinion','wheel'):
        model=context.models[member]
        if safe(model.request)!=manufacturing[member]['request'] or safe(model.settings)!=manufacturing[member]['blank']:raise ValueError('context实际request/blank不同于原规范制造')
        if context.cad_identity.get(member)!=cad[member]:raise ValueError('context的CAD身份不同于原actual_cad_sha256')
        cad_file=Path(cad_files[member]).resolve();cad_sha=sha(cad_file)
        if cad_sha!=cad[member]:raise ValueError('原CAD文件实际字节SHA不符')
        bindings[str(cad_file)]=cad_sha
        mesh_file=Path(mesh_files[member]).resolve();mesh_sha=sha(mesh_file)
        if mesh_sha!=source['mesh_sha256'][member]:raise ValueError('材料域mesh文件不是本态原mesh')
        try:
            import orjson
            raw=orjson.loads(mesh_file.read_bytes());parser='orjson'
        except ImportError:
            import json
            with mesh_file.open(encoding='utf8') as stream:raw=json.load(stream)
            parser='json'
        if set(raw)!= {'nodes_m','tetrahedra','boundary_triangles','boundary_tags','metadata'}:raise ValueError('未知原TetraMesh文件布局')
        for name in ('nodes_m','tetrahedra','boundary_triangles','boundary_tags'):
            original=np.asarray(raw[name]);current=np.asarray(getattr(meshes[member],name))
            if original.shape!=current.shape or not np.array_equal(original,current):raise ValueError('内存mesh数值数组不同于原文件：'+name)
            del original,current
        metadata=raw['metadata'];del raw;gc.collect()
        if metadata.get('cad_binding_state')=='synthetic_FE_fixture_no_CAD':raise ValueError('数值小TET明确无CAD，不能认证真实CAD材料域')
        if (metadata.get('actual_manufacturing_request')!=manufacturing[member]['request'] or
                metadata.get('actual_blank_settings')!=manufacturing[member]['blank'] or
                metadata.get('actual_relief_cad_sha256')!=cad[member]):raise ValueError('原mesh三项CAD/制造/blank身份与workflow不同或缺失')
        if safe(meshes[member].metadata)!=metadata:raise ValueError('内存mesh metadata不同于原文件，不能替换声明')
        if sha(mesh_file)!=mesh_sha:raise ValueError('读取metadata期间原mesh改变')
        bindings[str(mesh_file)]=mesh_sha
        records[member]=dict(mesh_sha256=mesh_sha,cad_sha256=cad_sha,manufacturing=manufacturing[member],metadata_reader=parser,
            request_blank_CAD_mesh_context_all_equal=True)
    return dict(passed=True,identity_only_not_CAD_geometry_proof=True,members=records,source_files=bindings)


def bound_projection_request(source,supplied=None,*,domain_request=None):
    """自动从原规范制造取投影request；显式覆盖只能与原值完全一致。"""
    from .generated_surface import _request
    from .native_fe_operator import safe
    workflow=source['inputs'].get('workflow_inputs',{})
    actual=workflow.get('actual_manufacturing',{}).get('pinion')
    raw=actual['request'] if actual is not None else workflow.get('base_case',{}).get('manufacturing',{}).get('pinion',{}).get('request')
    if raw is None:raise ValueError('缺原小轮制造request，不能猜测材料投影')
    request=_request(raw);canonical=safe(request)
    if supplied is not None and safe(_request(supplied))!=canonical:raise ValueError('调用方pinion_request与原source不同')
    if domain_request is not None and safe(_request(domain_request))!=canonical:raise ValueError('材料域pinion_request与原source不同')
    return request


def validate_material_projection_frames(reference,expected):
    """只核小轮材料投影实际使用的Rk字段；不误把原wheel参考φ换成finalφ。"""
    import numpy as np
    if not isinstance(reference,dict):raise ValueError('缺少原物理齿材料Rk')
    a=np.asarray(reference.get('pinion_rotation_rows'),float);b=np.asarray(expected.get('pinion_rotation_rows'),float)
    if a.shape!=(3,3) or b.shape!=(3,3) or not np.isfinite(a).all() or not np.isfinite(b).all() or not np.allclose(a,b,atol=2e-12,rtol=0):raise ValueError('原物理齿材料Rk rotation与同源主动角不同')
    translation=np.asarray(reference.get('translation_m',[0.,0.,0.]),float);target=np.asarray(expected.get('translation_m',[0.,0.,0.]),float)
    if translation.shape!=(3,) or target.shape!=(3,) or not np.isfinite(translation).all() or not np.isfinite(target).all() or np.abs(translation-target).max()>1e-15:raise ValueError('原材料Rk用于投影的translation与同源生成姿态不同')
    return dict(passed=True,checked_projection_fields=['pinion_rotation_rows','translation_m'],wheel_reference_rotation_replaced=False)
