"""实际展成模型的自包含JSON/BRep档案，原单位m/rad，BRep内部mm。

加载仅读原齿轮／齿坯BRep与原五面grid，重建开放shell，不再Boolean建模。
独立槽腔slot_shapes未包含，加载后的空tuple不表示原模型没有槽腔。
模块导入仅需标准库；实际CAD读写／建shell才需要可选OCP。
"""
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil

MEMBERS=('pinion','wheel')
PATCHES=('left_flank','left_fillet','root_bottom','right_fillet','right_flank')
CAPABILITIES=dict(full_gear_shape_included=True,blank_shape_included=True,grid_included=True,
    slot_shapes_included=False,shell_reconstructed_from_grid=True,boolean_reconstruction_performed=False)
BINDING_KEYS={'descriptor_file','descriptor_sha256','gear_brep_file','gear_brep_sha256','blank_brep_file','blank_brep_sha256'}


def _sha(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block:=stream.read(2**20):value.update(block)
    return value.hexdigest()


def _safe(value):
    for name in ('tolist','item'):
        if callable(getattr(value,name,None)):return getattr(value,name)()
    raise TypeError('模型档案包含不能保存的数值类型：'+type(value).__name__)


def _plain(value):return json.loads(json.dumps(value,allow_nan=False,default=_safe))
def _same(a,b):return _plain(a)==_plain(b)


def _read(path):
    if Path(path).stat().st_size>64*2**20:raise ValueError('单模型描述超过64MiB显式上限')
    def unique(items):
        result={}
        for key,value in items:
            if key in result:raise ValueError('模型档案包含重复JSON字段')
            result[key]=value
        return result
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),object_pairs_hook=unique,
        parse_constant=lambda value:(_ for _ in ()).throw(ValueError('模型档案不接受非有限JSON常数')))


def _digest(value):
    if not isinstance(value,str) or re.fullmatch('[0-9a-f]{64}',value) is None:raise ValueError('模型档案SHA须为完整小写SHA256')
    return value


def _file(parent,name):
    if (not isinstance(name,str) or not name or name in ('.','..') or any(x in name for x in ('/','\\',':','\0'))):
        raise ValueError('模型档案引用必须是安全basename，不能路径逃逸')
    root=Path(parent).resolve();path=(root/name).resolve()
    if path.parent!=root:raise ValueError('模型档案引用越出其目录')
    return path


def _manufacturing(value):
    from .generated_surface import _request
    from .generated_solid import GeneratedBlankSettings
    if not isinstance(value,dict) or set(value)!={'request','blank'}:raise ValueError('模型档案须绑定request/blank制造身份')
    return dict(request=asdict(_request(value['request'])),blank=asdict(GeneratedBlankSettings(**value['blank'])))


def _positive(value):return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and value>0


def _groups(groups,count):
    if not isinstance(groups,dict) or not groups:raise ValueError('模型缺少原面号分组')
    for name,indices in groups.items():
        if (not isinstance(name,str) or not isinstance(indices,list) or any(isinstance(x,bool) or not isinstance(x,int) or not 0<=x<count for x in indices)
                or len(indices)!=len(set(indices))):raise ValueError('模型面号分组须为原BRep范围内不重复整数')


def _validated(descriptor_path,*,expected_descriptor_sha256,expected_member,expected_manufacturing,expected_gear_sha256):
    """只核协议／原文件字节；真实CAD拓扑在loader读原BRep后另核。"""
    path=Path(descriptor_path).resolve()
    if expected_member not in MEMBERS or _sha(path)!=_digest(expected_descriptor_sha256):raise ValueError('模型描述身份或member不符')
    record=_read(path)
    keys={'schema','member','units','angle_units','cad_internal_length_unit','request','settings','grid','diagnostics',
        'blank_diagnostics','blank_profile_m','face_groups','slot_face_groups','shell_sewing_tolerance_m','files','capabilities'}
    if (not isinstance(record,dict) or set(record)!=keys or record['schema']!='bevel_generated_model_archive_v1' or
            record['member']!=expected_member or record['units']!='m' or record['angle_units']!='rad' or record['cad_internal_length_unit']!='mm'):
        raise ValueError('模型档案schema/member/单位不支持')
    actual=_manufacturing(dict(request=record['request'],blank=record['settings']));expected=_manufacturing(expected_manufacturing)
    if not _same(actual,expected) or not _same(actual,dict(request=record['request'],blank=record['settings'])):
        raise ValueError('模型档案与明确规范制造／齿坯身份不同')
    if record['capabilities']!=CAPABILITIES:raise ValueError('模型档案能力声明不支持，slot_shapes必须明确未包含')
    if not _positive(record['shell_sewing_tolerance_m']):raise ValueError('原shell缝合容差须为正有限m')
    grid=record['grid']
    if (not isinstance(grid,dict) or grid.get('units')!='m' or grid.get('angle_units')!='rad' or
            not _same(grid.get('inputs'),record['request']) or set(grid.get('patches',{}))!=set(PATCHES)):
        raise ValueError('原grid的制造／单位／五面与档案不一致')
    count=0;dimensions=None
    for patch in grid['patches'].values():
        width,profile=patch.get('width_count'),patch.get('profile_count')
        if any(isinstance(x,bool) or not isinstance(x,int) or x<4 for x in (width,profile)):raise ValueError('原grid尺寸须支持三次曲面')
        if dimensions is not None and dimensions!=(width,profile):raise ValueError('原五面grid尺寸不同')
        dimensions=(width,profile);points=patch.get('points_m');count+=width*profile
        if (not isinstance(points,list) or len(points)!=width*profile or any(not isinstance(p,list) or len(p)!=3 or
                any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in p) for p in points)):
            raise ValueError('原grid点必须保持有限SI三坐标和原索引布局')
    if count>2_000_000:raise ValueError('原grid超过两百万点上限')
    diagnostics=record['diagnostics'];blank=record['blank_diagnostics']
    if not isinstance(diagnostics,dict) or not isinstance(blank,dict):raise ValueError('模型缺少原拓扑诊断')
    faces=diagnostics.get('faces')
    if isinstance(faces,bool) or not isinstance(faces,int) or faces<1:raise ValueError('模型缺少原面数诊断')
    _groups(record['face_groups'],faces)
    slots=record['slot_face_groups']
    if not isinstance(slots,dict) or set(slots)!={str(i) for i in range(actual['request']['tooth_count'])}:raise ValueError('完整齿轮档案须保留每一物理齿的面号')
    for groups in slots.values():_groups(groups,faces)
    profile=record['blank_profile_m']
    if (not isinstance(profile,list) or len(profile)<3 or any(not isinstance(p,list) or len(p)!=2 or
            any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in p) for p in profile)):
        raise ValueError('原齿坯剖面须保留有限SI(radius,z)布局')
    files=record['files']
    if not isinstance(files,dict) or set(files)!={'gear_brep','blank_brep'}:raise ValueError('模型档案须包含原整轮和齿坯BRep')
    resolved={}
    for name,binding in files.items():
        if not isinstance(binding,dict) or set(binding)!={'file','sha256'}:raise ValueError('BRep绑定字段不符')
        file=_file(path.parent,binding['file']);digest=_digest(binding['sha256'])
        if not file.is_file() or _sha(file)!=digest:raise ValueError('原BRep文件SHA改变或缺失')
        resolved[name]=file
    if files['gear_brep']['file']==files['blank_brep']['file'] or files['gear_brep']['sha256']!=_digest(expected_gear_sha256):raise ValueError('模型档案原CAD身份与workflow不同')
    binding=dict(descriptor_file=path.name,descriptor_sha256=_sha(path),
        gear_brep_file=files['gear_brep']['file'],gear_brep_sha256=files['gear_brep']['sha256'],
        blank_brep_file=files['blank_brep']['file'],blank_brep_sha256=files['blank_brep']['sha256'])
    return record,resolved,binding


def _copy_original(source,target,digest):
    source=Path(source).resolve();target=Path(target)
    if target.exists():
        if not target.is_file() or _sha(target)!=digest:raise FileExistsError('模型档案目标已有不同字节，拒绝覆盖')
    else:
        try:os.link(source,target)
        except FileExistsError:
            if _sha(target)!=digest:raise FileExistsError('模型档案目标争用，拒绝覆盖')
        except OSError:
            with source.open('rb') as src,target.open('xb') as dst:shutil.copyfileobj(src,dst,2**20)
    if _sha(source)!=digest or _sha(target)!=digest:raise ArithmeticError('模型原文件链接／复制后SHA不同')


def save_generated_model_archive(model,descriptor_path,*,member,gear_brep_path=None,blank_brep_path=None):
    """保存原完整模型描述及双BRep；显式原BRep参数用于保留本次已有导出字节。"""
    from .generated_solid import GeneratedSolidModel
    from .cad_import import ImportedCadModel,export_model_brep
    if not isinstance(model,GeneratedSolidModel) or member not in MEMBERS:raise TypeError('需要明确member和GeneratedSolidModel')
    path=Path(descriptor_path).resolve()
    if path.exists():raise FileExistsError('模型描述已存在，拒绝覆盖')
    path.parent.mkdir(parents=True,exist_ok=True);files={}
    for name,source,shape,diagnostics,target_name in (
            ('gear_brep',gear_brep_path,model.shape,model.diagnostics,member+'.brep'),
            ('blank_brep',blank_brep_path,model.blank.shape,model.blank.diagnostics,member+'_blank.brep')):
        target=path.parent/target_name
        if source is None:
            if target.exists():raise FileExistsError('BRep目标已存在，须显式提供原文件以复用其字节')
            export_model_brep(ImportedCadModel(shape,diagnostics),target)
        else:_copy_original(source,target,_sha(source))
        files[name]=dict(file=target.name,sha256=_sha(target))
    record=_plain(dict(schema='bevel_generated_model_archive_v1',member=member,units='m',angle_units='rad',cad_internal_length_unit='mm',
        request=asdict(model.request),settings=asdict(model.settings),grid=model.grid,diagnostics=model.diagnostics,
        blank_diagnostics=model.blank.diagnostics,blank_profile_m=model.blank.profile_m,face_groups=model.face_groups,
        slot_face_groups=model.slot_face_groups,shell_sewing_tolerance_m=model.shell.diagnostics['sewing_tolerance_m'],
        files=files,capabilities=CAPABILITIES))
    with path.open('x',encoding='utf8') as stream:json.dump(record,stream,ensure_ascii=False,allow_nan=False,separators=(',',':'));stream.write('\n')
    return _validated(path,expected_descriptor_sha256=_sha(path),expected_member=member,
        expected_manufacturing=dict(request=asdict(model.request),blank=asdict(model.settings)),expected_gear_sha256=files['gear_brep']['sha256'])[2]


def load_generated_model_archive(descriptor_path,*,expected_descriptor_sha256,expected_member,expected_manufacturing,expected_gear_sha256):
    """仅读原BRep和grid恢复完整工作面模型；slot_shapes=()明确未加载独立槽腔。"""
    record,files,_=_validated(descriptor_path,expected_descriptor_sha256=expected_descriptor_sha256,expected_member=expected_member,
        expected_manufacturing=expected_manufacturing,expected_gear_sha256=expected_gear_sha256)
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.TopoDS import TopoDS_Shape
    from OCP.TopExp import TopExp
    from OCP.TopAbs import TopAbs_FACE,TopAbs_SOLID
    from OCP.TopTools import TopTools_IndexedMapOfShape
    from .generated_surface import _request
    from .generated_solid import GeneratedBlankSettings,GeneratedBlank,GeneratedSolidModel
    from .cad import build_generated_shell
    def shape(name,diagnostics):
        value=TopoDS_Shape()
        if not BRepTools.Read_s(value,str(files[name]),BRep_Builder()) or value.IsNull() or not BRepCheck_Analyzer(value).IsValid():raise ValueError('原BRep读取或真实拓扑有效性失败')
        for label,kind in (('faces',TopAbs_FACE),('solids',TopAbs_SOLID)):
            mapping=TopTools_IndexedMapOfShape();TopExp.MapShapes_s(value,kind,mapping)
            if diagnostics.get(label)!=mapping.Extent():raise ValueError('原BRep真实拓扑计数与原描述不同')
        return value
    request=_request(record['request']);settings=GeneratedBlankSettings(**record['settings'])
    blank=GeneratedBlank(shape('blank_brep',record['blank_diagnostics']),record['blank_diagnostics'],
        tuple(tuple(p) for p in record['blank_profile_m']),request,settings)
    model=GeneratedSolidModel(shape('gear_brep',record['diagnostics']),record['diagnostics'],blank,
        build_generated_shell(record['grid'],sewing_tolerance_m=record['shell_sewing_tolerance_m']),record['grid'],(),
        record['face_groups'],{int(k):v for k,v in record['slot_face_groups'].items()})
    return model


def copy_generated_model_archive(descriptor_path,destination,*,expected_descriptor_sha256,expected_member,expected_manufacturing,expected_gear_sha256):
    """逐level链接／复制原描述和依赖字节，不重新读CAD或重建模型。"""
    _,files,binding=_validated(descriptor_path,expected_descriptor_sha256=expected_descriptor_sha256,expected_member=expected_member,
        expected_manufacturing=expected_manufacturing,expected_gear_sha256=expected_gear_sha256)
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    for kind,file in files.items():_copy_original(file,destination/binding[kind+'_file'],binding[kind+'_sha256'])
    _copy_original(descriptor_path,destination/binding['descriptor_file'],binding['descriptor_sha256'])
    return binding


def load_native_model_archives(workflow_inputs,workflow_root):
    """按同workflow六字段绑定恢复两模型，不回退旧缓存或再次Boolean。"""
    bindings=workflow_inputs.get('generated_model_archives');manufacturing=workflow_inputs.get('actual_manufacturing')
    cad=workflow_inputs.get('actual_cad_sha256');cad_files=workflow_inputs.get('cad_files')
    if any(not isinstance(x,dict) or set(x)!=set(MEMBERS) for x in (bindings,manufacturing,cad,cad_files)):
        raise ValueError('workflow缺少双模型自包含档案／制造／CAD绑定')
    expected={m:_manufacturing(manufacturing[m]) for m in MEMBERS}
    if 'base_case' in workflow_inputs:
        base=workflow_inputs['base_case']['manufacturing']
        if any(not _same(_manufacturing(base[m]),expected[m]) for m in MEMBERS):raise ValueError('raw base制造默认值与实际模型档案不同')
    paths={}
    for m in MEMBERS:
        binding=bindings[m]
        if not isinstance(binding,dict) or set(binding)!=BINDING_KEYS or binding['gear_brep_file']!=cad_files[m] or binding['gear_brep_sha256']!=cad[m]:
            raise ValueError('workflow模型档案六字段与原CAD绑定不一致')
        descriptor=_file(workflow_root,binding['descriptor_file'])
        _,_,actual=_validated(descriptor,expected_descriptor_sha256=binding['descriptor_sha256'],expected_member=m,
            expected_manufacturing=expected[m],expected_gear_sha256=cad[m])
        if not _same(actual,binding):raise ValueError('workflow模型档案依赖与原描述不同')
        paths[m]=descriptor
    models={}
    for m,descriptor in paths.items():
        binding=bindings[m]
        models[m]=load_generated_model_archive(descriptor,expected_descriptor_sha256=binding['descriptor_sha256'],expected_member=m,
            expected_manufacturing=expected[m],expected_gear_sha256=cad[m])
    return models
