"""跨批次的真实啮合周期结果集合；保留各源协议和原始文件字节。

集合不是一次伪造的新求解业务：每个相位仍属于自己的原summary、
input_fingerprint与state。只允许相位采样、选择索引和文字标签不同；
制造、材料、载荷、接触及FE控制必须一致，FE两轮网格另以SHA256一致
绑定。大压力状态按需读取，避免将整个半解析周期的JSON数组常驻内存。
"""
from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
import hashlib
import json
import math
from numbers import Integral
from pathlib import Path
import shutil


COLLECTION_SCHEMA = 'bevel_mesh_cycle_collection_v1'


def _read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def _sha(path):
    """分块计算，避免仅为了身份检查再次读入整份大状态。"""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
        allow_nan=False, separators=(',', ':')).encode('utf-8')).hexdigest()


def _relative(root, reference):
    from .refined_fe_postprocessing import _relative_file
    return _relative_file(Path(root), reference)


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(name+'须为非负整数')
    return value


def _labels_removed(value):
    """只忽略明确的文字标签；几何、边界和所有数值控制仍严格比较。"""
    if isinstance(value, dict):
        return {k: _labels_removed(v) for k, v in value.items()
            if k not in ('name', 'description', 'label', 'title')}
    if isinstance(value, list):
        return [_labels_removed(v) for v in value]
    return value


def _physics_contract(inputs, method):
    base = _labels_removed(deepcopy(inputs['base_case']))
    key = 'refined_fe_case' if method == 'full_gear_fe' else 'surface_case'
    case = _labels_removed(deepcopy(inputs[key]))
    # 引用路径只负责找到已完整展开的base_case，不属于新的物理输入。
    case.pop('base_case', None)
    case.get('analysis', {}).pop('driver_angles_rad', None)
    if method == 'full_gear_fe':
        effective = deepcopy(inputs['effective'])
        effective.pop('phases', None); effective.pop('selected', None)
        return dict(method=method, base_case=base, case=case, effective=effective)
    analysis = deepcopy(inputs['effective_analysis'])
    analysis.pop('driver_angles_rad', None)
    return dict(method=method, base_case=base, case=case,
        effective_analysis=analysis, effective_solver=inputs['effective_solver'],
        pinion_signed_torque_nm=inputs['pinion_signed_torque_nm'],
        materials=inputs['materials'],
        maximum_pinion_torque_relative_residual=inputs['maximum_pinion_torque_relative_residual'])


class ReadOnlyCycleResults(Sequence):
    """按需读取原state.result的只读序列，不缓存整个周期或重编号原JSON。

    返回的新字典可由调用方临时处理；修改它不会回写文件。下次读取仍
    从原文件加载。每次索引均复核该状态、所属summary和inputs的SHA256。
    """

    def __init__(self, entries):
        self._entries = tuple(deepcopy(entry) for entry in entries)

    def __len__(self):
        return len(self._entries)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        if isinstance(index, bool) or not isinstance(index, Integral):
            raise TypeError('周期结果索引须为整数或切片')
        index = int(index)
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError('周期结果索引越界')
        entry = self._entries[index]
        for name in ('summary', 'input'):
            if _sha(entry[name+'_path']) != entry[name+'_sha256']:
                raise ValueError('延迟读取前原'+name+'的SHA256已变化')
        # 摘要和解码使用同一次读出的原字节，避免大状态三次I/O及检查/读取竞态。
        raw = entry['state_path'].read_bytes()
        if hashlib.sha256(raw).hexdigest() != entry['state_sha256']:
            raise ValueError('延迟读取前原state的SHA256已变化')
        state = json.loads(raw.decode('utf-8-sig'))
        del raw
        if (state.get('accepted') is not True or state.get('source_identity_unchanged') is not True or
                state.get('input_fingerprint') != entry['input_fingerprint'] or
                state.get('phase_index') != entry['source_phase_index']):
            raise ValueError('延迟读取的原状态身份或接受标志不一致')
        return state['result']


def _selected(summary, phases, accepted_phase_indices=None):
    selected = summary.get('selected_phase_indices')
    if (not isinstance(selected, list) or not selected or
            any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(phases) for i in selected) or
            len(set(selected)) != len(selected)):
        raise ValueError('源业务selected_phase_indices须为存在且不重复的原相位索引')
    rows = summary.get('states')
    if accepted_phase_indices is not None:
        # 中断批次只能由调用者明确列出已落盘相位；不按accepted自动过滤。
        chosen = accepted_phase_indices
        if (not isinstance(chosen, list) or not chosen or
                any(isinstance(i, bool) or not isinstance(i, int) or i not in selected for i in chosen) or
                len(set(chosen)) != len(chosen)):
            raise ValueError('显式接受子集须为原selected中不重复的非空索引')
        if not isinstance(rows, list):
            raise ValueError('中断源须保留原states清单')
        indices = [_integer(row.get('phase_index'), '源phase_index') for row in rows]
        if len(set(indices)) != len(indices) or any(i not in selected for i in indices):
            raise ValueError('中断源states重复或包含原selected之外的相位')
        lookup = {row['phase_index']: row for row in rows}
        if any(i not in lookup or lookup[i].get('accepted') is not True for i in chosen):
            raise ValueError('显式接受子集包含缺失或未接受的原状态')
        return list(chosen)
    if not isinstance(rows, list) or len(rows) != len(selected):
        raise ValueError('源业务须完整保存所有选中相位，不能过滤失败或缺失状态')
    indices = []
    for row in rows:
        i = _integer(row.get('phase_index'), '源phase_index')
        if i not in selected or row.get('accepted') is not True:
            raise ValueError('源业务含未选中或未接受状态')
        indices.append(i)
    if len(set(indices)) != len(indices) or set(indices) != set(selected):
        raise ValueError('源业务选中相位重复、缺失或混入其他状态')
    return selected


def _required_files(path, summary):
    """返回复算所需的原相对引用；不递归复制弹性缓存或未知目录。"""
    refs = [summary['input_file']]
    refs.extend(row['state_file'] for row in summary['states'])
    refs.extend(summary.get('mesh_files', {}).values())
    if summary.get('prepared_geometry_file'):
        refs.append(summary['prepared_geometry_file'])
    for row in summary['states']:
        if row.get('contact_csv') and _relative(path.parent, row['contact_csv']).is_file():
            refs.append(row['contact_csv'])
    files = {path.resolve()}
    for reference in refs:
        target = _relative(path.parent, reference)
        if not target.is_file():
            raise ValueError('源业务相对引用文件不存在')
        files.add(target)
    return sorted(files)


def load_cycle_source(summary_path, method=None, *, require_complete_cycle=False,
        accepted_phase_indices=None):
    """逐源验证所有选中相位，返回可延迟读取的原业务；不要求源独自覆盖周期。

    require_complete_cycle=True时还须原输入全部相位选中并覆盖完整均匀
    一齿距；False用于多个已完成批次的合并，不把部分源称为完整周期。
    accepted_phase_indices只允许半解析中断批次的明确接受子集，每个状态
    仍须通过完整接受门控及SHA核查；原summary完成标志保持原样。
    本函数复用正式FE/半解析接受门控，不运行CAD或接触求解。
    """
    from .refined_fe_postprocessing import _accepted_states
    from .pressure_domain_workflow import _surface_state
    from .contact_cycle import validate_mesh_cycle_angles
    if not isinstance(require_complete_cycle, bool):
        raise ValueError('完整周期开关须为布尔值')
    path = Path(summary_path).resolve(); summary = _read(path)
    inferred = {'bevel_refined_fe_workflow_v1': 'full_gear_fe',
        'bevel_surface_workflow_v1': 'local_halfspace'}.get(summary.get('schema'))
    if inferred is None or method not in (None, inferred):
        raise ValueError('源summary方法与所请求的FE/半解析方法不一致')
    method = inferred
    subset = accepted_phase_indices is not None
    if subset and (method != 'local_halfspace' or require_complete_cycle):
        raise ValueError('显式接受子集仅用于半解析分批集合，不能冒充单源完整周期')
    if not subset and summary.get('analysis_completed') is not True:
        raise ValueError('源业务未完成所有选中状态，不能加入正式周期集合')
    before = summary.get('source_identity_before')
    if (not isinstance(before, dict) or not before or
            before != summary.get('source_identity_after') or summary.get('source_identity_unchanged') is not True):
        raise ValueError('源计算模块的前后身份未通过')
    summary_sha = _sha(path)
    input_path = _relative(path.parent, summary['input_file']); inputs = _read(input_path)
    input_sha = _sha(input_path)
    phases = (inputs['effective']['phases'] if method == 'full_gear_fe' else
        inputs['surface_case']['analysis']['driver_angles_rad'])
    if (not isinstance(phases, list) or not phases or any(isinstance(a, bool) or
            not isinstance(a, (int, float)) or not math.isfinite(a) for a in phases)):
        raise ValueError('源相位列表须为有限弧度数组')
    selected = _selected(summary, phases, accepted_phase_indices)
    rows = {row['phase_index']: row for row in summary['states']}
    entries = []; mesh_hashes = {}
    if method == 'full_gear_fe':
        if inputs['effective'].get('selected') != selected:
            raise ValueError('FE原输入与summary的选中相位列表不一致')
        _, checked_inputs, accepted, excluded, manufacturing, _ = _accepted_states(path)
        if excluded or len(accepted) != len(selected) or {i for i, _, _ in accepted} != set(selected):
            raise ValueError('源FE含未接受、缺失或多余状态')
        bindings = {i: binding for i, _, binding in accepted}
        # 不在源加载器返回原result列表；压力数组由只读序列逐状态加载。
        del accepted, checked_inputs
        sampling = inputs['effective']['cad_sampling']
        if set(summary.get('mesh_files', {})) != {'pinion', 'wheel'}:
            raise ValueError('FE源必须绑定小轮及大轮两份原网格')
        mesh_hashes = {m: _sha(_relative(path.parent, reference))
            for m, reference in summary['mesh_files'].items()}
    else:
        manufacturing = inputs['base_case']['manufacturing']
        sampling = inputs['base_case'].get('cad_sampling', {})
        bindings = {}
        for i in selected:
            _, _, result, state_path = _surface_state(path, i)
            angle = result['geometry']['driving_angle_rad']
            if not math.isclose(angle, phases[i], rel_tol=0., abs_tol=1e-12):
                raise ValueError('半解析原相位索引与有效输入角度不一致')
            del result
            bindings[i] = dict(file=rows[i]['state_file'], sha256=_sha(state_path))
    for i in selected:
        state_path = _relative(path.parent, bindings[i]['file'])
        # 单独读取并释放协议外壳，核对未被既有半解析入口检查的原索引。
        state = _read(state_path)
        if (state.get('phase_index') != i or isinstance(state.get('phase_index'), bool) or
                not isinstance(state.get('driving_angle_rad'), (int, float)) or
                isinstance(state.get('driving_angle_rad'), bool) or
                not math.isclose(state['driving_angle_rad'], phases[i], rel_tol=0., abs_tol=1e-12) or
                state.get('input_fingerprint') != summary['input_fingerprint']):
            raise ValueError('源状态原相位索引、角度或指纹不一致')
        if (state.get('source_identity_unchanged') is not True or
                state.get('source_identity_before') != before or state.get('source_identity_after') != before):
            raise ValueError('源状态的原计算模块身份与业务汇总不一致')
        del state
        entries.append(dict(summary_path=path, summary_sha256=summary_sha,
            input_path=input_path, input_sha256=input_sha, state_path=state_path,
            state_sha256=bindings[i]['sha256'], input_fingerprint=summary['input_fingerprint'],
            source_phase_index=i, angle_rad=float(phases[i])))
    if _sha(path) != summary_sha or _sha(input_path) != input_sha:
        raise ValueError('源业务在验证过程中发生变化')
    required = _required_files(path, summary)
    file_hashes = {p: _sha(p) for p in required}
    expected = {path: summary_sha, input_path: input_sha,
        **{e['state_path']: e['state_sha256'] for e in entries},
        **{_relative(path.parent, summary['mesh_files'][m]): digest for m, digest in mesh_hashes.items()}}
    if any(file_hashes[p] != digest for p, digest in expected.items()):
        raise ValueError('源状态或网格在验证过程中发生变化')
    angles = [float(phases[i]) for i in selected]
    cycle = None
    if require_complete_cycle:
        if selected != list(range(len(phases))):
            raise ValueError('单源完整周期必须选择原输入中的全部相位')
        cycle = validate_mesh_cycle_angles(angles, manufacturing['pinion']['request']['tooth_count'])
    identity = dict(method=method, source_summary_sha256=summary_sha, input_sha256=input_sha,
        source_input_fingerprint=summary['input_fingerprint'],
        source_analysis_completed=summary.get('analysis_completed'),
        accepted_phase_indices=list(selected) if subset else None,
        states=[dict(phase_index=e['source_phase_index'], file=rows[e['source_phase_index']]['state_file'],
            sha256=e['state_sha256'], input_fingerprint=e['input_fingerprint']) for e in entries])
    return dict(summary=summary, inputs=inputs, manufacturing=manufacturing, sampling=sampling,
        results=ReadOnlyCycleResults(entries), cycle=cycle, identity=identity,
        phase_angles_rad=angles, method=method, mesh_sha256=mesh_hashes,
        physics_contract=_physics_contract(inputs, method), _entries=entries,
        _source_path=path, _required_files=required, _file_sha256=file_hashes)


def _compatible(sources):
    reference = sources[0]
    for source in sources[1:]:
        if source['physics_contract'] != reference['physics_contract']:
            raise ValueError('集合源的制造、材料、载荷、接触或离散控制不一致')
        if source['mesh_sha256'] != reference['mesh_sha256']:
            raise ValueError('集合FE源的两轮原网格SHA256不一致')
        if source['summary']['source_identity_before'] != reference['summary']['source_identity_before']:
            raise ValueError('集合源的原计算模块身份不一致')
    return _fingerprint(dict(physics=reference['physics_contract'], meshes=reference['mesh_sha256']))


def _choose_states(sources, selections, *, explicit):
    """保留全部唯一角度；不同状态的重复角只允许显式选择端点来源。"""
    available = {(j, e['source_phase_index']): e for j, source in enumerate(sources)
        for e in source['_entries']}
    ordered = sorted(available, key=lambda key: available[key]['angle_rad'])
    groups = []
    for key in ordered:
        if not groups or not math.isclose(available[key]['angle_rad'],
                available[groups[-1][0]]['angle_rad'], rel_tol=0., abs_tol=1e-12):
            groups.append([key])
        else:
            groups[-1].append(key)
    if selections is None:
        if explicit:
            raise ValueError('正式集合必须显式保存states相位选择清单')
        selections = [dict(source_index=group[0][0], phase_index=group[0][1]) for group in groups]
    if not isinstance(selections, list) or not selections:
        raise ValueError('集合states须显式列出非空原source_index/phase_index')
    chosen = []
    for item in selections:
        if not isinstance(item, dict):
            raise ValueError('集合状态选择须为具名字典')
        key = (_integer(item.get('source_index'), 'source_index'),
            _integer(item.get('phase_index'), '原phase_index'))
        if key not in available or key in chosen:
            raise ValueError('集合选择了不存在、未接受或重复的源相位')
        chosen.append(key)
    selected = set(chosen)
    for index, group in enumerate(groups):
        if len(selected.intersection(group)) != 1:
            raise ValueError('每个真实相位角必须恰好选择一个源状态，不能静默丢相')
        hashes = {available[key]['state_sha256'] for key in group}
        if len(hashes) > 1 and (not explicit or index not in (0, len(groups)-1)):
            raise ValueError('重复角的原state不同；仅允许显式选择周期端点来源')
    actual = sorted(chosen, key=lambda key: available[key]['angle_rad'])
    if explicit and chosen != actual:
        raise ValueError('集合states须按实际主动角升序排列')
    entries = [dict(available[key], source_index=key[0]) for key in actual]
    return actual, entries


def load_mesh_cycle_collection(collection_path):
    """安全加载完整一齿距集合；各源原协议与SHA逐项验证，压力结果延迟读取。"""
    from .contact_cycle import validate_mesh_cycle_angles
    path = Path(collection_path).resolve(); manifest = _read(path)
    if manifest.get('schema') != COLLECTION_SCHEMA:
        raise ValueError('需要正式啮合周期集合schema')
    if manifest.get('method') not in ('full_gear_fe', 'local_halfspace'):
        raise ValueError('集合必须明确FE或半解析方法')
    source_rows = manifest.get('sources')
    if not isinstance(source_rows, list) or not source_rows:
        raise ValueError('周期集合需要非空源summary清单')
    sources = []; seen = set()
    for row in source_rows:
        if not isinstance(row, dict):
            raise ValueError('周期集合源须为具名字典')
        summary_path = _relative(path.parent, row['summary_file'])
        if summary_path in seen:
            raise ValueError('同一源summary不能重复列入集合')
        seen.add(summary_path)
        if _sha(summary_path) != row.get('sha256'):
            raise ValueError('集合原summary的SHA256不一致')
        source = load_cycle_source(summary_path, manifest.get('method'),
            accepted_phase_indices=row.get('accepted_phase_indices'))
        if ('source_analysis_completed' in row and
                row['source_analysis_completed'] is not source['summary'].get('analysis_completed')):
            raise ValueError('集合不得改变原业务完成标志')
        if row.get('input_fingerprint') != source['identity']['source_input_fingerprint']:
            raise ValueError('集合声明的原输入指纹不一致')
        files = row.get('files')
        if not isinstance(files, list) or not files:
            raise ValueError('集合源必须保存原inputs/states/mesh文件SHA清单')
        actual_files = {}
        for item in files:
            target = _relative(path.parent, item['file'])
            if target in actual_files or _sha(target) != item.get('sha256'):
                raise ValueError('集合原文件SHA256变化或清单重复')
            actual_files[target] = item['sha256']
        if set(actual_files) != set(source['_required_files']):
            raise ValueError('集合原文件SHA清单与源summary所需文件不一致')
        sources.append(source)
    physics = _compatible(sources)
    if manifest.get('physics_fingerprint') != physics:
        raise ValueError('集合共同物理输入指纹不一致')
    chosen, entries = _choose_states(sources, manifest.get('states'), explicit=True)
    angles = [entry['angle_rad'] for entry in entries]
    cycle = validate_mesh_cycle_angles(angles, sources[0]['manufacturing']['pinion']['request']['tooth_count'])
    if manifest.get('analysis_completed') is not True:
        raise ValueError('集合尚未确认完整相位覆盖')
    states = [dict(phase_index=i, source_index=key[0], source_phase_index=key[1],
        driving_angle_rad=entry['angle_rad'], file=entry['state_path'].relative_to(path.parent).as_posix(),
        sha256=entry['state_sha256'], input_fingerprint=entry['input_fingerprint'],
        input_sha256=entry['input_sha256'], source_summary_sha256=entry['summary_sha256'])
        for i, (key, entry) in enumerate(zip(chosen, entries))]
    first = sources[0]; summary = deepcopy(manifest)
    # 集合schema始终保留，不把原state冒充新workflow。此入口仅补合法相对网格引用。
    if first['mesh_sha256']:
        summary['mesh_files'] = {m: _relative(first['_source_path'].parent, ref).relative_to(path.parent).as_posix()
            for m, ref in first['summary']['mesh_files'].items()}
    identity = dict(method=first['method'], source_summary_sha256=_sha(path),
        input_sha256=first['identity']['input_sha256'], source_input_fingerprint=physics,
        source_input_fingerprints=[s['identity']['source_input_fingerprint'] for s in sources],
        sources=deepcopy(source_rows), states=states, collection=True,
        original_input_files_modified=False, original_state_files_modified=False)
    cycle.update(collection=True, source_count=len(sources), original_state_files_modified=False)
    return dict(summary=summary, inputs=first['inputs'], manufacturing=first['manufacturing'],
        sampling=first['sampling'], results=ReadOnlyCycleResults(entries), cycle=cycle,
        identity=identity, phase_angles_rad=angles)


def create_mesh_cycle_collection(summaries, output, *, method=None, selections=None,
        accepted_source_subsets=None):
    """精确复制各源必要文件并建立集合，返回collection.json路径。

    summaries是原summary路径列表。selections可显式给出按角度升序的
    source_index/phase_index选择；默认包含全部唯一角度。重复角状态
    不同必须显式选择端点，内部重复冲突拒绝。绝不修改原输入、指纹、
    state索引或summary；目标中已有集合/源目录时拒绝覆盖。
    accepted_source_subsets可将源索引映射到明确的已接受相位列表，用于
    复用意外中断半解析批次中已保存的结果，不宣称该原业务已经完成。
    """
    from .contact_cycle import validate_mesh_cycle_angles
    if not isinstance(summaries, (list, tuple)) or not summaries:
        raise ValueError('需要非空原summary路径列表')
    paths = [Path(p).resolve() for p in summaries]
    if len(set(paths)) != len(paths):
        raise ValueError('同一原summary不能重复加入')
    subsets = {} if accepted_source_subsets is None else accepted_source_subsets
    if (not isinstance(subsets, dict) or any(isinstance(i, bool) or not isinstance(i, int)
            or not 0 <= i < len(paths) for i in subsets)):
        raise ValueError('接受子集须以存在的非负源索引为键')
    sources = [load_cycle_source(p, method, accepted_phase_indices=subsets.get(i))
        for i, p in enumerate(paths)]
    physics = _compatible(sources)
    chosen, entries = _choose_states(sources, selections, explicit=selections is not None)
    validate_mesh_cycle_angles([e['angle_rad'] for e in entries],
        sources[0]['manufacturing']['pinion']['request']['tooth_count'])
    output = Path(output).resolve(); target = output/'collection.json'
    if target.exists() or (output/'sources').exists():
        raise ValueError('集合目标已含collection.json或sources目录，拒绝覆盖')
    output.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema=COLLECTION_SCHEMA, method=sources[0]['method'], analysis_completed=True,
        physics_fingerprint=physics, sources=[], states=[dict(source_index=j, phase_index=i) for j, i in chosen],
        original_inputs_and_states_modified=False,
        definition='跨原业务的真实分相集合；保留原summary、inputs、state、指纹、索引与网格，不冒充一次新求解')
    for index, source in enumerate(sources):
        directory = output/'sources'/f'source_{index:04d}'; records = []
        for original in source['_required_files']:
            relative = original.relative_to(source['_source_path'].parent)
            destination = _relative(directory, relative.as_posix())
            before = source['_file_sha256'][original]
            if _sha(original) != before:
                raise ValueError('原文件在源验证后变化，拒绝复制')
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, destination)
            if _sha(original) != before or _sha(destination) != before:
                raise ValueError('复制过程中原文件变化，未创建正式集合')
            records.append(dict(file=destination.relative_to(output).as_posix(), sha256=before))
        summary_file = directory/source['_source_path'].name
        if _sha(summary_file) != source['identity']['source_summary_sha256']:
            raise ValueError('源summary在复制前已变化')
        manifest['sources'].append(dict(summary_file=summary_file.relative_to(output).as_posix(),
            sha256=source['identity']['source_summary_sha256'],
            source_analysis_completed=source['summary'].get('analysis_completed'),
            accepted_phase_indices=source['identity']['accepted_phase_indices'],
            input_fingerprint=source['identity']['source_input_fingerprint'], files=records))
    temporary = target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False)+'\n',
        encoding='utf-8', newline='\n')
    temporary.replace(target)
    # 按正式读取入口复核复制后的完整集合，再返回可发布容器路径。
    load_mesh_cycle_collection(target)
    return target
