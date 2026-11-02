"""读取完整已接受周期、绑定原状态与覆盖证据，输出可独立重绘的压力包络。"""
from __future__ import annotations

from pathlib import Path
from collections.abc import Mapping


def _postprocessing_identity():
    """独立保存实际参与最大值与显示的源码身份，不覆盖原力学源码记录。"""
    import hashlib
    directory = Path(__file__).resolve().parent
    return {'bevel_gear/'+name+'.py': hashlib.sha256((directory/(name+'.py')).read_bytes()).hexdigest()
        for name in ('contact_cycle', 'contact_cycle_workflow', 'contact_cycle_plotting',
            'cycle_collection', 'pressure_field', 'pressure_background', 'publication_style')}


class ReadOnlyCycleCertificates(Mapping):
    """逐证书读取及核验SHA，不将数GB分片树同时常驻内存。"""
    def __init__(self, entries):
        self._entries = dict(entries)

    def __iter__(self):
        return iter(self._entries)

    def __len__(self):
        return len(self._entries)

    def __getitem__(self, key):
        import hashlib
        import json
        path, expected = self._entries[key]
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError('全周期覆盖证书SHA256不一致')
        return json.loads(raw.decode('utf-8-sig'))


def _load_cycle(summary_path, method):
    """复用原业务验收并强制完整相位；禁止过滤失败状态后冒充全周期。"""
    from .refined_fe_postprocessing import _read
    from .cycle_collection import COLLECTION_SCHEMA, load_cycle_source, load_mesh_cycle_collection
    path = Path(summary_path).resolve()
    summary = _read(path)
    if summary.get('analysis_completed') is not True:
        raise ValueError('完整周期业务未全部接受；不能只选择成功相位绘制全过程')
    if summary.get('schema') == COLLECTION_SCHEMA:
        loaded = load_mesh_cycle_collection(path)
        if loaded['identity']['method'] != method:
            raise ValueError('周期集合方法与所请求方法不一致')
        return loaded
    return load_cycle_source(path, method, require_complete_cycle=True)


def load_cycle_domain_manifest(manifest_path):
    """安全读取通用相位/齿对证书清单；完整证书逐项绑定原压力与CAD背景。"""
    from .refined_fe_postprocessing import _read, _relative_file, _sha
    path = Path(manifest_path).resolve(); manifest = _read(path)
    if manifest.get('schema') not in ('bevel_mesh_cycle_domain_manifest_v1', 'bevel_cycle_domain_batch_v1'):
        raise ValueError('需要全周期压力域证书清单schema')
    certificates = {}
    for item in manifest.get('certificates', []):
        phase = item['phase_index']; tooth = item['tooth_offset']; key = f'{phase}:{tooth}'
        if key in certificates:
            raise ValueError('全周期覆盖清单含重复相位/齿对')
        source = _relative_file(path.parent, item.get('file', item.get('certificate_file')))
        expected = item.get('sha256', item.get('certificate_sha256'))
        if (not source.is_file() or not isinstance(expected, str) or len(expected) != 64 or
                any(c not in '0123456789abcdef' for c in expected)):
            raise ValueError('周期覆盖证书文件或SHA256格式不完整')
        certificates[key] = (source, expected)
    return ReadOnlyCycleCertificates(certificates), dict(manifest_sha256=_sha(path), certificate_count=len(certificates),
        source_summary_sha256=manifest.get('source_summary_sha256'))


def audit_mesh_cycle_pressure_workflow(summary_path, output, *, method, cad_models=None,
        cad_cache=None, maximum_depth=2, maximum_regions=300000, geometry_engine='brep',
        verify_window=False, progress=None):
    """完整周期批量域审核；共用一次CAD/弹性准备，保留每相位每齿原证书。

    manifest.complete仅表示所有请求均有真实证书，不能当作全部证书passed。
    部分证书的未决区域在最大压力图中仍须保持未知。原状态只读。
    """
    from .refined_fe_postprocessing import _load_or_build_meshes, _build_models, _validate_models, _read, _relative_file, _sha
    from .refined_fe_workflow import _write
    from .generated_halfspace import prepare_generated_halfspace
    from .generated_ltca import prepare_generated_ltca
    from .pressure_domain_batch import PressureDomainBatchAudit
    from .elasticity import TetraMesh
    path = Path(summary_path).resolve(); output = Path(output).resolve()
    if geometry_engine not in ('brep', 'triangles'):
        raise ValueError('周期几何审核引擎须为brep或triangles；三角距离按需需要VTK')
    if not isinstance(verify_window, bool):
        raise ValueError('外邻齿窗口审核开关须为布尔值')
    loaded = _load_cycle(path, method); manufacturing = loaded['manufacturing']; sampling = loaded['sampling']
    models = _build_models(manufacturing, sampling) if cad_models is None else cad_models
    _validate_models(models, manufacturing)
    meshes, cache, _ = _load_or_build_meshes(output,
        Path(cad_cache).resolve() if cad_cache is not None else output,
        manufacturing, sampling, models, flank_linear_deflection_m=1e-6)
    common = dict(cad_models=models, **sampling,
        pinion_blank_settings=manufacturing['pinion']['blank'], wheel_blank_settings=manufacturing['wheel']['blank'])
    if method == 'full_gear_fe':
        volume = {m: TetraMesh(**_read(_relative_file(path.parent, ref))) for m, ref in loaded['summary']['mesh_files'].items()}
        material = loaded['inputs']['effective']['materials']
        prepared = prepare_generated_ltca(manufacturing['pinion']['request'], manufacturing['wheel']['request'],
            pinion_material=material['pinion'], wheel_material=material['wheel'], meshes=volume,
            directory=output/'elastic_preparation', **common)
    else:
        prepared = prepare_generated_halfspace(manufacturing['pinion']['request'], manufacturing['wheel']['request'], **common)
    first = loaded['results'][0]
    requested = (sorted(int(k) for k in first['pair_force_shares']) if method == 'full_gear_fe'
        else sorted(first['geometry'].get('build_parameters', {}).get('tooth_offsets',
            [p['tooth_offset'] for p in first['geometry']['patches']])))
    del first
    context = PressureDomainBatchAudit(prepared, meshes['pinion_flank'], manufacturing['pinion']['request'],
        method=method, engine=geometry_engine, tooth_offsets=requested)
    window_context = None; window_rows = []
    if verify_window:
        if requested != list(range(min(requested), max(requested)+1)):
            raise ValueError('完整窗口审核不能跳过中间齿对')
        if method == 'full_gear_fe':
            from .pressure_window_audit import ExternalFeWindowAudit
            window_context = ExternalFeWindowAudit(models, volume, requested_tooth_offsets=requested,
                external_tooth_offsets=(min(requested)-1, max(requested)+1),
                engine='triangles' if geometry_engine == 'triangles' else 'exact')
        else:
            from .halfspace_window_audit import ExternalHalfspaceWindowAudit
            window_context = ExternalHalfspaceWindowAudit(prepared, requested_tooth_offsets=requested)
    manifest = dict(schema='bevel_mesh_cycle_domain_manifest_v1', method=method, complete=False,
        source_summary_sha256=loaded['identity']['source_summary_sha256'], cad_cache=cache,
        source_identity=loaded['identity'], maximum_depth=maximum_depth, geometry_engine=geometry_engine,
        certificates=[], failures=[])
    for i, result in enumerate(loaded['results']):
        first_certificate = None
        offsets = (sorted(int(k) for k in result['pair_force_shares']) if method == 'full_gear_fe'
            else sorted(result['geometry'].get('build_parameters', {}).get('tooth_offsets',
                [p['tooth_offset'] for p in result['geometry']['patches']])))
        for offset in offsets:
            identity = dict(phase_index=i, tooth_offset=offset)
            if progress is not None: progress(dict(stage='cycle_domain', **identity))
            try:
                certificate = context.audit(result, tooth_offset=offset, maximum_depth=maximum_depth,
                    maximum_regions=maximum_regions)
                if first_certificate is None:
                    first_certificate = certificate
                relative = f'certificates/phase_{i:04d}/k_{offset}/certificate.json'
                _write(certificate, output/relative)
                manifest['certificates'].append(dict(**identity, file=relative, sha256=_sha(output/relative),
                    passed=certificate['passed'], unresolved_region_count=certificate['unresolved_region_count']))
            except (ArithmeticError, RuntimeError, ValueError, TypeError, OSError) as error:
                manifest['failures'].append(dict(**identity, error_type=type(error).__name__, message=str(error)))
            _write(manifest, output/'manifest.json')
        if window_context is not None:
            record = dict(phase_index=i, driving_angle_rad=loaded['phase_angles_rad'][i],
                source_state_sha256=loaded['identity']['states'][i]['sha256'])
            try:
                proof = (window_context.audit(result, first_certificate) if method == 'full_gear_fe'
                    else window_context.audit(result))
                record['outer_pairs'] = [dict(row,
                    evidence_kind=row.get('evidence_kind', '完整CAD有限面和全轮FE边界双重距离下界，扣已绑定原状态全轮位移界'))
                    for row in proof['rows']]
                record['passed'] = proof['passed']
                _write(proof, output/f'window_checks/phase_{i:04d}.json')
            except (ArithmeticError, RuntimeError, ValueError, TypeError, KeyError, OSError) as error:
                record.update(passed=False, outer_pairs=[], error_type=type(error).__name__, message=str(error))
            window_rows.append(record)
    for item in loaded['identity']['states']:
        if _sha(_relative_file(path.parent, item['file'])) != item['sha256']:
            raise ArithmeticError('批量审核期间原周期压力状态已改变')
    manifest['complete'] = not manifest['failures']
    manifest['all_pressure_domains_certified'] = manifest['complete'] and all(r['passed'] for r in manifest['certificates'])
    if verify_window:
        from .contact_cycle import cycle_state_sequence_fingerprint
        passed = bool(len(window_rows) == len(loaded['phase_angles_rad']) and all(r['passed'] for r in window_rows))
        evidence = dict(schema='bevel_mesh_cycle_window_coverage_v1', method=method,
            requested_tooth_offsets=requested, outer_tooth_offsets=[min(requested)-1, max(requested)+1],
            state_sequence_sha256=cycle_state_sequence_fingerprint(method, loaded['phase_angles_rad'],
                [r['sha256'] for r in loaded['identity']['states']], requested), records=window_rows, passed=passed,
            checks=dict(all_requested_states_accepted=True, outer_pair_search_covered=passed,
                remaining_gap_nonnegative=passed, source_states_unchanged=True))
        _write({method: evidence}, output/'window_coverage.json')
        manifest['window_coverage_file'] = 'window_coverage.json'
        manifest['window_coverage_passed'] = passed
    _write(manifest, output/'manifest.json')
    return manifest


def plot_mesh_cycle_workflows(output, *, fe_summary=None, surface_summary=None,
        fe_domain_manifest=None, surface_domain_manifest=None, cad_models=None,
        cad_cache=None, spacing_m=2e-5, window_coverage=None, shared_pressure_scale=False, progress=None):
    """公开周期后处理入口；仅读取实际解，输出双语图及最大值/相位SI数据。"""
    from .refined_fe_postprocessing import _load_or_build_meshes, _read
    from .refined_fe_workflow import _write
    from .contact_cycle import build_mesh_cycle_pressure_data, validate_cycle_window_coverage
    from .contact_cycle_plotting import export_mesh_cycle_pressure_data, plot_mesh_cycle_pressure
    if fe_summary is None and surface_summary is None:
        raise ValueError('至少提供一种方法的完整周期业务汇总')
    if (fe_domain_manifest is not None and fe_summary is None or
            surface_domain_manifest is not None and surface_summary is None):
        raise ValueError('覆盖清单必须对应已提供的原周期结果')
    processing_before = _postprocessing_identity()
    output = Path(output).resolve(); data_sets = []; identities = []; reference = None
    for method, path, domains in (('full_gear_fe', fe_summary, fe_domain_manifest),
            ('local_halfspace', surface_summary, surface_domain_manifest)):
        if path is None:
            continue
        if progress is not None: progress(dict(stage='cycle_source_validation', method=method, completed=False))
        loaded = _load_cycle(path, method)
        if progress is not None: progress(dict(stage='cycle_source_validation', method=method, completed=True))
        physical = dict(manufacturing=loaded['manufacturing'],
            materials=loaded['inputs']['base_case']['materials'],
            load=loaded['inputs']['base_case']['load'])
        if reference is not None and reference != physical:
            raise ValueError('同图的两种周期模型必须采用相同实际制造输入、材料与载荷')
        reference = physical
        cache_root = Path(cad_cache).resolve() if cad_cache is not None else output
        meshes, cache, cad_mode = _load_or_build_meshes(output, cache_root, loaded['manufacturing'],
            loaded['sampling'], cad_models)
        certificates = None; domain_identity = None
        if domains is not None:
            certificates, domain_identity = load_cycle_domain_manifest(domains)
            if domain_identity['source_summary_sha256'] != loaded['identity']['source_summary_sha256']:
                raise ValueError('周期覆盖清单绑定其他业务汇总，不能混用')
        coverage = window_coverage
        if isinstance(coverage, (str, Path)):
            coverage = _read(coverage)
        elif isinstance(coverage, (list, tuple)):
            combined = {}
            for reference_path in coverage:
                item = _read(reference_path)
                if set(item) & set(combined):
                    raise ValueError('重复提供同方法窗口证据，不能覆盖已有身份')
                combined.update(item)
            coverage = combined
        evidence = coverage.get(method) if isinstance(coverage, dict) else None
        angles = loaded['phase_angles_rad']
        first_result = loaded['results'][0]
        offsets = (sorted(int(k) for k in first_result['pair_force_shares']) if method == 'full_gear_fe'
            else sorted(first_result['geometry'].get('build_parameters', {}).get('tooth_offsets',
                [p['tooth_offset'] for p in first_result['geometry']['patches']])))
        evidence = validate_cycle_window_coverage(evidence, method=method, angles_rad=angles,
            state_sha256=[r['sha256'] for r in loaded['identity']['states']], tooth_offsets=offsets)
        del first_result
        data = build_mesh_cycle_pressure_data(loaded['results'], method=method,
            pinion_request=loaded['manufacturing']['pinion']['request'], flank_mesh=meshes['pinion_flank'],
            domain_certificates=certificates, spacing_m=spacing_m, window_coverage=evidence, phase_angles_rad=angles, progress=progress)
        data['source_identity'] = loaded['identity']
        data['postprocessing_source_identity'] = processing_before
        if method == 'full_gear_fe':
            effective = loaded['inputs']['effective']
            data['solver_discretization'] = {k: effective[k] for k in ('finite_element', 'refinement', 'contact')}
            data['force_definition'] = '原FE节点法向力；压力按真实三维三角形集总面积定义'
        else:
            data['solver_discretization'] = dict(analysis=loaded['inputs']['effective_analysis'],
                solver=loaded['inputs']['effective_solver'])
            data['force_definition'] = '原半空间参考接触平面法向力；压力按原参考矩形格面积定义'
        data['cad_cache_identity'] = cache
        data['domain_manifest_identity'] = domain_identity
        export_mesh_cycle_pressure_data(data, output/'data')
        data_sets.append(data)
        identities.append(dict(**loaded['identity'], cad_mode=cad_mode,
            cad_cache_identity=cache, domain_manifest_identity=domain_identity))
    figures = [plot_mesh_cycle_pressure(data_sets, output, lang=lang,
        shared_pressure_scale=shared_pressure_scale) for lang in ('zh', 'en')]
    processing_after = _postprocessing_identity()
    if processing_before != processing_after:
        raise RuntimeError('压力包络后处理期间源码改变；原解未修改，本次输出不能视为可复现完成')
    manifest = dict(schema='bevel_mesh_cycle_plot_manifest_v1', source_identities=identities,
        figures=figures, original_physics_states_modified=False,
        full_cycle_sampled=True, phase_resolution_converged=False, pressure_mesh_converged=False,
        postprocessing_source_identity=processing_after, postprocessing_source_unchanged=True,
        single_tooth_path_coverage_confirmed=all(d['cycle']['single_tooth_path_coverage_confirmed'] for d in data_sets),
        data_files=[d['method']+'_cycle_pressure.json' for d in data_sets])
    _write(manifest, output/'mesh_cycle_plot_manifest.json')
    return manifest
