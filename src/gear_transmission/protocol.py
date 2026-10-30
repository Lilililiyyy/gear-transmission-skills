"""新写输入适配：显式单位换算、字段检查与现有内核对象构造，不计算接触公式。"""
from copy import deepcopy
from dataclasses import asdict
import inspect
import json
import math
from pathlib import Path

from gear_contact import Gear, GearPair, Material, OperatingPoint, analyze_pair
from gear_contact.config import modification_from_config
from gear_contact.iso_contact import ISOContactSettings
from gear_contact.extended_contact import ExtendedContactSettings

SCHEMA = 'gear_ltca_request_v1'
UNITS = {'_mm': ('_m', 1e-3), '_um': ('_m', 1e-6), '_deg': ('_rad', math.pi / 180),
         '_gpa': ('_pa', 1e9), '_mpa': ('_pa', 1e6)}


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError('JSON字段重复：' + key)
        result[key] = value
    return result


def read_json(path):
    """拒绝重复键和非标准 NaN／Infinity，避免输入被无声覆盖。"""
    def invalid(value):
        raise ValueError('JSON不接受非有限数：' + value)
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),
                      object_pairs_hook=_pairs, parse_constant=invalid)


def _scale(value, factor, path):
    if isinstance(value, list):
        return [_scale(v, factor, path) for v in value]
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(path + '的单位换算输入必须为有限数或数值数组')
    result = value * factor
    if not math.isfinite(result):
        raise ValueError(path + '换算后超出有限数范围')
    return result


def _convert(value, conversions, path=''):
    if isinstance(value, list):
        return [_convert(v, conversions, f'{path}[{i}]') for i, v in enumerate(value)]
    if not isinstance(value, dict):
        if isinstance(value, bool) and not path.endswith(('.internal', '.local_contact')):
            raise ValueError(path + '不能用布尔值替代数值或配置对象')
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(path + '必须为有限数')
        return value
    result = {}
    for key, original in value.items():
        target, converted = key, original
        for suffix, (si, factor) in UNITS.items():
            if key.endswith(suffix):
                target = key[:-len(suffix)] + si
                converted = _scale(original, factor, path + '.' + key)
                conversions.append({'input': path + '.' + key, 'output': path + '.' + target, 'factor': factor})
                break
        if target in result:
            raise ValueError(path + '.' + target + '重复指定了同一物理量，不能混用两个单位字段')
        result[target] = _convert(converted, conversions, path + '.' + target)
    return result


def _object(value, name):
    if not isinstance(value, dict):
        raise ValueError(name + '必须为JSON对象')
    return value


def _keys(value, allowed, name):
    unknown = set(value) - set(allowed)
    if unknown:
        raise ValueError(name + '含未支持字段：' + ', '.join(sorted(unknown)))


def _number(value, name, *, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name + '必须为有限实数')
    if minimum is not None and value < minimum:
        raise ValueError(name + f'必须不小于{minimum}')


def prepare_case(config, *, study=False):
    """构造原数据类并展开默认值；返回原函数可接受的命名参数和默认字段清单。"""
    config = deepcopy(_object(config, 'case'))
    defaults = []
    parameters = inspect.signature(analyze_pair).parameters
    allowed = set(parameters) - {'compliance_provider', 'boundary_gap'}
    _keys(config, allowed, 'case')
    if 'pair' not in config:
        raise ValueError('case缺少pair')
    pair_data = dict(_object(config.pop('pair'), 'case.pair'))

    def construct(cls, data, path):
        _object(data, path)
        for key, value in data.items():
            if key == 'internal' and not isinstance(value, bool):
                raise ValueError(path + '.internal必须为布尔值')
            if isinstance(value, bool) and key != 'internal':
                raise ValueError(path + '.' + key + '不能用布尔值替代数值')
        result = cls(**data)
        defaults.extend(path + '.' + key for key in asdict(result) if key not in data)
        return result

    try:
        pinion = construct(Gear, pair_data.pop('pinion'), 'case.pair.pinion')
        wheel = construct(Gear, pair_data.pop('wheel'), 'case.pair.wheel')
    except KeyError as error:
        raise ValueError('case.pair必须同时给定pinion和wheel') from error
    pair = GearPair(pinion=pinion, wheel=wheel, **pair_data)
    defaults.extend('case.pair.' + key for key in asdict(pair)
                    if key not in {'pinion', 'wheel'} and key not in pair_data)
    if study:
        if set(config) & {'torque_n_m', 'misalignment_m', 'modification_pinion'}:
            raise ValueError('比较任务的载荷、失准和小轮修形由study定义，case不能重复指定')
    elif 'torque_n_m' not in config:
        raise ValueError('单工况case必须给定torque_n_m，单位N·m')
    for name, parameter in parameters.items():
        if name in {'pair', 'compliance_provider', 'boundary_gap'}:
            continue
        if study and name in {'torque_n_m', 'misalignment_m', 'modification_pinion'}:
            continue
        if name not in config and parameter.default is not inspect.Parameter.empty:
            config[name] = parameter.default
            defaults.append('case.' + name)
    for name in ('material_pinion', 'material_wheel'):
        if isinstance(config[name], dict):
            config[name] = construct(Material, config[name], 'case.' + name)
        if not isinstance(config[name], Material):
            raise ValueError('case.' + name + '必须为材料配置对象')
    for name in ('modification_pinion', 'modification_wheel'):
        if name in config and isinstance(config[name], dict):
            config[name] = modification_from_config(config[name])
        if name in config and not hasattr(config[name], 'evaluate'):
            raise ValueError('case.' + name + '必须为修形配置对象')
    for name in ('phases', 'strips'):
        if isinstance(config[name], bool) or not isinstance(config[name], int) or config[name] < 2:
            raise ValueError(name + '必须为至少2的整数')
    for name in ('torque_n_m', 'misalignment_m', 'phase_offset'):
        if name in config:
            _number(config[name], name, minimum=0 if name == 'torque_n_m' else None)
    if not isinstance(config['local_contact'], bool):
        raise ValueError('local_contact必须为布尔值')
    if config['elastic_mode'] not in {'coupled', 'independent', 'iso', 'iso_love'}:
        raise ValueError('elastic_mode只支持coupled、independent、iso、iso_love')
    iso = config['elastic_mode'] in {'iso', 'iso_love'}
    if iso and (config['local_contact'] or config['iso_settings'] is None):
        raise ValueError('ISO模式须给定iso_settings并显式设置local_contact=false，避免重复计入压陷')
    if not iso and (config['iso_settings'] is not None or config['extended_contact_settings'] is not None):
        raise ValueError('ISO与延伸接触设置不能用于当前非ISO模式')
    for name, cls in (('iso_settings', ISOContactSettings), ('extended_contact_settings', ExtendedContactSettings)):
        if config[name] is not None:
            config[name] = construct(cls, _object(config[name], name), 'case.' + name)
    return pair, config, defaults


def prepare_request(request):
    """校验任务并展开比较对象；名称只用作标签，输出路径使用内部序号。"""
    request = deepcopy(_object(request, '请求'))
    _keys(request, {'schema', 'case', 'study', 'assumptions'}, '请求')
    if request.get('schema') != SCHEMA:
        raise ValueError('请求schema须为' + SCHEMA)
    assumptions = request.get('assumptions', [])
    if not isinstance(assumptions, list) or not all(isinstance(v, str) for v in assumptions):
        raise ValueError('assumptions须为字符串数组，记录用户授权的假设')
    conversions = []
    normalized = _convert(request, conversions)
    study = normalized.get('study')
    pair, kwargs, defaults = prepare_case(normalized.get('case'), study=study is not None)
    designs, points, pressure_limit = None, None, None
    if study is not None:
        _object(study, 'study')
        _keys(study, {'designs', 'operating_points', 'pressure_limit_pa'}, 'study')
        designs, points = {}, []
        for key in ('designs', 'operating_points'):
            if not isinstance(study.get(key), list) or not study[key]:
                raise ValueError('study.' + key + '须为非空数组')
        for record in study['designs']:
            _object(record, '修形方案')
            _keys(record, {'name', 'modification_pinion'}, '修形方案')
            name = record.get('name')
            if not isinstance(name, str) or not name.strip() or name in designs:
                raise ValueError('修形方案名称须非空且不能重复')
            if 'modification_pinion' not in record:
                raise ValueError('每个修形方案须显式给定modification_pinion')
            designs[name] = modification_from_config(record['modification_pinion'])
        for record in study['operating_points']:
            _object(record, '工况')
            _keys(record, {'name', 'torque_n_m', 'misalignment_m'}, '工况')
            if 'torque_n_m' not in record:
                raise ValueError('每个工况须显式给定torque_n_m')
            if not isinstance(record.get('name'), str) or not record['name'].strip():
                raise ValueError('工况名称须为非空字符串')
            _number(record['torque_n_m'], '工况扭矩', minimum=0)
            _number(record.get('misalignment_m', 0.0), '工况失准')
            points.append(OperatingPoint(**record))
        if len({p.name for p in points}) != len(points):
            raise ValueError('工况名称不能重复')
        pressure_limit = study.get('pressure_limit_pa')
        if pressure_limit is not None:
            _number(pressure_limit, 'pressure_limit_pa', minimum=0)
            if pressure_limit == 0:
                raise ValueError('压力上限须为正数')
    return {'normalized': normalized, 'pair': pair, 'kwargs': kwargs, 'defaults': defaults,
            'conversions': conversions, 'designs': designs, 'points': points,
            'pressure_limit_pa': pressure_limit, 'solve_count': len(designs) * len(points) if designs else 1}
