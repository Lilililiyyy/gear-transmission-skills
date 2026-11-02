"""将显式展成面与零偏置圆弧齿轮名义几何绑定，长度m、角度rad。

齿厚来自节锥上两槽面的方位角差；角度来自工作面的局部微分。
设计名义β/α与实际两侧值分别保存，不把名义近似门限称为制造公差。
"""
from dataclasses import asdict, dataclass
import math

from .geometry import GeometryInput, GeometryResult, calculate_geometry
from .generated_cutter import evaluate_generated_point
from .generated_surface import _request


@dataclass(frozen=True)
class BindingTolerances:
    """对象绑定门限；角度是名义面近似接受范围，不是制造验收公差。"""
    length_absolute_m: float = 2e-8
    length_relative: float = 1e-7
    angle_absolute_rad: float = 1e-9
    nominal_spiral_angle_rad: float = math.pi / 180
    nominal_pressure_angle_rad: float = math.pi / 180

    def __post_init__(self):
        for key, value in asdict(self).items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f'{key}须为非负有限数')


def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _derivative(a, b, step):
    return tuple((x-y)/(2*step) for x, y in zip(a, b))


def _pitch_section(request, fraction, nominal_spiral_angle_rad):
    """在真实节锥截线上求点与切平面；中央差分长度随锥距缩放。"""
    machine = request.machine
    delta = machine.pitch_cone_angle_rad
    station = machine.reference_cone_distance_m-request.face_width_m*fraction
    radius = station*math.sin(delta)
    step = machine.reference_cone_distance_m*1e-5
    tolerance = machine.reference_cone_distance_m*1e-12
    def value(r, s):
        return evaluate_generated_point(machine, 2*r, s, height_tolerance_m=tolerance)
    center = value(radius, station)
    radial_plus, radial_minus = value(radius+step, station), value(radius-step, station)
    station_plus, station_minus = value(radius, station+step), value(radius, station-step)
    flanks = {}
    for side in ('left', 'right'):
        data = center['flanks'][side]
        if data['region'] != 'straight_blade' or data['search_state'] != 'converged':
            raise ValueError('节锥测量点须处于已收敛的直刃工作面，不能用刀尖区外推名义角')
        point = data['point_m']
        azimuth = math.atan2(point[1], point[0])
        radial = (math.cos(azimuth), math.sin(azimuth), 0.)
        circumferential = (-radial[1], radial[0], 0.)
        generator = (radial[0]*math.sin(delta), radial[1]*math.sin(delta), math.cos(delta))
        cone_normal = (radial[0]*math.cos(delta), radial[1]*math.cos(delta), -math.sin(delta))
        dr = _derivative(radial_plus['flanks'][side]['point_m'], radial_minus['flanks'][side]['point_m'], step)
        ds = _derivative(station_plus['flanks'][side]['point_m'], station_minus['flanks'][side]['point_m'], step)
        trace = tuple(ds[i]+math.sin(delta)*dr[i] for i in range(3))
        normal = _cross(dr, ds)
        if not all(math.isfinite(x) for x in trace+normal) or _dot(normal, normal) <= 0 or _dot(trace, generator) <= 0:
            raise ValueError('节锥微分退化或宽向参数方向异常，不能返回默认零角')
        beta = math.atan2(_dot(trace, circumferential), _dot(trace, generator))
        alpha = math.atan2(abs(_dot(normal, cone_normal)), math.hypot(_dot(normal, generator), _dot(normal, circumferential)))
        flanks[side] = {'point_m': list(point), 'azimuth_rad': azimuth,
                        'signed_spiral_angle_rad': beta, 'spiral_angle_rad': abs(beta),
                        'normal_pressure_angle_rad': alpha}
    angular_slot = abs(math.remainder(flanks['right']['azimuth_rad']-flanks['left']['azimuth_rad'], 2*math.pi))
    pitch = 2*math.pi/request.tooth_count
    if not 0 < angular_slot < pitch:
        raise ValueError('节锥槽宽须在零与齿距之间，当前制造面不能绑定为正齿厚')
    transverse = radius*(pitch-angular_slot)
    measured_beta = sum(row['spiral_angle_rad'] for row in flanks.values())/2
    return {'face_fraction': fraction, 'cone_distance_m': station, 'pitch_radius_m': radius,
            'angular_slot_width_rad': angular_slot, 'transverse_tooth_thickness_m': transverse,
            'nominal_normal_tooth_thickness_m': transverse*math.cos(nominal_spiral_angle_rad),
            'actual_normal_tooth_thickness_m': transverse*math.cos(measured_beta),
            'flanks': flanks}


def measure_manufacturing_geometry(request, *, nominal_spiral_angle_rad, fractions=(0., .5, 1.)):
    """测量外/中/内节锥截面；法向名义齿厚投影使用显式设计β。"""
    request = _request(request)
    if isinstance(nominal_spiral_angle_rad, bool) or not math.isfinite(nominal_spiral_angle_rad) or not 0 <= nominal_spiral_angle_rad < math.pi/2:
        raise ValueError('设计螺旋角须为[0,π/2)中的有限弧度')
    if not fractions or any(isinstance(f, bool) or not math.isfinite(f) or not 0 <= f <= 1 for f in fractions):
        raise ValueError('宽向测量因子须在[0,1]内且不能为空')
    delta = request.machine.pitch_cone_angle_rad
    top_taper = request.tip_cone_angle_rad-delta
    root_taper = delta-request.machine.machine_root_angle_rad
    half_width = request.face_width_m/2
    return {'tooth_count': request.tooth_count, 'outer_pitch_diameter_m': request.tooth_count*request.module_m,
            'face_width_m': request.face_width_m, 'pitch_cone_angle_rad': delta,
            'outer_cone_distance_m': request.machine.reference_cone_distance_m,
            'cutter_radius_m': request.machine.cutter_diameter_m/2,
            'addendum_taper_angle_rad': top_taper, 'dedendum_taper_angle_rad': root_taper,
            'mean_addendum_m': request.addendum_m-half_width*math.tan(top_taper),
            'mean_dedendum_m': request.machine.dedendum_m-half_width*math.tan(root_taper),
            'sections': [_pitch_section(request, f, nominal_spiral_angle_rad) for f in fractions],
            'measurement': '节锥交线中央差分；齿厚为材料角宽的弧长，不是两侧点间直线弦长'}


def derive_nominal_geometry(pinion_request, wheel_request, *, spiral_angle_rad, pressure_angle_rad):
    """从制造尺寸反算名义几何；β/α为显式设计值，不由接触结果拟合。"""
    if isinstance(pressure_angle_rad, bool) or not math.isfinite(pressure_angle_rad) or not 0 < pressure_angle_rad < math.pi/2:
        raise ValueError('设计法向压力角须为(0,π/2)中的有限弧度')
    requests = tuple(_request(r) for r in (pinion_request, wheel_request))
    data = tuple(measure_manufacturing_geometry(r, nominal_spiral_angle_rad=spiral_angle_rad, fractions=(.5,)) for r in requests)
    p, w = data
    for field in ('outer_cone_distance_m', 'face_width_m', 'cutter_radius_m'):
        if not math.isclose(p[field], w[field], rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError(f'两轮{field}不一致，不能反算当前共同锥顶名义副')
    if not math.isclose(p['addendum_taper_angle_rad'], w['dedendum_taper_angle_rad'], abs_tol=1e-10) or not math.isclose(p['dedendum_taper_angle_rad'], w['addendum_taper_angle_rad'], abs_tol=1e-10):
        raise ValueError('两轮齿高锥度不互补，当前名义副不能表达该制造请求')
    mean = p['outer_cone_distance_m']-p['face_width_m']/2
    mn = 2*mean*math.sin(w['pitch_cone_angle_rad'])/w['tooth_count']*math.cos(spiral_angle_rad)
    ha = (p['mean_addendum_m']+w['mean_addendum_m'])/(2*mn)
    hf = (p['mean_dedendum_m']+w['mean_dedendum_m'])/(2*mn)
    shift = (p['mean_addendum_m']-w['mean_addendum_m'])/(2*mn)
    root_shift = (w['mean_dedendum_m']-p['mean_dedendum_m'])/(2*mn)
    if not math.isclose(shift, root_shift, rel_tol=1e-9, abs_tol=1e-10):
        raise ValueError('两轮齿根高与齿顶高不能由同一变位表达，请显式扩展模型而非平均抹平')
    sp, sw = (row['sections'][0]['nominal_normal_tooth_thickness_m'] for row in data)
    backlash = (math.pi*mn-sp-sw)*math.cos(pressure_angle_rad)
    if backlash < -1e-12:
        raise ValueError('制造齿厚合计超过名义齿距，不能反算为非负总侧隙')
    xs = (sp-sw)/(4*mn)-shift*math.tan(pressure_angle_rad)
    return GeometryInput(p['tooth_count'], w['tooth_count'], w['outer_pitch_diameter_m'], p['face_width_m'],
        shaft_angle_rad=p['pitch_cone_angle_rad']+w['pitch_cone_angle_rad'], spiral_angle_rad=spiral_angle_rad,
        pressure_angle_rad=pressure_angle_rad, addendum_coefficient=ha, dedendum_coefficient=hf,
        pinion_profile_shift_coefficient=shift, tooth_thickness_modification_coefficient=xs,
        backlash_m=max(0., backlash), cutter_radius_m=p['cutter_radius_m'], depth_taper='custom', tooth_trace='circular_arc',
        wheel_addendum_taper_angle_rad=w['addendum_taper_angle_rad'], wheel_dedendum_taper_angle_rad=w['dedendum_taper_angle_rad'])


def bind_manufacturing_geometry(geometry, pinion_request, wheel_request, *, tolerances=None):
    """比较名义尺寸/弧齿厚与制造测量；不满足门限时返回拒绝原因。"""
    if isinstance(geometry, GeometryInput):
        geometry = calculate_geometry(geometry)
    if not isinstance(geometry, GeometryResult):
        raise TypeError('绑定须给GeometryInput或GeometryResult')
    if tolerances is None:
        tolerances = BindingTolerances()
    if isinstance(tolerances, dict):
        tolerances = BindingTolerances(**tolerances)
    if not isinstance(tolerances, BindingTolerances):
        raise TypeError('tolerances须为BindingTolerances或具名字典')
    if geometry.diagnostics['tooth_trace'] != 'circular_arc':
        raise ValueError('当前制造绑定只支持零偏置圆弧齿线')
    measured, checks, angular = {}, [], []
    for member, request in (('pinion', pinion_request), ('wheel', wheel_request)):
        row = measure_manufacturing_geometry(request, nominal_spiral_angle_rad=geometry.inputs.spiral_angle_rad)
        measured[member] = row
        dimensions = getattr(geometry, member)
        pairs = {'tooth_count': dimensions.teeth, 'outer_pitch_diameter_m': dimensions.outer_pitch_diameter_m,
                 'face_width_m': geometry.inputs.face_width_m, 'pitch_cone_angle_rad': dimensions.pitch_cone_angle_rad,
                 'outer_cone_distance_m': geometry.outer_cone_distance_m, 'cutter_radius_m': geometry.inputs.cutter_radius_m,
                 'addendum_taper_angle_rad': dimensions.addendum_taper_angle_rad,
                 'dedendum_taper_angle_rad': dimensions.dedendum_taper_angle_rad,
                 'mean_addendum_m': dimensions.mean_addendum_m, 'mean_dedendum_m': dimensions.mean_dedendum_m}
        pairs['mean_normal_tooth_thickness_m'] = dimensions.mean_normal_tooth_thickness_m
        row['mean_normal_tooth_thickness_m'] = row['sections'][1]['nominal_normal_tooth_thickness_m']
        for key, expected in pairs.items():
            actual = row[key]
            if key == 'tooth_count':
                limit = 0.
            elif key.endswith('_rad'):
                limit = tolerances.angle_absolute_rad
            else:
                limit = tolerances.length_absolute_m+tolerances.length_relative*abs(expected)
            residual = actual-expected
            checks.append({'member': member, 'field': key, 'nominal': expected, 'actual': actual,
                           'residual': residual, 'limit': limit, 'accepted': abs(residual) <= limit})
        # 全宽β变化属于名义圆弧几何；逐侧α与β的差保留，不能当严格相同面。
        nominal_betas = (dimensions.outer_spiral_angle_rad, dimensions.mean_spiral_angle_rad, dimensions.inner_spiral_angle_rad)
        for section, nominal_beta in zip(row['sections'], nominal_betas):
            for side, flank in section['flanks'].items():
                for key, expected, limit in (('spiral_angle_rad', abs(nominal_beta), tolerances.nominal_spiral_angle_rad),
                                             ('normal_pressure_angle_rad', geometry.inputs.pressure_angle_rad, tolerances.nominal_pressure_angle_rad)):
                    residual = flank[key]-expected
                    angular.append({'member': member, 'side': side, 'face_fraction': section['face_fraction'],
                        'field': key, 'nominal': expected, 'actual': flank[key], 'residual': residual,
                        'limit': limit, 'accepted': abs(residual) <= limit})
    # 非零螺旋角的相交轴外啮合副须为相反旋向；不能只比较β幅值。
    if geometry.inputs.spiral_angle_rad > tolerances.angle_absolute_rad:
        signed = [sum(s['signed_spiral_angle_rad'] for s in measured[m]['sections'][1]['flanks'].values())/2
                  for m in ('pinion', 'wheel')]
        product = signed[0]*signed[1]
        opposite = product < 0
        checks.append({'member': 'pair', 'field': 'opposite_spiral_hands', 'nominal': -1,
                       'actual': -1 if opposite else 1, 'residual': 0 if opposite else 2,
                       'limit': 0, 'accepted': opposite})
    failures = [row for row in checks+angular if not row['accepted']]
    return {'accepted': not failures, 'tolerances': asdict(tolerances), 'dimensional_checks': checks,
            'nominal_angle_checks': angular, 'measured': measured, 'failures': failures,
            'scope': '同尺寸/中点齿高与弧齿厚绑定；β/α为显式名义面近似，不宣称实际双侧共轭或制造验收',
            'unverified': ['成对全周期共轭', '实际槽底完整刀具包络', '制造/试验公差认证']}
