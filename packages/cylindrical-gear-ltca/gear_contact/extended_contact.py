"""新写的圆柱齿轮齿顶延伸接触几何与附加间隙模型。

几何由渐开线、齿顶圆和转角相切关系计算，不扩大渐开线边界冒充实际
延伸曲线。长度m、角度rad；修正量采用法向正间隙，N个中心与N+1个
端点分别保存。当前为几何模型，不能替代有限齿体弹性或物理验证。
"""
from __future__ import annotations
from dataclasses import dataclass
from math import atan, cos, hypot, isfinite, isnan, pi, sin, sqrt, tan
import numpy as np
from .geometry import ContactGrid, GearPair


@dataclass
class _Line:
    source_line_index: int
    face_centres_m: tuple
    face_endpoints_m: tuple
    gear_a_roll_centres_m: tuple
    gear_b_roll_centres_m: tuple
    gear_a_roll_endpoints_m: tuple
    gear_b_roll_endpoints_m: tuple
    widths_m: tuple
    eligible: tuple
    gear_a_face_centres_m: tuple
    gear_b_face_centres_m: tuple
    contact_length_m: float
    source_gear_a_tooth_index: int
    source_gear_b_tooth_index: int
    partner_roll_mode: str
    diagnostics: tuple = ()
    correction_transverse_m: tuple | None = None
    gear_a_extended_centres: tuple | None = None
    gear_b_extended_centres: tuple | None = None
    gear_a_extended_endpoints: tuple | None = None
    gear_b_extended_endpoints: tuple | None = None


class _Gear:
    def __init__(self, gear):
        self.gear = gear
        self.sign = -1 if gear.internal else 1
        self.internal = gear.internal
        self.signed_teeth = self.sign*gear.teeth
        self.helix_angle_rad = abs(gear.helix_angle_rad)
        self.hand = 'left' if gear.helix_angle_rad < 0 else 'right'
        self.face_width_m = gear.face_width_m
        self.base_radius_m = self.sign*gear.base_radius_m
        self.base_diameter_m = 2*self.base_radius_m
        self.reference_diameter_m = self.sign*2*gear.reference_radius_m
        self.base_helix_angle_rad = abs(gear.base_helix_angle_rad)
        self.transverse_base_pitch_m = 2*pi*gear.base_radius_m/gear.teeth
        self.root_form_diameter_m = self.sign*2*gear.root_form_radius_m
        self.tip_form_diameter_m = self.sign*2*gear.tip_form_radius_m
        self.source_is_spur = self.helix_angle_rad <= 1e-8

    def roll_from_diameter(self, diameter):
        return .5*sqrt(max(0., diameter*diameter-self.base_diameter_m**2))

    def diameter_from_roll(self, roll, *, signed=False):
        diameter=sqrt(4*roll*roll+self.base_diameter_m**2)
        return self.sign*diameter if signed else diameter

    def below_root_form(self, roll):
        return self.diameter_from_roll(roll, signed=True) < self.root_form_diameter_m


class _Mesh:
    def __init__(self, pair, settings):
        self.pair = pair
        self.gear_a, self.gear_b = _Gear(pair.pinion), _Gear(pair.wheel)
        self.flank = 'right'
        self.center_distance_m = pair.center_distance_m
        self.base_radius_sum_m = self.gear_a.base_radius_m+self.gear_b.base_radius_m
        self.mesh_sign = -1 if pair.wheel.internal else 1
        self.working_pressure_angle_rad = pair.working_pressure_angle_rad
        self.effective_face_width_m = pair.face_overlap_m
        self.extension_roll_factor = settings.roll_extension_factor
        self.root_tolerance_policy = settings.root_tolerance_policy
        self.resolved_extended_face_width_factor = (1. if settings.face_extension_factor is None else settings.face_extension_factor)
        self.transverse_contact_ratio = pair.transverse_contact_ratio
        self.overlap_contact_ratio = pair.overlap_ratio

    def partner_roll_from_osep(self, roll):
        return sqrt(self.center_distance_m**2-self.base_radius_sum_m**2)-self.mesh_sign*roll

    def active_roll_bounds(self, gear='a'):
        low,high=self.pair.contact_roll_bounds_m
        if gear=='a':return low,high
        return self.partner_roll_from_osep(high),self.partner_roll_from_osep(low)

    def local_face_coordinate(self, gear, face):
        axial=self.pair.face_bounds_m[0]+face
        if gear=='a':return axial+self.pair.pinion.face_width_m/2
        return axial-self.pair.axial_offset_m+self.pair.wheel.face_width_m/2

    def geometric_root(self, function, minimum, maximum, seed):
        if self.root_tolerance_policy == "dimensionless_refined":
            scale = self.pair.pinion.reference_radius_m**2
            return _geometric_root(lambda value: function(value)/scale,
                minimum, maximum, seed, x_tolerance=1e-12, function_tolerance=1e-14)
        return _geometric_root(function, minimum, maximum, seed)



def _geometric_root(function, minimum, maximum, seed, x_tolerance=1e-10, function_tolerance=1e-10, evaluation_limit=100):
    points = {}
    evaluation_count = 0
    last_x = None

    def evaluate(x):
        nonlocal evaluation_count, last_x
        if x in points:
            return (x, points[x])
        if not isfinite(x):
            return (x, float('nan'))
        y = float(function(x))
        evaluation_count += 1
        last_x = x
        if not isnan(y):
            points[x] = y
        return (x, y)

    def bracket():
        ordered = sorted(points.items())
        for (a, b) in zip(ordered, ordered[1:]):
            if (a[1] > 0) - (a[1] < 0) != (b[1] > 0) - (b[1] < 0):
                return (a, b)
        return None

    def good(point):
        return point is not None and abs(point[1]) < function_tolerance

    def verge(x):
        return abs(x - minimum) <= x_tolerance or abs(x - maximum) <= x_tolerance

    def clip(x):
        return min(maximum, max(minimum, x))

    def secant_interval(a, b):
        (xa, ya) = a
        (xb, yb) = b
        remaining = evaluation_limit - evaluation_count
        for _ in range(max(0, remaining) + 1):
            if abs(ya) <= abs(yb) and abs(ya) < function_tolerance:
                return xa
            if abs(yb) < function_tolerance:
                return xb
            if abs(xb - xa) < x_tolerance:
                return (xa + xb) / 2
            if ya == yb:
                return float('nan')
            (x, y) = evaluate((xb * ya - xa * yb) / (ya - yb))
            if y * ya >= 0:
                (xa, ya) = (x, y)
            elif y * ya < 0:
                (xb, yb) = (x, y)
            else:
                return float('nan')
        return float('nan')

    def dip_search(a, b, current):
        sign = 1 if a[1] > 0 else -1 if a[1] < 0 else 0
        c = evaluate(a[0] + 0.3819660112501052 * (b[0] - a[0]))
        d = evaluate(a[0] + 0.6180339887498948 * (b[0] - a[0]))
        if isnan(c[1]) or isnan(d[1]):
            return (False, current)
        while evaluation_count < evaluation_limit and all(((p[1] > 0) - (p[1] < 0) == sign for p in (a, b, c, d))):
            if good(current) or abs(c[0] - d[0]) < x_tolerance:
                return (True, current)
            if a[1] * sign > c[1] * sign and d[1] * sign < b[1] * sign:
                if c[1] * sign < d[1] * sign:
                    (b, d) = (d, c)
                    c = evaluate(a[0] + 0.3819660112501052 * (b[0] - a[0]))
                    if abs(c[1]) < function_tolerance:
                        current = c
                else:
                    (a, c) = (c, d)
                    d = evaluate(a[0] + 0.6180339887498948 * (b[0] - a[0]))
                    if abs(d[1]) < function_tolerance:
                        current = d
                if isnan(c[1]) or isnan(d[1]):
                    return (False, current)
                continue
            return (True, current)
        return (True, current)
    if not all((isfinite(x) for x in (minimum, maximum, seed))) or maximum < minimum:
        return float('nan')
    if maximum - minimum <= x_tolerance:
        p = evaluate((minimum + maximum) / 2)
        return p[0] if good(p) else float('nan')
    initial_step = (maximum - minimum) / 20
    evaluate(seed)
    first = next(iter(sorted(points.items())), None)
    if first is None:
        first = evaluate((minimum + maximum) / 2)
    (current, far) = (first, first)
    step = initial_step
    reverse = False
    ordinal = 0
    while len(points) < 2 and evaluation_count < evaluation_limit:
        ordinal += 1
        x = clip(far[0] + step * ordinal)
        if x != current[0]:
            current = evaluate(x)
        if verge(current[0]):
            if reverse:
                break
            step = -step
            current = far
            reverse = True
    if not points:
        return float('nan')
    if len(points) == 1:
        current = far = next(iter(points.items()))
        previous = None
        step = initial_step
    else:
        by_y = sorted(points.items(), key=lambda p: p[1])
        first_y = next(iter(sorted(points.items())))[1]
        if first_y > 0:
            (current, previous, far) = (by_y[0], by_y[1], by_y[-1])
        else:
            (current, previous, far) = (by_y[-1], by_y[-2], by_y[0])
        step = current[0] - previous[0]
    if good(current):
        return current[0]

    def travel(current, previous, step):
        while True:
            interval = bracket()
            if interval:
                return (secant_interval(*interval), current)
            if evaluation_count > evaluation_limit or step == 0:
                return (float('nan'), current)
            x = clip(current[0] + step)
            step = x - current[0]
            if previous is not None and previous[1] > current[1] and (previous[1] != current[1]):
                extrapolated = previous[0] - previous[1] * (current[0] - previous[0]) / (current[1] - previous[1])
                if minimum <= extrapolated <= maximum and (step > 0 and extrapolated > x or (step < 0 and extrapolated < x)):
                    p = evaluate(extrapolated)
                    if not isnan(p[1]):
                        (previous, current) = (current, p)
                        continue
            p = evaluate(x)
            if isnan(p[1]):
                step *= 0.3
            else:
                (previous, current) = (current, p)
                step *= 4
            ordered = sorted(points.items())
            for i in range(len(ordered) - 3):
                values = [ordered[j][1] for j in range(i, i + 4)]
                if values[0] > values[1] and values[2] < values[3] or (values[0] < values[1] and values[2] > values[3]):
                    (ok, current) = dip_search(ordered[i], ordered[i + 3], current)
                    if not ok:
                        return (float('nan'), current)
                    break
            if bracket():
                continue
            if good(current):
                return (current[0], current)
            if verge(current[0]) or previous is None or abs(current[0] - previous[0]) <= x_tolerance:
                return (float('nan'), current)
    initial_step = step
    (value, current) = travel(current, previous, step)
    if isnan(value):
        (value, current) = travel(far, previous, -initial_step)
    if not isnan(value) and value != last_x:
        evaluate(value)
        interval = bracket()
        candidates = sorted(points.items())
        if interval:
            candidates = [p for p in candidates if p[0] <= interval[1][0]]
        return min(candidates, key=lambda p: abs(p[1]))[0]
    return value

def _periodic_wrap(value: float, minimum: float, maximum: float) -> float:
    spread = maximum - minimum
    if spread <= 0:
        raise ValueError('周期区间长度必须为正')
    if value > maximum:
        return minimum + (value - maximum) % spread
    if value < minimum:
        return maximum - (minimum - value) % spread
    return value

def _clip_line(start: tuple[float, float], end: tuple[float, float], minimum_roll: float, maximum_roll: float) -> tuple[tuple[float, float], tuple[float, float]] | None:
    (x0, q0) = start
    (x1, q1) = end
    if q0 > maximum_roll and q1 > maximum_roll or (q0 < minimum_roll and q1 < minimum_roll):
        return None
    if q0 == q1:
        return (start, end)
    (ta, tb) = sorted(((minimum_roll - q0) / (q1 - q0), (maximum_roll - q0) / (q1 - q0)))
    (ta, tb) = (max(0.0, ta), min(1.0, tb))
    if tb < ta:
        return None
    return ((x0 + ta * (x1 - x0), q0 + ta * (q1 - q0)), (x0 + tb * (x1 - x0), q0 + tb * (q1 - q0)))

def _acos(x):
    with np.errstate(invalid='ignore'):
        return float(np.arccos(x))

def _asin(x):
    with np.errstate(invalid='ignore'):
        return float(np.arcsin(x))

def _involute(x):
    return tan(x) - x

class _TipPair:

    def __init__(self, mesh, diameter_a_m, diameter_b_m, angle):
        self.mesh = mesh
        (self.r, self.s) = (mesh.gear_b, mesh.gear_a)
        self.xr = abs(diameter_b_m) / 2
        self.fr = abs(diameter_a_m) / 2
        self.angle = angle
        c0 = sqrt(mesh.center_distance_m ** 2 - mesh.base_radius_sum_m ** 2)
        self.ir = self.r.diameter_from_roll(self.s.sign * c0 - mesh.mesh_sign * self.s.roll_from_diameter(diameter_a_m)) / 2
        self.jr = self.s.diameter_from_roll(self.r.sign * c0 - mesh.mesh_sign * self.r.roll_from_diameter(diameter_b_m)) / 2

    def yt(self, which):
        return _asin(abs(self.s.base_radius_m) / self.fr) if which else _asin(abs(self.r.base_radius_m) / self.xr)

    def kt(self, which):
        if which:
            return self.s.sign * (pi / 2 - self.yt(True) - self.mesh.working_pressure_angle_rad)
        return _asin(self.xr * sin(self.ct(False)) / self.jr)

    def ct(self, which):
        if which:
            return _asin(self.fr * sin(self.kt(True)) / self.ir)
        return self.r.sign * (pi / 2 - self.yt(False) - self.mesh.working_pressure_angle_rad)

    def dt(self):
        return self.s.sign * abs(self.s.signed_teeth / self.r.signed_teeth) * self.angle

    def bt(self, which):
        return self.fr * sin(self.kt(True) + abs(self.angle)) if which else self.xr * sin(self.ct(False) + self.dt())

    def zt(self, which):
        return self.fr * cos(self.kt(True) + abs(self.angle)) if which else self.xr * cos(self.ct(False) + self.dt())

    def at(self, which):
        z = self.zt(which)
        return self.s.sign * (self.mesh.center_distance_m - self.r.sign * z) if which else self.r.sign * (self.mesh.center_distance_m - self.s.sign * z)

    def lt(self, which):
        return atan(self.bt(which) / self.at(which))

    def ut(self, roll_angle):
        radius = abs(self.r.base_radius_m) * sqrt(1 + roll_angle * roll_angle)
        return _involute(_acos(self.r.base_radius_m / radius)) - _involute(_acos(self.r.base_radius_m / self.ir))

    def qt(self, roll_angle):
        radius = abs(self.r.base_radius_m) * sqrt(1 + roll_angle * roll_angle)
        angle = self.ct(True) + self.dt() + self.ut(roll_angle)
        horizontal = radius * sin(angle)
        projected = radius * cos(angle)
        other = self.mesh.center_distance_m - self.s.sign * projected
        return self.fr * self.fr - other * other - horizontal * horizontal

    def mt(self):
        root = self.mesh.geometric_root(self.qt, 0, self.r.roll_from_diameter(2 * self.xr) / abs(self.r.base_radius_m), self.mesh.active_roll_bounds('b')[0] / abs(self.r.base_radius_m))
        if isnan(root):
            return float('nan')
        d = 2 * abs(self.r.base_radius_m) * sqrt(1 + root * root)
        (lo, hi) = sorted((abs(self.r.root_form_diameter_m), abs(self.r.tip_form_diameter_m)))
        return self.r.roll_from_diameter(min(hi, max(lo, d)))

    def wt(self, returned_roll):
        if isnan(returned_roll):
            return 1.0
        value = self.ct(True) + self.dt() + self.ut(returned_roll / abs(self.r.base_radius_m)) - self.lt(True)
        return value * abs(self.r.base_radius_m)

    def vt(self):
        d = 2 * hypot(self.bt(False), self.at(False))
        (lo, hi) = sorted((abs(self.s.root_form_diameter_m), abs(self.s.tip_form_diameter_m)))
        return self.s.roll_from_diameter(min(hi, max(lo, d)))

    def ot(self, returned_roll):
        radius = self.s.diameter_from_roll(returned_roll) / 2
        delta_inv = _involute(_acos(abs(self.s.base_radius_m) / radius)) - _involute(_acos(abs(self.s.base_radius_m) / self.jr))
        value = self.kt(False) + self.s.sign * self.angle + self.s.sign * delta_inv - self.lt(False)
        return value * abs(self.s.base_radius_m)

def _extend_pinion_tip(mesh, angle):
    p = _TipPair(mesh, mesh.gear_a.tip_form_diameter_m, mesh.gear_b.tip_form_diameter_m, angle)
    roll = p.mt()
    gap = p.wt(roll)
    if isnan(roll):
        roll = mesh.gear_b.roll_from_diameter(mesh.gear_b.tip_form_diameter_m)
    return (roll, gap)

def _extend_wheel_tip(mesh, angle):
    p = _TipPair(mesh, mesh.gear_a.tip_form_diameter_m, mesh.gear_b.tip_form_diameter_m, angle)
    roll = p.vt()
    return (roll, p.ot(roll))

@dataclass
class _Coverage:
    root: float
    tip: float
    contact_ratio: float

def _coverage(mesh, extend_root=False, extend_tip=False):
    (root, tip) = mesh.active_roll_bounds('a')
    ratio = mesh.transverse_contact_ratio + mesh.overlap_contact_ratio
    extension = mesh.extension_roll_factor
    original_length = abs(tip - root)
    pitch = mesh.gear_a.transverse_base_pitch_m
    total = ratio
    state = _Coverage(root, tip, ratio)

    def apply(factor):
        nonlocal total
        d = factor * original_length / pitch
        state.contact_ratio += d * (int(extend_root) + int(extend_tip))
        total += 2 * d
        if extend_tip:
            state.tip += mesh.gear_a.sign * factor * original_length
        if extend_root:
            state.root -= mesh.gear_a.sign * factor * original_length
    apply(extension)
    if abs(total - round(total)) <= 0.05 and (extend_root or extend_tip) and (extension != 0):
        apply(0.25 * pitch / original_length)
    return state

def _line_count(mesh, roll, state):
    pitch = mesh.gear_a.transverse_base_pitch_m
    n_min = int(np.floor(state.contact_ratio))
    n_max = int(np.ceil(state.contact_ratio))
    (lo, hi) = sorted((state.root, state.root + mesh.gear_a.sign * pitch))
    wrapped = _periodic_wrap(roll, lo, hi)
    shifted = roll - pitch * (roll - wrapped) / pitch + mesh.gear_a.sign * n_min * pitch
    shifted -= mesh.effective_face_width_m * tan(mesh.gear_a.base_helix_angle_rad) * mesh.gear_a.sign
    return n_max if (shifted < state.tip) != mesh.gear_a.internal else n_min

def _build_lines(mesh, roll, n, partner_roll_mode, tooth_a, tooth_b):
    if partner_roll_mode != 'osep':
        raise NotImplementedError('延伸接触仅支持当前接触坐标映射，未实现其他映射分支')
    (a, b) = (mesh.gear_a, mesh.gear_b)
    face_factor = mesh.resolved_extended_face_width_factor
    pitch = a.transverse_base_pitch_m
    count0 = _line_count(mesh, roll, _coverage(mesh))
    count_all = _line_count(mesh, roll, _coverage(mesh, True, True))
    count_tip = _line_count(mesh, roll, _coverage(mesh, False, True))
    count_root = _line_count(mesh, roll, _coverage(mesh, True, False))
    offset_root = count_root - count0 if count_all > count_tip else 0
    offset_tip = count_tip - count0 if count_all > count_root else 0
    adjusted_roll = roll - offset_root * a.sign * pitch
    width = mesh.effective_face_width_m
    starts_at_width = (a.hand == 'left') == (mesh.flank == 'left')
    flag = not starts_at_width
    (x0, x1) = (width, 0.0) if starts_at_width else (0.0, width)
    dq = -width * tan(a.base_helix_angle_rad) * a.sign
    (root, tip) = mesh.active_roll_bounds('a')
    s_tip = mesh.active_roll_bounds('b')[1]
    source_c0 = sqrt(mesh.center_distance_m ** 2 - mesh.base_radius_sum_m ** 2)

    def b_to_a(q):
        return b.sign * source_c0 - mesh.mesh_sign * q

    def vsep(face, which_a, reference_face, source_virtual_start=None):
        if which_a:
            base_angle = tip / abs(a.base_radius_m)
            virtual = source_virtual_start if source_virtual_start is not None else tip + a.sign * abs(a.base_radius_m) * abs(tan(a.helix_angle_rad) / (a.reference_diameter_m / 2) * (face - reference_face))
            angle = virtual / abs(a.base_radius_m) - base_angle
            (q_b, ct) = _extend_pinion_tip(mesh, angle)
            return (tip, q_b, ct, True, False)
        base_angle = s_tip / abs(b.base_radius_m)
        virtual = source_virtual_start if source_virtual_start is not None else s_tip + b.sign * abs(b.base_radius_m) * abs(tan(b.helix_angle_rad) / (b.reference_diameter_m / 2) * (face - reference_face))
        angle = b_to_a(s_tip) / abs(a.base_radius_m) - b_to_a(virtual) / abs(a.base_radius_m)
        (q_a, ct) = _extend_wheel_tip(mesh, angle)
        return (q_a, s_tip, ct, False, True)
    lines = []
    for index in range(count_all):
        q0 = adjusted_roll + index * a.sign * pitch
        segment = _clip_line((x0, q0), (x1, q0 + dq), min(root, tip), max(root, tip))
        seed = lambda x: q0 + dq * (x - x0) / (x1 - x0)
        ext_a = ext_b = 0.0
        border_a = border_b = None
        ordinary_length = 0.0
        only_extension = False
        which_a = False
        if segment is not None:
            ((xa, qa), (xb, qb)) = segment
            ordinary_length = hypot(xb - xa, qb - qa)
            if abs(qa - tip) <= 1e-09:
                border_a = xa
                ext_a = min(xa if flag else width - xa, face_factor * width)
            if abs(mesh.partner_roll_from_osep(qb) - s_tip) <= 1e-09:
                border_b = xb
                ext_b = min(width - xb if flag else xb, face_factor * width)
            extent = abs(xb - xa) + ext_a + ext_b
            start = (min(xa, xb) if border_a is None else border_a - ext_a) if flag else max(xa, xb) if border_a is None else border_a + ext_a
        else:
            only_extension = True
            expanded = _coverage(mesh, True, True)
            full = _clip_line((x0, q0), (x1, q0 + dq), min(expanded.root, expanded.tip), max(expanded.root, expanded.tip))
            if full is None:
                continue
            if a.source_is_spur:
                which_a = full[0][1] < tip if a.internal else full[0][1] > tip
                reference_face = 0.0
                extent = width
            else:
                at_a = x0 + (tip - q0) * (x1 - x0) / dq
                which_a = not flag and at_a < 0.0 or (flag and at_a > width)
                reference_face = at_a if which_a else x0 + (b_to_a(s_tip) - q0) * (x1 - x0) / dq
                extent = min(face_factor, 1.0) * width
            start = (width - extent if flag else extent) if which_a else 0.0 if flag else width
        if extent <= 0:
            continue
        direction = 1 if flag else -1
        xe = tuple((start + direction * extent * i / n for i in range(n + 1)))
        xc = tuple((start + direction * extent * (i + 0.5) / n for i in range(n)))

        def sample(x):
            if only_extension:
                virtual = full[0][1] if a.source_is_spur and which_a else mesh.partner_roll_from_osep(full[0][1]) if a.source_is_spur else None
                return vsep(x, which_a, reference_face, virtual)
            if border_a is not None and (not flag and x > border_a or (flag and x < border_a)):
                return vsep(x, True, border_a)
            if min(xa, xb) - 1e-12 <= x <= max(xa, xb) + 1e-12:
                q_a = seed(x)
                return (q_a, mesh.partner_roll_from_osep(q_a), 0.0, False, False)
            return vsep(x, False, border_b)
        centre = tuple((sample(x) for x in xc))
        endpoint = tuple((sample(x) for x in xe))
        ra = tuple((p[0] for p in centre))
        rb = tuple((p[1] for p in centre))
        rae = tuple((p[0] for p in endpoint))
        rbe = tuple((p[1] for p in endpoint))
        lengths = tuple((0.5 * (hypot(xe[i + 1] - xe[i], rae[i + 1] - rae[i]) + hypot(xe[i + 1] - xe[i], rbe[i + 1] - rbe[i])) for i in range(n)))
        eligible = tuple((not a.below_root_form(q_a) and (not b.below_root_form(q_b)) for (q_a, q_b) in zip(ra, rb)))
        flank_sign = -1 if mesh.flank == 'left' else 1
        warnings = ()
        lines.append(_Line(index, xc, xe, ra, rb, rae, rbe, lengths, eligible, tuple((mesh.local_face_coordinate('a', x) for x in xc)), tuple((mesh.local_face_coordinate('b', x) for x in xc)), ordinary_length, tooth_a + a.sign * flank_sign * (index - offset_root), tooth_b + b.sign * flank_sign * (count_all - 1 - index - offset_tip), partner_roll_mode, warnings, tuple((p[2] for p in centre)), tuple((p[3] for p in centre)), tuple((p[4] for p in centre)), tuple((p[3] for p in endpoint)), tuple((p[4] for p in endpoint))))
    return tuple(lines)


@dataclass(frozen=True)
class ExtendedContactSettings:
    """延伸候选范围；无量纲范围因子须与所用齿顶接触模型一致。

    face_extension_factor限制超出普通接触段的轴向候选长度/重叠齿宽；
    斜齿必须显式提供，直齿省略时取1。roll_extension_factor限制候选
    展长范围，不按求解结果调整。tip_radius_*_m为齿顶延伸段的曲率
    参照半径，单位m，省略时使用相应齿顶展长。

    root_tolerance_policy默认compatible保留固定绝对残差退出策略；
    dimensionless_refined是显式工程策略，对长度平方残差归一化后加严
    求根容差，使几何相似缩放不受单位尺度的绝对阈值影响。
    """
    face_extension_factor: float | None = None
    roll_extension_factor: float = 0.4
    tip_radius_pinion_m: float | None = None
    tip_radius_wheel_m: float | None = None
    root_tolerance_policy: str = "compatible"

    def __post_init__(self):
        if self.root_tolerance_policy not in ("compatible", "dimensionless_refined"):
            raise ValueError("root_tolerance_policy须为compatible或dimensionless_refined")
        for name in ("face_extension_factor", "tip_radius_pinion_m", "tip_radius_wheel_m"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not isfinite(value) or value <= 0):
                raise ValueError(f"{name}必须为正有限数")
        if (isinstance(self.roll_extension_factor, bool) or not isfinite(self.roll_extension_factor)
                or self.roll_extension_factor < 0):
            raise ValueError("延伸展长因子必须为非负有限数")


@dataclass(frozen=True)
class ExtendedEndpoints:
    """各接触线N+1端点，line_offsets逐线定位；长度m、曲率1/m。"""
    line_offsets: np.ndarray
    face_coordinate_m: np.ndarray
    roll_pinion_m: np.ndarray
    roll_wheel_m: np.ndarray
    tooth_pair_index: np.ndarray
    curvature_per_m: np.ndarray
    radius_pinion_m: np.ndarray
    radius_wheel_m: np.ndarray
    extended_pinion: np.ndarray
    extended_wheel: np.ndarray


@dataclass(frozen=True)
class ExtendedContactResult:
    """完整中心网格及端点，附加法向几何间隙应加入初始正间隙。

    eligible只判根形可用范围，不代表接触活动集合。求解器须对不可用
    单元回填零力，不得删除索引或将未求解域解释为真实零压区域。
    """
    grid: ContactGrid
    line_offsets: np.ndarray
    gap_correction_m: np.ndarray
    eligible: np.ndarray
    centre_extended_pinion: np.ndarray
    centre_extended_wheel: np.ndarray
    endpoints: ExtendedEndpoints


def _normal_radii(pair, qa, qb, a_tip, b_tip, settings):
    """延伸曲线几何与Hertz曲率分别定义；保留内啮合的曲率差号。"""
    cosine = cos(pair.base_helix_angle_rad)
    effective_a, effective_b = np.asarray(qa, float).copy(), np.asarray(qb, float).copy()
    tip_a = pair.pinion.involute_roll_bounds_m[1]
    tip_b = (pair.wheel.involute_roll_bounds_m[0] if pair.wheel.internal
             else pair.wheel.involute_roll_bounds_m[1])
    if settings.tip_radius_pinion_m is not None:
        mask = np.asarray(a_tip) & (np.abs(effective_a-tip_a) < 1e-6)
        effective_a[mask] = settings.tip_radius_pinion_m
    if settings.tip_radius_wheel_m is not None:
        mask = np.asarray(b_tip) & (np.abs(effective_b-tip_b) < 1e-6)
        effective_b[mask] = settings.tip_radius_wheel_m
    radius_a, radius_b = effective_a/cosine, effective_b/cosine
    with np.errstate(divide="ignore", invalid="ignore"):
        curvature = 1/radius_a + (-1 if pair.wheel.internal else 1)/radius_b
    return radius_a, radius_b, curvature


def extended_contact_grid(pair: GearPair, phase: float, strips: int = 21, *,
                          settings: ExtendedContactSettings) -> ExtendedContactResult:
    """由当前齿轮参数生成普通与延伸接触线，不读取对标响应数组。

    沿用contact_grid的基节相位、居中齿宽坐标及齿对index。公开齿轮
    的螺旋角带符号；内部只选择一个规范工作面，不额外反转材料坐标。
    几何根使用受限区间搜索，失败分支保留有限的非接触候选间隙；任何
    非有限几何均明确拒绝，不用微小曲率代替。候选参数域不等于精度界。
    """
    if not isinstance(settings, ExtendedContactSettings):
        raise TypeError("settings须为ExtendedContactSettings")
    if isinstance(strips, bool) or not isinstance(strips, (int, np.integer)) or strips < 2:
        raise ValueError("每条接触线至少需要2个条带")
    if not isfinite(phase):
        raise ValueError("接触相位必须为有限数")
    if abs(pair.pinion.helix_angle_rad) > 1e-8 and settings.face_extension_factor is None:
        raise ValueError("斜齿延伸接触须显式给出face_extension_factor")
    mesh = _Mesh(pair, settings)
    low, _ = pair.contact_roll_bounds_m
    pitch = pair.base_pitch_m
    slope = tan(pair.base_helix_angle_rad)
    origin = low + ((phase % 1.0)+0.5)*pitch
    left, right = pair.face_bounds_m
    start_z = left if slope < 0 else right
    roll = _periodic_wrap(origin+slope*start_z, low, low+pitch)
    lines = _build_lines(mesh, roll, int(strips), 'osep', 0, 0)
    if not lines:
        raise ValueError("当前相位没有可用延伸接触候选")
    ordinary_count = _line_count(mesh, roll, _coverage(mesh))
    all_count = _line_count(mesh, roll, _coverage(mesh, True, True))
    tip_count = _line_count(mesh, roll, _coverage(mesh, False, True))
    root_count = _line_count(mesh, roll, _coverage(mesh, True, False))
    root_offset = root_count-ordinary_count if all_count > tip_count else 0
    effective_origin = roll-root_offset*pitch-slope*start_z
    centre, endpoint = {}, {}
    centre_keys = ("qa", "qb", "face", "width", "id", "gap", "eligible", "a_tip", "b_tip")
    endpoint_keys = ("qa", "qb", "face", "id", "a_tip", "b_tip")
    for key in centre_keys:
        centre[key] = []
    for key in endpoint_keys:
        endpoint[key] = []
    offsets, end_offsets = [0], [0]
    for line in lines:
        identity = int(round((effective_origin+line.source_line_index*pitch-origin)/pitch))
        x = np.asarray(line.face_centres_m)
        xe = np.asarray(line.face_endpoints_m)
        order = np.argsort(x)
        end_order = np.argsort(xe)
        centre["qa"].extend(np.asarray(line.gear_a_roll_centres_m)[order])
        centre["qb"].extend(np.asarray(line.gear_b_roll_centres_m)[order])
        centre["face"].extend((x+left)[order])
        centre["width"].extend(np.asarray(line.widths_m)[order])
        centre["id"].extend([identity]*len(x))
        centre["gap"].extend(np.asarray(line.correction_transverse_m)[order]*cos(pair.base_helix_angle_rad))
        centre["eligible"].extend(np.asarray(line.eligible)[order])
        centre["a_tip"].extend(np.asarray(line.gear_a_extended_centres)[order])
        centre["b_tip"].extend(np.asarray(line.gear_b_extended_centres)[order])
        endpoint["qa"].extend(np.asarray(line.gear_a_roll_endpoints_m)[end_order])
        endpoint["qb"].extend(np.asarray(line.gear_b_roll_endpoints_m)[end_order])
        endpoint["face"].extend((xe+left)[end_order])
        endpoint["id"].extend([identity]*len(xe))
        endpoint["a_tip"].extend(np.asarray(line.gear_a_extended_endpoints)[end_order])
        endpoint["b_tip"].extend(np.asarray(line.gear_b_extended_endpoints)[end_order])
        offsets.append(offsets[-1]+len(x))
        end_offsets.append(end_offsets[-1]+len(xe))
    for values in (centre, endpoint):
        for key in values:
            dtype = bool if key in ("eligible", "a_tip", "b_tip") else int if key == "id" else float
            values[key] = np.asarray(values[key], dtype=dtype)
        for key in ("qa", "qb", "face"):
            if not np.all(np.isfinite(values[key])):
                raise ValueError("延伸接触几何根未得到有限坐标；当前参数分支不适用")
    if not np.all(np.isfinite(centre["gap"])) or np.any(centre["width"] <= 0):
        raise ValueError("延伸接触修正或条带长度无效")
    qa, qb, face = centre["qa"], centre["qb"], centre["face"]
    pl, ph = pair.pinion.involute_roll_bounds_m
    wl, wh = pair.wheel.involute_roll_bounds_m
    fp = (qa-pl)/(ph-pl)
    fw = ((wh-qb) if pair.wheel.internal else (qb-wl))/(wh-wl)
    lp = (face+pair.pinion.face_width_m/2)/pair.pinion.face_width_m
    lw = (face-pair.axial_offset_m+pair.wheel.face_width_m/2)/pair.wheel.face_width_m
    for factor in (fp, fw, lp, lw):
        if np.any(factor < -1e-9) or np.any(factor > 1+1e-9):
            raise ValueError("延伸接触点不在完整形径／齿宽范围，不能静默裁剪修形输入")
    ra, rb, curvature = _normal_radii(pair, qa, qb, centre["a_tip"], centre["b_tip"], settings)
    grid = ContactGrid(qa, qb, face, centre["width"], centre["id"],
        np.clip(fp, 0, 1), np.clip(fw, 0, 1), np.clip(lp, 0, 1), np.clip(lw, 0, 1), curvature, ra, rb)
    era, erb, ec = _normal_radii(pair, endpoint["qa"], endpoint["qb"], endpoint["a_tip"], endpoint["b_tip"], settings)
    ends = ExtendedEndpoints(np.asarray(end_offsets, int), endpoint["face"], endpoint["qa"], endpoint["qb"],
        endpoint["id"], ec, era, erb, endpoint["a_tip"], endpoint["b_tip"])
    return ExtendedContactResult(grid, np.asarray(offsets, int), centre["gap"], centre["eligible"],
        centre["a_tip"], centre["b_tip"], ends)
