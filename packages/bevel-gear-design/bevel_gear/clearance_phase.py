"""正间隙目标的刚性相位求解；角度为rad，长度为m，无载荷或材料模型。

距离回调负责真实几何、相邻齿及初始安装。返回角是安装后的从动角增量，
不能直接作为整轮绝对转角。本算法没有有符号间隙或实体侵入判据。
"""
from dataclasses import dataclass
import math
import struct


@dataclass(frozen=True)
class ClearancePhaseSettings:
    """显式保留启发修正、五点细化、扩窗和二分的参数与历史分支。"""
    target_gap_m: float = 3e-6
    residual_tolerance_m: float = 5e-7
    adaptive_residual_m: float = 2e-7
    overshoot_gap_m: float = 5e-6
    probe_angle_rad: float = math.radians(.0001)
    coarse_angle_tolerance_rad: float = math.radians(.001)
    refinement_radius_rad: float = math.radians(.002)
    adaptive_refinement_radius_rad: float = math.radians(.001)
    refinement_iterations: int = 1
    max_step_iterations: int = 100
    global_initial_radius_rad: float = math.radians(.01)
    global_expansions: int = 6
    global_intervals: int = 16
    bisection_iterations: int = 24
    bisection_residual_m: float = 1e-10
    max_evaluations: int = 1000

    def __post_init__(self):
        integer_fields = {"refinement_iterations", "max_step_iterations", "global_expansions",
                          "global_intervals", "bisection_iterations", "max_evaluations"}
        for name, value in self.__dict__.items():
            if name in integer_fields:
                if isinstance(value, bool) or not isinstance(value, int) or value < (0 if name == "max_step_iterations" else 1):
                    raise ValueError(f"{name}须为有效整数限次")
            elif isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name}须为正有限实数")


@dataclass(frozen=True)
class ClearancePhaseResult:
    """失败时角度和距离为最后一次有效评价；初始安装角始终由调用者保存。"""
    driving_angle_rad: float
    follower_increment_angle_rad: float | None
    nominal_follower_increment_angle_rad: float
    correction_from_ratio_rad: float | None
    distance_m: float | None
    residual_m: float | None
    evaluation_count: int
    step_iterations: int
    adaptive_refinements: int
    global_windows: int
    bisection_iterations: int
    step_iteration_limit_reached: bool
    converged: bool
    status: str
    failure: str | None
    penetration_checked: bool = False
    initial_installation_applied_by_solver: bool = False


class _Stopped(Exception):
    def __init__(self, status, message):
        self.status = status
        super().__init__(message)


def solve_clearance_phase(distance_callback, *, driving_angle_rad, driving_teeth, follower_teeth,
                          initial_correction_rad=0., flank="left", contact_radius_m=None,
                          contact_radius_callback=None, settings=None):
    """按最小距离正目标修正从动增量角，所有回调输入角均为rad。

    distance_callback(phi_driver, phi_follower)返回m；调用者负责当前及相邻齿的距离聚合。
    接触半径用于启发角步长，须在正常数contact_radius_m与逐点评价的
    contact_radius_callback之间二选一。若半径来自同一次最近点查询，可在适配器缓存。
    正半径近似值可以运行，但不能宣称与变半径原几何策略严格等价。
    设置的默认3微米目标是正距离策略，不是零间隙接触条件。
    """
    if settings is None:
        settings = ClearancePhaseSettings()
    elif isinstance(settings, dict):
        settings = ClearancePhaseSettings(**settings)
    if not isinstance(settings, ClearancePhaseSettings):
        raise TypeError("settings须为ClearancePhaseSettings或具名字典")
    if not callable(distance_callback):
        raise TypeError("distance_callback须可调用")
    for name, value in (("driving_angle_rad", driving_angle_rad), ("initial_correction_rad", initial_correction_rad)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{name}须为有限实数")
    for name, value in (("driving_teeth", driving_teeth), ("follower_teeth", follower_teeth)):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"{name}须为正整数")
    if flank not in ("left", "right"):
        raise ValueError("flank须为left或right")
    if (contact_radius_m is None) == (contact_radius_callback is None):
        raise ValueError("须且只能提供常数接触半径或接触半径回调之一")
    if contact_radius_callback is not None and not callable(contact_radius_callback):
        raise TypeError("contact_radius_callback须可调用")
    if contact_radius_m is not None and (isinstance(contact_radius_m, bool) or
            not isinstance(contact_radius_m, (int, float)) or not math.isfinite(contact_radius_m) or contact_radius_m <= 0):
        raise ValueError("contact_radius_m须为正有限实数")
    s = settings
    driver = math.degrees(driving_angle_rad)
    nominal = -driver * driving_teeth / follower_teeth + math.degrees(initial_correction_rad)
    counts = dict(evaluation_count=0, step_iterations=0, adaptive_refinements=0,
                  global_windows=0, bisection_iterations=0)
    last_angle = last_distance = None
    reached_limit = False

    def measure(angle):
        nonlocal last_angle, last_distance
        if counts["evaluation_count"] >= s.max_evaluations:
            raise _Stopped("evaluation_limit", "距离评价达到显式上限")
        counts["evaluation_count"] += 1
        try:
            value = distance_callback(driving_angle_rad, math.radians(angle))
        except Exception as error:
            raise _Stopped("callback_error", "距离回调失败：" + str(error)) from error
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise _Stopped("invalid_measurement", "距离回调须返回非负有限长度m")
        last_angle, last_distance = angle, float(value)
        return float(value)

    def radius(angle):
        if contact_radius_callback is None:
            return float(contact_radius_m)
        try:
            value = contact_radius_callback(driving_angle_rad, math.radians(angle))
        except Exception as error:
            raise _Stopped("callback_error", "接触半径回调失败：" + str(error)) from error
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise _Stopped("invalid_measurement", "接触半径回调须返回正有限长度m")
        return float(value)

    def refine(angle, span):
        # 只在本次五点评价器内按double位模式缓存，保留重复评价语义。
        cache = {}
        def evaluate(value):
            key = struct.pack("!d", value)
            if key not in cache:
                cache[key] = measure(value)
            return cache[key]
        best = angle
        residual = abs(evaluate(best) - s.target_gap_m)
        lower, upper = angle - span, angle + span
        for _ in range(s.refinement_iterations):
            step = (upper - lower) * .25
            for candidate in (lower, lower + step, (lower + upper) * .5, upper - step, upper):
                error = abs(evaluate(candidate) - s.target_gap_m)
                if error < residual:
                    best, residual = candidate, error
            lower, upper = best - step, best + step
        return best

    def global_refine(preferred):
        best = preferred
        residual = abs(measure(preferred) - s.target_gap_m)
        span = math.degrees(s.global_initial_radius_rad)
        for _ in range(s.global_expansions):
            counts["global_windows"] += 1
            lower = nominal - span
            step = 2 * span / s.global_intervals
            previous = lower
            previous_value = measure(previous) - s.target_gap_m
            if abs(previous_value) < residual:
                best, residual = previous, abs(previous_value)
            bracket = None
            bracket_distance = math.inf
            for index in range(1, s.global_intervals + 1):
                candidate = lower + index * step
                value = measure(candidate) - s.target_gap_m
                if abs(value) < residual:
                    best, residual = candidate, abs(value)
                if previous_value == 0:
                    return previous
                if value == 0:
                    return candidate
                if (previous_value > 0) != (value > 0):
                    proximity = abs(.5 * (previous + candidate) - preferred)
                    if proximity < bracket_distance:
                        bracket_distance = proximity
                        bracket = (previous, previous_value, candidate, value)
                previous, previous_value = candidate, value
            if bracket is not None:
                lo, flo, hi, fhi = bracket
                best = lo if abs(flo) <= abs(fhi) else hi
                residual = min(abs(flo), abs(fhi))
                for _ in range(s.bisection_iterations):
                    counts["bisection_iterations"] += 1
                    midpoint = .5 * (lo + hi)
                    value = measure(midpoint) - s.target_gap_m
                    if abs(value) < residual:
                        best, residual = midpoint, abs(value)
                    if residual <= s.bisection_residual_m:
                        break
                    if (flo > 0) == (value > 0):
                        lo, flo = midpoint, value
                    else:
                        hi, fhi = midpoint, value
                return best
            span *= 2
        return best

    def result(angle, distance, status, failure=None):
        error = None if distance is None else abs(distance - s.target_gap_m)
        follower = None if angle is None else math.radians(angle)
        return ClearancePhaseResult(driving_angle_rad, follower, math.radians(nominal),
            None if follower is None else follower + driving_angle_rad * driving_teeth / follower_teeth,
            distance, error, **counts, step_iteration_limit_reached=reached_limit,
            converged=status == "converged", status=status, failure=failure)

    try:
        angle = nominal
        distance = measure(angle)
        contact_radius = radius(angle)
        probe = measure(angle - math.degrees(s.probe_angle_rad))
        tolerance = min(contact_radius * s.coarse_angle_tolerance_rad, s.residual_tolerance_m)
        if (distance < probe) if flank == "left" else (distance > probe):
            step = (distance + s.overshoot_gap_m) * 2 / contact_radius * 180 / math.pi
            angle += step if flank == "left" else -step
        else:
            step = (distance - s.target_gap_m) / contact_radius * 180 / math.pi
            angle += -step if flank == "left" else step
        distance = measure(angle)
        contact_radius = radius(angle)
        while abs(distance - s.target_gap_m) > tolerance and counts["step_iterations"] < s.max_step_iterations:
            step = (distance - s.target_gap_m) / contact_radius * 180 / math.pi
            angle += -step if flank == "left" else step
            distance = measure(angle)
            contact_radius = radius(angle)
            probe = measure(angle - math.degrees(s.probe_angle_rad))
            if (distance < probe) if flank == "left" else (distance > probe):
                step = (distance + s.overshoot_gap_m) * 2 / contact_radius * 180 / math.pi
                angle += step if flank == "left" else -step
            counts["step_iterations"] += 1
        reached_limit = counts["step_iterations"] >= s.max_step_iterations and abs(distance - s.target_gap_m) > tolerance
        angle = refine(angle, math.degrees(s.refinement_radius_rad))
        distance = measure(angle)
        if abs(distance - s.target_gap_m) > s.adaptive_residual_m:
            counts["adaptive_refinements"] += 1
            angle = refine(angle, math.degrees(s.adaptive_refinement_radius_rad))
            distance = measure(angle)
        if abs(distance - s.target_gap_m) > s.residual_tolerance_m:
            counts["adaptive_refinements"] += 1
            angle = global_refine(angle)
            distance = measure(angle)
        if abs(distance - s.target_gap_m) <= s.residual_tolerance_m:
            return result(angle, distance, "converged")
        return result(angle, distance, "residual_not_met", "有限搜索返回最优候选，但未满足正间隙目标残差")
    except _Stopped as error:
        return result(last_angle, last_distance, error.status, str(error))
