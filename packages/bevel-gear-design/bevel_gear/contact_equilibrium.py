"""给定间隙与线性柔度的互补接触平衡，按明确物理方程新写。

本模块不生成齿面、齿体柔度或Hertz柔度，也不宣称构成完整锥齿轮LTCA。
f单位N，g及剩余间隙单位m，C单位m/N；力矩模式使用力臂m和转角rad。
求解凸二次规划：min(0.5*fᵀ*C*f+gᵀ*f)，f>=0，ellᵀ*f=T。
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ContactEquilibriumResult:
    """剩余间隙w=g+C*f−ell*theta；接触力是离散节点合力，不是压力。"""
    forces_n: list
    remaining_gaps_m: list
    elastic_deflections_m: list
    theta_rad: float | None
    normal_approach_m: float | None
    target_load: float
    resultant_load: float
    load_unit: str
    balance_residual: float
    relative_balance_residual: float
    complementarity_residual_nm: float
    contact_stationarity_residual_m: float
    minimum_force_n: float
    minimum_remaining_gap_m: float
    strain_energy_j: float
    active_indices: list
    iterations: int
    linear_solve_count: int
    converged: bool
    status: str
    failure: str | None
    diagnostics: dict


def _real(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name}须为有限实数")
    return float(value)


def _cholesky(matrix, minimum_pivot):
    n = len(matrix)
    # 大接触网格在已安装数值依赖时使用同一Cholesky方程的专业线性代数实现。
    # 未安装时保留原标准库算法；不增加对角正则、不更改正定或主元接受门槛。
    if n >= 96:
        try:
            import numpy as np
            from scipy.linalg import cholesky
        except ImportError:
            pass
        else:
            try:
                lower = cholesky(np.asarray(matrix, dtype=float), lower=True, check_finite=False)
            except np.linalg.LinAlgError as error:
                raise ValueError("柔度矩阵不是正定矩阵，或归一化Cholesky主元低于显式病态门槛") from error
            smallest = float(np.min(np.diag(lower) ** 2))
            if not math.isfinite(smallest) or smallest <= minimum_pivot:
                raise ValueError("柔度矩阵不是正定矩阵，或归一化Cholesky主元低于显式病态门槛")
            return lower, smallest
    lower = [[0.] * n for _ in range(n)]
    smallest = math.inf
    for i in range(n):
        for j in range(i + 1):
            value = matrix[i][j] - math.fsum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if not math.isfinite(value) or value <= minimum_pivot:
                    raise ValueError("柔度矩阵不是正定矩阵，或归一化Cholesky主元低于显式病态门槛")
                smallest = min(smallest, value)
                lower[i][j] = math.sqrt(value)
            else:
                lower[i][j] = value / lower[j][j]
    return lower, smallest


def _solve_cholesky(lower, rhs):
    n = len(rhs)
    if n >= 96 and not isinstance(lower, list):
        from scipy.linalg import solve_triangular
        y = solve_triangular(lower, rhs, lower=True, check_finite=False)
        return solve_triangular(lower.T, y, lower=False, check_finite=False).tolist()
    y = [0.] * n
    for i in range(n):
        y[i] = (rhs[i] - math.fsum(lower[i][j] * y[j] for j in range(i))) / lower[i][i]
    x = [0.] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - math.fsum(lower[j][i] * x[j] for j in range(i + 1, n))) / lower[i][i]
    return x


def solve_contact_equilibrium(gaps_m, compliance_m_per_n, *, lever_arms_m, torque_nm,
                              allow_initial_penetration=False, max_iterations=1000,
                              force_tolerance_n=1e-8, gap_tolerance_m=1e-10,
                              balance_tolerance_nm=1e-8, relative_tolerance=1e-9,
                              symmetry_tolerance=1e-10, minimum_normalized_pivot=1e-12):
    """力矩约束下解非负节点力与相对法向闭合转角；仅用标准库。

    ell为法向位移对该相对转角的导数，单位m/rad（rad按SI无量纲），
    必须与T采用同一轴和符号。允许正/负力臂，但接触力始终非负。
    负初始gap默认拒绝；显式允许时保留原值，theta可为负以表示初装打开修正。
    C须对称正定；舍入级不对称在门槛内取均值并报告，不修复非正定矩阵。
    """
    return _solve(gaps_m, compliance_m_per_n, lever_arms_m, torque_nm, "torque",
        allow_initial_penetration, max_iterations, force_tolerance_n, gap_tolerance_m,
        balance_tolerance_nm, relative_tolerance, symmetry_tolerance, minimum_normalized_pivot)


def solve_normal_contact_equilibrium(gaps_m, compliance_m_per_n, *, total_normal_force_n,
                                     allow_initial_penetration=False, max_iterations=1000,
                                     force_tolerance_n=1e-8, gap_tolerance_m=1e-10,
                                     balance_tolerance_n=1e-8, relative_tolerance=1e-9,
                                     symmetry_tolerance=1e-10, minimum_normalized_pivot=1e-12):
    """总法向力模式：ell=1，闭合量单位m，平衡残差单位N，不将其解释为转角。"""
    gaps = list(gaps_m)
    return _solve(gaps, compliance_m_per_n, [1.] * len(gaps), total_normal_force_n, "normal_force",
        allow_initial_penetration, max_iterations, force_tolerance_n, gap_tolerance_m,
        balance_tolerance_n, relative_tolerance, symmetry_tolerance, minimum_normalized_pivot)


def _solve(gaps, compliance, levers, target, mode, allow_penetration, max_iterations,
           force_tolerance, gap_tolerance, balance_tolerance, relative_tolerance,
           symmetry_tolerance, minimum_pivot):
    gaps = [_real(v, "gap") for v in gaps]
    levers = [_real(v, "lever") for v in levers]
    target = _real(target, "目标载荷")
    if not gaps or len(levers) != len(gaps):
        raise ValueError("间隙与力臂须为相同非零长度向量")
    n = len(gaps)
    matrix = [[_real(v, "C") for v in row] for row in compliance]
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("柔度矩阵须为与节点数相同的方阵")
    if not isinstance(allow_penetration, bool):
        raise ValueError("allow_initial_penetration须为布尔值")
    if min(gaps) < 0 and not allow_penetration:
        raise ValueError("检测到负初始间隙；须显式允许初始穿透，不能静默夹零")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations须为正整数")
    for name, value in (("force_tolerance",force_tolerance),("gap_tolerance",gap_tolerance),
                        ("balance_tolerance",balance_tolerance),("relative_tolerance",relative_tolerance),
                        ("symmetry_tolerance",symmetry_tolerance),("minimum_normalized_pivot",minimum_pivot)):
        if _real(value, name) <= 0:
            raise ValueError(f"{name}须为正值")
    compliance_scale = max(matrix[i][i] for i in range(n))
    if compliance_scale <= 0 or any(matrix[i][i] <= 0 for i in range(n)):
        raise ValueError("柔度矩阵对角项须为正")
    symmetry_error = max(abs(matrix[i][j] - matrix[j][i]) for i in range(n) for j in range(n))
    if symmetry_error > symmetry_tolerance * compliance_scale:
        raise ValueError("柔度矩阵不对称，未修改为另一个物理模型")
    matrix = [[.5 * (matrix[i][j] + matrix[j][i]) for j in range(n)] for i in range(n)]
    normalized = [[v / compliance_scale for v in row] for row in matrix]
    full_factor, smallest_pivot = _cholesky(normalized, minimum_pivot)
    lever_scale = max(abs(v) for v in levers)
    if lever_scale == 0:
        raise ValueError("至少一个广义力臂须非零")
    a = [v / lever_scale for v in levers]
    target_force = target / lever_scale
    force_scale = max(abs(target_force), 1.)
    displacement_scale = compliance_scale * force_scale
    if not math.isfinite(target_force) or not math.isfinite(force_scale) or not math.isfinite(displacement_scale) or displacement_scale <= 0:
        raise ValueError("载荷、力臂或柔度尺度导致数值溢出/下溢，未进入接触迭代")
    h = [v / displacement_scale for v in gaps]
    if not all(math.isfinite(value) for value in h):
        raise ValueError("间隙相对柔度/载荷尺度超出有限数范围")
    b = target_force / force_scale
    feasible = [i for i in range(n) if a[i] * b > 0] if b != 0 else [i for i in range(n) if a[i] != 0]
    if not feasible:
        raise ValueError("目标载荷不在非负接触力可形成的广义载荷锥内")
    initial = max(feasible, key=lambda i: abs(a[i]))
    x = [0.] * n
    x[initial] = b / a[initial]
    free = {initial}
    multiplier = 0.
    linear_solves = 0
    factors = {tuple(range(n)): full_factor}
    force_limit = force_tolerance + relative_tolerance * force_scale
    gap_limit = gap_tolerance + relative_tolerance * max(max(abs(v) for v in gaps), displacement_scale)
    balance_limit = balance_tolerance + relative_tolerance * max(abs(target), lever_scale * force_scale)
    complementarity_limit = force_scale * gap_limit
    diagnostics = {"mode":mode,"compliance_scale_m_per_n":compliance_scale,"force_scale_n":force_scale,
        "maximum_input_symmetry_error_m_per_n":symmetry_error,"minimum_normalized_cholesky_pivot":smallest_pivot,
        "initial_penetration_allowed":allow_penetration,"minimum_input_gap_m":min(gaps),
        "force_acceptance_n":force_limit,"gap_acceptance_m":gap_limit,
        "balance_acceptance":balance_limit,"complementarity_acceptance_nm":complementarity_limit,
        "model":"给定线性柔度的离散互补接触；未生成齿面或柔度，未计算压力"}

    def finish(iteration, status, failure=None):
        forces = [v * force_scale for v in x]
        elastic = [math.fsum(matrix[i][j] * forces[j] for j in range(n)) for i in range(n)]
        approach = multiplier * displacement_scale
        remaining = [gaps[i] + elastic[i] - a[i] * approach for i in range(n)]
        resultant = math.fsum(levers[i] * forces[i] for i in range(n))
        balance = resultant - target
        complementarity = max(abs(forces[i] * remaining[i]) for i in range(n))
        active = [i for i in range(n) if forces[i] > force_limit]
        stationarity = max((abs(remaining[i]) for i in active), default=0.)
        certified = (min(forces) >= -force_limit and min(remaining) >= -gap_limit
                     and abs(balance) <= balance_limit and complementarity <= complementarity_limit
                     and stationarity <= gap_limit)
        converged = status == "converged" and certified
        if status == "converged" and not certified:
            status, failure = "residual_not_met", "活动集候选未同时满足力平衡与互补残差"
        return ContactEquilibriumResult(forces, remaining, elastic,
            approach / lever_scale if mode == "torque" else None,
            approach if mode == "normal_force" else None, target, resultant,
            "N·m" if mode == "torque" else "N", balance,
            abs(balance) / max(abs(target), balance_limit), complementarity, stationarity,
            min(forces), min(remaining), .5 * math.fsum(forces[i]*elastic[i] for i in range(n)),
            active, iteration, linear_solves, converged, status, failure, diagnostics)

    if target == 0 and min(gaps) >= 0:
        return finish(0, "converged")
    for iteration in range(1, max_iterations + 1):
        indices = tuple(sorted(free))
        if indices not in factors:
            factors[indices] = _cholesky([[normalized[i][j] for j in indices] for i in indices], minimum_pivot)[0]
        factor = factors[indices]
        gradient = [h[i] + math.fsum(normalized[i][j] * x[j] for j in range(n)) for i in range(n)]
        u = _solve_cholesky(factor, [gradient[i] for i in indices])
        v = _solve_cholesky(factor, [a[i] for i in indices])
        linear_solves += 2
        denominator = math.fsum(a[i] * value for i, value in zip(indices, v))
        if denominator <= 0 or not math.isfinite(denominator):
            return finish(iteration, "singular_active_constraint", "当前接触活动集没有有效广义力臂")
        multiplier = math.fsum(a[i] * value for i, value in zip(indices, u)) / denominator
        direction = [0.] * n
        for index, ui, vi in zip(indices, u, v):
            direction[index] = -ui + vi * multiplier
        if max(abs(value) for value in direction) <= force_limit / force_scale:
            reduced = [gradient[i] - a[i] * multiplier for i in range(n)]
            candidates = [i for i in range(n) if i not in free and reduced[i] < -gap_limit / displacement_scale]
            if not candidates:
                return finish(iteration, "converged")
            free.add(min(candidates, key=lambda i: (reduced[i], i)))
            continue
        alpha = 1.
        blocker = None
        for i in indices:
            if direction[i] < -force_limit / force_scale:
                candidate = -x[i] / direction[i]
                if candidate < alpha:
                    alpha, blocker = max(0., candidate), i
        for i in indices:
            x[i] += alpha * direction[i]
        if blocker is not None:
            x[blocker] = 0.
            free.remove(blocker)
    return finish(max_iterations, "iteration_limit", "接触活动集达到显式迭代上限")
