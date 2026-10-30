"""由互补能最小化求解无摩擦单边接触；全部长度 m、力 N。

线性部分的目标为 fᵀCf/2+gᵀf，约束 f>=0、sum(f)=F。
局部压陷 h(f) 存在时加入逐单元积分 ∫h(t)dt。拉格朗日乘子
是公共法向接近量。主动集合允许节点退出和重新进入，避免只删负力
节点可能漏解的问题；停止条件同时检查接触协调和载荷平衡。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.linalg import cho_factor, cho_solve


class ContactConvergenceError(RuntimeError):
    """接触求解未满足声明的力与间隙容差。"""


@dataclass
class ContactSolution:
    """离散单元合力及法向接近量；刚度为 dF/dδ，单位 N/m。"""

    forces_n: np.ndarray
    approach_m: float
    residual_gaps_m: np.ndarray
    active_mask: np.ndarray
    tangent_stiffness_n_per_m: float
    diagnostics: dict


def _local_values(law, force):
    if law is None:
        return np.zeros_like(force), np.zeros_like(force)
    h, derivative = (np.asarray(v, dtype=float) for v in law(force))
    if h.shape != force.shape or derivative.shape != force.shape:
        raise ValueError("局部接触规律必须返回与接触力同形状的压陷和导数")
    if not np.all(np.isfinite(h)) or np.any(h < 0) or np.any(derivative < 0) or np.any(np.isnan(derivative)):
        raise ValueError("局部压陷应非负有限，导数应非负且不能为 NaN")
    if np.any(~np.isfinite(derivative[force > 0])):
        raise ValueError("正接触力处的局部压陷导数必须有限")
    return h, derivative


def _local_energy(law, force):
    if law is None:
        return 0.0
    if hasattr(law, "energy"):
        values = np.asarray(law.energy(force), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("局部互补能必须有限")
        return float(np.sum(values))
    # 自定义局部规律可只提供压陷和导数；此时用正权 Gauss 积分补互补能。
    nodes, weights = np.polynomial.legendre.leggauss(16)
    total = 0.0
    for node, weight in zip(nodes, weights):
        h, _ = _local_values(law, force * ((node + 1) / 2))
        total += weight * float(np.dot(force, h)) / 2
    return total


def solve_contact(gaps_m, compliance_m_per_n, total_force_n: float, *,
                  local_law: Callable | None = None, initial_forces_n=None,
                  relative_tolerance: float = 2e-8,
                  absolute_gap_tolerance_m: float = 2e-12,
                  max_iterations: int = 400) -> ContactSolution:
    """解总法向力控制的接触问题。

    C 是两轮结构柔度之和，必须对称正定，不能再包含 local_law 的压陷。
    局部规律返回 (h,dh/df)，可附带 energy(f) 返回逐单元互补能。
    输入间隙保留绝对基准，不按每个相位减去最小值。
    """
    g = np.asarray(gaps_m, dtype=float)
    c = np.asarray(compliance_m_per_n, dtype=float)
    force = float(total_force_n)
    if g.ndim != 1 or not g.size or c.shape != (g.size, g.size):
        raise ValueError("间隙必须为非空一维数组，柔度矩阵应为对应方阵")
    if not np.all(np.isfinite(g)) or not np.all(np.isfinite(c)) or not np.isfinite(force) or force < 0:
        raise ValueError("间隙、柔度与非负总力必须有限")
    if (not np.isfinite(relative_tolerance) or relative_tolerance <= 0
            or not np.isfinite(absolute_gap_tolerance_m) or absolute_gap_tolerance_m <= 0
            or isinstance(max_iterations, bool) or not isinstance(max_iterations, (int, np.integer))
            or max_iterations < 1):
        raise ValueError("接触求解容差必须正且有限，迭代上限必须为正整数")
    cscale = float(np.max(np.abs(c)))
    if cscale <= 0 or not np.allclose(c, c.T, rtol=1e-10, atol=1e-13*cscale):
        raise ValueError("柔度矩阵必须满足弹性互易性，即对称")
    c = (c+c.T)/2
    try:
        cho_factor(c/cscale, lower=True, check_finite=False)
    except np.linalg.LinAlgError as error:
        raise ValueError("柔度矩阵必须正定；请先消除刚体模态并检查边界") from error
    if force == 0:
        approach = float(np.min(g))
        mask = g == approach
        stiffness = 0.0
        local_zero, derivative_zero = _local_values(local_law, np.zeros_like(g))
        if np.any(local_zero != 0):
            raise ValueError("局部压陷规律在零载荷处必须为零")
        finite_mask = mask & np.isfinite(derivative_zero)
        if np.any(finite_mask):
            # 右切线同样受非负力约束。最小间隙节点不一定都能承受正力；
            # 以单位力解其线性化问题，避免直接C⁻¹1出现负权重。
            block = c[np.ix_(finite_mask, finite_mask)]+np.diag(derivative_zero[finite_mask])
            right_limit = solve_contact(np.zeros(np.count_nonzero(finite_mask)), block, 1.0,
                                        relative_tolerance=relative_tolerance,
                                        absolute_gap_tolerance_m=absolute_gap_tolerance_m,
                                        max_iterations=max_iterations)
            stiffness = right_limit.tangent_stiffness_n_per_m
        return ContactSolution(np.zeros_like(g), approach, g-approach, np.zeros_like(mask), stiffness,
                               {"converged": True, "iterations": 0, "force_balance_n": 0.0,
                                "maximum_penetration_m": 0.0, "active_gap_error_m": 0.0,
                                "complementarity_n_m": 0.0, "zero_load": True})

    # x=f/F：单纯形约束 sum(x)=1；归一化使牛顿方向不受 SI 数量级影响。
    deformation_scale = max(float(np.median(np.diag(c)))*force/g.size, 1e-9)
    length_scale = max(deformation_scale, float(np.ptp(g)))
    hessian_linear = c*force/length_scale
    q = (g-float(np.min(g)))/length_scale
    tolerance = max(relative_tolerance*deformation_scale, absolute_gap_tolerance_m)/length_scale
    if initial_forces_n is None:
        # 先用对角柔度的单调标量平衡得到可行初值，随后求完整耦合问题。
        diagonal = np.diag(hessian_linear)
        lo, hi = float(np.min(q)), float(np.max(q)+np.max(diagonal))
        for _ in range(70):
            mid = (lo+hi)/2
            if float(np.sum(np.maximum(mid-q, 0)/diagonal)) < 1:
                lo = mid
            else:
                hi = mid
        x = np.maximum((lo+hi)/2-q, 0)/diagonal
        x /= np.sum(x)
    else:
        x = np.asarray(initial_forces_n, dtype=float)/force
        if x.shape != g.shape or not np.all(np.isfinite(x)) or np.any(x < 0) or np.sum(x) <= 0:
            raise ValueError("接触力初值必须同形、非负、有限且总和为正")
        x = x/np.sum(x)
    free = x > 1e-14
    line_search_steps = 0
    added = removed = 0

    def objective(value):
        return (float(value @ hessian_linear @ value)/2 + float(q @ value)
                + _local_energy(local_law, value*force)/(force*length_scale))

    for iteration in range(1, max_iterations+1):
        local, _ = _local_values(local_law, x*force)
        gradient = hessian_linear@x + q + local/length_scale
        positive = x > 1e-13
        multiplier = float(np.mean(gradient[positive]))
        reduced = gradient-multiplier
        active_error = float(np.max(np.abs(reduced[positive])))
        inactive_error = float(max(0.0, -np.min(reduced[~positive]))) if np.any(~positive) else 0.0
        if max(active_error, inactive_error) <= tolerance:
            break
        inactive_indices = np.flatnonzero(~free)
        if inactive_indices.size:
            candidate = inactive_indices[np.argmin(reduced[inactive_indices])]
            if reduced[candidate] < -tolerance:
                free[candidate] = True
                added += 1

        # 零力处 Hertz 导数趋于无穷；只在形成方向时使用严格正小力，
        # 目标函数、协调残差和最终载荷始终使用真实接触力。
        probe = np.maximum(x*force, force*1e-12)
        _, derivative = _local_values(local_law, probe)
        full_hessian = hessian_linear + np.diag(derivative*force/length_scale)
        while True:
            indices = np.flatnonzero(free)
            block = full_hessian[np.ix_(indices, indices)]
            factor = cho_factor(block, lower=True, check_finite=False)
            inverse_one = cho_solve(factor, np.ones(indices.size), check_finite=False)
            inverse_gradient = cho_solve(factor, gradient[indices], check_finite=False)
            lam = float(np.sum(inverse_gradient)/np.sum(inverse_one))
            direction = np.zeros_like(x)
            direction[indices] = lam*inverse_one-inverse_gradient
            blocked = free & (x <= 1e-14) & (direction < -1e-14)
            if not np.any(blocked):
                break
            free[blocked] = False
            removed += int(np.count_nonzero(blocked))
        slope = float(gradient@direction)
        if slope >= 0:
            raise ContactConvergenceError("接触牛顿方向未下降；请检查局部压陷导数与互补能的一致性")
        negative = direction < -1e-16
        alpha = min(1.0, float(np.min(-x[negative]/direction[negative]))) if np.any(negative) else 1.0
        old_objective = objective(x)
        accepted = False
        for _ in range(55):
            trial = np.maximum(x+alpha*direction, 0)
            trial /= np.sum(trial)
            try:
                energy = objective(trial)
            except ValueError:
                energy = np.inf
            if energy <= old_objective+1e-4*alpha*slope+1e-15:
                accepted = True
                break
            alpha *= 0.5
            line_search_steps += 1
        if not accepted or alpha <= 1e-16:
            raise ContactConvergenceError("接触能量线搜索未收敛；当前结果不作为有效承载解返回")
        hit = free & (trial <= 1e-14) & (direction < 0)
        free[hit] = False
        removed += int(np.count_nonzero(hit))
        trial[hit] = 0
        trial /= np.sum(trial)
        x = trial
    else:
        raise ContactConvergenceError(f"接触求解达到 {max_iterations} 次上限，尚未满足互补条件")

    forces = force*x
    local, derivative = _local_values(local_law, forces)
    positive = forces > force*1e-13
    approach = float(np.mean((g+c@forces+local)[positive]))
    residual = g+c@forces+local-approach
    tangent_matrix = c[np.ix_(positive, positive)]+np.diag(derivative[positive])
    tangent = float(np.sum(np.linalg.solve(tangent_matrix, np.ones(np.count_nonzero(positive)))))
    diagnostics = {
        "converged": True, "iterations": iteration,
        "force_balance_n": float(np.sum(forces)-force),
        "maximum_penetration_m": float(max(0, -np.min(residual))),
        "active_gap_error_m": float(np.max(np.abs(residual[positive]))),
        "complementarity_n_m": float(np.max(np.abs(forces*residual))),
        "active_cells": int(np.count_nonzero(positive)),
        "energy_j": objective(x)*force*length_scale+force*float(np.min(g)),
        "added_cells": added, "removed_cells": removed,
        "line_search_reductions": line_search_steps,
        "gap_tolerance_m": tolerance*length_scale,
    }
    return ContactSolution(forces, approach, residual, positive, tangent, diagnostics)
