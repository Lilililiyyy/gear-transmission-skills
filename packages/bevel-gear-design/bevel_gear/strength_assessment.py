"""以显式设计应力限值校核标准名义应力，不臆造材料容量或寿命。

输入限值应已由使用者针对材料、热处理、寿命、表面及润滑条件确定。
安全系数S=设计应力限值/计算应力；许用值=限值/最低安全系数。
本层不将名义应力与其他标准版本的S-N容量暗中拼接。
"""
from dataclasses import asdict, dataclass
import math


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(name + '须为正有限实数')
    return float(value)


@dataclass(frozen=True)
class DesignStressLimits:
    """显式设计限值，Pa；不是由本模块自动推断的完整标准疲劳容量。"""
    contact_stress_limit_pa: float
    pinion_root_stress_limit_pa: float
    wheel_root_stress_limit_pa: float
    minimum_contact_safety: float
    minimum_root_safety: float
    basis: str

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name != 'basis':
                _positive(value, name)
        if not isinstance(self.basis, str) or not self.basis.strip():
            raise ValueError('设计限值必须说明材料/寿命或假设依据')
        if min(self.minimum_contact_safety, self.minimum_root_safety) < 1:
            raise ValueError('最低安全系数不得小于1')


def assess_nominal_strength(nominal_result, limits):
    """校核H、F1、F2；返回安全系数、目标、许用应力与利用率。"""
    if isinstance(limits, dict):
        limits = DesignStressLimits(**limits)
    if not isinstance(limits, DesignStressLimits):
        raise TypeError('需要DesignStressLimits或同字段字典')
    stresses = nominal_result['stress_pa']
    definitions = (('H', limits.contact_stress_limit_pa, limits.minimum_contact_safety),
                   ('F1', limits.pinion_root_stress_limit_pa, limits.minimum_root_safety),
                   ('F2', limits.wheel_root_stress_limit_pa, limits.minimum_root_safety))
    rows = []
    for quantity, resistance, required in definitions:
        stress = _positive(stresses[quantity], quantity + '应力Pa')
        allowable = resistance / required
        safety = resistance / stress
        rows.append(dict(quantity=quantity, calculated_stress_pa=stress, design_stress_limit_pa=resistance,
            minimum_safety=required, allowable_stress_pa=allowable, safety_factor=safety,
            utilization=stress/allowable, meets_target=safety >= required))
    return dict(rows=rows, meets_all_targets=all(row['meets_target'] for row in rows),
        limits=asdict(limits), automatic_material_capacity=False, fatigue_life_computed=False,
        definition='S=显式设计应力限值/计算应力；utilization=计算应力/(限值/目标安全系数)',
        scope='独立设计限值校核；不宣称自动完成标准材料容量、S-N寿命或LTCA应力场验证')
