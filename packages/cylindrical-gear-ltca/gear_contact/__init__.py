"""渐开线圆柱齿轮半解析准静态承载齿面接触分析（LTCA）。

源码和函数接口可直接复制使用；内部m、N、Pa、rad，扭矩N·m。
图形依赖只在明确调用plotting时加载。模型假设与验证范围见随附文档。
"""
from .geometry import Gear, GearPair, ContactGrid, contact_grid
from .surfaces import SurfaceModification, MeasuredCurve, MeasuredSurface, CompositeModification
from .elasticity import Material, build_structural_compliance
from .hertz import HertzLineResult, effective_modulus_pa, hertz_line_contact
from .solver import ContactSolution, ContactConvergenceError, solve_contact
from .analysis import AnalysisResult, RotationResult, analyze_pair, transmission_error_spectrum
from .design import OperatingPoint, DesignStudy, compare_designs, sweep_parameter, robust_design_search
from .io import export_result
from .config import analyze_config, modification_from_config
from .strength import MeshRatingFactors, GearRatingFactors, StressLimits, evaluate_iso6336_explicit_factors
from .iso_contact import ISOContactSettings, iso_single_tooth_stiffness, build_iso_compliance
from .extended_contact import ExtendedContactSettings, extended_contact_grid

__version__ = "0.2.0"
