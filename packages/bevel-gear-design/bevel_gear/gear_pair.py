"""按齿轮副对象组织现有计算；几何、标准评级与接触模型保持各自边界。

配置、材料与载荷不隐式补齐；CAD、数值求解和有限元依赖仅在实际
调用时加载。未实现的对象计算域明确拒绝，不用扭转直齿代替制造面。
"""
from __future__ import annotations

from copy import deepcopy
import math

from .geometry import GeometryInput, calculate_geometry


class _BevelPair:
    def __init__(self, config):
        if not isinstance(config, dict):
            raise TypeError("齿轮副配置须为具名字典")
        self._config = deepcopy(config)
        if "geometry" not in self._config:
            raise ValueError("齿轮副配置须包含geometry")
        self._geometry = calculate_geometry(GeometryInput(**self._config["geometry"]))

    @property
    def config(self):
        """返回输入副本，调用方修改不会改变已建立对象。"""
        return deepcopy(self._config)

    def geometry(self):
        """返回宏观几何副本；包含诊断字典，避免外部修改内部状态。"""
        return deepcopy(self._geometry)

    def analyze(self):
        """使用显式材料、载荷与标准配置运行现有工程评级。"""
        from .api import analyze_design
        return analyze_design(deepcopy(self._config))

    def nominal_strength_2014(self, **inputs):
        """显式2014 B1名义应力分支；不混入既有容量／安全系数或寿命。

        材料密度、实际单齿距偏差及刀具等输入须明确提供；未指定的扭矩
        和角速度可由本对象load读取，rpm仅在此边界换算为rad/s。
        """
        from .iso2014_nominal import evaluate_iso2014_nominal
        values=dict(inputs)
        load=self._config.get('load',{})
        if 'torque_nm' not in values:
            if 'pinion_torque_n_m' not in load:raise ValueError('2014名义链须显式给扭矩或配置load')
            values['torque_nm']=load['pinion_torque_n_m']
        if 'pinion_speed_rad_s' not in values:
            if 'pinion_speed_rpm' not in load:raise ValueError('2014名义链须显式给角速度或配置load')
            values['pinion_speed_rad_s']=load['pinion_speed_rpm']*2*math.pi/60
        return evaluate_iso2014_nominal(self.geometry(),**values)


class StraightBevelPair(_BevelPair):
    """零偏置理想直齿副；CAD工程根、孔壁约束及FE限制沿用原接口。"""

    def __init__(self, config):
        super().__init__(config)
        if self._geometry.inputs.spiral_angle_rad != 0 or self._geometry.diagnostics["tooth_trace"] != "straight":
            raise ValueError("StraightBevelPair要求螺旋角为0且实际齿线为straight")

    def build_models(self, pinion_settings, wheel_settings, *, full_gears=True):
        """按显式孔径和工程根参数生成两轮局部闭实体，不隐式安装。"""
        from .straight_cad import StraightCadSettings, build_straight_gear, build_straight_tooth_sector
        if not isinstance(full_gears, bool):
            raise ValueError("full_gears须为布尔值")
        builder = build_straight_gear if full_gears else build_straight_tooth_sector
        models = {}
        for member, settings in (("pinion", pinion_settings), ("wheel", wheel_settings)):
            if isinstance(settings, dict):
                settings = StraightCadSettings(**settings)
            models[member] = builder(self.geometry(), settings, member=member)
        return models

    def tca(self, driver, **kwargs):
        """主动增量角rad的理想相切分析，不以标准重合度替代接触求解。"""
        from .straight_tca import solve_straight_tca
        return solve_straight_tca(self.geometry(), driver, **kwargs)

    def prepare_ltca(self, **kwargs):
        """显式传入材料、CAD和网格参数；支承条件沿用孔壁固定模型。"""
        from .straight_ltca import prepare_straight_ltca
        return prepare_straight_ltca(self.geometry(), **kwargs)

    def solve_ltca(self, prepared, driver, torque, **kwargs):
        """复用本对象几何的FE准备结果；角度rad、输入扭矩N·m。"""
        from .straight_ltca import PreparedStraightLtca, solve_straight_fe_ltca
        if not isinstance(prepared, PreparedStraightLtca):
            raise TypeError("需要本几何对应的PreparedStraightLtca")
        if prepared.geometry != self._geometry:
            raise ValueError("已准备LTCA的几何与当前齿轮副不一致，不能跨对象复用")
        return solve_straight_fe_ltca(prepared, driver, torque, **kwargs)


class SpiralBevelPair(_BevelPair):
    """零偏置圆弧齿线副的宏观几何与标准评级入口。

    零螺旋角圆弧齿线仍可作为Zerol对象。实际展成齿面须按独立的
    机床/刀具输入入口构建；本门面不推定缺失的制造与支承参数。
    """

    def __init__(self, config, *, binding_tolerances=None):
        super().__init__(config)
        from .manufacturing_binding import BindingTolerances
        if self._geometry.diagnostics["tooth_trace"] != "circular_arc":
            raise ValueError("SpiralBevelPair要求圆弧齿线；直齿请使用StraightBevelPair")
        if isinstance(binding_tolerances, dict):
            binding_tolerances = BindingTolerances(**deepcopy(binding_tolerances))
        if binding_tolerances is None:
            binding_tolerances = BindingTolerances()
        if not isinstance(binding_tolerances, BindingTolerances):
            raise TypeError("binding_tolerances须为BindingTolerances或具名字典")
        self._binding_tolerances = deepcopy(binding_tolerances)

    def manufacturing_binding(self, pinion_request, wheel_request):
        """返回尺寸/齿高/齿厚和逐侧角度对照；角度门限不是制造公差。"""
        from .manufacturing_binding import bind_manufacturing_geometry
        return bind_manufacturing_geometry(self.geometry(), pinion_request, wheel_request,
                                           tolerances=self._binding_tolerances)

    def build_generated_models(self, pinion_request, wheel_request, pinion_settings, wheel_settings,
                                *, width_count=25, profile_count=25):
        """显式双轮机床请求与齿坯，构建实际展成工作面的有孔整轮。

        检查尺寸、中点齿高/弧齿厚和锥度，实际两侧角度按显式名义面
        接受门限核对；不将该门限称为制造公差或配对共轭证明。
        """
        from .generated_solid import build_generated_gear
        requests = self._manufacturing_requests(pinion_request, wheel_request)
        settings = {"pinion": pinion_settings, "wheel": wheel_settings}
        models = {member: build_generated_gear(request, settings[member], width_count=width_count,
                                              profile_count=profile_count) for member, request in requests.items()}
        binding = self.manufacturing_binding(pinion_request, wheel_request)
        for model in models.values():
            model.diagnostics["object_binding"] = deepcopy(binding)
        return models

    def _manufacturing_requests(self, pinion_request, wheel_request):
        """所有制造入口复用完整绑定，错误在CAD或FE计算前明确拒绝。"""
        from .generated_surface import _request
        requests = {"pinion": _request(pinion_request), "wheel": _request(wheel_request)}
        for member, request in requests.items():
            dimensions = getattr(self._geometry, member)
            if request.tooth_count != dimensions.teeth:
                raise ValueError(f"{member}机床请求齿数与宏观对象不一致")
            checks = ((request.tooth_count*request.module_m, dimensions.outer_pitch_diameter_m, "外端节径"),
                      (request.face_width_m, self._geometry.inputs.face_width_m, "齿宽"),
                      (request.machine.pitch_cone_angle_rad, dimensions.pitch_cone_angle_rad, "节锥角"))
            for actual, expected, label in checks:
                if not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12):
                    raise ValueError(f"{member}机床请求{label}与宏观对象不一致")
        binding = self.manufacturing_binding(requests['pinion'], requests['wheel'])
        if not binding['accepted']:
            first = binding['failures'][0]
            raise ValueError(f"{first['member']}制造绑定{first['field']}超出显式接受门限；残差{first['residual']}、门限{first['limit']}")
        return requests

    def build_models(self, *args, **kwargs):
        raise NotImplementedError("弧齿实体须使用build_generated_models并显式给定机床与齿坯，不能由宏观几何暗补")

    def prepare_tca(self, pinion_request, wheel_request, **kwargs):
        """显式制造工作面的共同锥顶TCA准备；可复用同请求的闭实体。"""
        from .generated_tca import prepare_generated_tca
        requests = self._manufacturing_requests(pinion_request, wheel_request)
        return prepare_generated_tca(requests['pinion'], requests['wheel'], **kwargs)

    def tca(self, driver=None, *, prepared=None, **kwargs):
        """实际展成面局部相切；必须提供本宏观尺寸对应的显式制造准备。"""
        if prepared is None:
            raise NotImplementedError("弧齿相切须先显式prepare_tca；宏观几何不能自动生成机床工作面")
        from .generated_tca import PreparedGeneratedTca, solve_generated_tca
        if not isinstance(prepared, PreparedGeneratedTca):
            raise TypeError("需要实际制造工作面的PreparedGeneratedTca")
        self._manufacturing_requests(prepared.requests['pinion'], prepared.requests['wheel'])
        return solve_generated_tca(prepared, driver, **kwargs)

    def prepare_ltca(self, pinion_request=None, wheel_request=None, **kwargs):
        """显式制造模型的整轮FE准备；不能用宏观齿线生成缺失工艺面。"""
        if pinion_request is None or wheel_request is None:
            raise NotImplementedError("弧齿承载需显式双轮机床请求、齿坯、材料和支承输入")
        from .generated_ltca import prepare_generated_ltca
        requests=self._manufacturing_requests(pinion_request,wheel_request)
        return prepare_generated_ltca(requests['pinion'],requests['wheel'],**kwargs)

    def solve_ltca(self, prepared=None, driver=None, torque=None, **kwargs):
        """实际展成整轮无摩擦承载；扭矩与弹性角相对参数轴保留符号。"""
        if prepared is None:
            raise NotImplementedError("弧齿承载需先准备同制造请求的整轮FE，不能使用直齿柔度")
        from .generated_ltca import PreparedGeneratedLtca,solve_generated_fe_ltca
        if not isinstance(prepared,PreparedGeneratedLtca):raise TypeError('需要PreparedGeneratedLtca')
        self._manufacturing_requests(prepared.tca.requests['pinion'],prepared.tca.requests['wheel'])
        return solve_generated_fe_ltca(prepared,driver,torque,**kwargs)
