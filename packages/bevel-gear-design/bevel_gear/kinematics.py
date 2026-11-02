"""新写理想齿数比运动学；不求齿面接触、安装误差或传动误差。"""
import math


def ideal_pair_pose(pinion_teeth: int, wheel_teeth: int, pinion_angle_rad: float, *,
                    pinion_speed_rpm: float = 0., initial_wheel_angle_rad: float | None = None) -> dict:
    """两轴都从锥顶指向外端时，外啮合角速度比为−z1/z2。

    初始大轮相位默认半齿距，使小轮齿中心对大轮齿槽中心。该相位用于
    规则阵列展示，不替代实际齿面接触对齐；角度不折回0…2π以保留圈数。
    """
    for name,value in (("小轮齿数",pinion_teeth),("大轮齿数",wheel_teeth)):
        if isinstance(value,bool) or not isinstance(value,int) or value < 3:
            raise ValueError(f"{name}必须为不小于3的整数")
    for name,value in (("小轮角度rad",pinion_angle_rad),("小轮转速rpm",pinion_speed_rpm)):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError(f"{name}必须为有限数")
    phase=math.pi/wheel_teeth if initial_wheel_angle_rad is None else initial_wheel_angle_rad
    if isinstance(phase,bool) or not isinstance(phase,(int,float)) or not math.isfinite(phase):
        raise ValueError("大轮初相位必须为有限角度rad")
    multiplier=-pinion_teeth/wheel_teeth
    return {"pinion_angle_rad":pinion_angle_rad,"wheel_angle_rad":phase+multiplier*pinion_angle_rad,
            "initial_wheel_angle_rad":phase,"wheel_relative_angle_rad":multiplier*pinion_angle_rad,
            "pinion_speed_rpm":pinion_speed_rpm,"wheel_speed_rpm":multiplier*pinion_speed_rpm,
            "angular_velocity_ratio":multiplier,
            "model":"ideal_external_ratio","limitations":["刚性理想齿数比，不含传动误差、弹性变形、碰撞或接触斑"]}
