"""新写三维装配视图：旋转包络、共同锥顶及明确装配基向量下的力。

长度m、力N。将已计算的齿顶／齿根边界绕轴旋转形成封闭齿高包络，
不构造实际齿面、孔、轮毂或加工齿坯。仅依赖标准库，不改动计算输入。
"""
from __future__ import annotations

from dataclasses import asdict, is_dataclass
from html import escape
from importlib.resources import files
import json
import math
from pathlib import Path


def _mapping(value):
    """接受计算数据类或其具名JSON字典。"""
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return value
    raise TypeError("三维视图输入须为计算数据类或具名字典")


def _add(a, b):
    return [x + y for x, y in zip(a, b)]


def _scale(vector, value):
    return [value * x for x in vector]


def _norm(vector):
    return math.sqrt(sum(x * x for x in vector))


def _point(axis, radial, tangent, axial, radius, angle):
    """共同锥顶为原点；轴向从锥顶向齿轮外端为正。"""
    return _add(_scale(axis, axial),
                _add(_scale(radial, radius * math.cos(angle)),
                     _scale(tangent, radius * math.sin(angle))))


def _envelope(dimensions, axis, radial, tangent, segments):
    """按内根→外根→外顶→内顶旋转，形成朝外法向的封闭四环网格。"""
    profile = [(dimensions[f"{end}_{kind}_axial_coordinate_m"],
                dimensions[f"{end}_{kind}_diameter_m"] / 2)
               for end, kind in (("inner", "root"), ("outer", "root"),
                                  ("outer", "tip"), ("inner", "tip"))]
    if any(not math.isfinite(v) for pair in profile for v in pair) or min(r for _, r in profile) <= 0:
        raise ValueError("齿顶／齿根边界须为有限坐标和正半径")
    vertices = [_point(axis, radial, tangent, axial, radius, 2 * math.pi * i / segments)
                for axial, radius in profile for i in range(segments)]
    triangles = []
    for ring in range(4):
        next_ring = (ring + 1) % 4
        for i in range(segments):
            j = (i + 1) % segments
            a, b, c, d = ring * segments + i, next_ring * segments + i, next_ring * segments + j, ring * segments + j
            triangles.extend(([a, b, c], [a, c, d]))
    return {"vertices_m": vertices, "triangles": triangles,
            "ring_order": ["inner_root", "outer_root", "outer_tip", "inner_tip"]}


def build_assembly_scene(geometry, forces=None, *, name="锥齿轮装配", segments=96) -> dict:
    """从宏观几何生成可复用的SI三维网格与力场景。

    全局Z为小轮轴，另一轴位于XZ平面；共同节锥母线同样位于XZ。
    每轮径向指向共同母线，轴向力正基向量指向锥顶，切向基向量取
    轴线×径向。该规范装配补足原局部力的方向约定，不推断实际转向。
    两轮作用点相同，结果力按各自基向量转换，保留全局力平衡残差。
    """
    if isinstance(segments, bool) or not isinstance(segments, int) or not 12 <= segments <= 256:
        raise ValueError("圆周分段数须为12至256的整数")
    g = _mapping(geometry)
    f = None if forces is None else _mapping(forces)
    delta1, delta2 = (g[role]["pitch_cone_angle_rad"] for role in ("pinion", "wheel"))
    if not all(math.isfinite(d) and 0 < d < math.pi / 2 for d in (delta1, delta2)):
        raise ValueError("三维装配要求两个节锥角均为锐角")
    shaft_angle = delta1 + delta2
    radius = g["outer_cone_distance_m"]
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("外锥距须为正有限长度")
    common = [math.sin(delta1), 0., math.cos(delta1)]
    contact = _scale(common, g["mean_cone_distance_m"])
    bases = (([0., 0., 1.], [1., 0., 0.], [0., 1., 0.]),
             ([math.sin(shaft_angle), 0., math.cos(shaft_angle)],
              [-math.cos(shaft_angle), 0., math.sin(shaft_angle)], [0., -1., 0.]))
    gears, vectors = [], []
    for role, label, color, basis in zip(("pinion", "wheel"), ("小轮", "大轮"),
                                        ("#2166AC", "#B2182B"), bases):
        axis, radial, tangent = basis
        dimensions = g[role]
        mesh = _envelope(dimensions, axis, radial, tangent, segments)
        rings = {}
        for end in ("inner", "mean", "outer"):
            axial, r = dimensions[f"{end}_axial_pitch_coordinate_m"], dimensions[f"{end}_pitch_diameter_m"] / 2
            rings[end] = [_point(axis, radial, tangent, axial, r, 2 * math.pi * i / segments)
                          for i in range(segments)]
        gear = {"id": role, "label": label, "color": color, "teeth": dimensions["teeth"],
                "pitch_angle_rad": dimensions["pitch_cone_angle_rad"],
                "outer_pitch_diameter_m": dimensions["outer_pitch_diameter_m"],
                "mean_pitch_diameter_m": dimensions["mean_pitch_diameter_m"],
                "frame": {"shaft_outward": axis, "radial": radial, "tangential": tangent,
                          "axial_toward_apex": _scale(axis, -1.)},
                "mesh": mesh, "pitch_rings_m": rings}
        if f is not None:
            components = {}
            for component, direction in (("tangential", tangent), ("radial", radial),
                                          ("axial", _scale(axis, -1.))):
                value = f[f"{role}_{component}_n"]
                if not math.isfinite(value):
                    raise ValueError("显示力分量须为有限值，单位N")
                components[component] = {"value_n": value, "vector_n": _scale(direction, value)}
            resultant = [sum(c["vector_n"][i] for c in components.values()) for i in range(3)]
            gear["forces"] = {"origin_m": contact, "components": components,
                              "resultant_n": resultant, "magnitude_n": _norm(resultant)}
            vectors.append(resultant)
        gears.append(gear)
    all_vertices = [[0., 0., 0.]] + [v for gear in gears for v in gear["mesh"]["vertices_m"]]
    bounds = {"min_m": [min(v[i] for v in all_vertices) for i in range(3)],
              "max_m": [max(v[i] for v in all_vertices) for i in range(3)]}
    scene = {"schema_version": "1.0", "name": str(name), "length_unit": "m", "force_unit": "N",
             "geometry_level": "axisymmetric_tooth_height_envelope", "segments": segments,
             "shaft_angle_rad": shaft_angle, "face_width_m": g["inputs"]["face_width_m"],
             "spiral_angle_rad": g["inputs"]["spiral_angle_rad"],
             "origin_m": [0., 0., 0.], "mean_contact_point_m": contact,
             "common_generator_outer_m": _scale(common, radius), "gears": gears, "bounds": bounds,
             "limitations": ["显示节锥及齿顶／齿根旋转包络，未生成实际齿形、轮毂和孔",
                             "两轮齿高包络可能交叠，不据此判断实体干涉或接触斑",
                             "力方向采用明示规范装配；观察旋转不表示齿轮运动仿真"]}
    if vectors:
        residual = _norm(_add(*vectors))
        normal = max(_norm(v) for v in vectors)
        scene["force_display_scale_m_per_n"] = .52 * radius / normal if normal else 0.
        scene["force_balance_residual_n"] = residual
        scene["load"] = {"flank": f["flank"], "spiral_hand": f["spiral_hand"],
                         "pinion_torque_n_m": f["pinion_torque_n_m"]}
    return scene


def _surface_transform(point, frame, phase):
    """将同一解析齿面的局部坐标旋转到规范装配；正交变换保留法向。"""
    x,y,z=point
    x,y=x*math.cos(phase)-y*math.sin(phase),x*math.sin(phase)+y*math.cos(phase)
    return [x*frame["radial"][i]+y*frame["tangential"][i]+z*frame["shaft_outward"][i] for i in range(3)]


def build_tooth_assembly_scene(geometry, forces=None, *, name="理想直齿齿面", width_count=5, profile_count=13) -> dict:
    """把有效解析球面渐开线及其周期阵列放入装配；不封造根部和轮体。"""
    from .geometry import GeometryInput,GeometryResult,calculate_geometry
    from .tooth_surface import build_spherical_tooth_surface
    g=geometry if isinstance(geometry,GeometryResult) else calculate_geometry(GeometryInput(**_mapping(geometry)["inputs"]))
    scene=build_assembly_scene(g,forces,name=name)
    for gear in scene["gears"]:
        mesh=asdict(build_spherical_tooth_surface(g,gear["id"]).mesh_full_gear(width_count,profile_count))
        phase=0. if gear["id"]=="pinion" else math.pi/gear["teeth"]
        mesh["vertices_m"]=[_surface_transform(p,gear["frame"],phase) for p in mesh["vertices_m"]]
        mesh["normals"]=[_surface_transform(n,gear["frame"],phase) for n in mesh["normals"]]
        mesh["metadata"].update(coordinate_frame="规范装配全局坐标；锥顶原点，Z沿小轮轴，另一轴位于XZ平面",
                                initial_assembly_phase_rad=phase)
        gear["surface_mesh"]=mesh
        gear["rotation_ratio"]=1. if gear["id"]=="pinion" else -g.inputs.pinion_teeth/g.inputs.wheel_teeth
    scene["geometry_level"]="ideal_straight_spherical_involute_flanks"
    scene["limitations"]=["只含理想直齿有效球面渐开线，不含真实刀具包络、基圆以下根过渡、孔和轮体",
                           "转角按理想齿数比，不求接触对齐、传动误差或干涉",
                           "啮合力仍标在平均节点，不追踪实际随动接触点"]
    return scene


def export_assembly_view(result: dict, output_dir, *, segments=96, with_tooth_surfaces=False,
                         width_count=5, profile_count=13) -> list[str]:
    """导出自包含离线HTML及具名SI场景JSON；不需要绘图库或服务器。"""
    if with_tooth_surfaces:
        scene=build_tooth_assembly_scene(result["geometry"],result.get("forces"),name=result.get("name","理想直齿齿面"),
                                         width_count=width_count,profile_count=profile_count)
    else:
        scene = build_assembly_scene(result["geometry"], result.get("forces"),
                                     name=result.get("name", "锥齿轮装配"), segments=segments)
    raw = json.dumps(scene, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    # 转义HTML脚本结束标记；标题作为纯文本，设计名称不能注入标签或脚本。
    embedded = raw.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    template = files("bevel_gear").joinpath("assets/viewer.html").read_text(encoding="utf-8")
    page = template.replace("__TITLE__", escape(scene["name"], quote=True)).replace("__SCENE_JSON__", embedded)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    data_path, html_path = output / "assembly_scene.json", output / "assembly.html"
    data_path.write_text(json.dumps(scene, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    html_path.write_text(page, encoding="utf-8")
    paths=[str(html_path.resolve()), str(data_path.resolve())]
    if with_tooth_surfaces:
        from .tooth_surface import ToothSurfaceMesh,export_surface_obj,export_surface_json
        for gear in scene["gears"]:
            mesh=ToothSurfaceMesh(**gear["surface_mesh"])
            paths.extend(str(p.resolve()) for p in (export_surface_obj(mesh,output/f"{gear['id']}_flanks.obj"),
                                                     export_surface_json(mesh,output/f"{gear['id']}_flanks.json")))
    return paths
