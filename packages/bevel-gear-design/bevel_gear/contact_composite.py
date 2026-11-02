"""在原P1三角域直接构造压力等值面及逐点最大值，不经过笛卡尔掩码。

采用恒等式 {max_i p_i >= L} = union_i {p_i >= L}。每个源三角内
压力为原仿射函数，Matplotlib原生三角等值线给出阈值交线，GEOS保留
孔洞并完成并集与实际轮廓裁切。色阶只量化显示颜色，不修改源压力。
瞬时状态和离散全过程分别提供覆盖域，不能把全过程证书套用为瞬时证书。
"""
from __future__ import annotations

import hashlib
import gzip
import json
import math
import zlib
from pathlib import Path


_SCHEMA = 'bevel_triangle_pressure_map_v1'
_METHODS = ('full_gear_fe', 'local_halfspace')
_KINDS = ('instantaneous', 'sampled_cycle_maximum')


def _area_geometry(value, name, *, allow_empty=True):
    """统一保留派生面域的正面积原Polygon；不修复或改变二维区域。

    零面积Polygon须满足64ε×自身包围盒面积的凸包舍入界；有效点/线直接丢弃。
    原多边形坐标/顺序、孔洞与压力不改；outline仍必须有正面积。
    """
    import numpy as np
    import shapely
    from shapely.geometry import shape, MultiPolygon
    geometry = shape(value) if isinstance(value, dict) else value
    if geometry is None:
        if not allow_empty:
            raise ValueError(name+'不得为空')
        return shapely.GeometryCollection()
    if (not hasattr(geometry, 'is_valid') or not geometry.is_valid or geometry.has_z or
            not np.isfinite(shapely.get_coordinates(geometry)).all()):
        raise ValueError(name+'须为有效有限二维面域')
    polygons = []
    def collect(part):
        if part.is_empty:
            return
        if part.geom_type in ('MultiPolygon', 'GeometryCollection', 'MultiLineString', 'MultiPoint'):
            for child in part.geoms:
                collect(child)
            return
        if not part.is_valid:
            raise ValueError(name+'含无效分量，不修复或删除自交面')
        if part.geom_type == 'Polygon':
            if part.area == 0:
                x0,y0,x1,y1 = part.bounds
                bbox_area = (x1-x0)*(y1-y0)
                if part.convex_hull.area > 64*np.finfo(float).eps*bbox_area:
                    raise ValueError(name+'的零净面积多边形仍包围正面域，拒绝删除')
            else:
                polygons.append(part)
        elif part.area != 0:
            raise ValueError(name+'的非多边形分量含非零面积')
        # 有效点/线本来为0D/1D；弯折线和闭合ring不以凸包面积冒充面域。
    collect(geometry)
    # 保留原Polygon及其顺序，不union；只容许同一组面积求和的机器舍入。
    result = (polygons[0] if len(polygons) == 1 else MultiPolygon(polygons)) if polygons else shapely.GeometryCollection()
    area_roundoff = 64*np.finfo(float).eps*max(abs(geometry.area),abs(result.area),1e-30)
    if (not result.is_valid or abs(result.area-geometry.area)>area_roundoff or
            result.symmetric_difference(geometry).area != 0):
        raise ArithmeticError(name+'规范化改变正面积区域，拒绝输出')
    if result.is_empty and not allow_empty:
        raise ValueError(name+'须为正面积Polygon或MultiPolygon')
    return result


def _contour_geometry(vertices, codes):
    """按MOVETO分出真实封闭环；GEOS构面保留内孔和孔内孤岛。"""
    import numpy as np
    import shapely
    from shapely.geometry import MultiLineString
    rings = []
    for points, kinds in zip(vertices, codes):
        starts = np.flatnonzero(kinds == 1)
        for start, end in zip(starts, np.r_[starts[1:], len(points)]):
            ring = points[start:end]
            if len(ring) < 3:
                continue
            if not np.array_equal(ring[0], ring[-1]):
                ring = np.vstack((ring, ring[0]))
            if len(np.unique(ring[:-1], axis=0)) >= 3:
                rings.append(ring)
    if not rings:
        return shapely.GeometryCollection()
    result = shapely.build_area(MultiLineString(rings))
    if not result.is_valid:
        raise ArithmeticError('原P1等值环未构成有效面域，停止而非修补拓扑')
    return result


class OriginalP1Contours:
    """隔离原生三角等值线适配；输入m、Pa，保留原连接和标量。

    允许相邻单元使用各自的重合顶点，以保留原P0或不连续P1压力。
    不把单元常压力中心值主动平均成另一种连续模型。
    """

    def __init__(self, field):
        import numpy as np
        import shapely
        from matplotlib.tri import Triangulation
        from matplotlib import _tri
        xy = np.asarray(field['points_width_height_m'], dtype=float)
        triangles = np.asarray(field['triangles'])
        pressure = np.asarray(field['pressure_pa'], dtype=float)
        if (xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 3 or
                pressure.shape != (len(xy),) or not np.isfinite(xy).all() or
                not np.isfinite(pressure).all() or np.any(pressure < 0)):
            raise ValueError('原P1坐标须为有限二维m坐标，压力须为同长非负Pa数组')
        if (triangles.ndim != 2 or triangles.shape[1] != 3 or not len(triangles) or
                triangles.dtype.kind not in 'iu' or np.any(triangles < 0) or
                np.any(triangles >= len(xy))):
            raise ValueError('原P1连接须为有效三角形整数索引')
        triangles = triangles.astype(int, copy=True)
        self.original_triangles = triangles.copy()
        signed = np.cross(xy[triangles[:, 1]]-xy[triangles[:, 0]],
            xy[triangles[:, 2]]-xy[triangles[:, 0]])
        if np.any(signed == 0):
            raise ValueError('原P1场包含退化投影三角形')
        triangles[signed < 0] = triangles[signed < 0][:, [0, 2, 1]]
        self.xy = xy.copy()
        self.triangles = triangles
        self.pressure = pressure.copy()
        polygons = shapely.polygons(xy[triangles])
        full_domain = shapely.union_all(polygons)
        area_sum = float(np.abs(signed).sum()*.5)
        if area_sum-full_domain.area > 64*np.finfo(float).eps*area_sum:
            raise ValueError('同一源场的真实三角形存在正面积交叠，压力函数不唯一')
        positive_triangles = np.any(pressure[triangles] > 0, axis=1)
        # 合法P0场可以有坐标重合、标量不同的独立单元顶点，原生定位器
        # 不接受这种拓扑；等值线遍历使用原连接，并由GEOS严格检查面积交叠。
        self.mesh = Triangulation(xy[:, 0], xy[:, 1], triangles, mask=~positive_triangles)
        try:
            self.generator = _tri.TriContourGenerator(self.mesh.get_cpp_triangulation(), self.pressure)
        except AttributeError as error:
            raise RuntimeError('当前Matplotlib缺少已验证的原生三角等值线接口') from error
        self.maximum = float(pressure.max())
        self.upper = float(np.nextafter(max(self.maximum, 1.), math.inf))
        self.positive_support = full_domain if positive_triangles.all() else shapely.union_all(polygons[positive_triangles])
        self.fingerprint = _field_fingerprint(self.xy, self.original_triangles, self.pressure)

    def superlevel(self, level_pa):
        """正阈值返回p≥L；零阈值专指正压支撑闭包，不冒充零压证明。"""
        import shapely
        if isinstance(level_pa, bool) or not isinstance(level_pa, (int, float)) or not math.isfinite(level_pa) or level_pa < 0:
            raise ValueError('压力阈值须为有限非负Pa数值')
        if level_pa == 0:
            return self.positive_support
        if level_pa > self.maximum or self.maximum == 0:
            return shapely.GeometryCollection()
        return _contour_geometry(*self.generator.create_filled_contour(float(level_pa), self.upper))


def _field_fingerprint(xy, triangles, pressure):
    """固定小端数值布局的内容摘要；不把标签当作数值身份。"""
    import numpy as np
    digest = hashlib.sha256()
    for name, array, dtype in (('xy', xy, '<f8'), ('triangles', triangles, '<i8'), ('pressure', pressure, '<f8')):
        value = np.ascontiguousarray(array, dtype=dtype)
        digest.update(json.dumps([name, list(value.shape), dtype], separators=(',', ':')).encode('ascii'))
        digest.update(value.tobytes())
    return digest.hexdigest()


def triangle_pressure_superlevel(fields, level_pa, *, outline_geometry):
    """直接求任意拓扑原P1场的最大值超水平集；不规则采样、不填孔。"""
    import shapely
    if isinstance(level_pa, bool) or not isinstance(level_pa, (int, float)) or not math.isfinite(level_pa) or level_pa < 0:
        raise ValueError('压力阈值须为有限非负Pa数值')
    outline = _area_geometry(outline_geometry, '实际工作齿面轮廓', allow_empty=False)
    engines = [field if isinstance(field, OriginalP1Contours) else OriginalP1Contours(field) for field in fields]
    return shapely.union_all([engine.superlevel(level_pa) for engine in engines]).intersection(outline)


def build_triangle_pressure_map(fields, *, method, kind, outline_geometry, levels_pa,
        unknown_geometry=None, certified_zero_geometry=None, source_identity=None, progress=None):
    """构造一个瞬时或离散最大压力图的数据，完整保留原源场及域语义。

    fields每项含points_width_height_m[N,2]、triangles[M,3]、pressure_pa[N]。
    phase_index、tooth_offset和source_state_sha256为可选身份字段。瞬时图
    只允许同一相位、同一物理齿；周期图可以合并周期等价的不同相位/齿。
    unknown_geometry必须来自本图对应的审核；未给时整齿面保持未认证。
    certified_zero_geometry只能传已有明确零压证明的域，不能按支撑外推零。
    """
    import numpy as np
    import shapely
    if method not in _METHODS or kind not in _KINDS:
        raise ValueError('需要明确的FE/半解析方法及瞬时/离散全过程图种类')
    if not isinstance(fields, (list, tuple)) or not fields:
        raise ValueError('需要非空原P1场列表，不能以占位场代替真实解')
    if progress is not None and not callable(progress):
        raise ValueError('进度回调须可调用')
    levels = np.asarray(levels_pa, dtype=float)
    if levels.ndim != 1 or len(levels) < 2 or not np.isfinite(levels).all() or levels[0] != 0 or np.any(np.diff(levels) <= 0):
        raise ValueError('显示压力色阶须从0开始、严格升序且为有限Pa数组')
    spacing = np.diff(levels)
    if not np.allclose(spacing, spacing[0], rtol=1e-12, atol=0):
        raise ValueError('本接口采用线性绝对压力色阶，不接受百分位或非线性截断')
    outline = _area_geometry(outline_geometry, '实际工作齿面轮廓', allow_empty=False)
    unknown = outline if unknown_geometry is None else _area_geometry(unknown_geometry, '本图未认证域').intersection(outline)
    zero = _area_geometry(certified_zero_geometry, '有证据零压域').intersection(outline)
    # GEOS的差集/交集可能留下浮点尺度共边薄片；只容许64ε面积量级，
    # 不以接触斑大小、压力百分位或图形分辨率放宽覆盖证据。
    coverage_tolerance = 64*np.finfo(float).eps*outline.area
    zero_unknown_overlap = zero.intersection(unknown).area
    if zero_unknown_overlap > coverage_tolerance:
        raise ValueError('零压证明域与本图未认证域存在正面积冲突')
    identities = [(field.get('phase_index'), field.get('tooth_offset')) for field in fields]
    if kind == 'instantaneous' and (len(set(identities)) != 1 or
            len(fields) > 1 and any(phase is None or tooth is None for phase, tooth in identities)):
        raise ValueError('瞬时图必须限定同一实际相位和同一物理齿，不能叠加不同齿的斑点')
    engines = [OriginalP1Contours(field) for field in fields]
    maximum = max(engine.maximum for engine in engines)
    if levels[-1] < maximum:
        raise ValueError('线性色标不得低于原源场峰值，不能隐去狭窄离散峰')
    support = shapely.union_all([engine.positive_support for engine in engines]).intersection(outline)
    if zero.intersection(support).area > coverage_tolerance:
        raise ValueError('零压证明域与原正压三角支撑存在冲突')
    geometries = [support]
    previous = support
    for index, level in enumerate(levels[1:-1], start=1):
        geometry = shapely.union_all([engine.superlevel(float(level)) for engine in engines
            if engine.maximum >= level]).intersection(outline)
        # 此检查仅接受GEOS舍入量级边界误差，不将相交域改成凸包或平滑斑点。
        excess = geometry.difference(previous).area
        if excess > outline.area*1e-12:
            raise ArithmeticError('最大压力超水平集不满足嵌套性')
        geometries.append(geometry)
        previous = geometry
        if progress is not None:
            progress(dict(stage='triangle_pressure_superlevels', method=method, kind=kind,
                completed_level_count=index, total_level_count=len(levels)-1))
    sources = []
    for field, engine in zip(fields, engines):
        identity = {key: field[key] for key in ('phase_index', 'tooth_offset', 'driver_angle_rad', 'source_state_sha256') if key in field}
        sources.append(dict(identity, points_width_height_m=engine.xy.copy(), triangles=engine.original_triangles.copy(),
            pressure_pa=engine.pressure.copy(), numeric_sha256=engine.fingerprint))
    unresolved_support = support.intersection(unknown)
    return dict(schema=_SCHEMA, method=method, kind=kind, levels_pa=levels.copy(),
        superlevel_geometries=geometries, outline_geometry=outline, unknown_geometry=unknown,
        certified_zero_geometry=zero, positive_support_geometry=support,
        unclassified_geometry=outline.difference(support.union(zero).union(unknown)),
        source_fields=sources, source_identity={} if source_identity is None else source_identity,
        source_peak_pressure_pa=maximum, source_field_count=len(fields),
        positive_support_projected_area_m2=float(support.area),
        unknown_projected_area_fraction=float(unknown.area/outline.area),
        unknown_positive_support_area_m2=float(unresolved_support.area),
        zero_unknown_roundoff_overlap_area_m2=float(zero_unknown_overlap),
        coverage_geometry_roundoff_tolerance_m2=float(coverage_tolerance),
        numerical_definition='P1超水平集逐场并集：{max_i p_i >= L}=union_i{p_i >= L}',
        zero_threshold_definition='L=0时仅表示正压支撑闭包；有证据零压域单独给定',
        query_grid_resampling=False, original_pressure_values_modified=False,
        physical_validation=False, temporal_convergence_certified=False,
        force_integral_defined=False,
        integration_note='子午投影不是等面积展开；法向力须沿原求解三维面积与原积分规则计算')


def export_triangle_pressure_map(data, path, *, include_source_fields=True):
    """保存可独立重绘的SI JSON；.json.gz仅无损压缩同一UTF-8载荷。

    gzip采用固定时间戳且不存文件名，不携带工作目录或额外对象协议。
    原三角场可随文件保存供完全重算；默认plain文件兼容原接口。
    """
    import numpy as np
    from shapely.geometry import mapping
    if data.get('schema') != _SCHEMA:
        raise ValueError('需要原三角域压力图数据schema')
    output = {}
    for key, value in data.items():
        if key == 'source_fields':
            if include_source_fields:
                output[key] = [{k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in field.items()}
                    for field in value]
        elif key == 'superlevel_geometries':
            output[key] = [mapping(geometry) for geometry in value]
        elif hasattr(value, 'geom_type'):
            output[key] = mapping(value)
        else:
            output[key] = value.tolist() if isinstance(value, np.ndarray) else value
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload=(json.dumps(output, ensure_ascii=False, allow_nan=False, separators=(',', ':'))+'\n').encode('utf-8')
    compressed=target.name.endswith('.json.gz')
    if target.suffix=='.gz' and not compressed:
        raise ValueError('压力JSON压缩文件须使用.json.gz后缀')
    target.write_bytes(gzip.compress(payload,compresslevel=6,mtime=0) if compressed else payload)
    return dict(file=target.name, sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        original_source_fields_included=include_source_fields,
        compression='gzip' if compressed else 'plain',decoded_bytes=len(payload),
        decoded_sha256=hashlib.sha256(payload).hexdigest())


def load_triangle_pressure_map(path, *, recompute=False):
    """透明读取SI JSON/.json.gz；解压后执行相同schema和数值摘要核验。"""
    import numpy as np
    from shapely.geometry import shape
    def invalid_constant(value):
        raise ValueError('压力数据不得含NaN或Infinity')
    target=Path(path)
    try:
        payload=target.read_bytes()
        if target.name.endswith('.json.gz'):
            payload=gzip.decompress(payload)
        elif target.suffix=='.gz':
            raise ValueError('压力JSON压缩文件须使用.json.gz后缀')
        data=json.loads(payload.decode('utf-8-sig'),parse_constant=invalid_constant)
    except (OSError,EOFError,UnicodeDecodeError,zlib.error) as error:
        raise ValueError('压力文件不是有效的UTF-8 JSON或标准gzip JSON流') from error
    if not isinstance(data,dict):
        raise ValueError('压力JSON顶层必须是相同schema的对象')
    if data.get('schema') != _SCHEMA:
        raise ValueError('需要原三角域压力图数据schema')
    if data.get('method') not in _METHODS or data.get('kind') not in _KINDS:
        raise ValueError('原三角域文件的方法或图种类无效')
    levels = np.asarray(data['levels_pa'], dtype=float)
    peak = data['source_peak_pressure_pa']
    if (levels.ndim != 1 or len(levels) < 2 or not np.isfinite(levels).all() or
            levels[0] != 0 or np.any(np.diff(levels) <= 0) or
            not np.allclose(np.diff(levels), np.diff(levels)[0], rtol=1e-12, atol=0) or
            isinstance(peak, bool) or not isinstance(peak, (int, float)) or
            not math.isfinite(peak) or peak < 0 or levels[-1] < peak or
            len(data['superlevel_geometries']) != len(levels)-1):
        raise ValueError('原三角域文件的线性色阶、原峰或超水平域数量不一致')
    for key in ('outline_geometry', 'unknown_geometry', 'certified_zero_geometry',
            'positive_support_geometry', 'unclassified_geometry'):
        data[key] = _area_geometry(shape(data[key]), key, allow_empty=key != 'outline_geometry')
    data['superlevel_geometries'] = [_area_geometry(shape(g), '压力超水平域') for g in data['superlevel_geometries']]
    data['levels_pa'] = levels
    for field in data.get('source_fields', []):
        field['points_width_height_m'] = np.asarray(field['points_width_height_m'], dtype=float)
        field['triangles'] = np.asarray(field['triangles'], dtype=int)
        field['pressure_pa'] = np.asarray(field['pressure_pa'], dtype=float)
        if _field_fingerprint(field['points_width_height_m'], field['triangles'], field['pressure_pa']) != field['numeric_sha256']:
            raise ValueError('所附原P1坐标/连接/压力数值摘要已变化')
    if recompute:
        if not data.get('source_fields'):
            raise ValueError('该文件未包含原场，不能声称已完全重算')
        return build_triangle_pressure_map(data['source_fields'], method=data['method'], kind=data['kind'],
            outline_geometry=data['outline_geometry'], levels_pa=data['levels_pa'],
            unknown_geometry=data['unknown_geometry'], certified_zero_geometry=data['certified_zero_geometry'],
            source_identity=data['source_identity'])
    return data


def _draw_geometry(ax, geometry, **style):
    """按原环朝向绘制含孔面域；不把孔洞作为独立实心面片。"""
    import numpy as np
    from matplotlib.path import Path as PlotPath
    from matplotlib.patches import PathPatch
    from shapely.geometry.polygon import orient
    polygons = [geometry] if geometry.geom_type == 'Polygon' else [g for g in getattr(geometry, 'geoms', []) if g.geom_type == 'Polygon']
    for polygon in polygons:
        polygon = orient(polygon, sign=1.)
        vertices = []; codes = []
        for ring in (polygon.exterior, *polygon.interiors):
            xy = np.asarray(ring.coords)*1000.
            vertices.extend(xy.tolist())
            codes.extend([PlotPath.MOVETO]+[PlotPath.LINETO]*(len(xy)-2)+[PlotPath.CLOSEPOLY])
        ax.add_patch(PathPatch(PlotPath(vertices, codes), **style))


def plot_contact_composite(rows, output, *, lang='zh', stem='17_contact_composite',
        display_scope='certified', comparison_evidence=None, export_sources=True,source_compression='plain',
        instantaneous_zoom=False):
    """两行方法×两列瞬时/离散最大压力；四图共用同一绝对线性色标。

    rows每项含method、instantaneous和cycle_maximum压力图数据。
    certified模式用原未知面域遮罩；computed模式只在真实已计算原三角域
    保留颜色并加透明网纹，不把未认证压力声称为物理真值或压力上界。
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Patch
    from .publication_style import contact_pressure_colormap
    from .workflow_plotting import _paper, _save, _text
    if display_scope not in ('certified', 'computed') or not isinstance(rows, (list, tuple)) or len(rows) != 2:
        raise ValueError('需要两种方法及明确的certified/computed显示范围')
    if source_compression not in ('plain','gzip'):
        raise ValueError('源压力图数据压缩方式须为plain或gzip')
    if not isinstance(instantaneous_zoom,bool):
        raise ValueError('瞬时压力区细节开关须为布尔值')
    if {row.get('method') for row in rows} != set(_METHODS):
        raise ValueError('四面板必须同时包含真实FE和半解析方法')
    rows = sorted(rows, key=lambda row: _METHODS.index(row['method']))
    maps = []
    for row in rows:
        for key, kind in (('instantaneous', 'instantaneous'), ('cycle_maximum', 'sampled_cycle_maximum')):
            data = row[key]
            if data.get('schema') != _SCHEMA or data['method'] != row['method'] or data['kind'] != kind:
                raise ValueError('四面板的数据方法或瞬时/全过程语义不一致')
            maps.append(data)
    levels = np.asarray(maps[0]['levels_pa'])
    if any(not np.array_equal(levels, data['levels_pa']) for data in maps):
        raise ValueError('四面板必须使用同一绝对Pa色阶，不能按各图峰值单独变色')
    if any(not maps[0]['outline_geometry'].equals(data['outline_geometry']) for data in maps[1:]):
        raise ValueError('四面板必须表示同一参考材料工作齿面')
    norm = Normalize(0., float(levels[-1])*1e-6)
    cmap = contact_pressure_colormap()
    source_rows = []; source_files = []; summaries = [];zoom_records=[]
    with _paper(lang) as plt:
        fig, axes = plt.subplots(2, 2)
        if instantaneous_zoom:
            # 细节窗放在行间/图底留白，主图完整齿面和真实长宽比例均保留。
            fig.set_size_inches(190/25.4,148/25.4)
            fig.subplots_adjust(left=.075,right=.865,bottom=.30,top=.86,wspace=.32,hspace=1.5)
        else:
            fig.subplots_adjust(left=.075, right=.865, bottom=.20, top=.86, wspace=.32, hspace=.85)
        for index, data in enumerate(maps):
            ax = axes.ravel()[index]
            outline = data['outline_geometry']; unknown = data['unknown_geometry']
            _draw_geometry(ax, outline, facecolor='#F3F4F5', edgecolor='#626B72', linewidth=.5, zorder=0)
            _draw_geometry(ax, data['certified_zero_geometry'], facecolor=cmap(0.), edgecolor='none', zorder=1, antialiased=False)
            for level_index, geometry in enumerate(data['superlevel_geometries']):
                colour = cmap(norm((levels[level_index]+levels[level_index+1])*.5e-6))
                _draw_geometry(ax, geometry, facecolor=colour, edgecolor='none', zorder=2, rasterized=True, antialiased=False)
            _draw_geometry(ax, unknown, facecolor='#D5DADD' if display_scope == 'certified' else 'none',
                edgecolor='#89959C', hatch='////', linewidth=0., zorder=3)
            xmin, ymin, xmax, ymax = outline.bounds
            ax.set(xlim=(xmin*1000-.1, xmax*1000+.1), ylim=(ymin*1000-.1, ymax*1000+.1))
            ax.set_aspect('equal')
            method_name = _text(lang, '整轮有限元', 'Full-gear FE') if data['method'] == 'full_gear_fe' else _text(lang, '局部半解析', 'Local semi-analytical')
            kind_name = _text(lang, '瞬时压力', 'Instantaneous pressure') if data['kind'] == 'instantaneous' else _text(lang, '采样周期最大压力', 'Sampled cycle maximum')
            # 周期标题直接统计真实源场相位；没有角度元数据时不猜测采样数。
            if data['kind'] == 'sampled_cycle_maximum':
                angles = [field.get('driver_angle_rad') for field in data['source_fields']]
                if angles and all(isinstance(value, (int, float)) and not isinstance(value, bool)
                        and math.isfinite(value) for value in angles):
                    kind_name += f' (N={len(set(angles))})'
            ax.set_title(method_name+'\n'+kind_name, loc='left', fontsize=8, pad=8)
            ax.annotate('abcd'[index], (0., 1.), xycoords='axes fraction', xytext=(-8, 10), textcoords='offset points',
                ha='right', weight='bold', fontsize=9, annotation_clip=False)
            ax.set_xlabel(_text(lang, r'齿宽位置 $w_\mathrm{f}$ (mm)', r'Face-width position, $w_\mathrm{f}$ (mm)'))
            ax.set_ylabel(_text(lang, '法向高度 h (mm)', 'Normal height, h (mm)'))
            ax.text(.98, 1.02, rf'$p_\mathrm{{max}}={data["source_peak_pressure_pa"]*1e-6:.4g}$ MPa',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=6)
            name=data['method']+'_'+data['kind']+'_triangle_pressure.json'
            if source_compression=='gzip':
                name+='.gz'
            source_files.append(export_triangle_pressure_map(data, Path(output)/'data'/name, include_source_fields=export_sources))
            summary={key:data[key] for key in ('method', 'kind', 'source_field_count', 'source_peak_pressure_pa',
                'positive_support_projected_area_m2', 'unknown_projected_area_fraction', 'unknown_positive_support_area_m2')}
            source_rows.append(summary);summaries.append(summary)
        if instantaneous_zoom:
            from matplotlib.patches import Rectangle
            from matplotlib.ticker import MaxNLocator,ScalarFormatter
            fig.canvas.draw()
            for row_index in range(2):
                data=maps[row_index*2];support=data['positive_support_geometry']
                if support.is_empty:
                    zoom_records.append(dict(method=data['method'],available=False,reason='没有原正压支撑，不生成拟合斑点'))
                    continue
                xmin,ymin,xmax,ymax=support.bounds
                padx=(xmax-xmin)*.10;pady=(ymax-ymin)*.20
                bounds=[xmin-padx,ymin-pady,xmax+padx,ymax+pady]
                parent=axes[row_index,0];position=parent.get_position()
                detail=fig.add_axes([position.x0+.14*position.width,.54 if row_index==0 else .12,
                    position.width*.72,.081])
                _draw_geometry(detail,data['outline_geometry'],facecolor='#F3F4F5',edgecolor='#626B72',linewidth=.4,zorder=0)
                _draw_geometry(detail,data['certified_zero_geometry'],facecolor=cmap(0.),edgecolor='none',antialiased=False,zorder=1)
                for j,geometry in enumerate(data['superlevel_geometries']):
                    _draw_geometry(detail,geometry,facecolor=cmap(norm((levels[j]+levels[j+1])*.5e-6)),
                        edgecolor='none',antialiased=False,rasterized=True,zorder=2)
                _draw_geometry(detail,data['unknown_geometry'],facecolor='#D5DADD' if display_scope=='certified' else 'none',
                    edgecolor='#89959C',hatch='////',linewidth=0.,zorder=3)
                detail.set(xlim=(bounds[0]*1000,bounds[2]*1000),ylim=(bounds[1]*1000,bounds[3]*1000))
                detail.set_aspect('equal',adjustable='box');detail.tick_params(labelsize=5.5)
                detail.set_xlabel(r'$w_\mathrm{f}$ (mm)',fontsize=6,labelpad=1)
                detail.set_ylabel('h (mm)',fontsize=6,labelpad=1)
                detail.set_title(_text(lang,'压力区细节（真实比例）','Pressure-region detail (equal aspect)'),fontsize=6.5,pad=3)
                for axis in (detail.xaxis,detail.yaxis):
                    axis.set_major_locator(MaxNLocator(nbins=3));axis.set_major_formatter(ScalarFormatter(useOffset=False))
                parent.add_patch(Rectangle((bounds[0]*1000,bounds[1]*1000),(bounds[2]-bounds[0])*1000,
                    (bounds[3]-bounds[1])*1000,fill=False,edgecolor='#40484D',linestyle='--',linewidth=.5,zorder=4))
                zoom_records.append(dict(method=data['method'],available=True,bounds_width_height_m=bounds,
                    padding_fraction=[.10,.20],geometry_aspect='equal',colour_scale_unchanged=True,
                    pressure_values_modified=False,location='完整主齿面之外的图内留白'))
        bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap), cax=fig.add_axes([.91,.27,.018,.46]))
        bar.set_label(_text(lang, '接触压力 p (MPa)', 'Contact pressure p (MPa)'))
        # 图例仅说明四图中实际具有正面积的类别，不为不存在的未知域占位。
        handles=[];legend_categories=[]
        if any(data['certified_zero_geometry'].area>0 for data in maps):
            handles.append(Patch(facecolor=cmap(0.),label=_text(lang,r'零压区 $p=0$',r'Zero-pressure region $p=0$')))
            legend_categories.append('zero_pressure')
        if any(data['unknown_geometry'].area>0 for data in maps):
            handles.append(Patch(facecolor='#D5DADD' if display_scope=='certified' else 'white',edgecolor='#89959C',hatch='////',
                label=_text(lang,'未认证域','Uncertified domain')))
            legend_categories.append('uncertified')
        if any(data['unclassified_geometry'].area > max(data.get('coverage_geometry_roundoff_tolerance_m2', 0.), 256*math.ulp(1.)*data['outline_geometry'].area) for data in maps):
            handles.append(Patch(facecolor='#F3F4F5',edgecolor='#A5ADB2',label=_text(lang,'未给定压力/零压证据','No pressure/zero evidence')))
            legend_categories.append('unclassified')
        if handles:
            fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.47,.025),ncol=min(3,len(handles)),fontsize=6)
        if not instantaneous_zoom:
            fig.text(.075, .105, _text(lang,
                'N为实际求解相位数；周期图取离散采样最大值，未进行相位插值。',
                'N counts solved phases; cycle maxima use discrete samples, without phase interpolation.'),
                fontsize=6.1, color='#40474C', va='bottom')
        title=_text(lang,'同一材料齿面：瞬时与采样周期最大压力',
            'One material flank: instantaneous and sampled maximum contact pressure')
        conditions=[]
        for key in ('instantaneous_phase_rad','signed_pinion_torque_nm'):
            if comparison_evidence is not None and key in comparison_evidence:
                value=comparison_evidence[key]
                if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                    raise ValueError('比较工况的相位和有符号扭矩须为有限数值')
                conditions.append(rf'$\theta_1={value*180/math.pi:.3g}^\circ$' if key=='instantaneous_phase_rad'
                    else rf'$T_1={value:g}$ N·m')
        fig.suptitle(title+('\n'+', '.join(conditions) if conditions else ''),fontsize=9,y=.98)
        caption=_text(lang,
            '行分别为整轮有限元与半解析，列分别为同一相位同一物理齿的瞬时压力和离散全过程逐点最大压力。'
            '四图共用绝对线性MPa色标，保持真实齿宽/齿高比例。原P1三角域按压力阈值精确裁切后并集，'
            '不经笛卡尔角点掩码、不填孔、不模糊、不拟合斑点。零压蓝仅来自明确零压域，未认证域加网纹；'
            '未提供数值或零压证据的区域为中性色。源三角连接和压力保留；最大值不相加，也不是同时受载面积。',
            'Rows compare full-gear FE and semi-analytical results; columns show one phase/tooth and the '
            'pointwise maximum over sampled phases. All four panels share an absolute linear MPa scale '
            'and true coordinate aspect. Original P1 triangle superlevel sets are unioned without Cartesian '
            'corner masking, hole filling, smoothing or footprint fitting. Certified zero is blue; '
            'uncertified domains are hatched; regions without pressure/zero evidence are neutral. '
            'Original triangle connectivity and pressures are retained. The maximum is neither a pressure sum nor simultaneous loaded area.')
        if display_scope=='computed':
            caption+=_text(lang,' 本图显式显示已计算压力；网纹处不声称压力最大域已认证，更不构成物理验证。',
                ' Computed pressures remain visible beneath hatching; these regions are not certified maxima or physical validation.')
        if instantaneous_zoom:
            caption+=_text(lang,' 瞬时图另附原正压支撑包围盒加边距的等比例细节窗，坐标为真实mm；放大不增加求解分辨率、不改变压力域或色标。',
                ' Instantaneous details use the original positive-support bounding box with padding and true mm coordinates; zoom adds no solver resolution and changes neither pressure domain nor color scale.')
        return _save(fig,output,stem,lang,source_rows,caption,width_mm=190,height_mm=148 if instantaneous_zoom else 108,
            extra=dict(schema='bevel_contact_composite_v1',display_scope=display_scope,models=summaries,
                pressure_limit_mpa=float(levels[-1])*1e-6,colour_mapping='linear',colour_scale_shared_by_all_panels=True,
                original_pressure_values_modified=False,query_grid_resampling=False,geometry_aspect='equal',
                instantaneous_zoom=instantaneous_zoom,instantaneous_zoom_windows=zoom_records,
                legend_categories=legend_categories,
                source_pressure_maps=source_files,source_compression=source_compression,
                comparison_evidence={} if comparison_evidence is None else comparison_evidence))
