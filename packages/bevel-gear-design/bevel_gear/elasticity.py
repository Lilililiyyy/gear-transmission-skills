"""新写三维小变形线弹性计算：四节点四面体、稀疏求解与接触柔度。

本模块从线弹性弱式实现，不以经验啮合刚度替代实体柔度。位移m、
节点力N、材料Pa，工程剪应变顺序xx/yy/zz/xy/yz/zx。仅实际调用时
加载NumPy/SciPy；中性CAD网格划分另外按需加载Gmsh。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
import uuid


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name}须为正有限数")
    return float(value)


def _numeric():
    try:
        import numpy as np
        from scipy import sparse
        from scipy.sparse.linalg import splu
    except ImportError as error:
        raise ImportError("三维弹性计算需要可选NumPy/SciPy依赖") from error
    return np, sparse, splu


@dataclass(frozen=True)
class ElasticMaterial:
    """各向同性线弹性；近不可压缩域明确拒绝，避免低阶单元体积锁死。"""
    young_modulus_pa: float
    poisson_ratio: float

    def __post_init__(self):
        _positive(self.young_modulus_pa, "弹性模量Pa")
        v = self.poisson_ratio
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not -1 < v < .49:
            raise ValueError("泊松比须位于(-1,0.49)；本低阶位移单元不处理近不可压缩材料")


@dataclass
class TetraMesh:
    """SI节点、0起始四面体与边界三角形；boundary_tags为网格面实体编号。"""
    nodes_m: object
    tetrahedra: object
    boundary_triangles: object = field(default_factory=tuple)
    boundary_tags: object = field(default_factory=tuple)
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ContactLoadChannel:
    """单位法向总力的节点分配；权重和为1，使位移平均与载荷功共轭。"""
    node_indices: tuple
    weights: tuple
    normal: tuple


def contact_elasticity_factor(pinion_material, wheel_material):
    """由同一对E/ν求sqrt(E*/π)，供标准接触应力与实体材料一致使用。

    返回√Pa，1/E*=(1−ν1²)/E1+(1−ν2²)/E2；不覆盖调用方明确给定
    的标准系数，也不改变既有评级接口的兼容默认值。
    """
    values=[]
    for material in (pinion_material,wheel_material):
        if isinstance(material,dict):material=ElasticMaterial(**material)
        if not isinstance(material,ElasticMaterial):raise TypeError("需要两轮ElasticMaterial或具名字典")
        values.append((1-material.poisson_ratio**2)/material.young_modulus_pa)
    return math.sqrt(1/(math.pi*sum(values)))


def _constitutive(material, np):
    if isinstance(material, dict):
        material = ElasticMaterial(**material)
    if not isinstance(material, ElasticMaterial):
        raise TypeError("材料须为ElasticMaterial或具名字典")
    e, v = material.young_modulus_pa, material.poisson_ratio
    mu = e / (2 * (1 + v))
    lam = e * v / ((1 + v) * (1 - 2 * v))
    d = np.zeros((6, 6))
    d[:3, :3] = lam
    for i in range(3):
        d[i, i] += 2 * mu
        d[i + 3, i + 3] = mu
    return material, d


def _element_kinematics(points, np):
    # 在平移后的局部坐标计算梯度，避免大型装配坐标污染条件数。
    edges = points[:, 1:] - points[:, :1]
    matrices = np.transpose(edges, (0, 2, 1))
    determinants = np.linalg.det(matrices)
    size = np.max(np.linalg.norm(edges, axis=2), axis=1)
    if np.any(np.abs(determinants) <= 1e-12 * size**3):
        raise ValueError("存在退化或极薄四面体，须重新划分实体网格")
    inverse = np.linalg.inv(matrices)
    gradients = np.empty((len(points), 4, 3))
    gradients[:, 1:, :] = inverse
    gradients[:, 0, :] = -inverse.sum(axis=1)
    b = np.zeros((len(points), 6, 12))
    for i in range(4):
        x, y, z = gradients[:, i, :].T
        j = 3 * i
        b[:, 0, j] = x
        b[:, 1, j + 1] = y
        b[:, 2, j + 2] = z
        b[:, 3, j] = y
        b[:, 3, j + 1] = x
        b[:, 4, j + 1] = z
        b[:, 4, j + 2] = y
        b[:, 5, j] = z
        b[:, 5, j + 2] = x
    return b, np.abs(determinants) / 6


class ElasticSystem:
    """已组装并固定边界的系统；重复载荷复用同一稀疏分解。"""

    def __init__(self, mesh, material, fixed_nodes, *, maximum_dofs=120000,
            factorization_ordering='symmetric',factorization_factory=None):
        np, sparse, splu = _numeric()
        if factorization_ordering not in ('symmetric','colamd'):
            raise ValueError('稀疏排序须为symmetric或colamd；不改变弹性刚度矩阵')
        if factorization_factory is not None and not callable(factorization_factory):
            raise ValueError('自定义分解工厂须为接收Kff的可调用对象')
        if not isinstance(mesh, TetraMesh):
            raise TypeError("实体网格须为TetraMesh")
        nodes = np.asarray(mesh.nodes_m, dtype=float)
        raw = np.asarray(mesh.tetrahedra)
        if nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) < 4 or not np.isfinite(nodes).all():
            raise ValueError("节点须为至少四个有限三维SI坐标")
        if raw.ndim != 2 or raw.shape[1] != 4 or len(raw) == 0 or raw.dtype.kind not in "iu":
            raise ValueError("四面体须为非空N×4整数索引，不能截断浮点索引")
        tetra = raw.astype(int)
        if tetra.min() < 0 or tetra.max() >= len(nodes):
            raise ValueError("四面体节点索引越界")
        if len(np.unique(tetra)) != len(nodes):
            raise ValueError("存在未被实体单元使用的孤立节点")
        fixed_raw = np.asarray(fixed_nodes)
        if fixed_raw.ndim != 1 or fixed_raw.dtype.kind not in "iu":
            raise ValueError("固定节点须为一维整数索引")
        fixed = np.unique(fixed_raw.astype(int))
        if len(fixed) < 3 or fixed.min() < 0 or fixed.max() >= len(nodes):
            raise ValueError("至少须固定三个合法且不共线的节点")
        if np.linalg.matrix_rank(nodes[fixed] - nodes[fixed[0]]) < 2:
            raise ValueError("固定节点共线，仍存在刚体转动自由度")
        ndof = 3 * len(nodes)
        if isinstance(maximum_dofs, bool) or not isinstance(maximum_dofs, int) or not 12 <= maximum_dofs or ndof > maximum_dofs:
            raise ValueError("自由度超过显式规模限制")
        self.mesh = mesh
        self.material, self.constitutive = _constitutive(material, np)
        self.nodes, self.tetrahedra = nodes, tetra
        self.b, self.volumes = _element_kinematics(nodes[tetra], np)
        self.dofs = (3 * tetra[:, :, None] + np.arange(3)).reshape((-1, 12))
        local = np.einsum("eai,ab,ebj,e->eij", self.b, self.constitutive, self.b, self.volumes, optimize=True)
        rows = np.repeat(self.dofs, 12, axis=1).ravel()
        cols = np.tile(self.dofs, (1, 12)).ravel()
        self.stiffness = sparse.coo_matrix((local.ravel(), (rows, cols)), shape=(ndof, ndof)).tocsc()
        # 全局K已独立组装；单元K及重复COO索引不再使用。先释放，避免
        # 百万单元的临时数组与大型稀疏因子同时占用数十GB内存。
        del local,rows,cols
        self.fixed_dofs = (3 * fixed[:, None] + np.arange(3)).ravel()
        self.free_dofs = np.setdiff1d(np.arange(ndof), self.fixed_dofs)
        if len(self.free_dofs) == 0:
            raise ValueError("全部自由度均被固定，不能计算接触柔度")
        try:
            options=(dict(permc_spec='MMD_AT_PLUS_A',options={'SymmetricMode':True})
                if factorization_ordering=='symmetric' else dict(permc_spec='COLAMD'))
            # 线弹性Kff对称；只优化稀疏消元排序，不改变材料、单元、边界或载荷。
            constrained=self.stiffness[self.free_dofs][:, self.free_dofs]
            self.factorization = (splu(constrained,**options) if factorization_factory is None
                else factorization_factory(constrained))
        except RuntimeError as error:
            raise ArithmeticError("受约束弹性矩阵奇异；检查断开实体及每个实体的边界约束") from error
        relative_pivot=None;factor_diagnostics={}
        if factorization_factory is None:
            pivots = np.abs(self.factorization.U.diagonal())
            if pivots.min() <= pivots.max() * 1e-13:
                raise ArithmeticError("弹性矩阵病态或仍存在未约束刚体模态")
            relative_pivot=float(pivots.min()/pivots.max())
        else:
            # Cholesky没有LU的U主元；要求真实完成SPD分解，不伪造LU诊断。
            factor_diagnostics=getattr(self.factorization,'analysis_diagnostics',None)
            if (not callable(getattr(self.factorization,'solve',None)) or
                    getattr(self.factorization,'positive_definite',None) is not True or
                    not isinstance(factor_diagnostics,dict) or not factor_diagnostics):
                raise ValueError('自定义分解须提供solve、完成SPD分解的positive_definite=True和analysis_diagnostics')
        self.diagnostics = dict(nodes=len(nodes), tetrahedra=len(tetra), dofs=ndof,
            free_dofs=len(self.free_dofs), fixed_nodes=len(fixed), volume_m3=float(self.volumes.sum()),
            minimum_element_volume_m3=float(self.volumes.min()),
            model="三维各向同性小变形线弹性；四节点常应变四面体",
            boundary="所列节点的三个平移自由度均固定；不包含轴/轴承柔度",
            units="m/N/Pa", minimum_relative_lu_pivot=relative_pivot,
            sparse_factorization_ordering=factorization_ordering if factorization_factory is None else 'external_factory',
            sparse_factorization_diagnostics=dict(factor_diagnostics))

    def solve(self, nodal_forces_n):
        """求位移、单元应力和支反力；不将单元常应力当齿根收敛峰值。"""
        np, _, _ = _numeric()
        force = np.asarray(nodal_forces_n, dtype=float)
        if force.size != self.stiffness.shape[0] or not np.isfinite(force).all():
            raise ValueError("节点力须覆盖所有三平移自由度且均有限，单位N")
        force = force.reshape(-1)
        displacement = np.zeros_like(force)
        displacement[self.free_dofs] = self.factorization.solve(force[self.free_dofs])
        residual = self.stiffness @ displacement - force
        relative_residual=float(np.linalg.norm(residual[self.free_dofs])/max(np.linalg.norm(force[self.free_dofs]),1e-30))
        if not np.isfinite(displacement).all() or relative_residual>1e-8:
            raise ArithmeticError('弹性方程响应非有限或自由度相对残差超限')
        stress = np.einsum("ab,ebi,ei->ea", self.constitutive, self.b, displacement[self.dofs])
        sx, sy, sz, xy, yz, zx = stress.T
        mises = np.sqrt(.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2) + 3*(xy**2+yz**2+zx**2))
        return dict(displacements_m=displacement.reshape((-1,3)).tolist(),
            element_stress_pa=stress.tolist(), element_von_mises_pa=mises.tolist(),
            reaction_forces_n=residual.reshape((-1,3)).tolist(),
            maximum_free_residual_n=float(np.abs(residual[self.free_dofs]).max()),
            relative_free_residual=relative_residual,
            total_force_balance_n=(residual.reshape((-1,3)).sum(axis=0)+force.reshape((-1,3)).sum(axis=0)).tolist(),
            strain_energy_j=float(.5 * displacement @ force),
            stress_definition="四节点单元内常应力；未作节点外推、热点拟合或网格收敛认证")


def contact_load_matrix(system, channels, *, matrix_format='dense',maximum_sparse_channels=2048):
    """构造单位总法向力列矩阵W；节点位移由Wᵀ投影，保持功共轭。"""
    np, _, _ = _numeric()
    channels = list(channels)
    if matrix_format not in ('dense', 'csc'):
        raise ValueError("单位力矩阵格式须为dense或csc")
    if isinstance(maximum_sparse_channels,bool) or not isinstance(maximum_sparse_channels,int) or not 1<=maximum_sparse_channels<=100000:
        raise ValueError('稀疏载荷通道上限须为1至100000整数，需由调用方显式控制')
    maximum = 512 if matrix_format == 'dense' else maximum_sparse_channels
    if not channels or len(channels) > maximum:
        raise ValueError(f"所选单位力矩阵格式的通道数须为1至{maximum}")
    row_indices, column_indices, values = [], [], []
    for index, channel in enumerate(channels):
        if isinstance(channel, dict):
            channel = ContactLoadChannel(**channel)
        if not isinstance(channel, ContactLoadChannel):
            raise TypeError("接触通道须为ContactLoadChannel")
        ids = np.asarray(channel.node_indices)
        weights = np.asarray(channel.weights, dtype=float)
        normal = np.asarray(channel.normal, dtype=float)
        if ids.ndim != 1 or ids.dtype.kind not in "iu" or len(ids) == 0 or len(np.unique(ids)) != len(ids) or ids.min() < 0 or ids.max() >= len(system.nodes):
            raise ValueError("接触通道节点索引非法或重复")
        if weights.shape != ids.shape or not np.isfinite(weights).all() or np.any(weights < 0) or abs(weights.sum()-1) > 1e-12:
            raise ValueError("接触节点权重须非负且总和为1")
        if normal.shape != (3,) or not np.isfinite(normal).all() or abs(np.linalg.norm(normal)-1) > 1e-10:
            raise ValueError("接触通道法向须为单位三维向量")
        for node, weight in zip(ids, weights):
            for component in range(3):
                row_indices.append(3 * int(node) + component)
                column_indices.append(index)
                values.append(float(weight * normal[component]))
    from scipy.sparse import csc_matrix
    sparse = csc_matrix((values, (row_indices, column_indices)),
                        shape=(system.stiffness.shape[0], len(channels)))
    if np.any(np.asarray(sparse[system.free_dofs].power(2).sum(axis=0)).ravel() == 0):
        raise ValueError("接触载荷全部施加在固定自由度，不能形成有效柔度")
    return sparse if matrix_format == 'csc' else sparse.toarray()


def condense_contact_flexibility(system, channels, *, require_positive_definite=True):
    """C=WᵀKff⁻¹W，单位m/N；保留耦合项，不增加经验对角刚度。"""
    np, _, _ = _numeric()
    w = contact_load_matrix(system, channels, matrix_format='csc')
    free = w[system.free_dofs]
    # 分批解Kff U=W，避免把自由度数×全部接触通道同时展开到内存。
    # 分批只改变线性代数存储，C仍保留每一对接触通道的完整耦合。
    count = free.shape[1]
    c = np.empty((count, count))
    batch_size = 32
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        displacements = system.factorization.solve(free[:, start:end].toarray())
        c[:, start:end] = free.T @ displacements
    asymmetry = float(np.abs(c-c.T).max())
    scale = float(np.abs(c).max())
    if not np.isfinite(c).all() or scale <= 0 or asymmetry > 1e-9 * scale:
        raise ArithmeticError("凝聚柔度无效或不满足Maxwell互等关系")
    # 稀疏分解的最后位误差显式对称化，不改负特征值或增加正则项。
    c = .5 * (c+c.T)
    eigenvalues = np.linalg.eigvalsh(c)
    if not isinstance(require_positive_definite,bool):
        raise ValueError("require_positive_definite须为布尔值")
    if eigenvalues[0] < -eigenvalues[-1]*1e-10:
        raise ArithmeticError("凝聚柔度出现超出舍入尺度的负特征值")
    if require_positive_definite and eigenvalues[0] <= eigenvalues[-1] * 1e-12:
        raise ArithmeticError("接触通道线性相关或柔度非正定，请减少重复通道")
    return dict(compliance_m_per_n=c.tolist(), channels=len(c), linear_solve_batch_size=batch_size,
        symmetry_absolute_error_m_per_n=asymmetry,
        eigenvalue_min_m_per_n=float(eigenvalues[0]), eigenvalue_max_m_per_n=float(eigenvalues[-1]),
        positive_definite=bool(eigenvalues[0]>eigenvalues[-1]*1e-12),
        definition="WᵀKff⁻¹W；单位力与平均位移功共轭，含非对角实体耦合",
        limitations=["单体W可能列相关，显式允许半正定时应对两体组合柔度另检验；不添加正则刚度", "离散网格与固定边界决定柔度；未作齿轮物理试验标定", "低阶四面体须加密检查，不能用粗网格应力峰值声明Hertz精度"])


def mesh_step_tetra(path, *, element_size_m, maximum_nodes=40000, output_mesh=None,
                    refinement_points_m=(), refined_size_m=None, refinement_radius_m=None,
                    refinement_polyline=False, refinement_polylines_m=(),refinement_core_radius_m=None,
                    refinement_regions=()):
    """用专业Gmsh读STEP并按文件单位转换为m；仅生成一次四面体网格。

    保留已有Gmsh会话与数值选项，不启动GUI。此入口要求闭合体，开放
    曲面壳明确失败。若需硬超时，应在有界子进程中调用。
    refinement_core_radius_m显式控制保持细单元的实体近场厚度；省略时
    保留原一个refined_size厚度。只缩小表面尺寸不代表接触深层同步收敛。
    refinement_regions可额外指定独立折线/尺寸/core/影响半径。各Distance+
    Threshold场与原背景经Gmsh Min组合；空列表保持原控制场和默认行为。
    """
    h = _positive(element_size_m, "网格目标尺寸m")
    if isinstance(maximum_nodes, bool) or not isinstance(maximum_nodes, int) or maximum_nodes < 4:
        raise ValueError("maximum_nodes须为不小于4的整数")
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError("待划分的STEP文件不存在")
    np, _, _ = _numeric()
    refine=np.asarray(refinement_points_m,dtype=float)
    tracks=[]
    for track in refinement_polylines_m:
        points=np.asarray(track,dtype=float)
        if points.ndim!=2 or points.shape[1]!=3 or len(points)<2 or not np.isfinite(points).all():
            raise ValueError("每条连续加密轨迹须为至少两个有限三维SI点")
        if np.any(np.linalg.norm(np.diff(points,axis=0),axis=1)<=1e-12):
            raise ValueError("连续加密轨迹相邻点不得重合")
        tracks.append(points)
    if tracks and refine.size:
        raise ValueError("独立加密点与多折线轨迹不能同时指定；避免混淆控制场")
    if tracks:
        refine=np.concatenate(tracks,axis=0)
    if not isinstance(refinement_polyline,bool):raise ValueError("refinement_polyline须为布尔值")
    if refine.size:
        if refine.ndim!=2 or refine.shape[1]!=3 or not np.isfinite(refine).all():
            raise ValueError("局部加密控制点须为有限N×3 SI坐标")
        fine=_positive(refined_size_m,"局部单元尺寸m")
        reach=_positive(refinement_radius_m,"局部加密影响半径m")
        core=fine if refinement_core_radius_m is None else _positive(refinement_core_radius_m,'近场细网格核心半径m')
        if fine>=h:raise ValueError("局部加密尺寸须小于主体单元尺寸")
        if reach<=fine:raise ValueError("局部影响半径须大于加密单元尺寸")
        if not fine<=core<reach:raise ValueError('细网格核心半径须不小于细单元尺寸且小于影响半径')
        if refinement_polyline and len(refine)<2:raise ValueError("连续加密折线至少需要两个点")
    elif refined_size_m is not None or refinement_radius_m is not None or refinement_core_radius_m is not None:
        raise ValueError("局部加密尺寸/半径必须同时提供非空控制点")
    regions=[]
    if not isinstance(refinement_regions,(list,tuple)):raise ValueError('附加加密区域须为列表或元组')
    for region in refinement_regions:
        allowed={'refinement_polylines_m','refined_size_m','refinement_core_radius_m','refinement_radius_m'}
        if not isinstance(region,dict) or set(region)-allowed:raise ValueError('附加加密区域含未知字段或不是字典')
        lines=[]
        for track in region.get('refinement_polylines_m',()):
            points=np.asarray(track,dtype=float)
            if points.ndim!=2 or points.shape[1]!=3 or len(points)<2 or not np.isfinite(points).all():
                raise ValueError('附加加密轨迹须有至少两个有限三维SI点')
            if np.any(np.linalg.norm(np.diff(points,axis=0),axis=1)<=1e-12):raise ValueError('附加加密轨迹相邻点不得重合')
            lines.append(points)
        if not lines:raise ValueError('附加加密区域必须有非空折线')
        size=_positive(region.get('refined_size_m'),'附加局部单元尺寸m')
        radius=_positive(region.get('refinement_radius_m'),'附加局部影响半径m')
        thickness=region.get('refinement_core_radius_m')
        thickness=size if thickness is None else _positive(thickness,'附加近场细网格核心半径m')
        if not size<h or not size<=thickness<radius:raise ValueError('附加区域须满足fine<bulk且fine≤core<radius')
        regions.append(dict(refinement_polylines_m=lines,refined_size_m=size,
            refinement_core_radius_m=thickness,refinement_radius_m=radius))
    try:
        import gmsh
    except ImportError as error:
        raise ImportError("实体四面体划分需要可选Gmsh依赖") from error
    owned = not gmsh.isInitialized()
    if owned:
        gmsh.initialize([], readConfigFiles=False, run=False)
    previous_model = gmsh.model.getCurrent()
    label = "elastic_" + uuid.uuid4().hex
    settings = {"General.Terminal":0., "Mesh.MeshSizeMin":min(h*.5,fine) if refine.size else h*.5, "Mesh.MeshSizeMax":h,
        "Mesh.ElementOrder":1., "Mesh.Algorithm3D":1., "Mesh.Optimize":1., "Mesh.MeshSizeFromCurvature":12.}
    if regions:settings['Mesh.MeshSizeMin']=min(settings['Mesh.MeshSizeMin'],*(r['refined_size_m'] for r in regions))
    previous_options = {k:gmsh.option.getNumber(k) for k in settings}
    previous_unit = gmsh.option.getString("Geometry.OCCTargetUnit")
    try:
        for key, value in settings.items():
            gmsh.option.setNumber(key, value)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.add(label)
        gmsh.model.occ.importShapes(str(source.resolve()))
        control_tags=[]
        control_curves=[]
        if refine.size:
            if tracks:
                for track in tracks:
                    tags=[gmsh.model.occ.addPoint(float(x),float(y),float(z)) for x,y,z in track]
                    control_tags.extend(tags)
                    control_curves.extend(gmsh.model.occ.addLine(a,b) for a,b in zip(tags[:-1],tags[1:]))
            else:
                control_tags=[gmsh.model.occ.addPoint(float(x),float(y),float(z)) for x,y,z in refine]
            if refinement_polyline and not tracks:
                if np.any(np.linalg.norm(np.diff(refine,axis=0),axis=1)<=1e-12):
                    raise ValueError("加密折线相邻控制点不得重合")
                control_curves=[gmsh.model.occ.addLine(a,b) for a,b in zip(control_tags[:-1],control_tags[1:])]
        extra_controls=[]
        for region in regions:
            curves=[]
            for track in region['refinement_polylines_m']:
                tags=[gmsh.model.occ.addPoint(float(x),float(y),float(z)) for x,y,z in track]
                curves.extend(gmsh.model.occ.addLine(a,b) for a,b in zip(tags[:-1],tags[1:]))
            extra_controls.append(curves)
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(3)
        if not volumes:
            raise ValueError("STEP没有可划分的闭合实体；开放壳不能用于实体弹性")
        thresholds=[]
        if control_tags:
            distance_field=gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(distance_field,"CurvesList" if control_curves else "PointsList",control_curves or control_tags)
            if control_curves:gmsh.model.mesh.field.setNumber(distance_field,"Sampling",100)
            threshold=gmsh.model.mesh.field.add("Threshold")
            for key,value in (("InField",distance_field),("SizeMin",fine),("SizeMax",h),("DistMin",core),("DistMax",reach)):
                gmsh.model.mesh.field.setNumber(threshold,key,value)
            gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
            thresholds.append(threshold)
        for region,curves in zip(regions,extra_controls):
            distance=gmsh.model.mesh.field.add('Distance')
            gmsh.model.mesh.field.setNumbers(distance,'CurvesList',curves)
            gmsh.model.mesh.field.setNumber(distance,'Sampling',100)
            threshold=gmsh.model.mesh.field.add('Threshold')
            for key,value in (('InField',distance),('SizeMin',region['refined_size_m']),('SizeMax',h),
                    ('DistMin',region['refinement_core_radius_m']),('DistMax',region['refinement_radius_m'])):
                gmsh.model.mesh.field.setNumber(threshold,key,value)
            thresholds.append(threshold)
        if regions:
            combined=thresholds[0]
            if len(thresholds)>1:
                combined=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(combined,'FieldsList',thresholds)
            gmsh.model.mesh.field.setAsBackgroundMesh(combined)
        gmsh.model.mesh.generate(3)
        tags, coordinates, _ = gmsh.model.mesh.getNodes()
        raw_node_count = len(tags)
        points = np.asarray(coordinates).reshape((-1,3))
        tag_to_index = {int(tag):i for i,tag in enumerate(tags)}
        _, tetra_nodes = gmsh.model.mesh.getElementsByType(4)
        if len(tetra_nodes) == 0:
            raise ArithmeticError("专业网格器未生成线性四面体")
        tetra = np.array([tag_to_index[int(t)] for t in tetra_nodes], dtype=int).reshape((-1,4))
        triangles, face_ids = [], []
        for _, tag in gmsh.model.getEntities(2):
            _, face_nodes = gmsh.model.mesh.getElementsByType(2, tag)
            mapped = [tag_to_index[int(t)] for t in face_nodes]
            triangles.extend([mapped[i:i+3] for i in range(0,len(mapped),3)])
            face_ids.extend([tag] * (len(mapped)//3))
        used = np.unique(tetra)
        if len(used) != len(points):
            # STEP孤立顶点不作为弹性自由度；保留实体单元覆盖的节点并重编索引。
            lookup = np.full(len(points), -1, dtype=int)
            lookup[used] = np.arange(len(used))
            tetra = lookup[tetra]
            mapped_triangles = lookup[np.asarray(triangles, dtype=int)]
            valid = np.all(mapped_triangles >= 0, axis=1)
            triangles = mapped_triangles[valid].tolist()
            face_ids = np.asarray(face_ids)[valid].tolist()
            points = points[used]
        # 距离场控制线也会被Gmsh离散，但未参与四面体的控制节点不是
        # 弹性自由度。上限约束实际实体节点，不能在剔除孤立控制节点前误拒绝。
        if len(points) > maximum_nodes:
            raise ValueError(f"实体网格生成{len(points)}节点，超过显式上限{maximum_nodes}；须检查小几何特征与局部网格控制")
        if output_mesh is not None:
            target = Path(output_mesh)
            target.parent.mkdir(parents=True, exist_ok=True)
            gmsh.write(str(target.resolve()))
        metadata=dict(units="m", mesh_engine="Gmsh", mesh_engine_version=gmsh.__version__,
                element="4-node tetrahedron", requested_element_size_m=h, volume_entities=len(volumes),
                raw_mesh_node_count=raw_node_count, removed_unused_node_count=raw_node_count-len(points),
                refinement_points_m=refine.tolist() if refine.size else [],
                refinement_polylines_m=[track.tolist() for track in tracks],
                refined_size_m=refined_size_m, refinement_radius_m=refinement_radius_m,
                refinement_core_radius_m=core if refine.size else None,
                refinement_geometry="continuous_polyline" if control_curves else "points",
                bounds_m=[points.min(axis=0).tolist(),points.max(axis=0).tolist()],
                scale="STEP文件单位由Geometry.OCCTargetUnit=M显式转换；不按坐标量级猜单位")
        if regions:
            metadata.update(refinement_regions=[dict(region,refinement_polylines_m=[track.tolist() for track in region['refinement_polylines_m']]) for region in regions],
                refinement_field_combination='Min原背景及全部独立Distance/Threshold；不放粗其它区域')
        return TetraMesh(points.tolist(),tetra.tolist(),triangles,face_ids,metadata)
    finally:
        if gmsh.model.getCurrent() == label:
            gmsh.model.remove()
        for key, value in previous_options.items():
            gmsh.option.setNumber(key, value)
        gmsh.option.setString("Geometry.OCCTargetUnit", previous_unit)
        if previous_model:
            gmsh.model.setCurrent(previous_model)
        if owned:
            gmsh.finalize()
