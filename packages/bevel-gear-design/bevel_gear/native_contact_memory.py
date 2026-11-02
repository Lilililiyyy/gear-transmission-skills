"""接触专用内存策略：消费已prepare系统的组装缓存，保留原K／原因子。

不提供FE单元应力；几何、W、QP、节点响应和标准强度业务不变。
"""
from __future__ import annotations
from copy import deepcopy
import gc
from types import MethodType
import weakref

MEMBERS = ('pinion', 'wheel')


RETAINED = ('mesh', 'material', 'constitutive', 'nodes', 'tetrahedra', 'stiffness',
    'fixed_dofs', 'free_dofs', 'factorization')


CACHES = ('b', 'dofs', 'volumes')


def _resource():
    try:import psutil
    except ImportError as error:raise ImportError('接触专用内存计量需要可选psutil依赖') from error
    process = psutil.Process(); memory = process.memory_info(); io = process.io_counters()
    return dict(rss_bytes=memory.rss, peak_rss_bytes=getattr(memory, 'peak_wset', memory.rss),
        host_available_bytes=psutil.virtual_memory().available,
        process_read_bytes=io.read_bytes, process_write_bytes=io.write_bytes)


def _stress_solve_unavailable(self, nodal_forces_n):
    raise NotImplementedError('此对象已显式释放b/dofs/volumes，仅供原生接触；不支持单元应力求解。'
        '如需应力须另建完整ElasticSystem或另行审核按原公式重算，不能调用本对象冒充完整应力输出')


class LeanContactElasticSystem:
    """接触专用对象，明确不是ElasticSystem；保留K、网格、材料和原因子引用。"""
    def __init__(self, full):
        for name in RETAINED:
            setattr(self, name, getattr(full, name))
        self.diagnostics = deepcopy(full.diagnostics)
        self.diagnostics['lean_contact_contract'] = dict(element_stress_available=False,
            released_assembly_cache_names=list(CACHES), original_K_preserved=True)

    solve = _stress_solve_unavailable

    @property
    def b(self):
        raise AttributeError('接触精简对象已释放单元B矩阵，不支持单元应力接口')

    @property
    def dofs(self):
        raise AttributeError('接触精简对象已释放单元DOF映射；free_dofs/fixed_dofs仍保留')

    @property
    def volumes(self):
        raise AttributeError('接触精简对象已释放逐单元体积；已核验总体积诊断仍保留')


def release_prepared_assembly_caches(prepared, *, journal=None, require_released=True):
    """消费prepared内完整系统，释放b/dofs/volumes；保留原K/因子/mesh等引用。

    只在原ElasticSystem初始化及diagnostics完成之后调用。旧系统外部别名
    也被标记为已消费且solve明确拒绝应力；外部若另持缓存数组，则记录并拒绝
    宣称完成释放。不得同时把旧系统当完整ElasticSystem继续使用。
    """
    import numpy as np
    from bevel_gear.elasticity import ElasticSystem
    if set(prepared.systems) != set(MEMBERS): raise ValueError('须提供完整pinion/wheel准备对象')
    if not isinstance(require_released, bool): raise ValueError('require_released须为布尔值')
    before = _resource(); replacements = {}; proofs = {}; roots = []
    for member in MEMBERS:
        full = prepared.systems[member]
        if not isinstance(full, ElasticSystem) or getattr(full, '_lean_caches_consumed', False):
            raise TypeError('须显式消费完整且未被精简的ElasticSystem：'+member)
        elements = len(full.tetrahedra); info = {}
        for name, shape in (('b', (elements, 6, 12)), ('dofs', (elements, 12)), ('volumes', (elements,))):
            array = getattr(full, name)
            if not isinstance(array, np.ndarray) or array.shape != shape:
                raise ValueError('原组装缓存布局异常：'+member+'.'+name)
            root = array
            while isinstance(root.base, np.ndarray): root = root.base
            roots.append((member, name, weakref.ref(root)))
            info[name] = dict(shape=list(array.shape), dtype=array.dtype.name, bytes=int(array.nbytes))
        del array, root
        replacements[member] = LeanContactElasticSystem(full)
        proofs[member] = dict(caches=info, payload_bytes=sum(row['bytes'] for row in info.values()),
            retained_original_K_object=replacements[member].stiffness is full.stiffness,
            retained_original_factor_object=replacements[member].factorization is full.factorization,
            element_stress_available=False)
    # 两体均通过布局预检后才提交，避免检查第二体时已先删除第一体。
    for member in MEMBERS:
        full = prepared.systems[member]
        for name in CACHES: delattr(full, name)
        full._lean_caches_consumed = True
        full.solve = MethodType(_stress_solve_unavailable, full)
        prepared.systems[member] = replacements[member]
    del full; gc.collect()
    retained = [dict(member=member, cache=name) for member, name, ref in roots if ref() is not None]
    receipt = dict(schema='bevel_native_contact_memory_release_v1', before=before, after=_resource(),
        members=proofs, removed_array_payload_bytes=sum(p['payload_bytes'] for p in proofs.values()),
        retained_external_cache_buffers=retained, complete_cache_release=not retained,
        K_or_factor_values_modified=False, W_or_QP_or_stress_formula_modified=False,
        original_systems_consumed=True, element_stress_available=False,
        note='payload为已移除数组字节；实际RSS/文件缓存变化须实测，不宣称等额释放或加速')
    if journal: journal.emit('assembly_cache_release', receipt=receipt)
    if require_released and retained:
        error = RuntimeError('组装缓存仍有外部强引用，不能宣称精简完成；准备对象已消费，必须保留证据并按原流程释放因子')
        error.lean_release_receipt = receipt
        raise error
    return receipt

