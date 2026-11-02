"""可选oneMKL PARDISO_64实对称正定求解器，含逐列原K残差验收。

新写实现；只依赖调用者显式选择的数值运行库，不附带二进制。
导入本模块只需标准库；构造时才需要NumPy、SciPy和psutil。
每个handle必须提供未存在的独立缓存目录。OOC core限制内部工作内存，
不是进程RSS硬限制；单次C调用的峰值仍需调用者独立进程监控。
要求实际oneMKL至少2025.2；当前验证版本为2025.2。
"""
from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
import ctypes as ct
import hashlib
import math
from numbers import Integral,Real
import os
from pathlib import Path
import re
import threading
import time

_ENV_LOCK=threading.RLock()
_INTEGER=ct.c_longlong
_INTEGER_POINTER=ct.POINTER(_INTEGER)


def array_identity(array):
    """记录实际数值字节、dtype与形状，不序列化或修改原数组。"""
    import numpy as np
    value=np.ascontiguousarray(array)
    return dict(dtype=value.dtype.name,shape=list(value.shape),sha256=hashlib.sha256(memoryview(value).cast('B')).hexdigest())


def _positive_integer(value,name,*,maximum=2**63-1):
    if isinstance(value,bool) or not isinstance(value,Integral) or not 1<=value<=maximum:
        raise ValueError(name+'须为范围内的正整数')
    return int(value)


def _positive_real(value,name):
    if isinstance(value,bool) or not isinstance(value,Real) or not math.isfinite(value) or value<=0:
        raise ValueError(name+'须为正有限实数')
    return float(value)


def _explicit_path(value,name):
    if isinstance(value,(bool,bytes)) or not isinstance(value,(str,os.PathLike)) or not str(value).strip():
        raise ValueError(name+'须显式指定非空路径')
    return Path(value).expanduser().resolve()


def _load_library(path):
    """只加载明确指定的库；要求线程局部设置，拒绝隐式全局fallback。"""
    if not path.is_file():raise FileNotFoundError('显式数值运行库文件不存在')
    directory=os.add_dll_directory(str(path.parent)) if os.name=='nt' else None
    try:
        library=ct.CDLL(str(path))
        if not hasattr(library,'pardiso_64'):raise RuntimeError('指定库不导出pardiso_64，不能降为32位接口')
        if not hasattr(library,'MKL_Set_Num_Threads_Local') or not hasattr(library,'MKL_Get_Max_Threads'):
            raise RuntimeError('指定库缺少线程局部设置／查询符号，拒绝隐式全局线程设置')
        if not hasattr(library,'MKL_Get_Version_String'):raise RuntimeError('指定库缺少实际版本查询符号')
        library.MKL_Get_Version_String.argtypes=[ct.c_void_p,ct.c_int]
        library.MKL_Get_Version_String.restype=None
        library.MKL_Set_Num_Threads_Local.argtypes=[ct.c_int]
        library.MKL_Set_Num_Threads_Local.restype=ct.c_int
        library.MKL_Get_Max_Threads.argtypes=[];library.MKL_Get_Max_Threads.restype=ct.c_int
        library.pardiso_64.argtypes=[ct.POINTER(ct.c_void_p),_INTEGER_POINTER,_INTEGER_POINTER,
            _INTEGER_POINTER,_INTEGER_POINTER,_INTEGER_POINTER,ct.c_void_p,_INTEGER_POINTER,
            _INTEGER_POINTER,_INTEGER_POINTER,_INTEGER_POINTER,_INTEGER_POINTER,_INTEGER_POINTER,
            ct.c_void_p,ct.c_void_p,_INTEGER_POINTER]
        library.pardiso_64.restype=None
        return library,directory
    except BaseException:
        if directory is not None:directory.close()
        raise


class PardisoSPD64:
    """上三角CSR、零基索引、float64及全int64的SPD Cholesky。

    实际原生RHS批次在phase11前固定；默认24。不足列严格补零，回原形状。
    auto_factor=False只执行phase11和预算检查，调用者审核后显式factor。
    不切换进程全局MKL接口层；pardiso_64每一个整数实参均显式64位。
    """
    def __init__(self,matrix,*,library_path,cache_directory,max_core_mib,mode='OOC',
            maximum_IC_gib=1.,threads=1,maximum_rhs_columns=24,minimum_available_gib=1.,
            maximum_process_rss_gib=2.,minimum_free_disk_gib=1.,auto_factor=True):
        self._released=True;self._factored=False;self._analysis_approved=False;self._native_touched=False
        self._dll_directory=None;self.calls=[];self._lock=threading.RLock()
        self.last_solve_certificate=None;self.last_failure=None
        self._fixed_rhs_context=None;self.last_fixed_rhs_metadata=None
        if mode not in ('IC','OOC'):raise ValueError('模式仅允许IC或OOC')
        if not isinstance(auto_factor,bool):raise ValueError('auto_factor须为布尔值')
        self.mode=mode;self.max_core_mib=_positive_integer(max_core_mib,'max_core_mib')
        self.threads=_positive_integer(threads,'threads',maximum=2**31-1)
        self.maximum_rhs_columns=_positive_integer(maximum_rhs_columns,'maximum_rhs_columns')
        for name,value in dict(maximum_IC_gib=maximum_IC_gib,minimum_available_gib=minimum_available_gib,
                maximum_process_rss_gib=maximum_process_rss_gib,minimum_free_disk_gib=minimum_free_disk_gib).items():
            setattr(self,name,_positive_real(value,name))
        self.library_path=_explicit_path(library_path,'library_path')
        self.cache_directory=_explicit_path(cache_directory,'cache_directory')
        if self.cache_directory.exists():raise FileExistsError('缓存目录已经存在，每个handle须使用独立新目录')
        if ct.sizeof(_INTEGER)!=8:raise RuntimeError('当前C long long不是8字节')
        import numpy as np
        from scipy import sparse
        import psutil
        self.np=np;self._psutil=psutil
        if not sparse.issparse(matrix) or len(matrix.shape)!=2 or matrix.dtype.kind not in 'fiu':
            raise ValueError('须传入实数数值型SciPy稀疏矩阵，不能丢弃复数虚部')
        if matrix.shape[0]!=matrix.shape[1] or not 1<=matrix.shape[0]<2**63:
            raise ValueError('须传入64位范围内的非空方阵')
        try:
            # 先复制为实际求解精度，再合并重复项，避免整数CSR先求和溢出。
            canonical=matrix.astype(np.float64,copy=True).tocsr();canonical.check_format(full_check=True)
            if not np.isfinite(canonical.data).all():raise ValueError('稀疏矩阵包含非有限值')
            canonical.sum_duplicates();canonical.eliminate_zeros();canonical.sort_indices()
            if canonical.nnz==0 or not np.isfinite(canonical.data).all():raise ValueError('矩阵为空或重复项求和后非有限')
        except (TypeError,OverflowError) as error:raise ValueError('不能规范化为有限实CSR矩阵') from error
        if np.any(canonical.diagonal()<=0):raise ValueError('SPD矩阵必须具有严格正对角线')
        asymmetry=canonical-canonical.T;scale=max(float(np.max(np.abs(canonical.data))),1.)
        symmetry_error=float(np.max(np.abs(asymmetry.data)))/scale if asymmetry.nnz else 0.
        if symmetry_error>1e-12:raise ValueError('拒绝非对称矩阵，不自动对称化')
        upper=sparse.triu(canonical,format='csr');upper.sum_duplicates();upper.sort_indices()
        if upper.nnz>=2**63:raise ValueError('CSR非零数超出有符号64位范围')
        self.values=np.array(upper.data,dtype=np.float64,copy=True,order='C')
        self.columns=np.array(upper.indices,dtype=np.int64,copy=True,order='C')
        self.pointers=np.array(upper.indptr,dtype=np.int64,copy=True,order='C')
        if not np.isfinite(self.values).all():raise ValueError('float64转换产生非有限矩阵值')
        self.input_identity=self._matrix_identity()
        for array in (self.values,self.columns,self.pointers):array.flags.writeable=False
        self.size=_INTEGER(upper.shape[0]);self.maxfct=_INTEGER(1);self.mnum=_INTEGER(1)
        self.mtype=_INTEGER(2);self.msglvl=_INTEGER(0)
        self.pt=(ct.c_void_p*64)();self.perm=np.zeros(upper.shape[0],dtype=np.int64)
        self.iparm=np.zeros(64,dtype=np.int64)
        # 64位手工初始化，不把int64控制数组传给LP64 pardisoinit。
        self.iparm[0]=1;self.iparm[1]=2;self.iparm[7]=2;self.iparm[17]=-1;self.iparm[18]=-1
        self.iparm[26]=1;self.iparm[34]=1;self.iparm[59]=2 if mode=='OOC' else 0
        self.environment=dict(MKL_PARDISO_OOC_PATH=str(self.cache_directory)+os.sep,
            MKL_PARDISO_OOC_FILE_NAME=str(self.cache_directory/'factor64'),
            MKL_PARDISO_OOC_MAX_CORE_SIZE=str(self.max_core_mib),MKL_PARDISO_OOC_MAX_SWAP_SIZE='0',
            MKL_PARDISO_OOC_KEEP_FILE='0',MKL_PARDISO_OOC_CFG_PATH=str(self.cache_directory)+os.sep,
            MKL_PARDISO_OOC_CFG_FILE_NAME='environment_only_not_present.cfg')
        try:
            self.library,self._dll_directory=_load_library(self.library_path)
            # mkdir的exist_ok=False是实际目录所有权争用门控，不删除旧目录。
            self.cache_directory.mkdir(parents=True,exist_ok=False)
            self._released=False
            self._verify_runtime_version()
            self._call(11)
            if any(int(self.iparm[i])<0 for i in (14,15,16,62)):
                raise ArithmeticError('phase11返回非法负内存估算，拒绝分解')
            ic_kib=max(int(self.iparm[14]),int(self.iparm[15])+int(self.iparm[16]))
            ooc_kib=max(int(self.iparm[14]),int(self.iparm[15])+int(self.iparm[62]))
            self.analysis_diagnostics=dict(backend='oneMKL_PARDISO_64',mode=mode,matrix_type=2,
                factorization='SPD_Cholesky',integer_width_bits=64,float_width_bits=64,index_base=0,
                threads=self.threads,maximum_rhs_columns=self.maximum_rhs_columns,
                analysis_peak_KiB=int(self.iparm[14]),symbolic_permanent_KiB=int(self.iparm[15]),
                IC_numerical_factor_KiB=int(self.iparm[16]),extra_minimum_OOC_KiB=int(self.iparm[62]),
                estimated_IC_GiB=ic_kib/2**20,minimum_OOC_GiB=ooc_kib/2**20,
                max_core_MiB=self.max_core_mib,factor_nonzeros=int(self.iparm[17]),
                source_symmetry_max_relative_error=symmetry_error,input_identity=self.input_identity,
                defaults_initialization='手工int64控制数组，不调用LP64初始化函数',
                global_interface_layer_changed=False,thread_setting='逐调用线程局部设置并恢复',
                runtime_version=self.runtime_version,runtime_minimum_version='2025.2',
                cache_directory=str(self.cache_directory),scalar_integer_bytes=ct.sizeof(_INTEGER),
                integer_array_dtypes={key:getattr(self,key).dtype.name for key in ('columns','pointers','perm','iparm')})
            self.analysis_diagnostics.update(native_rhs_batch_policy='fixed_maximum_zero_padded',
                fixed_native_rhs_columns=self.maximum_rhs_columns,original_maximum_iterative_refinement_steps=2,
                all_actual_and_zero_columns_use_original_Ku_guard=True,
                original_rhs_identity_is_separate_from_native_padded_identity=True)
            if mode=='OOC' and ooc_kib>=self.max_core_mib*1024:
                raise MemoryError('phase11最低OOC内存不小于显式core上限，拒绝phase22')
            if mode=='IC' and ic_kib/2**20>self.maximum_IC_gib:
                raise MemoryError('phase11的IC内存估算超出显式上限，拒绝phase22')
            self._prepare_residual_operator()
            self._analysis_approved=True
            if auto_factor:self.factor()
        except BaseException as error:
            self._fail_and_close(error,'construction_or_phase11')
            if self._dll_directory is not None:self._dll_directory.close();self._dll_directory=None
            raise

    def _verify_runtime_version(self):
        """版本查询也在受锁DLL／PATH／线程上下文内；早于任何phase11。"""
        text=ct.create_string_buffer(512)
        with self._native_context():self.library.MKL_Get_Version_String(text,ct.c_int(512))
        value=text.value.decode('ascii','replace')
        parsed=re.search(r'\bVersion\s+(\d+)\.(\d+)(?:\.(\d+))?',value)
        self.runtime_version=dict(text=value,major=None,minor=None,update=None)
        if parsed is None:raise RuntimeError('无法解析实际oneMKL版本，要求至少2025.2')
        major,minor=int(parsed[1]),int(parsed[2]);update=int(parsed[3]) if parsed[3] else 0
        self.runtime_version.update(major=major,minor=minor,update=update)
        if (major,minor)<(2025,2):raise RuntimeError('实际oneMKL版本低于最低2025.2：'+value)

    def _prepare_residual_operator(self):
        """共享原64位upper CSR并缓存对角；先预留可能的构造副本。"""
        from scipy import sparse
        np=self.np;n=self.size.value;nnz=self.values.size
        # 实际要求共享三组数组；仍预留构造器可能产生的整份数据／索引
        # 副本，避免版本差异在共享检查之前造成未纳入预算的瞬时峰。
        reserve=16*nnz+8*(n+1)+24*n
        self._resource_check(reserve)
        operator=sparse.csr_array((self.values,self.columns,self.pointers),shape=(n,n),copy=False)
        shared={key:bool(np.shares_memory(actual,original)) for key,actual,original in
            (('values',operator.data,self.values),('columns',operator.indices,self.columns),('pointers',operator.indptr,self.pointers))}
        self.analysis_diagnostics['residual_guard']=dict(relative_tolerance=1e-8,column_block=1,
            representation='原upper CSR的A+A.T-diag(A)，不改原矩阵',shared_arrays=shared,
            indices_dtype=operator.indices.dtype.name,pointers_dtype=operator.indptr.dtype.name,
            construction_reserve_bytes=reserve,single_column_reserve_bytes=24*n)
        if not all(shared.values()) or operator.indices.dtype!=np.dtype('int64') or operator.indptr.dtype!=np.dtype('int64'):
            raise RuntimeError('当前SciPy不能按契约共享三组原64位upper CSR数组，拒绝隐藏副本')
        self._residual_upper=operator
        self._residual_diagonal=np.asarray(operator.diagonal(),dtype=np.float64)
        self._residual_diagonal.flags.writeable=False
        self.analysis_diagnostics['residual_guard']['cached_diagonal_bytes']=self._residual_diagonal.nbytes

    def _certify_response(self,rhs,result):
        """逐列验收||K_upper u-b||₂/||b||₂；不改U，不生成N×RHS的Ku。"""
        np=self.np;upper=self._residual_upper;count=rhs.shape[1]
        certificate=dict(passed=False,relative_tolerance=1e-8,column_block=1,requested_columns=count,columns=[])
        self.last_solve_certificate=certificate
        if self.calls and self.calls[-1]['phase']==33:self.calls[-1]['original_K_residual_certificate']=certificate
        for column in range(count):
            self._resource_check(24*self.size.value)
            with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
                r=upper@result[:,column]
                r+=upper.T@result[:,column]
                r-=self._residual_diagonal*result[:,column]
                r-=rhs[:,column]
                norm_b=float(np.linalg.norm(rhs[:,column]));norm_r=float(np.linalg.norm(r))
                maximum_abs=float(np.max(np.abs(r)))
            finite=math.isfinite(norm_b) and math.isfinite(norm_r)
            relative=(norm_r/norm_b if norm_b else (0. if norm_r==0 else float('inf'))) if finite else float('inf')
            passed=math.isfinite(relative) and relative<=1e-8
            record=dict(local_column=column,norm_B=norm_b if math.isfinite(norm_b) else None,
                norm_R=norm_r if math.isfinite(norm_r) else None,relative_residual=relative if math.isfinite(relative) else None,
                nonfinite_residual_count=int(np.count_nonzero(~np.isfinite(r))),
                maximum_abs_residual=maximum_abs if math.isfinite(maximum_abs) else None,passed=passed)
            certificate['columns'].append(record)
            if not passed:
                certificate.update(first_failed_column=column,checked_columns=column+1,
                    unprocessed_columns=count-column-1,maximum_relative_residual=max(
                        (row['relative_residual'] for row in certificate['columns'] if row['relative_residual'] is not None),default=None))
                raise ArithmeticError(f'原upper K残差未通过：局部列{column}，相对值{relative}，门限1e-8')
            del r
        certificate.update(passed=True,checked_columns=count,unprocessed_columns=0,
            maximum_relative_residual=max(row['relative_residual'] for row in certificate['columns']))

    def _fail_and_close(self,error,stage):
        """内部求解／资源／验收失败保留量值并释放，不把有限坏U返回。"""
        if self._released and self.last_failure is not None and getattr(error,'spd64_failure',None) is self.last_failure:return
        failure=dict(stage=stage,exception_type=type(error).__name__,message=str(error),
            runtime_version=getattr(self,'runtime_version',None),solve_certificate=self.last_solve_certificate,
            analysis_diagnostics=getattr(self,'analysis_diagnostics',None))
        self.last_failure=failure
        try:self.close()
        except BaseException as close_error:failure['close_failure']=dict(exception_type=type(close_error).__name__,message=str(close_error))
        failure.update(handle_released=self._released,calls=list(self.calls))
        try:error.spd64_failure=failure
        except (AttributeError,TypeError):pass

    def _matrix_identity(self):
        return {key:array_identity(getattr(self,key)) for key in ('values','columns','pointers')}

    def _resource_check(self,extra_bytes=0):
        psutil=self._psutil;available=psutil.virtual_memory().available
        if available-extra_bytes<self.minimum_available_gib*2**30:raise MemoryError('调用前系统可用内存不足，包括显式RHS预留')
        if psutil.Process().memory_info().rss+extra_bytes>self.maximum_process_rss_gib*2**30:
            raise MemoryError('调用前进程RSS与RHS预留之和超过预算')
        if psutil.disk_usage(str(self.cache_directory)).free<self.minimum_free_disk_gib*2**30:
            raise OSError('独立缓存卷剩余空间低于显式门限')

    @contextmanager
    def _native_context(self):
        """本模块所有handle串行设置／恢复环境及当前调用线程的MKL线程数。"""
        with _ENV_LOCK:
            managed=dict(self.environment)
            # 部分Windows数值库内部仍按PATH加载线程DLL；只在本次受锁
            # 调用中加入调用者指定目录，完成后恢复原值，不永久改PATH。
            if os.name=='nt':managed['PATH']=str(self.library_path.parent)+os.pathsep+os.environ.get('PATH','')
            previous={key:os.environ.get(key) for key in managed}
            setter=self.library.MKL_Set_Num_Threads_Local;getter=self.library.MKL_Get_Max_Threads
            old_local=None;info={}
            try:
                os.environ.update(managed)
                info['effective_threads_before']=int(getter())
                old_local=int(setter(self.threads))
                info.update(previous_local_threads=old_local,requested_threads=self.threads,effective_threads=int(getter()))
                if info['effective_threads']!=self.threads:raise RuntimeError('实际oneMKL线程数与显式请求不符')
                yield info
            finally:
                try:
                    if old_local is not None:
                        info['restore_setter_return']=int(setter(old_local))
                        info['effective_threads_after']=int(getter())
                        if info['effective_threads_after']!=info['effective_threads_before']:
                            raise RuntimeError('oneMKL线程局部设置未恢复到调用前有效值')
                finally:
                    for key,value in previous.items():
                        if value is None:os.environ.pop(key,None)
                        else:os.environ[key]=value

    def cache_inventory(self):
        """只列自己目录的文件；不递归、删除或搬移任何缓存。"""
        return [dict(file=p.name,bytes=p.stat().st_size) for p in sorted(self.cache_directory.iterdir()) if p.is_file()]

    def _read_process_bytes(self):
        """进程级只读I/O计量；缺接口／权限时明确unavailable。"""
        try:return int(self._psutil.Process().io_counters().read_bytes)
        except Exception:return None

    def _call(self,phase,rhs=None):
        """只附到本次新记录，保留原细化控制与实际nrhs/原请求分列。"""
        if phase!=-1 and int(self.iparm[7])!=2:raise ArithmeticError('原最大细化控制须保持2，拒绝隐式改写')
        if phase==33 and (rhs is None or rhs.ndim!=2 or rhs.shape[1]!=self.maximum_rhs_columns):
            raise ValueError('phase33实际RHS列数必须等于phase11固定批次')
        start=len(self.calls);before=self._read_process_bytes()
        try:return self._call_native(phase,rhs)
        finally:
            if len(self.calls)==start+1 and self.calls[start]['phase']==phase:
                after=self._read_process_bytes();row=self.calls[start]
                row['fixed_rhs_telemetry']=dict(actual_nrhs=row['nrhs'],
                    original_request=deepcopy(self._fixed_rhs_context) if phase==33 else None,
                    actual_native_rhs_identity=deepcopy(getattr(self,'last_rhs_identity',None)) if phase==33 else None,
                    actual_iterative_refinement_steps=int(self.iparm[6]),maximum_iterative_refinement_steps=int(self.iparm[7]),
                    backward_error_tolerance_parameter=int(self.iparm[8]),perturbed_pivots=int(self.iparm[13]),
                    process_read_bytes_before=before if before is not None else 'unavailable',
                    process_read_bytes_after=after if after is not None else 'unavailable',
                    process_read_bytes_delta=after-before if before is not None and after is not None else 'unavailable',
                    read_bytes_scope='整个当前进程，不等同于本handle缓存读取量')

    def _call_native(self,phase,rhs=None):
        np=self.np
        if phase!=-1 and self._matrix_identity()!=self.input_identity:raise ArithmeticError('实际CSR输入字节已改变，拒绝调用')
        if rhs is None:
            count=self.maximum_rhs_columns if phase in (11,22) else 1
            length=self.size.value if phase in (11,22) else 1
            # input/output各8B，C-order身份hash可能8B，finite临时bool 1B。
            if phase!=-1:self._resource_check(length*count*25)
            rhs=np.zeros((length,count),dtype=np.float64,order='F')
        else:self._resource_check(rhs.size*17)  # 输入已存在，预留output/hash/finite。
        output=np.zeros_like(rhs,order='F');nrhs=_INTEGER(rhs.shape[1]);error=_INTEGER(0);p=_INTEGER(phase)
        before=array_identity(rhs);started=time.perf_counter();thread_info={}
        with self._native_context() as thread_info:
            self._native_touched=True
            self.library.pardiso_64(self.pt,ct.byref(self.maxfct),ct.byref(self.mnum),ct.byref(self.mtype),
                ct.byref(p),ct.byref(self.size),self.values.ctypes.data_as(ct.c_void_p),
                self.pointers.ctypes.data_as(_INTEGER_POINTER),self.columns.ctypes.data_as(_INTEGER_POINTER),
                self.perm.ctypes.data_as(_INTEGER_POINTER),ct.byref(nrhs),self.iparm.ctypes.data_as(_INTEGER_POINTER),
                ct.byref(self.msglvl),rhs.ctypes.data_as(ct.c_void_p),output.ctypes.data_as(ct.c_void_p),ct.byref(error))
        memory=self._psutil.Process().memory_info()
        self.calls.append(dict(phase=phase,nrhs=nrhs.value,error=error.value,seconds=time.perf_counter()-started,
            rss_bytes=memory.rss,thread_context=thread_info,cache=self.cache_inventory()))
        if error.value:raise ArithmeticError(f'PARDISO_64阶段{phase}返回错误码{error.value}')
        if before!=array_identity(rhs):raise ArithmeticError('数值运行库修改了RHS输入字节')
        if phase!=-1:
            if self._matrix_identity()!=self.input_identity:raise ArithmeticError('数值运行库修改了CSR输入字节')
            if not np.isfinite(output).all():raise ArithmeticError('数值运行库返回非有限响应')
            if memory.rss>self.maximum_process_rss_gib*2**30:raise MemoryError('调用后进程RSS超限；C调用内硬峰值仍须外部监控')
        return output

    def factor(self):
        """只在phase11及预算已通过时执行一次phase22。"""
        with self._lock:
            if self._released or not self._analysis_approved or self._factored:
                raise RuntimeError('因子已释放、分析未通过或已经分解')
            try:self._call(22)
            except BaseException as error:
                self._fail_and_close(error,'phase22')
                raise
            self._factored=True;self.positive_definite=True
            self.analysis_diagnostics['factor_nonzeros']=int(self.iparm[17])

    def _solve_native_rhs(self,rhs):
        """同阶有限实向量或二维列矩阵；复制RHS，保持调用者原数组不变。"""
        with self._lock:
            if self._released or not self._factored:raise RuntimeError('尚无可用数值因子，须先分解且未释放')
            np=self.np;values=np.asarray(rhs)
            if values.dtype.kind not in 'fiu':raise ValueError('RHS须为实数数值型，不能丢弃复数虚部')
            vector=values.ndim==1
            if vector:values=values[:,None]
            if values.ndim!=2 or values.shape[0]!=self.size.value or not 1<=values.shape[1]<=self.maximum_rhs_columns:
                raise ValueError('RHS须同阶且列数在phase11约定范围内，不允许零列或截断')
            self._resource_check(values.size)
            if not np.isfinite(values).all():raise ValueError('RHS包含非有限值')
            self.last_solve_certificate=None
            try:
                self._resource_check(values.size*25+24*self.size.value)
                values=np.array(values,dtype=np.float64,copy=True,order='F')
                if not np.isfinite(values).all():raise ValueError('float64转换产生非有限RHS')
                self.last_rhs_identity=array_identity(values)
                result=self._call(33,values)
                self._certify_response(values,result)
            except BaseException as error:
                self._fail_and_close(error,'phase33_or_original_K_residual')
                raise
            return result[:,0] if vector else result

    def solve(self,rhs):
        """原输入有限实数→固定max列零补→全部实际列原Ku→owned原形状。

        last_rhs_identity/last_solve_certificate描述实际原生批次；
        last_fixed_rhs_metadata单列调用方原输入，二者不能混称。
        """
        with self._lock:
            if self._released or not self._factored:raise RuntimeError('尚无可用数值因子，须先分解且未释放')
            np=self.np;original=np.asarray(rhs)
            if original.dtype.kind not in 'fiu':raise ValueError('RHS须为实数数值型，不能丢弃复数虚部')
            vector=original.ndim==1;shaped=original[:,None] if vector else original
            if shaped.ndim!=2 or shaped.shape[0]!=self.size.value or not 1<=shaped.shape[1]<=self.maximum_rhs_columns:
                raise ValueError('RHS须同阶且列数在phase11约定范围内，不允许零列或截断')
            # 无效输入不销毁有效因子；资源预留失败则保留失败证据并关闭。
            try:self._resource_check(shaped.size)
            except BaseException as error:
                self._fail_and_close(error,'fixed_rhs_input_validation_budget');raise
            if not np.isfinite(shaped).all():raise ValueError('RHS包含非有限值')
            count=int(shaped.shape[1]);n=self.size.value;actual_count=self.maximum_rhs_columns
            self.last_solve_certificate=None;self.last_rhs_identity=None
            metadata=dict(original_input_shape=list(original.shape),original_input_dtype=original.dtype.name,
                original_nrhs=count,original_is_vector=vector,actual_native_nrhs=actual_count,
                zero_padding_columns=actual_count-count,copied_output=False,status='准备零补；尚未通过全部实际列原Ku')
            self.last_fixed_rhs_metadata=metadata
            try:
                padding=n*actual_count*8 if count<actual_count else 0
                reserve=padding+int(original.nbytes)+(n*count if count<actual_count else 0)
                metadata.update(padding_preallocation_reserve_bytes=reserve,
                    padding_reserve_breakdown=dict(float64_padding_bytes=padding,original_hash_scratch_bytes=int(original.nbytes),
                        converted_prefix_finite_bytes=n*count if count<actual_count else 0))
                self._resource_check(reserve)
                metadata['original_rhs_identity']=array_identity(original)
                self._fixed_rhs_context=deepcopy(metadata)
                if count<actual_count:
                    actual=np.zeros((n,actual_count),dtype=np.float64,order='F');actual[:,:count]=shaped
                    if not np.isfinite(actual[:,:count]).all():raise ValueError('float64零补转换产生非有限值')
                else:actual=shaped
                full=self._solve_native_rhs(actual);certificate=self.last_solve_certificate
                if (not isinstance(certificate,dict) or certificate.get('passed') is not True or
                        certificate.get('requested_columns')!=actual_count or len(certificate.get('columns',[]))!=actual_count):
                    raise ArithmeticError('原后端没有对全部实际列含零列完成Ku验收')
                if count<actual_count:
                    self._resource_check(n*(actual_count-count))
                    if np.any(full[:,count:]!=0):raise ArithmeticError('原零RHS列产生非零响应，拒绝回切')
                self._resource_check(int(original.nbytes));unchanged=array_identity(original)==metadata['original_rhs_identity']
                if not unchanged:raise ArithmeticError('调用方原RHS字节改变')
                if vector:
                    self._resource_check(n*8);answer=np.array(full[:,0],dtype=np.float64,copy=True)
                elif count<actual_count:
                    self._resource_check(n*count*8);answer=np.array(full[:,:count],dtype=np.float64,copy=True,order='F')
                else:answer=full
                if not answer.flags.owndata:raise ArithmeticError('返回响应必须拥有原形状数据，拒绝隐藏整批次view')
                metadata.update(status='全部实际列含零列通过原Ku，已返回原形状',copied_output=vector or count<actual_count,
                    returned_shape=list(answer.shape),returned_array_owns_data=True,
                    actual_native_rhs_identity=deepcopy(self.last_rhs_identity),all_actual_original_Ku_columns_passed=True,
                    original_K_relative_residual_limit=1e-8,zero_column_maximum_abs_U=0.,
                    original_rhs_identity_unchanged=unchanged,original_rhs_identity_checked_before_trim=True)
                return answer
            except BaseException as error:
                metadata.update(status='失败，不返回响应',failure_type=type(error).__name__,failure_message=str(error))
                self._fail_and_close(error,'fixed_rhs_padding_or_original_Ku')
                if isinstance(getattr(error,'spd64_failure',None),dict):error.spd64_failure['fixed_rhs_padding']=deepcopy(metadata)
                raise
            finally:self._fixed_rhs_context=None

    def close(self):
        """幂等释放本handle；始终关闭本实例的DLL搜索目录句柄。"""
        with self._lock:
            try:
                if not self._released and self._native_touched:self._call(-1)
            finally:
                self._released=True;self._factored=False;self._analysis_approved=False
                self._residual_upper=None;self._residual_diagonal=None
                if self._dll_directory is not None:self._dll_directory.close();self._dll_directory=None

    def __enter__(self):
        if self._released:raise RuntimeError('不能进入已经释放的因子上下文')
        return self

    def __exit__(self,*exc):self.close()

    def __del__(self):
        try:self.close()
        except Exception:pass
