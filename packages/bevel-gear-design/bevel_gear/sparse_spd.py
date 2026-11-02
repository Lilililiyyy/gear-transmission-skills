"""可选专业稀疏SPD求解器：oneMKL PARDISO Cholesky。

按oneMKL官方C协议使用LP64、上三角CSR及零基索引；库路径由调用者
明确传入；本项目不附带二进制运行库。分解前单独分析内存，超限即释放拒绝。
线程和LP64接口设置作用于调用进程中的oneMKL；Windows还将明确选择的
运行库目录加入本进程PATH，供oneMKL加载其依赖。宜在独立计算进程使用。
工厂接口可直接传给prepare_generated_ltca的factorization_factory。
"""
from pathlib import Path
import ctypes as ct
import os
import math


class PardisoSPD:
    def __init__(self,matrix,*,library_path,maximum_factor_memory_gib=20.,threads=1):
        import numpy as np
        from scipy import sparse
        self.np=np;self._released=True
        if isinstance(threads,bool) or not isinstance(threads,int) or not 1<=threads<2**31:
            raise ValueError('稀疏求解线程数须为正整数')
        if (isinstance(maximum_factor_memory_gib,bool) or not isinstance(maximum_factor_memory_gib,(int,float))
                or not math.isfinite(maximum_factor_memory_gib) or maximum_factor_memory_gib<=0):
            raise ValueError('分解分析内存上限须为正有限GiB')
        if not sparse.issparse(matrix) or len(matrix.shape)!=2 or np.iscomplexobj(matrix.data):
            raise ValueError('SPD后端要求实数SciPy稀疏矩阵')
        if matrix.shape[0]!=matrix.shape[1] or matrix.shape[0]<1:
            raise ValueError('刚度须为非空方阵')
        if not np.isfinite(matrix.data).all():raise ValueError('刚度含非有限值')
        if matrix.nnz==0:raise ValueError('SPD后端拒绝空刚度')
        asymmetry=matrix-matrix.T
        scale=max(float(np.max(np.abs(matrix.data))),1.)
        if asymmetry.nnz and np.max(np.abs(asymmetry.data))>scale*1e-12:
            raise ValueError('SPD后端拒绝非对称刚度')
        upper=sparse.triu(matrix,format='csr');upper.sum_duplicates();upper.sort_indices()
        if upper.nnz>=2**31 or upper.shape[0]>=2**31:raise ValueError('本适配采用LP64，索引超出32位')
        self.values=np.ascontiguousarray(upper.data,dtype=np.float64)
        self.columns=np.ascontiguousarray(upper.indices,dtype=np.int32)
        self.pointers=np.ascontiguousarray(upper.indptr,dtype=np.int32)
        self.size=ct.c_int(upper.shape[0]);self.iparm=np.zeros(64,dtype=np.int32)
        self.pt=(ct.c_void_p*64)();self.perm=np.zeros(upper.shape[0],dtype=np.int32)
        self.mtype=ct.c_int(2);self.maxfct=ct.c_int(1);self.mnum=ct.c_int(1)
        self.msglvl=ct.c_int(0);self._directory=None
        library_path=Path(library_path).resolve()
        if os.name=='nt':
            self._directory=os.add_dll_directory(str(library_path.parent))
            os.environ['PATH']=str(library_path.parent)+os.pathsep+os.environ.get('PATH','')
        try:
            self.library=ct.CDLL(str(library_path))
        except BaseException:
            if self._directory is not None:self._directory.close();self._directory=None
            raise
        self.library.MKL_Set_Interface_Layer.argtypes=[ct.c_int]
        self.library.MKL_Set_Interface_Layer.restype=ct.c_int
        if self.library.MKL_Set_Interface_Layer(0)!=0:raise RuntimeError('无法设为oneMKL LP64接口')
        self.library.MKL_Set_Num_Threads.argtypes=[ct.c_int]
        self.library.MKL_Set_Num_Threads(ct.c_int(threads))
        self.library.pardisoinit.argtypes=[ct.c_void_p,ct.POINTER(ct.c_int),ct.c_void_p]
        self.library.pardisoinit(self.pt,ct.byref(self.mtype),self.iparm.ctypes.data_as(ct.c_void_p))
        self.iparm[34]=1;self.iparm[26]=1;self.iparm[7]=2
        self.library.pardiso.argtypes=[ct.c_void_p]+[ct.c_void_p]*15
        self.library.pardiso.restype=None
        self._released=False
        try:
            self._call(11)
            # 官方iparm[14]为分析峰值KB，[15]为永久符号KB，[16]为数值因子KB。
            peak_kib=max(int(self.iparm[14]),int(self.iparm[15])+int(self.iparm[16]))
            self.analysis_diagnostics=dict(backend='oneMKL_PARDISO',factorization='SPD_Cholesky',
                matrix_type=2,index_base=0,integer_width_bits=32,threads=threads,
                analysis_peak_kib=int(self.iparm[14]),symbolic_permanent_kib=int(self.iparm[15]),
                numerical_factor_kib=int(self.iparm[16]),estimated_peak_gib=peak_kib/2**20,
                factor_nonzeros=int(self.iparm[17]),source_symmetry_max_relative_error=
                float(np.max(np.abs(asymmetry.data))/scale) if asymmetry.nnz else 0.)
            if peak_kib/2**20>maximum_factor_memory_gib:
                raise MemoryError('SPD分解分析内存超出显式上限：'+str(self.analysis_diagnostics))
            self._call(22)
            self.positive_definite=True
            self.analysis_diagnostics['factor_nonzeros']=int(self.iparm[17])
        except BaseException:
            self.close();raise

    def _call(self,phase,rhs=None):
        np=self.np;phase_value=ct.c_int(phase);error=ct.c_int(0)
        if rhs is None:rhs=np.zeros((self.size.value,1),dtype=np.float64,order='F')
        output=np.zeros_like(rhs,order='F');nrhs=ct.c_int(rhs.shape[1])
        self.library.pardiso(self.pt,ct.byref(self.maxfct),ct.byref(self.mnum),ct.byref(self.mtype),
            ct.byref(phase_value),ct.byref(self.size),self.values.ctypes.data_as(ct.c_void_p),
            self.pointers.ctypes.data_as(ct.c_void_p),self.columns.ctypes.data_as(ct.c_void_p),
            self.perm.ctypes.data_as(ct.c_void_p),ct.byref(nrhs),self.iparm.ctypes.data_as(ct.c_void_p),
            ct.byref(self.msglvl),rhs.ctypes.data_as(ct.c_void_p),output.ctypes.data_as(ct.c_void_p),ct.byref(error))
        if error.value:raise ArithmeticError(f'PARDISO阶段{phase}返回错误码{error.value}')
        return output

    def solve(self,rhs):
        np=self.np
        if self._released:raise RuntimeError('Cholesky因子已经释放')
        if np.iscomplexobj(rhs):raise ValueError('载荷须为实数，不能丢弃虚部')
        values=np.asarray(rhs,dtype=np.float64);vector=values.ndim==1
        if vector:values=values[:,None]
        if values.ndim!=2 or values.shape[0]!=self.size.value or not np.isfinite(values).all():
            raise ValueError('载荷须为同阶有限实向量或列向量矩阵')
        output=self._call(33,np.asfortranarray(values))
        return output[:,0] if vector else output

    def close(self):
        if not self._released:
            try:self._call(-1)
            finally:
                self._released=True
                if self._directory is not None:self._directory.close();self._directory=None

    def __enter__(self):
        return self

    def __exit__(self,*exc):
        self.close()

    def __del__(self):
        try:self.close()
        except Exception:pass
        directory=getattr(self,'_directory',None)
        if directory is not None:
            try:directory.close()
            except Exception:pass
