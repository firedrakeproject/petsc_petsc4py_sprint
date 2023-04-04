# --------------------------------------------------------------------

class VecType(object):
    SEQ        = S_(VECSEQ)
    MPI        = S_(VECMPI)
    STANDARD   = S_(VECSTANDARD)
    SHARED     = S_(VECSHARED)
    SEQVIENNACL= S_(VECSEQVIENNACL)
    MPIVIENNACL= S_(VECMPIVIENNACL)
    VIENNACL   = S_(VECVIENNACL)
    SEQCUDA    = S_(VECSEQCUDA)
    MPICUDA    = S_(VECMPICUDA)
    CUDA       = S_(VECCUDA)
    SEQHIP     = S_(VECSEQHIP)
    MPIHIP     = S_(VECMPIHIP)
    HIP        = S_(VECHIP)
    NEST       = S_(VECNEST)
    SEQKOKKOS  = S_(VECSEQKOKKOS)
    MPIKOKKOS  = S_(VECMPIKOKKOS)
    KOKKOS     = S_(VECKOKKOS)

class VecOption(object):
    IGNORE_OFF_PROC_ENTRIES = VEC_IGNORE_OFF_PROC_ENTRIES
    IGNORE_NEGATIVE_INDICES = VEC_IGNORE_NEGATIVE_INDICES

# --------------------------------------------------------------------

cdef class Vec(Object):
    """A vector object.

    Vectors are further documented in
    `the PETSc manual <petsc:chapter_vectors>`

    See Also
    --------
    petsc.Vec

    """

    Type = VecType
    Option = VecOption

    #

    def __cinit__(self):
        self.obj = <PetscObject*> &self.vec
        self.vec = NULL

    # unary operations

    def __pos__(self):
        return vec_pos(self)

    def __neg__(self):
        return vec_neg(self)

    def __abs__(self):
        return vec_abs(self)

    # inplace binary operations

    def __iadd__(self, other):
        return vec_iadd(self, other)

    def __isub__(self, other):
        return vec_isub(self, other)

    def __imul__(self, other):
        return vec_imul(self, other)

    def __idiv__(self, other):
        return vec_idiv(self, other)

    def __itruediv__(self, other):
        return vec_idiv(self, other)

    # binary operations

    def __add__(self, other):
        if isinstance(self, Vec):
            return vec_add(self, other)
        else:
            return vec_radd(other, self)

    def __sub__(self, other):
        if isinstance(self, Vec):
            return vec_sub(self, other)
        else:
            return vec_rsub(other, self)

    def __mul__(self, other):
        if isinstance(self, Vec):
            return vec_mul(self, other)
        else:
            return vec_rmul(other, self)

    def __div__(self, other):
        if isinstance(self, Vec):
            return vec_div(self, other)
        else:
            return vec_rdiv(other, self)

    def __truediv__(self, other):
        if isinstance(self, Vec):
            return vec_div(self, other)
        else:
            return vec_rdiv(other, self)

    #

    #def __len__(self):
    #    cdef PetscInt size = 0
    #    CHKERR( VecGetSize(self.vec, &size) )
    #    return <Py_ssize_t>size

    def __getitem__(self, i):
        return vec_getitem(self, i)

    def __setitem__(self, i, v):
        vec_setitem(self, i, v)

    # buffer interface (PEP 3118)

    def __getbuffer__(self, Py_buffer *view, int flags):
        cdef _Vec_buffer buf = _Vec_buffer(self)
        buf.acquirebuffer(view, flags)

    def __releasebuffer__(self, Py_buffer *view):
        cdef _Vec_buffer buf = <_Vec_buffer>(view.obj)
        buf.releasebuffer(view)
        <void>self # unused

    # 'with' statement (PEP 343)

    def __enter__(self):
        cdef _Vec_buffer buf = _Vec_buffer(self)
        self.set_attr('__buffer__', buf)
        return buf.enter()

    def __exit__(self, *exc):
        cdef _Vec_buffer buf = self.get_attr('__buffer__')
        self.set_attr('__buffer__', None)
        return buf.exit()

    #

    def view(self, Viewer viewer=None) -> None:
        """Display the vector.

        Collective.

        Parameters
        ----------
        viewer
            The viewer instance, defaults to printing vector contents.

        See Also
        --------
        petsc.VecView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( VecView(self.vec, vwr) )

    def destroy(self) -> Self:
        """Destroy the vector.

        Collective.

        See Also
        --------
        Vec.create, petsc.VecDestroy

        """

        CHKERR( VecDestroy(&self.vec) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create an empty vector object.

        After creation the vector type can then be set with `Vec.setType`
        or `Vec.setFromOptions`.

        Collective.

        Parameters
        ----------
        comm
            Communicator for the vector, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        Vec.destroy, petsc.VecCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscVec newvec = NULL
        CHKERR( VecCreate(ccomm, &newvec) )
        PetscCLEAR(self.obj); self.vec = newvec
        return self

    def setType(self, vec_type: Vec.Type | str) -> None:
        """Set the type of the vector.

        Parameters
        ----------
        vec_type
            The type of the vector.

        See Also
        --------
        petsc_options, Vec.setFromOptions, petsc.VecSetType

        """
        cdef PetscVecType cval = NULL
        vec_type = str2bytes(vec_type, &cval)
        CHKERR( VecSetType(self.vec, cval) )

    def setSizes(
        self,
        size: tuple[int, int] | int,
        bsize: int | None = None,
    ) -> None:
        """Set the local and global sizes of the vector.

        Collective.

        Parameters
        ----------
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. For more information see `Sys.splitOwnership`.
        bsize
            Block size, defaults to ``1``.

        See Also
        --------
        Vec.getSizes, petsc.VecSetSizes

        """
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        CHKERR( VecSetSizes(self.vec, n, N) )
        if bs != PETSC_DECIDE:
            CHKERR( VecSetBlockSize(self.vec, bs) )

    #

    def createSeq(
        self,
        size: tuple[int, int] | int,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a standard, sequential vector.

        Collective.

        Parameters
        ----------
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. For more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        Vec.createMPI, petsc.VecCreateSeq

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_SELF)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        cdef PetscVec newvec = NULL
        CHKERR( VecCreate(ccomm,&newvec) )
        CHKERR( VecSetSizes(newvec, n, N) )
        CHKERR( VecSetBlockSize(newvec, bs) )
        CHKERR( VecSetType(newvec, VECSEQ) )
        PetscCLEAR(self.obj); self.vec = newvec
        return self

    def createMPI(
        self,
        size: tuple[int, int] | int,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a parallel vector.

        Collective.

        Parameters
        ----------
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. For more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        Vec.createSeq, petsc.VecCreateMPI

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        cdef PetscVec newvec = NULL
        CHKERR( VecCreate(ccomm, &newvec) )
        CHKERR( VecSetSizes(newvec, n, N) )
        CHKERR( VecSetBlockSize(newvec, bs) )
        CHKERR( VecSetType(newvec, VECMPI) )
        PetscCLEAR(self.obj); self.vec = newvec
        return self

    def createWithArray(
        self,
        array: Sequence[Scalar],
        size: tuple[int, int] | int = None,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a vector using a provided array.

        This method will create either a `Vec.Type.SEQ` or `Vec.Type.MPI`
        depending on whether the communicator has more than one rank.

        Collective.

        Parameters
        ----------
        array
            Array to store the vector values. Must be at least as large as
            the local size of the vector.
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. If `None` then defaults to the size of ``array``. For
            more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        petsc.VecCreateSeqWithArray, petsc.VecCreateMPIWithArray

        """
        cdef PetscInt na=0
        cdef PetscScalar *sa=NULL
        array = iarray_s(array, &na, &sa)
        if size is None: size = (toInt(na), toInt(PETSC_DECIDE))
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        if na < n:  raise ValueError(
            "array size %d and vector local size %d block size %d" %
            (toInt(na), toInt(n), toInt(bs)))
        cdef PetscVec newvec = NULL
        if comm_size(ccomm) == 1:
            CHKERR( VecCreateSeqWithArray(ccomm,bs,N,sa,&newvec) )
        else:
            CHKERR( VecCreateMPIWithArray(ccomm,bs,n,N,sa,&newvec) )
        PetscCLEAR(self.obj); self.vec = newvec
        self.set_attr('__array__', array)
        return self

    def createCUDAWithArrays(
        self,
        cpuarray: Sequence[Scalar] | None = None,
        cudahandle: Any | None = None,  # FIXME What type is appropriate here?
        size: tuple[int, int] | int = None,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a vector with type `Vec.Type.CUDA` with optional arrays.

        Collective.

        Parameters
        ----------
        cpuarray
            Host array. Will be lazily allocated if not provided.
        cudahandle
            Address of the array on the GPU. Will be lazily allocated if
            not provided.
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. If `None` then defaults to the size of ``cpuarray``. For
            more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        petsc.VecCreateSeqCUDAWithArrays, petsc.VecCreateMPICUDAWithArrays

        """
        cdef PetscInt na=0
        cdef PetscScalar *sa=NULL
        cdef PetscScalar *gpuarray = NULL
        if cudahandle:
            gpuarray = <PetscScalar*>(<Py_uintptr_t>cudahandle)
        if cpuarray is not None:
            cpuarray = iarray_s(cpuarray, &na, &sa)

        if size is None: size = (toInt(na), toInt(PETSC_DECIDE))
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        if na < n:  raise ValueError(
            "array size %d and vector local size %d block size %d" %
            (toInt(na), toInt(n), toInt(bs)))
        cdef PetscVec newvec = NULL
        if comm_size(ccomm) == 1:
            CHKERR( VecCreateSeqCUDAWithArrays(ccomm,bs,N,sa,gpuarray,&newvec) )
        else:
            CHKERR( VecCreateMPICUDAWithArrays(ccomm,bs,n,N,sa,gpuarray,&newvec) )
        PetscCLEAR(self.obj); self.vec = newvec

        if cpuarray is not None:
            self.set_attr('__array__', cpuarray)
        return self

    def createHIPWithArrays(
        self,
        cpuarray: Sequence[Scalar] | None = None,
        hiphandle: Any | None = None,  # FIXME What type is appropriate here?
        size: tuple[int, int] | int | None = None,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a vector with type `Vec.Type.HIP` with optional arrays.

        Collective.

        Parameters
        ----------
        cpuarray
            Host array. Will be lazily allocated if not provided.
        hiphandle
            Address of the array on the GPU. Will be lazily allocated if
            not provided.
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. If `None` then defaults to the size of ``cpuarray``. For
            more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        petsc.VecCreateSeqHIPWithArrays, petsc.VecCreateMPIHIPWithArrays

        """
        cdef PetscInt na=0
        cdef PetscScalar *sa=NULL
        cdef PetscScalar *gpuarray = NULL
        if hiphandle:
            gpuarray = <PetscScalar*>(<Py_uintptr_t>hiphandle)
        if cpuarray is not None:
            cpuarray = iarray_s(cpuarray, &na, &sa)

        if size is None: size = (toInt(na), toInt(PETSC_DECIDE))
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        if na < n:  raise ValueError(
            "array size %d and vector local size %d block size %d" %
            (toInt(na), toInt(n), toInt(bs)))
        cdef PetscVec newvec = NULL
        if comm_size(ccomm) == 1:
            CHKERR( VecCreateSeqHIPWithArrays(ccomm,bs,N,sa,gpuarray,&newvec) )
        else:
            CHKERR( VecCreateMPIHIPWithArrays(ccomm,bs,n,N,sa,gpuarray,&newvec) )
        PetscCLEAR(self.obj); self.vec = newvec

        if cpuarray is not None:
            self.set_attr('__array__', cpuarray)
        return self

    def createViennaCLWithArrays(
        self,
        cpuarray: Sequence[Scalar] | None = None,
        viennaclvechandle: Any | None = None,  # FIXME What type is appropriate here?
        size: tuple[int, int] | int | None = None,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a vector with type `Vec.Type.VIENNACL` with optional arrays.

        Collective.

        Parameters
        ----------
        cpuarray
            Host array. Will be lazily allocated if not provided.
        viennaclvechandle
            Address of the array on the GPU. Will be lazily allocated if
            not provided.
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. If `None` then defaults to the size of ``cpuarray``. For
            more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        petsc.VecCreateSeqViennaCLWithArrays
        petsc.VecCreateMPIViennaCLWithArrays

        """
        cdef PetscInt na=0
        cdef PetscScalar *sa=NULL
        cdef PetscScalar *vclvec = NULL
        if viennaclvechandle:
            vclvec = <PetscScalar*>(<Py_uintptr_t>viennaclvechandle)
        if cpuarray is not None:
            cpuarray = iarray_s(cpuarray, &na, &sa)

        if size is None: size = (toInt(na), toInt(PETSC_DECIDE))
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        if na < n:  raise ValueError( "array size %d and vector local size %d block size %d" % (toInt(na), toInt(n), toInt(bs)))
        cdef PetscVec newvec = NULL
        if comm_size(ccomm) == 1:
            CHKERR( VecCreateSeqViennaCLWithArrays(ccomm,bs,N,sa,vclvec,&newvec) )
        else:
            CHKERR( VecCreateMPIViennaCLWithArrays(ccomm,bs,n,N,sa,vclvec,&newvec) )
        PetscCLEAR(self.obj); self.vec = newvec

        if cpuarray is not None:
            self.set_attr('__array__', cpuarray)
        return self

    def createWithDLPack(
        self,
        object dltensor,
        size: tuple[int, int] | int | None = None,
        bsize: int | None = None,
        comm: Comm | None = None
    ) -> Self:
        """Create a vector wrapping a DLPack object, sharing the same memory.

        This operation does not modify the storage of the original tensor and
        should be used with contiguous tensors only. If the tensor is stored in
        row-major order (e.g. PyTorch tensors), the resulting vector will look
        like an unrolled tensor using row-major order.

        The resulting vector type will be one of `Vec.Type.SEQ`, `Vec.Type.MPI`,
        `Vec.Type.SEQCUDA`, `Vec.Type.MPICUDA`, `Vec.Type.SEQHIP` or
        `Vec.Type.MPIHIP` depending on the type of ``dltensor`` and the number
        of ranks in the communicator.

        Collective.

        Parameters
        ----------
        dltensor
            Either an object with a ``__dlpack__`` method or a DLPack tensor object.
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. If `None` then defaults to the flattened size of
            ``dltensor``. For more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        """
        cdef DLManagedTensor* ptr = NULL
        cdef int bits = 0
        cdef PetscInt nz = 1
        cdef int64_t ndim = 0
        cdef int64_t* shape = NULL
        cdef int64_t* strides = NULL
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs = 0,n = 0,N = 0
        cdef DLContext* ctx = NULL

        if not PyCapsule_CheckExact(dltensor):
            dltensor = dltensor.__dlpack__()

        if PyCapsule_IsValid(dltensor, 'dltensor'):
            ptr = <DLManagedTensor*>PyCapsule_GetPointer(dltensor, 'dltensor')
            bits = ptr.dl_tensor.dtype.bits
            if bits != 8*sizeof(PetscScalar):
                raise TypeError("Tensor dtype = {} does not match PETSc precision".format(ptr.dl_tensor.dtype))
            ndim = ptr.dl_tensor.ndim
            shape = ptr.dl_tensor.shape
            for s in shape[:ndim]:
                nz = nz*s
            strides = ptr.dl_tensor.strides
            PyCapsule_SetName(dltensor, 'used_dltensor')
        else:
            raise ValueError("Expect a dltensor field, pycapsule.PyCapsule can only be consumed once")
        if size is None: size = (toInt(nz), toInt(PETSC_DECIDE))
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if bs == PETSC_DECIDE: bs = 1
        if nz < n:  raise ValueError(
            "array size %d and vector local size %d block size %d" %
            (toInt(nz), toInt(n), toInt(bs)))
        cdef PetscVec newvec = NULL
        cdef PetscDLDeviceType dltype = ptr.dl_tensor.ctx.device_type
        if dltype in [kDLCUDA,kDLCUDAManaged]:
            if comm_size(ccomm) == 1:
                CHKERR( VecCreateSeqCUDAWithArray(ccomm,bs,N,<PetscScalar*>(ptr.dl_tensor.data),&newvec) )
            else:
                CHKERR( VecCreateMPICUDAWithArray(ccomm,bs,n,N,<PetscScalar*>(ptr.dl_tensor.data),&newvec) )
        elif dltype in [kDLCPU,kDLCUDAHost,kDLROCMHost]:
            if comm_size(ccomm) == 1:
                CHKERR( VecCreateSeqWithArray(ccomm,bs,N,<PetscScalar*>(ptr.dl_tensor.data),&newvec) )
            else:
                CHKERR( VecCreateMPIWithArray(ccomm,bs,n,N,<PetscScalar*>(ptr.dl_tensor.data),&newvec) )
        elif dltype == kDLROCM:
            if comm_size(ccomm) == 1:
                CHKERR( VecCreateSeqHIPWithArray(ccomm,bs,N,<PetscScalar*>(ptr.dl_tensor.data),&newvec) )
            else:
                CHKERR( VecCreateMPIHIPWithArray(ccomm,bs,n,N,<PetscScalar*>(ptr.dl_tensor.data),&newvec) )
        else:
            raise TypeError("Device type {} not supported".format(dltype))

        PetscCLEAR(self.obj); self.vec = newvec
        self.set_attr('__array__', dltensor)
        cdef int64_t* shape_arr = NULL
        cdef int64_t* strides_arr = NULL
        cdef object s1 = oarray_p(empty_p(ndim), NULL, <void**>&shape_arr)
        cdef object s2 = oarray_p(empty_p(ndim), NULL, <void**>&strides_arr)
        for i in range(ndim):
            shape_arr[i] = shape[i]
            strides_arr[i] = strides[i]
        self.set_attr('__dltensor_ctx__', (ptr.dl_tensor.ctx.device_type, ptr.dl_tensor.ctx.device_id, ndim, s1, s2))
        if ptr.manager_deleter != NULL:
            ptr.manager_deleter(ptr) # free the manager
        return self

    def attachDLPackInfo(
        self,
        Vec vec=None,
        object dltensor=None
    ) -> Self:
        """Attach tensor information from another vector or DLPack tensor.

        This tensor information is required when converting a `Vec` to a
        DLPack object.

        Logically collective.

        Parameters
        ----------
        vec
            Vector with attached tensor information. This is typically created
            by calling `Vec.createWithDLPack`.
        dltensor
            DLPack tensor. This will only be used if ``vec`` is `None`.

        Notes
        -----
        This operation does not copy any data from ``vec`` or ``dltensor``.

        See Also
        --------
        Vec.clearDLPackInfo, Vec.createWithDLPack

        """
        cdef object ctx0 = self.get_attr('__dltensor_ctx__'), ctx = None
        cdef DLManagedTensor* ptr = NULL
        cdef int64_t* shape_arr = NULL
        cdef int64_t* strides_arr = NULL
        cdef object s1 = None, s2 = None

        if vec is None and dltensor is None:
            raise ValueError('Missing input parameters')
        if vec is not None:
            ctx = (<Object>vec).get_attr('__dltensor_ctx__')
            if ctx is None:
                raise ValueError('Input vector has no tensor information')
            self.set_attr('__dltensor_ctx__', ctx)
        else:
            if PyCapsule_IsValid(dltensor, 'dltensor'):
                ptr = <DLManagedTensor*>PyCapsule_GetPointer(dltensor, 'dltensor')
            elif PyCapsule_IsValid(dltensor, 'used_dltensor'):
                ptr = <DLManagedTensor*>PyCapsule_GetPointer(dltensor, 'used_dltensor')
            else:
                raise ValueError("Expect a dltensor or used_dltensor field")
            bits = ptr.dl_tensor.dtype.bits
            if bits != 8*sizeof(PetscScalar):
                raise TypeError("Tensor dtype = {} does not match PETSc precision".format(ptr.dl_tensor.dtype))
            ndim = ptr.dl_tensor.ndim
            shape = ptr.dl_tensor.shape
            strides = ptr.dl_tensor.strides
            s1 = oarray_p(empty_p(ndim), NULL, <void**>&shape_arr)
            s2 = oarray_p(empty_p(ndim), NULL, <void**>&strides_arr)
            for i in range(ndim):
                shape_arr[i] = shape[i]
                strides_arr[i] = strides[i]
            self.set_attr('__dltensor_ctx__', (ptr.dl_tensor.ctx.device_type, ptr.dl_tensor.ctx.device_id, ndim, s1, s2))
        return self

    def clearDLPackInfo(self) -> Self:
        """Clear tensor information.

        Logically collective.

        See Also
        --------
        Vec.attachDLPackInfo, Vec.createWithDLPack

        """
        self.set_attr('__dltensor_ctx__', None)
        return self

    # TODO Stream
    def __dlpack__(self, stream=-1):
        return self.toDLPack('rw')

    def __dlpack_device__(self):
        (dltype, devId, _, _, _) = vec_get_dlpack_ctx(self)
        return (dltype, devId)

    # FIXME Not sure what the return type should be
    def toDLPack(self, mode: Literal["rw", "r", "w"] | None = "rw") -> Any:
        """Return a DLPack `PyCapsule` wrapping the vector data.

        Collective.

        Parameters
        ----------
        mode
            Access mode for the vector. Must be read-write (``"rw"``), read
            (``"r"``) or write (``"w"``). If `None` defaults to ``"rw"``.

        Returns
        -------
        `PyCapsule`
            Capsule of a DLPack tensor wrapping a `Vec`.

        Notes
        -----
        It is important that the access mode is respected by the consumer
        as this is not enforced internally.

        See Also
        --------
        Vec.createWithDLPack

        """
        if mode is None: mode = 'rw'
        if mode not in ['rw', 'r', 'w']:
            raise ValueError("Invalid mode: expected 'rw', 'r', or 'w'")

        cdef int64_t ndim = 0
        (device_type, device_id, ndim, shape, strides) = vec_get_dlpack_ctx(self)
        hostmem = (device_type == kDLCPU)

        cdef DLManagedTensor* dlm_tensor = <DLManagedTensor*>malloc(sizeof(DLManagedTensor))
        cdef DLTensor* dl_tensor = &dlm_tensor.dl_tensor
        cdef PetscScalar *a = NULL
        cdef int64_t* shape_strides = NULL
        dl_tensor.byte_offset = 0

        # DLPack does not currently play well with our get/restore model
        # Call restore right-away and hope that the consumer will do the right thing
        # and not modify memory requested with read access
        # By restoring now, we guarantee the sanity of the ObjectState
        if mode == 'w':
            if hostmem:
                CHKERR( VecGetArrayWrite(self.vec, <PetscScalar**>&a) )
                CHKERR( VecRestoreArrayWrite(self.vec, NULL) )
            else:
                CHKERR( VecGetArrayWriteAndMemType(self.vec, <PetscScalar**>&a, NULL) )
                CHKERR( VecRestoreArrayWriteAndMemType(self.vec, NULL) )
        elif mode == 'r':
            if hostmem:
                CHKERR( VecGetArrayRead(self.vec, <const PetscScalar**>&a) )
                CHKERR( VecRestoreArrayRead(self.vec, NULL) )
            else:
                CHKERR( VecGetArrayReadAndMemType(self.vec, <const PetscScalar**>&a, NULL) )
                CHKERR( VecRestoreArrayReadAndMemType(self.vec, NULL) )
        else:
            if hostmem:
                CHKERR( VecGetArray(self.vec, <PetscScalar**>&a) )
                CHKERR( VecRestoreArray(self.vec, NULL) )
            else:
                CHKERR( VecGetArrayAndMemType(self.vec, <PetscScalar**>&a, NULL) )
                CHKERR( VecRestoreArrayAndMemType(self.vec, NULL) )
        dl_tensor.data = <void *>a

        cdef DLContext* ctx = &dl_tensor.ctx
        ctx.device_type = device_type
        ctx.device_id = device_id
        shape_strides = <int64_t*>malloc(sizeof(int64_t)*2*ndim)
        for i in range(ndim):
            shape_strides[i] = shape[i]
        for i in range(ndim):
            shape_strides[i+ndim] = strides[i]
        dl_tensor.ndim = ndim
        dl_tensor.shape = shape_strides
        dl_tensor.strides = shape_strides + ndim

        cdef DLDataType* dtype = &dl_tensor.dtype
        dtype.code = <uint8_t>DLDataTypeCode.kDLFloat
        if sizeof(PetscScalar) == 8:
            dtype.bits = <uint8_t>64
        elif sizeof(PetscScalar) == 4:
            dtype.bits = <uint8_t>32
        else:
            raise ValueError('Unsupported PetscScalar type')
        dtype.lanes = <uint16_t>1
        dlm_tensor.manager_ctx = <void *>self.vec
        CHKERR( PetscObjectReference(<PetscObject>self.vec) )
        dlm_tensor.manager_deleter = manager_deleter
        dlm_tensor.del_obj = <dlpack_manager_del_obj>PetscDEALLOC
        return PyCapsule_New(dlm_tensor, 'dltensor', pycapsule_deleter)

    def createGhost(
        self,
        ghosts: Sequence[int],
        size: tuple[int, int] | int,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a parallel vector with ghost padding on each processor.

        Collective.

        Parameters
        ----------
        ghosts
            Global indices of ghost points. These do not need to be sorted.
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. For more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        Vec.createGhostWithArray, petsc.VecCreateGhost
        petsc.VecCreateGhostBlock

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt ng=0, *ig=NULL
        ghosts = iarray_i(ghosts, &ng, &ig)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        cdef PetscVec newvec = NULL
        if bs == PETSC_DECIDE:
            CHKERR( VecCreateGhost(
                    ccomm, n, N, ng, ig, &newvec) )
        else:
            CHKERR( VecCreateGhostBlock(
                    ccomm, bs, n, N, ng, ig, &newvec) )
        PetscCLEAR(self.obj); self.vec = newvec
        return self

    def createGhostWithArray(
        self,
        ghosts: Sequence[int],
        array: Sequence[Scalar],
        size: tuple[int, int] | int | None = None,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a parallel vector with ghost padding and provided arrays.

        Collective.

        Parameters
        ----------
        ghosts
            Global indices of ghost points. These do not need to be sorted.
        array
            Array to store the vector values. Must be at least as large as
            the local size of the vector (including ghost points).
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. For more information see `Sys.splitOwnership`.
        bsize
            The block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        Vec.createGhost, petsc.VecCreateGhostWithArray
        petsc.VecCreateGhostBlockWithArray

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt ng=0, *ig=NULL
        ghosts = iarray_i(ghosts, &ng, &ig)
        cdef PetscInt na=0
        cdef PetscScalar *sa=NULL
        array = oarray_s(array, &na, &sa)
        cdef PetscInt b = 1 if bsize is None else asInt(bsize)
        if size is None: size = (toInt(na-ng*b), toInt(PETSC_DECIDE))
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        if na < (n+ng*b): raise ValueError(
            "ghosts size %d, array size %d, and "
            "vector local size %d block size %d" %
            (toInt(ng), toInt(na), toInt(n), toInt(b)))
        cdef PetscVec newvec = NULL
        if bs == PETSC_DECIDE:
            CHKERR( VecCreateGhostWithArray(
                    ccomm, n, N, ng, ig, sa, &newvec) )
        else:
            CHKERR( VecCreateGhostBlockWithArray(
                    ccomm, bs, n, N, ng, ig, sa, &newvec) )
        PetscCLEAR(self.obj); self.vec = newvec
        self.set_attr('__array__', array)
        return self

    def createShared(
        self,
        size: tuple[int, int] | int,
        bsize: int | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a parallel vector that uses shared memory.

        Collective.

        Parameters
        ----------
        size
            Global size ``N`` or 2-tuple ``(n, N)`` with local and global
            sizes. For more information see `Sys.splitOwnership`.
        bsize
            Block size, defaults to 1.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        petsc.VecCreateShared

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs=0, n=0, N=0
        Vec_Sizes(size, bsize, &bs, &n, &N)
        Sys_Layout(ccomm, bs, &n, &N)
        cdef PetscVec newvec = NULL
        CHKERR( VecCreateShared(ccomm, n, N, &newvec) )
        PetscCLEAR(self.obj); self.vec = newvec
        if bs != PETSC_DECIDE:
            CHKERR( VecSetBlockSize(self.vec, bs) )
        return self

    # FIXME Is the description of isets correct?
    def createNest(
        self,
        vecs: Sequence[Vec],
        isets: Sequence[IS] = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a vector containing multiple nested subvectors.

        The subvectors are stored separately.

        Collective.

        Parameters
        ----------
        vecs
            Iterable of subvectors.
        isets
            Iterable of index sets describing a reordering for each of the
            nested vectors, defaults to no reordering.
        comm
            MPI communicator, defaults to `Sys.getDefaultComm`.

        See Also
        --------
        petsc.VecCreateNest

        """
        vecs = list(vecs)
        if isets:
            isets = list(isets)
            assert len(isets) == len(vecs)
        else:
            isets = None
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef Py_ssize_t i, m = len(vecs)
        cdef PetscInt n = <PetscInt>m
        cdef PetscVec *cvecs  = NULL
        cdef PetscIS  *cisets = NULL
        cdef object tmp1, tmp2
        tmp1 = oarray_p(empty_p(n), NULL, <void**>&cvecs)
        for i from 0 <= i < m: cvecs[i] = (<Vec?>vecs[i]).vec
        if isets is not None:
            tmp2 = oarray_p(empty_p(n), NULL, <void**>&cisets)
            for i from 0 <= i < m: cisets[i] = (<IS?>isets[i]).iset
        cdef PetscVec newvec = NULL
        CHKERR( VecCreateNest(ccomm, n, cisets, cvecs,&newvec) )
        PetscCLEAR(self.obj); self.vec = newvec
        return self

    #

    def setOptionsPrefix(self, prefix: str) -> None:
        """Set the prefix used to index the vector in the options database.

        Logically collective.

        Parameters
        ----------
        prefix
            Prefix prepended to all options.

        See Also
        --------
        petsc_options, Vec.getOptionsPrefix, petsc.VecSetOptionsPrefix

        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( VecSetOptionsPrefix(self.vec, cval) )

    def getOptionsPrefix(self) -> str:
        """Return the prefix used to index the vector in the options database.

        Not collective.

        See Also
        --------
        petsc_options, Vec.setOptionsPrefix, petsc.VecGetOptionsPrefix

        """
        cdef const char *cval = NULL
        CHKERR( VecGetOptionsPrefix(self.vec, &cval) )
        return bytes2str(cval)

    def appendOptionsPrefix(self, prefix: str) -> None:
        """Append to prefix used to index the vector in the options database.

        Logically collective.

        Parameters
        ----------
        prefix
            Prefix appended to existing options prefix.

        See Also
        --------
        petsc_options, Vec.getOptionsPrefix, Vec.setOptionsPrefix
        petsc.VecAppendOptionsPrefix

        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( VecAppendOptionsPrefix(self.vec, cval) )

    def setFromOptions(self) -> None:
        """Configure the vector using the options database.

        Collective.

        See Also
        --------
        petsc_options, petsc.VecSetFromOptions

        """
        CHKERR( VecSetFromOptions(self.vec) )

    def setUp(self) -> Self:
        """Prepare vector for use.

        Most users will not need to call this directly as it will be called
        automatically.

        Collective.

        See Also
        --------
        Vec.create, Vec.destroy, petsc.VecSetUp

        """
        CHKERR( VecSetUp(self.vec) )
        return self

    def setOption(self, option: Vec.Option | str, flag: bool) -> None:
        """Set a database option to control the vector's behaviour.

        Collective.

        Parameters
        ----------
        option
            The option to set.
        flag
            Whether to enable (`True`) or disable (`False`) the option.

        See Also
        --------
        petsc_options, petsc.VecSetOption

        """
        CHKERR( VecSetOption(self.vec, option, flag) )

    def getType(self) -> str:
        """Return the type of the vector.

        Not collective.

        See Also
        --------
        petsc.VecGetType

        """
        cdef PetscVecType cval = NULL
        CHKERR( VecGetType(self.vec, &cval) )
        return bytes2str(cval)

    def getSize(self) -> int:
        """Return the global number of elements in the vector.

        Not collective.

        See Also
        --------
        Vec.setSizes, Vec.getLocalSize, petsc.VecGetSize

        """
        cdef PetscInt N = 0
        CHKERR( VecGetSize(self.vec, &N) )
        return toInt(N)

    def getLocalSize(self) -> int:
        """Return the number of elements stored locally by the vector.

        Not collective.

        See Also
        --------
        Vec.setSizes, Vec.getSize, petsc.VecGetLocalSize

        """
        cdef PetscInt n = 0
        CHKERR( VecGetLocalSize(self.vec, &n) )
        return toInt(n)

    def getSizes(self) -> tuple[int, int]:
        """Return the local and global vector sizes.

        Not collective.

        Returns
        -------
        local_size : int
            Local number of vector elements.
        global_size : int
            Global number of vector elements.

        See Also
        --------
        Vec.getSize, Vec.getLocalSize, petsc.VecGetLocalSize, petsc.VecGetSize

        """
        cdef PetscInt n = 0, N = 0
        CHKERR( VecGetLocalSize(self.vec, &n) )
        CHKERR( VecGetSize(self.vec, &N) )
        return (toInt(n), toInt(N))

    def setBlockSize(self, bsize: int) -> None:
        """Set the block size.

        Logically collective.

        Parameters
        ----------
        bsize
            Block size.

        See Also
        --------
        petsc.VecSetBlockSize

        """
        cdef PetscInt bs = asInt(bsize)
        CHKERR( VecSetBlockSize(self.vec, bs) )

    def getBlockSize(self) -> int:
        """Return the block size of the vector.

        Not collective.

        See Also
        --------
        petsc.VecGetBlockSize

        """
        cdef PetscInt bs=0
        CHKERR( VecGetBlockSize(self.vec, &bs) )
        return toInt(bs)

    def getOwnershipRange(self) -> tuple[int, int]:
        """Return the locally owned range of indices.

        Not collective.

        Returns
        -------
        low : int
            The first local element.
        high : int
            One more than the last local element.

        Notes
        -----
        This assumes that vectors are laid out with the first ``n1`` elements
        stored on the first processor, then the next ``n2`` elements stored
        on the second, etc. For certain parallel layouts this may not be the
        case in which case this method is not well defined.

        See Also
        --------
        Vec.getOwnershipRanges, petsc.VecGetOwnershipRange

        """
        cdef PetscInt low=0, high=0
        CHKERR( VecGetOwnershipRange(self.vec, &low, &high) )
        return (toInt(low), toInt(high))

    def getOwnershipRanges(self) -> ArrayInt:
        """Return the range of indices owned by each process.

        Not collective.

        Returns
        -------
        ArrayInt
            Array with length one greater than the number of communicator
            ranks storing start and end + 1 indices for each process.

        See Also
        --------
        Vec.getOwnershipRange, petsc.VecGetOwnershipRanges

        """
        cdef const PetscInt *rng = NULL
        CHKERR( VecGetOwnershipRanges(self.vec, &rng) )
        cdef MPI_Comm comm = MPI_COMM_NULL
        CHKERR( PetscObjectGetComm(<PetscObject>self.vec, &comm) )
        cdef int size = -1
        CHKERR( <PetscErrorCode>MPI_Comm_size(comm, &size) )
        return array_i(size+1, rng)

    def createLocalVector(self) -> Vec:
        """Create a vector storing the local portion of the current vector.

        `Vec.destroy` must be called when this vector is no longer needed.

        Not collective.

        Returns
        -------
        Vec
            The local vector.

        See Also
        --------
        Vec.getLocalVector, petsc.VecCreateLocalVector

        """
        lvec = Vec()
        CHKERR( VecCreateLocalVector(self.vec, &lvec.vec) )
        return lvec

    def getLocalVector(self, Vec lvec, readonly: bool = False) -> None:
        """Load the local portion of the vector.

        `Vec.restoreLocalVector` (with the same ``readonly`` value) should
        be called once the local vector is no longer needed.

        Not collective if ``readonly`` is `True`, collective otherwise.

        Parameters
        ----------
        lvec
            Local vector.
        readonly
            Whether the local vector is only read. If `True` then the local
            vector will preserve cached information.

        See Also
        --------
        Vec.createLocalVector, Vec.restoreLocalVector
        petsc.VecGetLocalVectorRead, petsc.VecGetLocalVector

        """
        if readonly:
            CHKERR( VecGetLocalVectorRead(self.vec, lvec.vec) )
        else:
            CHKERR( VecGetLocalVector(self.vec, lvec.vec) )

    def restoreLocalVector(self, Vec lvec, readonly: bool = False) -> None:
        """Unmap a local vector mapping the local portion of the vector.

        Not collective if ``readonly`` is `True`, logically collective otherwise.

        Parameters
        ----------
        lvec
            Local vector. `Vec.getLocalVector` should already have been
            called with it.
        readonly
            Whether the local vector is only read. If `True` then the local
            vector will preserve cached information.

        See Also
        --------
        Vec.createLocalVector, Vec.getLocalVector
        petsc.VecRestoreLocalVectorRead, petsc.VecRestoreLocalVector

        """
        if readonly:
            CHKERR( VecRestoreLocalVectorRead(self.vec, lvec.vec) )
        else:
            CHKERR( VecRestoreLocalVector(self.vec, lvec.vec) )

    # FIXME Return type should be more specific
    def getBuffer(self, readonly: bool = False) -> Any:
        """Return a buffered view of the local portion of the vector.

        Not collective if ``readonly`` is `True`, logically collective
        otherwise.

        Parameters
        ----------
        readonly
            Whether the vector is only read. If `True` then the vector will
            preserve cached information.

        Returns
        -------
        typing.Any
            `Buffer object <python:c-api/buffer>` wrapping the local portion of
            the vector data. This can be used either as a context manager
            providing access as a numpy array or can be passed to array
            constructors accepting buffered objects such as `numpy.asarray`.

        Examples
        --------
        Accessing the data with a context manager:

        >>> vec = PETSc.Vec().createWithArray([1, 2, 3])
        >>> with vec.getBuffer() as arr:
        ...     arr
        array([1., 2., 3.])

        Converting the buffer to an `ndarray`:

        >>> buf = PETSc.Vec().createWithArray([1, 2, 3]).getBuffer()
        >>> np.asarray(buf)
        array([1., 2., 3.])

        See Also
        --------
        Vec.getArray

        """
        if readonly:
            return vec_getbuffer_r(self)
        else:
            return vec_getbuffer_w(self)

    def getArray(self, readonly: bool=False) -> ArrayScalar:
        """Return local portion of the vector as an `ndarray`.

        Not collective if ``readonly`` is `True`, logically collective
        otherwise.

        Parameters
        ----------
        readonly
            Whether the vector is only read. If `True` then the vector will
            preserve cached information.

        See Also
        --------
        Vec.setArray, Vec.getBuffer

        """
        if readonly:
            return vec_getarray_r(self)
        else:
            return vec_getarray_w(self)

    def setArray(self, array: Sequence[Scalar]) -> None:
        """Set the local portion of the vector.

        This will fail if ``array`` has a different size to the local portion
        of the vector.

        Logically collective.

        Parameters
        ----------
        array
            Values to set the local portion of the vector to. These will be
            copied.

        See Also
        --------
        Vec.placeArray

        """
        vec_setarray(self, array)

    def placeArray(self, array: Sequence[Scalar]) -> None:
        """Set the local portion of the vector to a provided array.

        This method can be used instead of `Vec.setArray` to avoid copying
        data.

        The original array can be returned to using `Vec.resetArray`. It is
        the user's responsibility to free the provided array.

        Not collective.

        Parameters
        ----------
        array
            Data to set as the local portion of the vector. An error will be
            raised if the size does not match `Vec.getLocalSize`.

        See Also
        --------
        Vec.resetArray, Vec.setArray, petsc.VecPlaceArray

        """
        cdef PetscInt nv=0
        cdef PetscInt na=0
        cdef PetscScalar *a = NULL
        CHKERR( VecGetLocalSize(self.vec, &nv) )
        array = oarray_s(array, &na, &a)
        if (na != nv): raise ValueError(
            "cannot place input array size %d, vector size %d" %
            (toInt(na), toInt(nv)))
        CHKERR( VecPlaceArray(self.vec, a) )
        self.set_attr('__placed_array__', array)

    def resetArray(self, force: bool = False) -> ArrayScalar | None:
        """Reset the vector to use its default array for its elements.

        Not collective.

        Parameters
        ----------
        force
            Force the calling of `petsc.VecResetArray` even if no user array
            has been placed with `Vec.placeArray`.

        Returns
        -------
        ArrayScalar
            The array previously provided by the user with `Vec.placeArray`.
            Can be `None` if ``force`` is `True` and no array was placed
            before.

        See Also
        --------
        Vec.placeArray, petsc.VecResetArray

        """
        cdef object array = None
        array = self.get_attr('__placed_array__')
        if array is None and not force: return None
        CHKERR( VecResetArray(self.vec) )
        self.set_attr('__placed_array__', None)
        return array

    def bindToCPU(self, flg: bool) -> None:
        """Indicate to the vector that it will only be accessed on the CPU.

        Logically collective.

        Parameters
        ----------
        flg
            If `True` then subsequent operations must be performed on the CPU.
            If `False` then they must all be performed on the device, assuming
            that the vector type (`Vec.Type`) is capable of offloading.

        See Also
        --------
        Vec.boundToCPU, petsc.VecBindToCPU

        """
        cdef PetscBool bindFlg = asBool(flg)
        CHKERR( VecBindToCPU(self.vec, bindFlg) )

    def boundToCPU(self) -> bool:
        """Return whether the vector has been bound to the CPU.

        See Also
        --------
        Vec.bindToCPU, petsc.VecBoundToCPU

        """
        cdef PetscBool flg = PETSC_TRUE
        CHKERR( VecBoundToCPU(self.vec, &flg) )
        return toBool(flg)

    def getCUDAHandle(
        self,
        mode: Literal["rw", "r", "w"] | None = "rw",
    ) -> Any:  # FIXME What is the right return type?
        """Return a pointer to the CUDA buffer inside the vector.

        The returned pointer should be released using `Vec.restoreCUDAHandle`
        with the same access mode.

        Not collective.

        Parameters
        ----------
        mode
            Access mode for the vector. Must be read-write (``"rw"``), read
            (``"r"``) or write (``"w"``). If `None` defaults to ``"rw"``.

        Returns
        -------
        typing.Any
            CUDA device pointer.

        Notes
        -----
        This method may incur a host-to-device copy if the device data is
        out of date and ``mode`` is ``"r"`` or ``"rw"``.

        See Also
        --------
        Vec.restoreCUDAHandle, petsc.VecCUDAGetArray
        petsc.VecCUDAGetArrayRead, petsc.VecCUDAGetArrayWrite

        """
        cdef PetscScalar *hdl = NULL
        cdef const char *m = NULL
        if mode is not None: mode = str2bytes(mode, &m)
        if m == NULL or (m[0] == c'r' and m[1] == c'w'):
            CHKERR( VecCUDAGetArray(self.vec, &hdl) )
        elif m[0] == c'r':
            CHKERR( VecCUDAGetArrayRead(self.vec, <const PetscScalar**>&hdl) )
        elif m[0] == c'w':
            CHKERR( VecCUDAGetArrayWrite(self.vec, &hdl) )
        else:
            raise ValueError("Invalid mode: expected 'rw', 'r', or 'w'")
        return <Py_uintptr_t>hdl

    def restoreCUDAHandle(
        self,
        handle: Any,  # FIXME What type hint is appropriate?
        mode: Literal["rw", "r", "w"] | None = "rw",
    ) -> None:
        """Restore a pointer to the CUDA buffer inside the vector.

        The pointer should have been obtained by calling `Vec.getCUDAHandle`
        with the same access mode.

        Not collective.

        Parameters
        ----------
        handle
            CUDA device pointer.
        mode
            Access mode for the vector. Must be read-write (``"rw"``), read
            (``"r"``) or write (``"w"``). If `None` defaults to ``"rw"``.

        Notes
        -----
        This method will mark the host data as out of date if ``mode`` is
        ``"w"`` or ``"rw"``, resulting in a host-to-device copy the next time
        data is accessed on the host.

        See Also
        --------
        Vec.getCUDAHandle, petsc.VecCUDARestoreArray
        petsc.VecCUDARestoreArrayRead, petsc.VecCUDARestoreArrayWrite

        """
        cdef PetscScalar *hdl = <PetscScalar*>(<Py_uintptr_t>handle)
        cdef const char *m = NULL
        if mode is not None: mode = str2bytes(mode, &m)
        if m == NULL or (m[0] == c'r' and m[1] == c'w'):
            CHKERR( VecCUDARestoreArray(self.vec, &hdl) )
        elif m[0] == c'r':
            CHKERR( VecCUDARestoreArrayRead(self.vec, <const PetscScalar**>&hdl) )
        elif m[0] == c'w':
            CHKERR( VecCUDARestoreArrayWrite(self.vec, &hdl) )
        else:
            raise ValueError("Invalid mode: expected 'rw', 'r', or 'w'")

    def getHIPHandle(
        self,
        mode: Literal["rw", "r", "w"] | None = "rw",
    ) -> Any:  # FIXME What is the right return type?
        """Return a pointer to the HIP buffer inside the vector.

        The returned pointer should be released using `Vec.restoreHIPHandle`
        with the same access mode.

        Not collective.

        Parameters
        ----------
        mode
            Access mode for the vector. Must be read-write (``"rw"``), read
            (``"r"``) or write (``"w"``). If `None` defaults to ``"rw"``.

        Returns
        -------
        typing.Any
            HIP device pointer.

        Notes
        -----
        This method may incur a host-to-device copy if the device data is
        out of date and ``mode`` is ``"r"`` or ``"rw"``.

        See Also
        --------
        Vec.restoreHIPHandle, petsc.VecHIPGetArray, petsc.VecHIPGetArrayRead
        petsc.VecHIPGetArrayWrite

        """
        cdef PetscScalar *hdl = NULL
        cdef const char *m = NULL
        if mode is not None: mode = str2bytes(mode, &m)
        if m == NULL or (m[0] == c'r' and m[1] == c'w'):
            CHKERR( VecHIPGetArray(self.vec, &hdl) )
        elif m[0] == c'r':
            CHKERR( VecHIPGetArrayRead(self.vec, <const PetscScalar**>&hdl) )
        elif m[0] == c'w':
            CHKERR( VecHIPGetArrayWrite(self.vec, &hdl) )
        else:
            raise ValueError("Invalid mode: expected 'rw', 'r', or 'w'")
        return <Py_uintptr_t>hdl

    def restoreHIPHandle(
        self,
        handle: Any,  # FIXME What type hint is appropriate?
        mode: Literal["rw", "r", "w"] | None = "rw",
    ) -> None:
        """Restore a pointer to the HIP buffer inside the vector.

        The pointer should have been obtained by calling `Vec.getHIPHandle`
        with the same access mode.

        Not collective.

        Parameters
        ----------
        handle
            HIP device pointer.
        mode
            Access mode for the vector. Must be read-write (``"rw"``), read
            (``"r"``) or write (``"w"``). If `None` defaults to ``"rw"``.

        Notes
        -----
        This method will mark the host data as out of date if ``mode`` is
        ``"w"`` or ``"rw"``, resulting in a host-to-device copy the next time
        data is accessed on the host.

        See Also
        --------
        Vec.getHIPHandle, petsc.VecHIPRestoreArray
        petsc.VecHIPRestoreArrayRead, petsc.VecHIPRestoreArrayWrite

        """
        cdef PetscScalar *hdl = <PetscScalar*>(<Py_uintptr_t>handle)
        cdef const char *m = NULL
        if mode is not None: mode = str2bytes(mode, &m)
        if m == NULL or (m[0] == c'r' and m[1] == c'w'):
            CHKERR( VecHIPRestoreArray(self.vec, &hdl) )
        elif m[0] == c'r':
            CHKERR( VecHIPRestoreArrayRead(self.vec, <const PetscScalar**>&hdl) )
        elif m[0] == c'w':
            CHKERR( VecHIPRestoreArrayWrite(self.vec, &hdl) )
        else:
            raise ValueError("Invalid mode: expected 'rw', 'r', or 'w'")

    def getOffloadMask(self) -> int:
        """Return the offloading status of the vector.

        Not collective.

        Common return values include:

        - 1: ``PETSC_OFFLOAD_CPU`` - CPU has valid entries
        - 2: ``PETSC_OFFLOAD_GPU`` - GPU has valid entries
        - 3: ``PETSC_OFFLOAD_BOTH`` - both CPU and GPU have valid entries

        Returns
        -------
        int
            Enum value from `petsc.PetscOffloadMask` describing the offloading
            status.

        See Also
        --------
        petsc.VecGetOffloadMask, petsc.PetscOffloadMask

        """
        cdef PetscOffloadMask mask
        CHKERR( VecGetOffloadMask(self.vec, &mask) )
        return mask

    def getCLContextHandle(self) -> int:
        """Return the OpenCL context associated with the vector.

        Not collective.

        Returns
        -------
        int
            Pointer to underlying CL context. This can be used with
            `pyopencl` through `pyopencl.Context.from_int_ptr`.

        See Also
        --------
        Vec.getCLQueueHandle, petsc.VecViennaCLGetCLContext

        """
        cdef Py_uintptr_t ctxhdl = 0
        CHKERR( VecViennaCLGetCLContext(self.vec, &ctxhdl) )
        return ctxhdl

    def getCLQueueHandle(self) -> int:
        """Return the OpenCL command queue associated with the vector.

        Not collective.

        Returns
        -------
        int
            Pointer to underlying CL command queue. This can be used with
            `pyopencl` through `pyopencl.Context.from_int_ptr`.

        See Also
        --------
        Vec.getCLContextHandle, petsc.VecViennaCLGetCLQueue

        """
        cdef Py_uintptr_t queuehdl = 0
        CHKERR( VecViennaCLGetCLQueue(self.vec, &queuehdl) )
        return queuehdl

    def getCLMemHandle(
        self,
        mode: Literal["rw", "r", "w"] = "rw",
    ) -> int:
        """Return the OpenCL buffer associated with the vector.

        Not collective.

        Parameters
        ----------
        mode
            Access mode for the vector. Must be read-write (``"rw"``), read
            (``"r"``) or write (``"w"``).

        Returns
        -------
        int
            Pointer to the device buffer. This can be used with
            `pyopencl` through `pyopencl.Context.from_int_ptr`.

        Notes
        -----
        This method may incur a host-to-device copy if the device data is
        out of date and ``mode`` is ``"r"`` or ``"rw"``.

        See Also
        --------
        Vec.restoreCLMemHandle, petsc.VecViennaCLGetCLMem
        petsc.VecViennaCLGetCLMemRead, petsc.VecViennaCLGetCLMemWrite

        """
        cdef Py_uintptr_t memhdl = 0
        cdef const char *m = NULL
        mode = str2bytes(mode, &m)
        if m == NULL or (m[0] == c'r' and m[1] == c'w'):
            CHKERR( VecViennaCLGetCLMem(self.vec, &memhdl) )
        elif m[0] == c'r':
            CHKERR( VecViennaCLGetCLMemRead(self.vec, &memhdl) )
        elif m[0] == c'w':
            CHKERR( VecViennaCLGetCLMemWrite(self.vec, &memhdl) )
        else:
            raise ValueError("Invalid mode: expected 'r', 'w' or 'rw'")
        return memhdl

    def restoreCLMemHandle(self) -> None:
        """Restore a pointer to the OpenCL buffer inside the vector.

        This method only needs to be called after accessing the buffer
        (with `Vec.getCLMemHandle`) with ``"w"`` or ``"rw"`` modes.

        Not collective.

        Notes
        -----
        This method will mark the host data as out of date and will cause
        a host-to-device copy the next time data is accessed on the host.

        See Also
        --------
        Vec.getCLMemHandle, petsc.VecViennaCLRestoreCLMemWrite

        """
        CHKERR( VecViennaCLRestoreCLMemWrite(self.vec) )

    def duplicate(self, array: Sequence[Scalar] | None = None) -> Vec:
        """Create a new vector with the same type, optionally with data.

        Collective.

        Parameters
        ----------
        array
            Values to store in the new vector. The size must match the local
            size of the current vector. If not provided then the new vector
            will be empty.

        Notes
        -----
        This method will *not* copy the vector entries to the new vector. If
        ``array`` is provided then these values will be *copied*.

        See Also
        --------
        Vec.copy, petsc.VecDuplicate

        """
        cdef Vec vec = type(self)()
        CHKERR( VecDuplicate(self.vec, &vec.vec) )
        # duplicate tensor context
        cdef object ctx0 = self.get_attr('__dltensor_ctx__')
        if ctx0 is not None:
            vec.set_attr('__dltensor_ctx__', ctx0)
        if array is not None:
            vec_setarray(vec, array)
        return vec

    def copy(self, Vec result=None) -> Vec:
        """Return a copy of the vector.

        This operation copies vector entries to the new vector.

        Logically collective.

        Parameters
        ----------
        result
            Target vector for the copy. If `None` then a new vector is
            allocated with `Vec.duplicate`.

        Returns
        -------
        Vec
            The new vector. If ``result`` is not `None` then it will be
            returned.

        See Also
        --------
        Vec.duplicate, petsc.VecCopy

        """
        if result is None:
            result = type(self)()
        if result.vec == NULL:
            CHKERR( VecDuplicate(self.vec, &result.vec) )
        CHKERR( VecCopy(self.vec, result.vec) )
        return result

    def chop(self, tol: float) -> None:
        """Set all vector entries less than some tolerance to zero.

        Parameters
        ----------
        tol
            The tolerance below which entries are set to zero.

        See Also
        --------
        petsc.VecChop

        """
        cdef PetscReal rval = asReal(tol)
        CHKERR( VecChop(self.vec, rval) )

    def load(self, Viewer viewer) -> Self:
        """Load a vector that has been stored in a binary format.

        The vector must have been stored with `Vec.view` with a viewer with
        type `Viewer.Type.BINARY` or `Viewer.Type.HDF5`.

        Collective

        Parameters
        ----------
        viewer
            Binary file viewer, either `Viewer.Type.BINARY` or
            `Viewer.Type.HDF5`.

        Notes
        -----
        Vector type defaults to either `Vec.Type.SEQ` or `Vec.Type.MPI`. To
        load other types `Vec.setType` or `Vec.setFromOptions` should be
        called in advance.

        See Also
        --------
        Vec.view, petsc.VecLoad

        """
        cdef MPI_Comm comm = MPI_COMM_NULL
        cdef PetscObject obj = <PetscObject>(viewer.vwr)
        if self.vec == NULL:
            CHKERR( PetscObjectGetComm(obj, &comm) )
            CHKERR( VecCreate(comm, &self.vec) )
        CHKERR( VecLoad(self.vec, viewer.vwr) )
        return self

    def equal(self, Vec vec) -> bool:
        """Return whether the vector is equal to another.

        Collective.

        Parameters
        ----------
        vec
            Vector to compare with.

        Returns
        -------
        bool
            Whether the vectors are equal or not. The vectors are considered
            equal if they point to the same memory buffer or are bitwise
            identical with the same local and global layouts. Rounding errors
            are not taken into account.

        See Also
        --------
        petsc.VecEqual

        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( VecEqual(self.vec, vec.vec, &flag) )
        return toBool(flag)

    def dot(self, Vec vec) -> Scalar:
        """Return the dot product with ``vec``.

        For complex numbers this computes yᴴ·x with ``self`` as x, ``vec``
        as y and where yᴴ denotes the conjugate transpose of y.

        Use `Vec.tDot` for the indefinite form yᵀ·x where yᵀ denotes the
        transpose of y.

        Collective.

        Parameters
        ----------
        vec
            Vector to compute the dot product with.

        See Also
        --------
        Vec.dotBegin, Vec.dotEnd, Vec.tDot, petsc.VecDot

        """
        cdef PetscScalar sval = 0
        CHKERR( VecDot(self.vec, vec.vec, &sval) )
        return toScalar(sval)

    def dotBegin(self, Vec vec) -> None:
        """Begin computing the dot product.

        This should be paired with a call to `Vec.dotEnd`.

        Parameters
        ----------
        vec
            Vector to compute the dot product with.

        See Also
        --------
        Vec.dotEnd, Vec.dot, petsc.VecDotBegin

        """
        cdef PetscScalar sval = 0
        CHKERR( VecDotBegin(self.vec, vec.vec, &sval) )

    def dotEnd(self, Vec vec) -> Scalar:
        """Finish computing the dot product.

        Parameters
        ----------
        vec
            Vector to compute the dot product with.

        Returns
        -------
        Scalar
            The dot product.

        See Also
        --------
        Vec.dotBegin, Vec.dot, petsc.VecDotEnd

        """
        cdef PetscScalar sval = 0
        CHKERR( VecDotEnd(self.vec, vec.vec, &sval) )
        return toScalar(sval)

    def tDot(self, Vec vec) -> Scalar:
        """Return the indefinite dot product with ``vec``.

        This computes yᵀ·x with ``self`` as x, ``vec``
        as y and where yᵀ denotes the transpose of y.

        Use `Vec.dot` for the inner product yᴴ·x where yᴴ denotes the
        conjugate transpose of y.

        Collective.

        Parameters
        ----------
        vec
            Vector to compute the indefinite dot product with.

        Returns
        -------
        Scalar
            The indefinite dot product.

        See Also
        --------
        Vec.tDotBegin, Vec.tDotEnd, Vec.dot, petsc.VecTDot

        """
        cdef PetscScalar sval = 0
        CHKERR( VecTDot(self.vec, vec.vec, &sval) )
        return toScalar(sval)

    def tDotBegin(self, Vec vec) -> None:
        """Begin computing the indefinite dot product.

        This should be paired with a call to `Vec.tDotEnd`.

        Parameters
        ----------
        vec
            Vector to compute the indefinite dot product with.

        See Also
        --------
        Vec.tDotEnd, Vec.tDot, petsc.VecTDotBegin

        """
        cdef PetscScalar sval = 0
        CHKERR( VecTDotBegin(self.vec, vec.vec, &sval) )

    def tDotEnd(self, Vec vec) -> Scalar:
        """Finish computing the indefinite dot product.

        Parameters
        ----------
        vec
            Vector to compute the indefinite dot product with.

        Returns
        -------
        Scalar
            The indefinite dot product.

        See Also
        --------
        Vec.tDotBegin, Vec.tDot, petsc.VecTDotEnd

        """
        cdef PetscScalar sval = 0
        CHKERR( VecTDotEnd(self.vec, vec.vec, &sval) )
        return toScalar(sval)

    def mDot(self, vecs, out=None):
        """Not implemented."""
        <void>self; <void>vecs; <void>out; # unused
        raise NotImplementedError

    def mDotBegin(self, vecs, out=None):
        """Not implemented."""
        <void>self; <void>vecs; <void>out; # unused
        raise NotImplementedError

    def mDotEnd(self, vecs, out=None):
        """Not implemented."""
        <void>self; <void>vecs; <void>out; # unused
        raise NotImplementedError

    def mtDot(self, vecs, out=None):
        """Not implemented."""
        <void>self; <void>vecs; <void>out; # unused
        raise NotImplementedError

    def mtDotBegin(self, vecs, out=None):
        """Not implemented."""
        <void>self; <void>vecs; <void>out; # unused
        raise NotImplementedError

    def mtDotEnd(self, vecs, out=None):
        """Not implemented."""
        <void>self; <void>vecs; <void>out; # unused
        raise NotImplementedError

    def norm(
        self,
        norm_type: NormType | int | None = None,
    ) -> float | tuple[float, float]:
        """Compute the vector norm.

        Collective.

        Parameters
        ----------
        norm_type
            The type of norm requested. Possible values (assuming ``self`` as
            x) include:

            - `NormType.NORM_1` Compute Σₙ abs(xₙ)

            - `NormType.NORM_2` Compute √(Σₙ abs(xₙ)²)

            - `NormType.NORM_INFINITY` Compute maxₙ abs(xₙ)

            - `NormType.NORM_1_AND_2` Compute both `NormType.NORM_1` and
              `NormType.NORM_2`.

            If `None`, defaults to `NormType.NORM_2`.

        Returns
        -------
        typing.Any
            The computed norm. A 2-tuple is returned if `NormType.NORM_1_AND_2`
            is specified.

        See Also
        --------
        NormType, petsc.VecNorm, petsc.NormType

        """
        cdef PetscNormType norm_1_2 = PETSC_NORM_1_AND_2
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal rval[2]
        CHKERR( VecNorm(self.vec, ntype, rval) )
        if ntype != norm_1_2: return toReal(rval[0])
        else: return (toReal(rval[0]), toReal(rval[1]))

    def normBegin(
        self,
        norm_type: NormType | int | None = None,
    ) -> None:
        """Begin computing the vector norm.

        This should be paired with a call to `Vec.normEnd`.

        Parameters
        ----------
        norm_type
            The type of norm to compute, defaults to `NormType.NORM_2`. See
            `Vec.norm` for more information.

        See Also
        --------
        Vec.normEnd, Vec.norm, petsc.VecNormBegin

        """
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal dummy[2]
        CHKERR( VecNormBegin(self.vec, ntype, dummy) )

    def normEnd(
        self,
        norm_type: NormType | int | None = None,
    ) -> float | tuple[float, float]:
        """Finish computing the vector norm.

        `Vec.normBegin` should already have been called.

        Parameters
        ----------
        norm_type
            The type of norm to compute, defaults to `NormType.NORM_2`. See
            `Vec.norm` for more information.

        Returns
        -------
        typing.Any
            The computed norm. A 2-tuple is returned if `NormType.NORM_1_AND_2`
            is specified.

        See Also
        --------
        Vec.normBegin, Vec.norm, petsc.VecNormEnd

        """
        cdef PetscNormType norm_1_2 = PETSC_NORM_1_AND_2
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal rval[2]
        CHKERR( VecNormEnd(self.vec, ntype, rval) )
        if ntype != norm_1_2: return toReal(rval[0])
        else: return (toReal(rval[0]), toReal(rval[1]))

    def sum(self) -> Scalar:
        """Compute the sum of all the entries of the vector.

        Collective.

        See Also
        --------
        petsc.VecSum

        """
        cdef PetscScalar sval = 0
        CHKERR( VecSum(self.vec, &sval) )
        return toScalar(sval)

    def min(self) -> tuple[int, Scalar]:
        """Return the entry in the vector with minimum real part.

        Collective.

        Returns
        -------
        p : int
            Location of the minimum value. If multiple entries exist with the
            same value then the smallest index will be returned.
        val : Scalar
            Minimum value.

        Notes
        -----
        Returns ``PETSC_MAX_REAL`` for ``val`` and negative ``p`` if the
        vector has length ``0``.

        See Also
        --------
        Vec.max, petsc.VecMin

        """
        cdef PetscInt  ival = 0
        cdef PetscReal rval = 0
        CHKERR( VecMin(self.vec, &ival, &rval) )
        return (toInt(ival), toReal(rval))

    def max(self) -> tuple[int, Scalar]:
        """Return the entry in the vector with maximum real part.

        Collective.

        Returns
        -------
        p : int
            Location of the maximum value. If multiple entries exist with the
            same value then the smallest index will be returned.
        val : Scalar
            Minimum value.

        Notes
        -----
        Returns ``PETSC_MIN_REAL`` for ``val`` and negative ``p`` if the
        vector has length ``0``.

        See Also
        --------
        Vec.min, petsc.VecMax

        """
        cdef PetscInt  ival = 0
        cdef PetscReal rval = 0
        CHKERR( VecMax(self.vec, &ival, &rval) )
        return (toInt(ival), toReal(rval))

    def normalize(self) -> float:
        """Normalize the vector by its 2-norm.

        Collective.

        Returns
        -------
        float
            The vector norm before normalization.

        See Also
        --------
        Vec.norm, petsc.VecNormalize

        """
        cdef PetscReal rval = 0
        CHKERR( VecNormalize(self.vec, &rval) )
        return toReal(rval)

    def reciprocal(self) -> None:
        """Replace each entry in the vector by its reciprocal.

        Logically collective.

        See Also
        --------
        petsc.VecReciprocal

        """
        CHKERR( VecReciprocal(self.vec) )

    def exp(self) -> None:
        """Replace each entry (xₙ) in the vector by exp(xₙ).

        Not collective.

        See Also
        --------
        Vec.log, petsc.VecExp

        """
        CHKERR( VecExp(self.vec) )

    def log(self) -> None:
        """Replace each entry in the vector by its natural logarithm.

        Not collective.

        See Also
        --------
        Vec.exp, petsc.VecLog

        """
        CHKERR( VecLog(self.vec) )

    def sqrtabs(self) -> None:
        """Replace each entry (xₙ) in the vector by √abs(xₙ).

        Not collective.

        See Also
        --------
        petsc.VecSqrtAbs

        """
        CHKERR( VecSqrtAbs(self.vec) )

    def abs(self) -> None:
        """Replace each entry (xₙ) in the vector by abs(xₙ).

        Logically collective.

        See Also
        --------
        petsc.VecAbs

        """
        CHKERR( VecAbs(self.vec) )

    def conjugate(self):
        """Conjugates the vector.

        Logically collective.

        See Also
        --------
        petsc.VecConjugate

        """
        CHKERR( VecConjugate(self.vec) )

    def setRandom(self, Random random=None) -> None:
        """Set all components of the vector to random numbers.

        Logically collective.

        Parameters
        ----------
        random
            Random number generator. If `None` then one will be created
            internally.

        See Also
        --------
        petsc.VecSetRandom

        """
        cdef PetscRandom rnd = NULL
        if random is not None: rnd = random.rnd
        CHKERR( VecSetRandom(self.vec, rnd) )

    def permute(self, IS order, invert: bool = False) -> None:
        """Permute the vector in-place with a provided ordering.

        Parameters
        ----------
        order
            Ordering for the permutation.
        invert
            Whether to invert the permutation.

        Notes
        -----
        Parallel index sets with non-local permutations are not currently
        supported.

        See Also
        --------
        petsc.VecPermute

        """
        cdef PetscBool cinvert = PETSC_FALSE
        if invert: cinvert = PETSC_TRUE
        CHKERR( VecPermute(self.vec, order.iset, cinvert) )

    def zeroEntries(self) -> None:
        """Set all entries in the vector to zero.

        Logically collective.

        See Also
        --------
        Vec.set, petsc.VecZeroEntries

        """
        CHKERR( VecZeroEntries(self.vec) )

    def set(self, alpha: Scalar) -> None:
        """Set all components of the vector to the same value.

        Logically collective.

        Parameters
        ----------
        alpha
            Value to set all vector entries to.

        Notes
        -----
        This method should not be called between `Vec.setValues` and
        `Vec.assemblyBegin`.

        See Also
        --------
        Vec.zeroEntries, Vec.isset, petsc.VecSet

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecSet(self.vec, sval) )

    def isset(self, IS idx, alpha: Scalar) -> None:
        """Set specific elements of the vector to the same value.

        Parameters
        ----------
        idx
            Index set specifying the vector entries to set.
        alpha
            Value to set the selected entries to.

        See Also
        --------
        Vec.set, Vec.zeroEntries, petsc.VecISSet

        """
        cdef PetscScalar aval = asScalar(alpha)
        CHKERR( VecISSet(self.vec, idx.iset, aval) )

    def scale(self, alpha: Scalar) -> None:
        """Scale all entries of the vector by some value.

        This method sets each entry (xₙ) in the vector to ɑ·xₙ.

        Not collective.

        Parameters
        ----------
        alpha
            The scaling factor.

        See Also
        --------
        Vec.shift, petsc.VecScale

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecScale(self.vec, sval) )

    def shift(self, alpha: Scalar) -> None:
        """Shift all entries in the vector.

        This method sets each entry (xₙ) in the vector to xₙ + ɑ.

        Logically collective.

        Parameters
        ----------
        alpha
            The shift to apply to the vector values.

        See Also
        --------
        Vec.scale, petsc.VecShift

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecShift(self.vec, sval) )

    def swap(self, Vec vec) -> None:
        """Swap the data stored by two vectors.

        Logically collective.

        Parameters
        ----------
        vec
            The vector to swap data with.

        See Also
        --------
        petsc.VecSwap

        """
        CHKERR( VecSwap(self.vec, vec.vec) )

    def axpy(self, alpha: Scalar, Vec x) -> None:
        """Compute and store y = ɑ·x + y.

        Logically collective.

        Parameters
        ----------
        alpha
            Scale factor.
        x
            Input vector, must not be the current vector.

        See Also
        --------
        Vec.isaxpy, petsc.VecAXPY

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecAXPY(self.vec, sval, x.vec) )

    def isaxpy(self, IS idx, alpha: Scalar, Vec x) -> None:
        """Add a scaled reduced-space vector to a subset of the vector.

        This is equivalent to ``y[idx[i]] += alpha*x[i]``.

        Parameters
        ----------
        idx
            Index set for the reduced space. Negative indices are skipped.
        alpha
            Scale factor.
        x
            Reduced-space vector.

        See Also
        --------
        Vec.axpy, Vec.aypx, Vec.axpby, petsc.VecISAXPY

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecISAXPY(self.vec, idx.iset, sval, x.vec) )

    def aypx(self, alpha: Scalar, Vec x) -> None:
        """Compute and store y = x + ɑ·y.

        Logically collective.

        Parameters
        ----------
        alpha
            Scale factor.
        x
            Input vector, must not be the current vector.

        See Also
        --------
        Vec.axpy, Vec.axpby, petsc.VecAYPX

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecAYPX(self.vec, sval, x.vec) )

    def axpby(self, alpha: Scalar, beta: Scalar, Vec x) -> None:
        """Compute and store y = ɑ·x + β·y.

        Logically collective.

        Parameters
        ----------
        alpha
            First scale factor.
        beta
            Second scale factor.
        x
            Input vector, must not be the current vector.

        See Also
        --------
        Vec.axpy, Vec.aypx, Vec.waxpy, petsc.VecAXPBY

        """
        cdef PetscScalar sval1 = asScalar(alpha)
        cdef PetscScalar sval2 = asScalar(beta)
        CHKERR( VecAXPBY(self.vec, sval1, sval2, x.vec) )

    def waxpy(self, alpha: Scalar, Vec x, Vec y) -> None:
        """Compute and store w = ɑ·x + y.

        Logically collective.

        Parameters
        ----------
        alpha
            Scale factor.
        x
            First input vector.
        y
            Second input vector.

        Notes
        -----
        The current vector (``w``) cannot be used for ``x`` or ``y`` but
        ``x`` and ``y`` can be the same.

        See Also
        --------
        Vec.axpy, Vec.aypx, Vec.axpby, Vec.maxpy, petsc.VecWAXPY

        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecWAXPY(self.vec, sval, x.vec, y.vec) )

    def maxpy(self, alphas: Sequence[Scalar], vecs: Sequence[Vec]) -> None:
        """Compute and store y = Σₙ(ɑₙ·Xₙ) + y with X an array of vectors.

        Equivalent to ``y[:] = alphas[i]*vecs[i, :] + y[:]``.

        Logically collective.

        Parameters
        ----------
        alphas
            Array of scale factors, one for each vector in ``vecs``.
        vecs
            Array of vectors.

        Notes
        -----
        The current vector cannot be used as any of ``vecs``, but ``vecs``
        can contain duplicates.

        See Also
        --------
        Vec.axpy, Vec.aypx, Vec.axpby, Vec.waxpy, petsc.VecMAXPY

        """
        cdef PetscInt n = 0
        cdef PetscScalar *a = NULL
        cdef PetscVec *v = NULL
        cdef object tmp1 = iarray_s(alphas, &n, &a)
        cdef object tmp2 = oarray_p(empty_p(n),NULL, <void**>&v)
        assert n == len(vecs)
        cdef Py_ssize_t i=0
        for i from 0 <= i < n:
            v[i] = (<Vec?>(vecs[i])).vec
        CHKERR( VecMAXPY(self.vec, n, a, v) )

    def pointwiseMult(self, Vec x, Vec y) -> None:
        """Compute and store the component-wise multiplication of two vectors.

        Equivalent to ``w[i] = x[i] * y[i]``.

        Logically collective.

        Parameters
        ----------
        x, y
            Input vectors to multiply component-wise.

        Notes
        -----
        Any subset of this vector, ``x`` or ``y`` may be the same vector.

        See Also
        --------
        Vec.pointwiseDivide, petsc.VecPointwiseMult

        """
        CHKERR( VecPointwiseMult(self.vec, x.vec, y.vec) )

    def pointwiseDivide(self, Vec x, Vec y) -> None:
        """Compute and store the component-wise division of two vectors.

        Equivalent to ``w[i] = x[i] / y[i]``.

        Logically collective.

        Parameters
        ----------
        x
            Numerator vector.
        y
            Denominator vector.

        Notes
        -----
        Any subset of this vector, ``x`` or ``y`` may be the same vector.

        See Also
        --------
        Vec.pointwiseMult, petsc.VecPointwiseDivide

        """
        CHKERR( VecPointwiseDivide(self.vec, x.vec, y.vec) )

    def pointwiseMin(self, Vec x, Vec y) -> None:
        """Compute and store the component-wise minimum of two vectors.

        Equivalent to ``w[i] = min(x[i], y[i])``.

        Logically collective.

        Parameters
        ----------
        x, y
            Input vectors to find the component-wise minima.

        Notes
        -----
        Any subset of this vector, ``x`` or ``y`` may be the same vector.

        For complex numbers only the real part is compared.

        See Also
        --------
        Vec.pointwiseMax, Vec.pointwiseMaxAbs, petsc.VecPointwiseMin

        """
        CHKERR( VecPointwiseMin(self.vec, x.vec, y.vec) )

    def pointwiseMax(self, Vec x, Vec y) -> None:
        """Compute and store the component-wise maximum of two vectors.

        Equivalent to ``w[i] = max(x[i], y[i])``.

        Logically collective.

        Parameters
        ----------
        x, y
            Input vectors to find the component-wise maxima.

        Notes
        -----
        Any subset of this vector, ``x`` or ``y`` may be the same vector.

        For complex numbers only the real part is compared.

        See Also
        --------
        Vec.pointwiseMin, Vec.pointwiseMaxAbs, petsc.VecPointwiseMax

        """
        CHKERR( VecPointwiseMax(self.vec, x.vec, y.vec) )

    def pointwiseMaxAbs(self, Vec x, Vec y) -> None:
        """Compute and store the component-wise maximum absolute values.

        Equivalent to ``w[i] = max(abs(x[i]), abs(y[i]))``.

        Logically collective.

        Parameters
        ----------
        x, y
            Input vectors to find the component-wise maxima.

        Notes
        -----
        Any subset of this vector, ``x`` or ``y`` may be the same vector.

        See Also
        --------
        Vec.pointwiseMin, Vec.pointwiseMax, petsc.VecPointwiseMaxAbs

        """
        CHKERR( VecPointwiseMaxAbs(self.vec, x.vec, y.vec) )

    def maxPointwiseDivide(self, Vec vec) -> float:
        """Return the maximum of the component-wise absolute value division.

        Equivalent to ``result = max_i abs(x[i] / y[i])``.

        Logically collective.

        Parameters
        ----------
        x
            Numerator vector.
        y
            Denominator vector.

        Notes
        -----
        ``x`` and ``y`` may be the same vector.

        If a particular ``y[i]`` is zero it will be treated as one.

        See Also
        --------
        Vec.pointwiseMin, Vec.pointwiseMax, Vec.pointwiseMaxAbs
        petsc.VecMaxPointwiseDivide

        """
        cdef PetscReal rval = 0
        CHKERR( VecMaxPointwiseDivide(self.vec, vec.vec, &rval) )
        return toReal(rval)

    def getValue(self, index: int) -> Scalar:
        """Return a single value from the vector.

        Not collective.

        Parameters
        ----------
        index
            Location of the value to read.

        Notes
        -----
        If `Vec.setValues` has been called then `Vec.assemblyBegin` and
        `Vec.assemblyEnd` must precede this method.

        See Also
        --------
        Vec.getValues, petsc.VecGetValues

        """
        cdef PetscInt    ival = asInt(index)
        cdef PetscScalar sval = 0
        CHKERR( VecGetValues(self.vec, 1, &ival, &sval) )
        return toScalar(sval)

    def getValues(
        self,
        indices: Sequence[int],
        values: Sequence[Scalar] | None = None,
    ) -> ArrayScalar:
        """Return values from certain locations in the vector.

        Only values on the same processor may be accessed.

        Not collective.

        Parameters
        ----------
        indices
            Locations of the values to read.
        values
            Location to store the collected values. If not provided then a new
            array will be allocated.

        Returns
        -------
        ArrayScalar
            Collected values. If ``values`` is provided then that is returned.

        Notes
        -----
        If `Vec.setValues` has been called then `Vec.assemblyBegin` and
        `Vec.assemblyEnd` must precede this method.

        See Also
        --------
        Vec.getValue, Vec.setValues, petsc.VecGetValues

        """
        return vecgetvalues(self.vec, indices, values)

    def getValuesStagStencil(self, indices, values=None):
        """Not implemented."""
        raise NotImplementedError('getValuesStagStencil not yet implemented in petsc4py')

    def setValue(
        self,
        index: int,
        value: Scalar,
        addv: InsertMode | int | None = None,
    ) -> None:
        """Insert or add a single value in the vector.

        Not collective.

        Parameters
        ----------
        index
            Location to write to.
        value
            Value to insert at ``index``.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entry with new
              value (default).

            - `InsertMode.ADD_VALUES` Add new value to existing one.

        Notes
        -----
        The values may be cached so `Vec.assemblyBegin` and `Vec.assemblyEnd`
        must be called after all calls of this method are completed.

        Multiple calls to `Vec.setValue` cannot be made with different values
        for ``addv`` without intermediate calls to `Vec.assemblyBegin` and
        `Vec.assemblyEnd`.

        See Also
        --------
        Vec.getValue, Vec.getValues, Vec.setValues, petsc.VecSetValues

        """
        cdef PetscInt    ival = asInt(index)
        cdef PetscScalar sval = asScalar(value)
        cdef PetscInsertMode caddv = insertmode(addv)
        CHKERR( VecSetValues(self.vec, 1, &ival, &sval, caddv) )

    def setValues(
        self,
        indices: Sequence[int],
        values: Sequence[Scalar],
        addv: InsertMode | int | None = None,
    ) -> None:
        """Insert or add multiple values in the vector.

        Not collective.

        Parameters
        ----------
        indices
            Locations to write to. Any negative indices are ignored.
        values
            Values to insert at ``indices``.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

        Notes
        -----
        The values may be cached so `Vec.assemblyBegin` and `Vec.assemblyEnd`
        must be called after all calls of this method are completed.

        Multiple calls to `Vec.setValues` cannot be made with different values
        for ``addv`` without intermediate calls to `Vec.assemblyBegin` and
        `Vec.assemblyEnd`.

        See Also
        --------
        Vec.getValue, Vec.getValues, Vec.setValue, petsc.VecSetValues

        """
        vecsetvalues(self.vec, indices, values, addv, 0, 0)

    def setValuesBlocked(
        self,
        indices: Sequence[int],
        values: Sequence[Scalar],
        addv: InsertMode | int | None = None,
    ) -> None:
        """Insert or add blocks of values in the vector.

        Equivalent to ``x[bs*indices[i]+j] = y[bs*i+j]`` for
        ``0 <= i < len(indices)``, ``0 <= j < bs`` and ``bs`` `Vec.block_size`.

        Not collective.

        Parameters
        ----------
        indices
            Block indices to write to. Any negative indices are ignored.
        values
            Values to insert at ``indices``. Should have length
            ``len(indices) * vec.block_size``.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

        Notes
        -----
        The values may be cached so `Vec.assemblyBegin` and `Vec.assemblyEnd`
        must be called after all calls of this method are completed.

        Multiple calls to `Vec.setValuesBlocked` cannot be made with different
        values for ``addv`` without intermediate calls to `Vec.assemblyBegin`
        and `Vec.assemblyEnd`.

        See Also
        --------
        Vec.getValues, Vec.setValues, petsc.VecSetValues

        """
        vecsetvalues(self.vec, indices, values, addv, 1, 0)

    def setValuesStagStencil(self, indices, values, addv=None):
        """Not implemented."""
        raise NotImplementedError('setValuesStagStencil not yet implemented in petsc4py')

    def setLGMap(self, LGMap lgmap) -> None:
        """Set the local-to-global numbering.

        This allows users to insert vector entries using a local numbering
        with `Vec.setValuesLocal`.

        Logically collective.

        Parameters
        ----------
        lgmap
            Local-to-global mapping.

        Notes
        -----
        Vectors created with `Vec.duplicate` inherit the same mapping.

        See Also
        --------
        Vec.setValues, Vec.setValuesLocal, Vec.getLGMap
        petsc.VecSetLocalToGlobalMapping

        """
        CHKERR( VecSetLocalToGlobalMapping(self.vec, lgmap.lgm) )

    def getLGMap(self) -> LGMap:
        """Return the local-to-global numbering.

        Not collective.

        See Also
        --------
        Vec.setLGMap, petsc.VecGetLocalToGlobalMapping

        """
        cdef LGMap cmap = LGMap()
        CHKERR( VecGetLocalToGlobalMapping(self.vec, &cmap.lgm) )
        PetscINCREF(cmap.obj)
        return cmap

    def setValueLocal(
        self,
        index: int,
        value: Scalar,
        addv: InsertMode | int | None = None,
    ):
        """Insert or add a single value in the vector using a local numbering.

        Not collective.

        Parameters
        ----------
        index
            Location to write to.
        value
            Value to insert at ``index``.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entry with new
              value (default).

            - `InsertMode.ADD_VALUES` Add new value to existing one.

        Notes
        -----
        The values may be cached so `Vec.assemblyBegin` and `Vec.assemblyEnd`
        must be called after all calls of this method are completed.

        Multiple calls to `Vec.setValueLocal` cannot be made with different
        values for ``addv`` without intermediate calls to `Vec.assemblyBegin`
        and `Vec.assemblyEnd`.

        See Also
        --------
        Vec.setValuesLocal, petsc.VecSetValuesLocal

        """
        cdef PetscInt    ival = asInt(index)
        cdef PetscScalar sval = asScalar(value)
        cdef PetscInsertMode caddv = insertmode(addv)
        CHKERR( VecSetValuesLocal(self.vec, 1, &ival, &sval, caddv) )

    def setValuesLocal(
        self,
        indices: Sequence[int],
        values: Sequence[Scalar],
        addv: InsertMode | int | None = None,
    ) -> None:
        """Insert or add multiple values in the vector with a local numbering.

        Not collective.

        Parameters
        ----------
        indices
            Locations to write to. Any negative indices are ignored.
        values
            Values to insert at ``indices``.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

        Notes
        -----
        The values may be cached so `Vec.assemblyBegin` and `Vec.assemblyEnd`
        must be called after all calls of this method are completed.

        Multiple calls to `Vec.setValuesLocal` cannot be made with different
        values for ``addv`` without intermediate calls to `Vec.assemblyBegin`
        and `Vec.assemblyEnd`.

        See Also
        --------
        Vec.setValueLocal, petsc.VecSetValuesLocal

        """
        vecsetvalues(self.vec, indices, values, addv, 0, 1)

    def setValuesBlockedLocal(
        self,
        indices: Sequence[int],
        values: Sequence[Scalar],
        addv: InsertMode | int | None = None,
    ) -> None:
        """Insert or add blocks of values in the vector with a local numbering.

        Equivalent to ``x[bs*indices[i]+j] = y[bs*i+j]`` for
        ``0 <= i < len(indices)``, ``0 <= j < bs`` and ``bs`` `Vec.block_size`.

        Not collective.

        Parameters
        ----------
        indices
            Block indices to write to. Any negative indices are ignored.
        values
            Values to insert at ``indices``. Should have length
            ``len(indices) * vec.block_size``.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

        Notes
        -----
        The values may be cached so `Vec.assemblyBegin` and `Vec.assemblyEnd`
        must be called after all calls of this method are completed.

        Multiple calls to `Vec.setValuesBlockedLocal` cannot be made with
        different values for ``addv`` without intermediate calls to
        `Vec.assemblyBegin` and `Vec.assemblyEnd`.

        See Also
        --------
        Vec.setValuesBlocked, Vec.setValuesLocal
        petsc.VecSetValuesBlockedLocal

        """
        vecsetvalues(self.vec, indices, values, addv, 1, 1)

    def assemblyBegin(self) -> None:
        """Begin assembling the vector.

        This routine should be called after completing all calls to
        `Vec.setValues`.

        Collective.

        See Also
        --------
        Vec.assemblyEnd, Vec.setValues, petsc.VecAssemblyBegin

        """
        CHKERR( VecAssemblyBegin(self.vec) )

    # FIXME I don't understand what is happening re viewing the vector
    # post assembly
    def assemblyEnd(self) -> None:
        """Finish assembling the vector.

        This routine should be called after `Vec.assemblyBegin`.

        Collective.

        Notes
        -----
        The vector can be viewed after assembly. See `petsc_options` and
        `petsc.VecAssemblyEnd` for information on the right keys for the
        options database.

        See Also
        --------
        Vec.assemblyBegin, petsc.VecAssemblyEnd

        """
        CHKERR( VecAssemblyEnd(self.vec) )

    def assemble(self) -> None:
        """Assemble the vector in one step.

        To interleave communication and computation `Vec.assemblyBegin` and
        `Vec.assemblyEnd` should be used instead.

        Collective.

        See Also
        --------
        Vec.assemblyBegin, Vec.assemblyEnd

        """
        CHKERR( VecAssemblyBegin(self.vec) )
        CHKERR( VecAssemblyEnd(self.vec) )

    # --- methods for strided vectors ---

    def strideScale(self, field: int, alpha: Scalar) -> None:
        """Scale a component of the vector.

        Logically collective.

        Parameters
        ----------
        field
            Component index. Must be between ``0`` and `Vec.block_size`.
        alpha
            Factor to multiple the component entries by.

        See Also
        --------
        Vec.strideSum, Vec.strideMin, Vec.strideMax, petsc.VecStrideScale

        """
        cdef PetscInt    ival = asInt(field)
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( VecStrideScale(self.vec, ival, sval) )

    def strideSum(self, field: int) -> Scalar:
        """Sum subvector entries.

        Equivalent to ``sum(x[field], x[field+bs], x[field+2*bs], ...)`` where
        ``bs`` is `Vec.block_size`.

        Collective.

        Parameters
        ----------
        field
            Component index. Must be between ``0`` and `Vec.block_size`.

        See Also
        --------
        Vec.strideScale, Vec.strideMin, Vec.strideMax, petsc.VecStrideSum

        """
        cdef PetscInt    ival = asInt(field)
        cdef PetscScalar sval = 0
        CHKERR( VecStrideSum(self.vec, ival, &sval) )
        return toScalar(sval)

    def strideMin(self, field: int) -> tuple[int, float]:
        """Return the minimum of entries in a subvector.

        Equivalent to ``min(x[field], x[field+bs], x[field+2*bs], ...)`` where
        ``bs`` is `Vec.block_size`.

        Collective.

        Parameters
        ----------
        field
            Component index. Must be between ``0`` and `Vec.block_size`.

        Returns
        -------
        int
            Location of minimum.
        float
            Minimum value.

        See Also
        --------
        Vec.strideScale, Vec.strideSum, Vec.strideMax, petsc.VecStrideMin

        """
        cdef PetscInt  ival1 = asInt(field)
        cdef PetscInt  ival2 = 0
        cdef PetscReal rval  = 0
        CHKERR( VecStrideMin(self.vec, ival1, &ival2, &rval) )
        return (toInt(ival2), toReal(rval))

    def strideMax(self, field: int) -> tuple[int, float]:
        """Return the maximum of entries in a subvector.

        Equivalent to ``max(x[field], x[field+bs], x[field+2*bs], ...)`` where
        ``bs`` is `Vec.block_size`.

        Collective.

        Parameters
        ----------
        field
            Component index. Must be between ``0`` and `Vec.block_size`.

        Returns
        -------
        int
            Location of maximum.
        float
            Maximum value.

        See Also
        --------
        Vec.strideScale, Vec.strideSum, Vec.strideMin, petsc.VecStrideMax

        """
        cdef PetscInt  ival1 = asInt(field)
        cdef PetscInt  ival2 = 0
        cdef PetscReal rval  = 0
        CHKERR( VecStrideMax(self.vec, ival1, &ival2, &rval) )
        return (toInt(ival2), toReal(rval))

    def strideNorm(
        self,
        field: int,
        norm_type: NormType | int | None = None,
    ) -> float | tuple[float, float]:
        """Return the norm of entries in a subvector.

        Equivalent to ``norm(x[field], x[field+bs], x[field+2*bs], ...)`` where
        ``bs`` is `Vec.block_size`.

        Collective.

        Parameters
        ----------
        field
            Component index. Must be between ``0`` and `Vec.block_size`.
        norm_type
            The norm type. See `Vec.norm` for more information.

        Returns
        -------
        typing.Any
            The computed norm. A 2-tuple is returned if `NormType.NORM_1_AND_2`
            is specified.

        See Also
        --------
        Vec.strideScale, Vec.strideSum, petsc.VecStrideNorm

        """
        cdef PetscInt ival = asInt(field)
        cdef PetscNormType norm_1_2 = PETSC_NORM_1_AND_2
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal rval[2]
        CHKERR( VecStrideNorm(self.vec, ival, ntype, rval) )
        if ntype != norm_1_2: return toReal(rval[0])
        else: return (toReal(rval[0]), toReal(rval[1]))

    def strideScatter(
        self,
        field: int,
        Vec vec,
        addv: InsertMode | int | None = None,
    ) -> None:
        """Scatter entries into a component of another vector.

        The current vector is expected to be single-component
        (`Vec.block_size` of ``1``) and the target vector is expected to be
        multi-component.

        Collective.

        Parameters
        ----------
        field
            Component index of ``vec``. Must be between ``0`` and
            `Vec.block_size`.
        vec
            Multi-component vector to be scattered into.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

            - `InsertMode.MAX_VALUES` Keep the maximum value in for each
              entry.

        Notes
        -----
        The parallel layouts of the vectors must match.

        See Also
        --------
        Vec.strideGather, petsc.VecStrideScatter

        """
        cdef PetscInt ival = asInt(field)
        cdef PetscInsertMode caddv = insertmode(addv)
        CHKERR( VecStrideScatter(self.vec, ival, vec.vec, caddv) )

    def strideGather(
        self,
        field: int,
        Vec vec,
        addv: InsertMode | int | None = None,
    ) -> None:
        """Insert component values into a single-component vector.

        The current vector is expected to be multi-component (`Vec.block_size`
        greater than ``1``) and the target vector is expected to be
        single-component.

        Collective.

        Parameters
        ----------
        field
            Component index of the current vector. Must be between ``0`` and
            `Vec.block_size`.
        vec
            Single-component vector to be inserted into.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

            - `InsertMode.MAX_VALUES` Keep the maximum value in for each
              entry.

        Notes
        -----
        The parallel layouts of the vectors must match.

        See Also
        --------
        Vec.strideScatter, petsc.VecStrideScatter

        """
        cdef PetscInt ival = asInt(field)
        cdef PetscInsertMode caddv = insertmode(addv)
        CHKERR( VecStrideGather(self.vec, ival, vec.vec, caddv) )

    # --- methods for vectors with ghost values ---

    def localForm(self) -> Any:
        """Return context manager for viewing ghost vectors in local form.

        Logically collective.

        Returns
        -------
        typing.Any
            Context manager yielding the vector in local (ghosted) form.

        Notes
        -----
        This operation does not perform a copy. To obtain up-to-date ghost
        values `Vec.ghostUpdateBegin` and `Vec.ghostUpdateEnd` must be called
        first.

        Non-ghost values can be found
        at ``values[0:nlocal]`` and ghost values at
        ``values[nlocal:nlocal+nghost]``.

        Examples
        --------
        >>> with vec.localForm() as lf:
        ...     # compute with lf

        See Also
        --------
        Vec.createGhost, Vec.ghostUpdateBegin, Vec.ghostUpdateEnd
        petsc.VecGhostGetLocalForm, petsc.VecGhostRestoreLocalForm

        """
        return _Vec_LocalForm(self)

    def ghostUpdateBegin(
        self,
        addv: InsertMode | int | None = None,
        mode: ScatterMode | int | None = None,
    ) -> None:
        """Start updating ghosted vector entries.

        See `Vec.ghostUpdate` for more information.

        Neighbour-wise collective.

        See Also
        --------
        Vec.ghostUpdateEnd, Vec.ghostUpdate, Vec.createGhost
        petsc.VecGhostUpdateBegin

        """
        cdef PetscInsertMode  caddv = insertmode(addv)
        cdef PetscScatterMode csctm = scattermode(mode)
        CHKERR( VecGhostUpdateBegin(self.vec, caddv, csctm) )

    def ghostUpdateEnd(
        self,
        addv: InsertMode | int | None = None,
        mode: ScatterMode | int | None = None,
    ) -> None:
        """Finish updating ghosted vector entries.

        See `Vec.ghostUpdate` for more information.

        Neighbour-wise collective.

        See Also
        --------
        Vec.ghostUpdateBegin, Vec.ghostUpdate, Vec.createGhost
        petsc.VecGhostUpdateEnd

        """
        cdef PetscInsertMode  caddv = insertmode(addv)
        cdef PetscScatterMode csctm = scattermode(mode)
        CHKERR( VecGhostUpdateEnd(self.vec, caddv, csctm) )

    # FIXME addv should also support InsertMode.MIN_VALUES but this isn't defined
    def ghostUpdate(
        self,
        addv: InsertMode | int | None = None,
        mode: ScatterMode | int | None = None,
    ) -> None:
        """Update ghosted vector entries.

        Neighbour-wise collective.

        Parameters
        ----------
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing ones.

            - `InsertMode.MAX_VALUES` Keep the maximum value in for each
              entry.
        mode
            Scatter mode. Possible values are:

            - `ScatterMode.FORWARD` Update ghost regions with correct values
              from the owning process.

            - `ScatterMode.REVERSE` Accumulate ghost region values onto the
              owning process.

        Notes
        -----
        This operation is blocking. To interleave computation and
        communication use `Vec.ghostUpdateBegin` and `Vec.ghostUpdateEnd`
        instead.

        Examples
        --------
        To accumulate ghost region values onto owning processes and then
        update ghost regions correctly one should do the following:

        >>> vec.ghostUpdate(InsertMode.ADD_VALUES, ScatterMode.REVERSE)
        >>> vec.ghostUpdate(InsertMode.INSERT_VALUES, ScatterMode.FORWARD)

        See Also
        --------
        Vec.ghostUpdateBegin, Vec.ghostUpdateEnd

        """
        cdef PetscInsertMode  caddv = insertmode(addv)
        cdef PetscScatterMode csctm = scattermode(mode)
        CHKERR( VecGhostUpdateBegin(self.vec, caddv, csctm) )
        CHKERR( VecGhostUpdateEnd(self.vec, caddv, csctm) )

    def setMPIGhost(self, ghosts: Sequence[int]) -> None:
        """Set the ghost points for a ghosted vector.

        This method is an alternative to calling `Vec.createGhost`.

        Collective.

        Parameters
        ----------
        ghosts
            Global indices of ghost points. These do not need to be sorted.

        See Also
        --------
        Vec.createGhost

        """
        cdef PetscInt ng=0, *ig=NULL
        ghosts = iarray_i(ghosts, &ng, &ig)
        CHKERR( VecMPISetGhost(self.vec, ng, ig) )

    #

    def getSubVector(self, IS iset, Vec subvec=None) -> Vec:
        """Return a subvector from given indices.

        Once finished with the subvector it should be returned with
        `Vec.restoreSubVector`.

        Collective.

        Parameters
        ----------
        iset
            Index set describing which indices to extract into the subvector.
        subvec
            Subvector to copy entries into. If `None` then a new `Vec` will
            be created.

        Returns
        -------
        Vec
            Subvector containing the extracted entries. If ``subvec`` is
            provided this is returned.

        Notes
        -----
        This function may return a subvector without making a copy, therefore
        it is not safe to use the original vector while modifying the
        subvector. Other non-overlapping subvectors can still be obtained
        using this function.

        The resulting subvector inherits the block size from ``iset`` if
        greater than one. Otherwise, the block size is guessed from the block
        size of the original vector.

        See Also
        --------
        Vec.restoreSubVector, petsc.VecGetSubVector

        """
        if subvec is None: subvec = Vec()
        else: CHKERR( VecDestroy(&subvec.vec) )
        CHKERR( VecGetSubVector(self.vec, iset.iset, &subvec.vec) )
        return subvec

    def restoreSubVector(self, IS iset, Vec subvec) -> None:
        """Restore a subvector extracted using `Vec.getSubVector`.

        Collective.

        Parameters
        ----------
        iset
            Index set describing the indices represented by the subvector.
        subvec
            Subvector to restore.

        See Also
        --------
        Vec.getSubVector, petsc.VecRestoreSubVector

        """
        CHKERR( VecRestoreSubVector(self.vec, iset.iset, &subvec.vec) )

    def getNestSubVecs(self) -> list[Vec]:
        """Return all the vectors contained in the nested vector.

        Not collective.

        See Also
        --------
        Vec.setNestSubVecs, petsc.VecNestGetSubVecs

        """
        cdef PetscInt N=0
        cdef PetscVec* sx=NULL
        CHKERR( VecNestGetSubVecs(self.vec, &N, &sx) )
        output = []
        for i in range(N):
          pyvec = Vec()
          pyvec.vec = sx[i]
          CHKERR( PetscObjectReference(<PetscObject> pyvec.vec) )
          output.append(pyvec)

        return output

    def setNestSubVecs(
        self,
        sx: Sequence[Vec],
        idxm: Sequence[int] | None = None,
    ) -> None:
        """Set the component vectors at specified indices in the nested vector.

        Not collective.

        Parameters
        ----------
        sx
            Array of component vectors.
        idxm
            Indices of the component vectors, defaults to ``range(len(sx))``.

        See Also
        --------
        Vec.getNestSubVecs, petsc.VecNestSetSubVecs

        """
        if idxm is None: idxm = range(len(sx))
        else: assert len(idxm) == len(sx)
        cdef PetscInt N = 0
        cdef PetscInt* cidxm = NULL
        idxm = iarray_i(idxm, &N, &cidxm)


        cdef PetscVec* csx = NULL
        tmp = oarray_p(empty_p(N), NULL, <void**>&csx)
        for i from 0 <= i < N: csx[i] = (<Vec?>sx[i]).vec

        CHKERR( VecNestSetSubVecs(self.vec, N, cidxm, csx) )

    #

    def setDM(self, DM dm) -> None:
        """Set the DM describing the data layout of the vector.

        Not collective.

        Notes
        -----
        This method is rarely needed as `DM.getLocalVec` or
        `DM.getGlobalVec` will set this appropriately.

        See Also
        --------
        Vec.getDM, petsc.VecSetDM

        """
        CHKERR( VecSetDM(self.vec, dm.dm) )

    def getDM(self) -> DM:
        """Return the DM describing the data layout of the vector.

        Not collective.

        See Also
        --------
        Vec.setDM, petsc.VecGetDM

        """
        cdef DM dm = DM()
        CHKERR( VecGetDM(self.vec, &dm.dm) )
        return dm

    #

    property sizes:
        """The local and global vector sizes.

        See Also
        --------
        Vec.getSizes, Vec.setSizes

        """
        def __get__(self) -> tuple[int, int]:
            return self.getSizes()
        def __set__(self, value):
            self.setSizes(value)

    property size:
        """The global vector size.

        See Also
        --------
        Vec.getSize

        """
        def __get__(self) -> int:
            return self.getSize()

    property local_size:
        """The local vector size.

        See Also
        --------
        Vec.getLocalSize

        """
        def __get__(self) -> int:
            return self.getLocalSize()

    property block_size:
        """The block size.

        See Also
        --------
        Vec.getBlockSize

        """
        def __get__(self) -> int:
            return self.getBlockSize()

    property owner_range:
        """The locally owned range of indices in the form ``[low, high)``.

        See Also
        --------
        Vec.getOwnershipRange

        """
        def __get__(self) -> tuple[int, int]:
            return self.getOwnershipRange()

    property owner_ranges:
        """The range of indices owned by each process.

        See Also
        --------
        Vec.getOwnershipRanges

        """
        def __get__(self) -> ArrayInt:
            return self.getOwnershipRanges()

    property buffer_w:
        """Writeable buffered view of the local portion of the vector.

        See Also
        --------
        Vec.getBuffer

        """
        def __get__(self) -> Any:
            return self.getBuffer()

    property buffer_r:
        """Read-only buffered view of the local portion of the vector.

        See Also
        --------
        Vec.getBuffer

        """
        def __get__(self) -> Any:
            return self.getBuffer(True)

    property array_w:
        """Writeable numpy array containing the local portion of the vector.

        See Also
        --------
        Vec.getArray

        """
        def __get__(self) -> ArrayScalar:
            return self.getArray()
        def __set__(self, value):
            cdef buf = self.getBuffer()
            with buf as array: array[:] = value

    property array_r:
        """Read-only numpy array containing the local portion of the vector.

        See Also
        --------
        Vec.getArray

        """
        def __get__(self) -> ArrayScalar:
            return self.getArray(True)

    property buffer:
        """Alias for `Vec.buffer_w`."""
        def __get__(self) -> Any:
            return self.buffer_w

    property array:
        """Alias for `Vec.array_w`."""
        def __get__(self) -> ArrayScalar:
            return self.array_w
        def __set__(self, value):
            self.array_w = value

    # --- NumPy array interface (legacy) ---

    property __array_interface__:
        def __get__(self):
            cdef buf = self.getBuffer()
            return buf.__array_interface__

# --------------------------------------------------------------------

del VecType
del VecOption

# --------------------------------------------------------------------
