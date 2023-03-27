# --------------------------------------------------------------------

class ISType(object):
    """TODO"""
    GENERAL = S_(ISGENERAL)
    BLOCK   = S_(ISBLOCK)
    STRIDE  = S_(ISSTRIDE)

# --------------------------------------------------------------------

cdef class IS(Object):
    """PETSc object used for efficient indexing.

    Attributes
    ----------
    permutation
        TODO
    identity
        TODO
    sorted
        TODO
    sizes
        TODO
    size
        TODO
    local_size
        TODO
    block_size
        TODO
    indices
        TODO
    array
        TODO
    """

    Type = ISType

    #

    def __cinit__(self):
        self.obj = <PetscObject*> &self.iset
        self.iset = NULL

    # buffer interface (PEP 3118)

    def __getbuffer__(self, Py_buffer *view, int flags):
        cdef _IS_buffer buf = _IS_buffer(self)
        buf.acquirebuffer(view, flags)

    def __releasebuffer__(self, Py_buffer *view):
        cdef _IS_buffer buf = <_IS_buffer>(view.obj)
        buf.releasebuffer(view)
        <void>self # unused


    # 'with' statement (PEP 343)

    def __enter__(self):
        cdef _IS_buffer buf = _IS_buffer(self)
        self.set_attr('__buffer__', buf)
        return buf.enter()

    def __exit__(self, *exc):
        cdef _IS_buffer buf = self.get_attr('__buffer__')
        self.set_attr('__buffer__', None)
        return buf.exit()
    #

    def view(self, Viewer viewer=None):
        """Display the IS.

        Parameters
        ----------
        viewer : Viewer, optional
            Viewer used to display the IS.

        See Also
        --------
        TODO
        """
        cdef PetscViewer cviewer = NULL
        if viewer is not None: cviewer = viewer.vwr
        CHKERR( ISView(self.iset, cviewer) )

    def destroy(self):
        """Destroy the IS.

        Returns
        -------
        self

        See Also
        --------
        TODO https://petsc.org/release/docs/manualpages/IS/ISDestroy/
        """
        CHKERR( ISDestroy(&self.iset) )
        return self

    def create(self, comm=None):
        """Create an IS.

        Parameters
        ----------
        comm : PetscComm?, optional

        Returns
        -------
        self

        See Also
        --------
        manualpages/IS/ISCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscIS newiset = NULL
        CHKERR( ISCreate(ccomm, &newiset) )
        PetscCLEAR(self.obj); self.iset = newiset
        return self

    def setType(self, is_type):
        """Build an IS for a particular :class:`ISType`.

        Parameters
        ----------
        is_type : str
            The name of the index set type.

        See Also
        --------
        TODO https://petsc.org/release/docs/manualpages/IS/ISSetType/
        """
        cdef PetscISType cval = NULL
        is_type = str2bytes(is_type, &cval)
        CHKERR( ISSetType(self.iset, cval) )

    def getType(self):
        """Return the index set type name associated with the IS.

        Returns
        -------
        str
            Index set type name.

        See Also
        --------
        petsc:ISGetType
        """
        cdef PetscISType cval = NULL
        CHKERR( ISGetType(self.iset, &cval) )
        return bytes2str(cval)

    def createGeneral(self, indices, comm=None):
        """Create an IS with indices.

        Parameters
        ----------
        indices
            Integer array.
        comm : optional
            MPI communicator. 

        Returns
        -------
        IS
            A new IS.

        See Also
        --------
        petsc:ISCreateGeneral
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt nidx = 0, *idx = NULL
        cdef PetscCopyMode cm = PETSC_COPY_VALUES
        cdef PetscIS newiset = NULL
        indices = iarray_i(indices, &nidx, &idx)
        CHKERR( ISCreateGeneral(ccomm, nidx, idx, cm, &newiset) )
        PetscCLEAR(self.obj); self.iset = newiset
        return self

    def createBlock(self, bsize, indices, comm=None):
        """Create an IS where each integer represents a fixed block of indices.

        Parameters
        ----------
        bsize : int
            The block size.
        indices : array_like
            Integer array of indices.
        comm : optional
            MPI communicator.

        Returns
        -------
        IS
            The new, blocked, IS.

        See Also
        --------
        petsc:ISCreateBlock
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs = asInt(bsize)
        cdef PetscInt nidx = 0, *idx = NULL
        cdef PetscCopyMode cm = PETSC_COPY_VALUES
        cdef PetscIS newiset = NULL
        indices = iarray_i(indices, &nidx, &idx)
        CHKERR( ISCreateBlock(ccomm, bs, nidx, idx, cm, &newiset) )
        PetscCLEAR(self.obj); self.iset = newiset
        return self

    def createStride(self, size, first=0, step=0, comm=None):
        """Create an IS consisting of evenly spaced integers.

        Parameters
        ----------
        size : int
            The length of the locally owned portion of the index set.
        first : int, default 0
            The first element of the index set.
        step : int, default 0
            The difference between adjacent indices.
        comm : optional
            The MPI communicator.

        Returns
        -------
        TODO or self?
        IS
            The new IS.

        See Also
        --------
        petsc:ISCreateStride
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt csize  = asInt(size)
        cdef PetscInt cfirst = asInt(first)
        cdef PetscInt cstep  = asInt(step)
        cdef PetscIS newiset = NULL
        CHKERR( ISCreateStride(ccomm, csize, cfirst, cstep, &newiset) )
        PetscCLEAR(self.obj); self.iset = newiset
        return self

    def duplicate(self):
        """Create a copy of the index set.

        Returns
        -------
        IS
            The new IS.

        See Also
        --------
        petsc:ISDuplicate
        """
        cdef IS iset = type(self)()
        CHKERR( ISDuplicate(self.iset, &iset.iset) )
        return iset

    def copy(self, IS result=None):
        """Copy the contents of the IS into another.

        Parameters
        ----------
        result : IS, optional
            The target IS. If ``None`` then `IS.duplicate` is called first.

        Returns
        -------
        IS
            The copied IS.

        See Also
        --------
        petsc:ISCopy
        """
        if result is None:
            result = type(self)()
        if result.iset == NULL:
            CHKERR( ISDuplicate(self.iset, &result.iset) )
        CHKERR( ISCopy(self.iset, result.iset) )
        return result

    def load(self, Viewer viewer):
        """Load a stored index set.

        Parameters
        ----------
        viewer : Viewer
            Binary (``BINARY`` or ``HDF5``) file viewer.

        Returns
        -------
        # TODO or self?
        IS
            The loaded index set.

        See Also
        --------
        petsc:ISLoad
        """
        cdef MPI_Comm comm = MPI_COMM_NULL
        cdef PetscObject obj = <PetscObject>(viewer.vwr)
        if self.iset == NULL:
            CHKERR( PetscObjectGetComm(obj, &comm) )
            CHKERR( ISCreate(comm, &self.iset) )
        CHKERR( ISLoad(self.iset, viewer.vwr) )
        return self

    def allGather(self):
        """Concatenate index sets stored across processors.

        Returns
        -------
        IS
            The concatenated index set (same on different processes).

        See Also
        --------
        petsc:ISAllGather
        """
        cdef IS iset = IS()
        CHKERR( ISAllGather(self.iset, &iset.iset) )
        return iset

    def toGeneral(self):
        """Convert the index set type to `ISType.GENERAL`.

        Returns
        -------
        IS
            self

        See Also
        --------
        petsc:ISToGeneral
        """
        CHKERR( ISToGeneral(self.iset) )
        return self

    def buildTwoSided(self, IS toindx=None):
        """
        TODO
        """
        cdef PetscIS ctoindx = NULL
        if toindx is not None: ctoindx = toindx.iset
        cdef IS result = IS()
        CHKERR( ISBuildTwoSided(self.iset, ctoindx, &result.iset) )
        return result

    def invertPermutation(self, nlocal=None):
        """Creates a new permutation that is the inverse of a given permutation.

        Parameters
        ----------
        nlocal : int, optional
            The number of indices on this processor in the result
            (default ``PETSC_DECIDE``).

        Returns
        -------
        IS
            The inverted permutation.

        See Also
        --------
        petsc:ISInvertPermutation
        """
        cdef PetscInt cnlocal = PETSC_DECIDE
        if nlocal is not None: cnlocal = asInt(nlocal)
        cdef IS iset = IS()
        CHKERR( ISInvertPermutation(self.iset, cnlocal, &iset.iset) )
        return iset

    def getSize(self):
        """Return the global length of an index set.

        Returns
        -------
        int
            The global size.

        See Also
        --------
        petsc:ISGetSize
        """
        cdef PetscInt N = 0
        CHKERR( ISGetSize(self.iset, &N) )
        return toInt(N)

    def getLocalSize(self):
        """Return the process-local length of an index set.

        Returns
        -------
        int
            The local size.

        See Also
        --------
        petsc:ISGetLocalSize
        """
        cdef PetscInt n = 0
        CHKERR( ISGetLocalSize(self.iset, &n) )
        return toInt(n)

    def getSizes(self):
        """Return the local and global sizes of an index set.

        Returns
        -------
        local_size : int
            The local size.
        global_size : int
            The global size.

        See Also
        --------
        IS.getLocalSize
        IS.getGlobalSize
        """
        cdef PetscInt n = 0, N = 0
        CHKERR( ISGetLocalSize(self.iset, &n) )
        CHKERR( ISGetSize(self.iset, &N) )
        return (toInt(n), toInt(N))

    def getBlockSize(self):
        """Return the number of elements in a block.

        Returns
        -------
        int
            The number of elements in a block.

        See Also
        --------
        petsc:ISGetBlockSize
        """
        cdef PetscInt bs = 1
        CHKERR( ISGetBlockSize(self.iset, &bs) )
        return toInt(bs)

    def setBlockSize(self, bs):
        cdef PetscInt cbs = asInt(bs)
        CHKERR( ISSetBlockSize(self.iset, cbs) )

    def sort(self):
        CHKERR( ISSort(self.iset) )
        return self

    def isSorted(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( ISSorted(self.iset, &flag) )
        return toBool(flag)

    def setPermutation(self):
        CHKERR( ISSetPermutation(self.iset) )
        return self

    def isPermutation(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( ISPermutation(self.iset, &flag) )
        return toBool(flag)

    def setIdentity(self):
        CHKERR( ISSetIdentity(self.iset) )
        return self

    def isIdentity(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( ISIdentity(self.iset, &flag) )
        return toBool(flag)

    def equal(self, IS iset):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( ISEqual(self.iset, iset.iset, &flag) )
        return toBool(flag)

    def sum(self, IS iset):
        cdef IS out = IS()
        CHKERR( ISSum(self.iset, iset.iset, &out.iset) )
        return out

    def expand(self, IS iset):
        cdef IS out = IS()
        CHKERR( ISExpand(self.iset, iset.iset, &out.iset) )
        return out

    def union(self, IS iset): # XXX review this
        cdef PetscBool flag1=PETSC_FALSE, flag2=PETSC_FALSE
        CHKERR( ISSorted(self.iset, &flag1) )
        CHKERR( ISSorted(iset.iset, &flag2) )
        cdef IS out = IS()
        if flag1==PETSC_TRUE and flag2==PETSC_TRUE:
            CHKERR( ISSum(self.iset, iset.iset, &out.iset) )
        else:
            CHKERR( ISExpand(self.iset, iset.iset, &out.iset) )
        return out

    def difference(self, IS iset):
        cdef IS out = IS()
        CHKERR( ISDifference(self.iset, iset.iset, &out.iset) )
        return out

    def complement(self, nmin, nmax):
        cdef PetscInt cnmin = asInt(nmin)
        cdef PetscInt cnmax = asInt(nmax)
        cdef IS out = IS()
        CHKERR( ISComplement(self.iset, cnmin, cnmax, &out.iset) )
        return out

    def embed(self, IS iset, drop):
        cdef PetscBool bval = drop
        cdef IS out = IS()
        CHKERR( ISEmbed(self.iset, iset.iset, bval, &out.iset) )
        return out

    def renumber(self, IS mult=None):
        cdef PetscIS mlt = NULL
        if mult is not None: mlt = mult.iset
        cdef IS out = IS()
        cdef PetscInt n = 0
        CHKERR( ISRenumber(self.iset, mlt, &n, &out.iset) )
        return (toInt(n), out)
    #

    def setIndices(self, indices):
        cdef PetscInt nidx = 0, *idx = NULL
        cdef PetscCopyMode cm = PETSC_COPY_VALUES
        indices = iarray_i(indices, &nidx, &idx)
        CHKERR( ISGeneralSetIndices(self.iset, nidx, idx, cm) )

    def getIndices(self):
        cdef PetscInt size = 0
        cdef const PetscInt *indices = NULL
        CHKERR( ISGetLocalSize(self.iset, &size) )
        CHKERR( ISGetIndices(self.iset, &indices) )
        cdef object oindices = None
        try:
            oindices = array_i(size, indices)
        finally:
            CHKERR( ISRestoreIndices(self.iset, &indices) )
        return oindices

    def setBlockIndices(self, bsize, indices):
        cdef PetscInt bs = asInt(bsize)
        cdef PetscInt nidx = 0, *idx = NULL
        cdef PetscCopyMode cm = PETSC_COPY_VALUES
        indices = iarray_i(indices, &nidx, &idx)
        CHKERR( ISBlockSetIndices(self.iset, bs, nidx, idx, cm) )

    def getBlockIndices(self):
        cdef PetscInt size = 0, bs = 1
        cdef const PetscInt *indices = NULL
        CHKERR( ISGetLocalSize(self.iset, &size) )
        CHKERR( ISGetBlockSize(self.iset, &bs) )
        CHKERR( ISBlockGetIndices(self.iset, &indices) )
        cdef object oindices = None
        try:
            oindices = array_i(size//bs, indices)
        finally:
            CHKERR( ISBlockRestoreIndices(self.iset, &indices) )
        return oindices

    def setStride(self, size, first=0, step=1):
        cdef PetscInt csize = asInt(size)
        cdef PetscInt cfirst = asInt(first)
        cdef PetscInt cstep = asInt(step)
        CHKERR( ISStrideSetStride(self.iset, csize, cfirst, cstep) )

    def getStride(self):
        cdef PetscInt size=0, first=0, step=0
        CHKERR( ISGetLocalSize(self.iset, &size) )
        CHKERR( ISStrideGetInfo(self.iset, &first, &step) )
        return (toInt(size), toInt(first), toInt(step))

    def getInfo(self):
        cdef PetscInt first = 0, step = 0
        CHKERR( ISStrideGetInfo(self.iset, &first, &step) )
        return (toInt(first), toInt(step))

    #

    property permutation:
        def __get__(self):
            return self.isPermutation()

    property identity:
        def __get__(self):
            return self.isIdentity()

    property sorted:
        def __get__(self):
            return self.isSorted()

    #

    property sizes:
        def __get__(self):
            return self.getSizes()

    property size:
        def __get__(self):
            return self.getSize()

    property local_size:
        def __get__(self):
            return self.getLocalSize()

    property block_size:
        def __get__(self):
            return self.getBlockSize()

    property indices:
        def __get__(self):
            return self.getIndices()

    property array:
        def __get__(self):
            return asarray(self)

    # --- NumPy array interface (legacy) ---

    property __array_interface__:
        def __get__(self):
            cdef _IS_buffer buf = _IS_buffer(self)
            return buf.__array_interface__

# --------------------------------------------------------------------


class GLMapMode(object):
    MASK = PETSC_IS_GTOLM_MASK
    DROP = PETSC_IS_GTOLM_DROP


class LGMapType(object):
    BASIC = S_(ISLOCALTOGLOBALMAPPINGBASIC)
    HASH  = S_(ISLOCALTOGLOBALMAPPINGHASH)


# --------------------------------------------------------------------

cdef class LGMap(Object):

    MapMode = GLMapMode

    Type = LGMapType
    #

    def __cinit__(self):
        self.obj = <PetscObject*> &self.lgm
        self.lgm = NULL

    def __call__(self, indices, result=None):
        self.apply(indices, result)

    #

    def setType(self, lgmap_type):
        cdef PetscISLocalToGlobalMappingType cval = NULL
        lgmap_type = str2bytes(lgmap_type, &cval)
        CHKERR( ISLocalToGlobalMappingSetType(self.lgm, cval) )

    def setFromOptions(self):
        CHKERR( ISLocalToGlobalMappingSetFromOptions(self.lgm) )

    def view(self, Viewer viewer=None):
        cdef PetscViewer cviewer = NULL
        if viewer is not None: cviewer = viewer.vwr
        CHKERR( ISLocalToGlobalMappingView(self.lgm, cviewer) )

    def destroy(self):
        CHKERR( ISLocalToGlobalMappingDestroy(&self.lgm) )
        return self

    def create(self, indices, bsize=None, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscInt bs = 1, nidx = 0, *idx = NULL
        cdef PetscCopyMode cm = PETSC_COPY_VALUES
        cdef PetscLGMap newlgm = NULL
        if bsize is not None: bs = asInt(bsize)
        if bs == PETSC_DECIDE: bs = 1
        indices = iarray_i(indices, &nidx, &idx)
        CHKERR( ISLocalToGlobalMappingCreate(
                ccomm, bs, nidx, idx, cm, &newlgm) )
        PetscCLEAR(self.obj); self.lgm = newlgm
        return self

    def createIS(self, IS iset):
        cdef PetscLGMap newlgm = NULL
        CHKERR( ISLocalToGlobalMappingCreateIS(
            iset.iset, &newlgm) )
        PetscCLEAR(self.obj); self.lgm = newlgm
        return self

    def createSF(self, SF sf, start):
        cdef PetscLGMap newlgm = NULL
        cdef PetscInt cstart = asInt(start)
        CHKERR( ISLocalToGlobalMappingCreateSF(sf.sf, cstart, &newlgm) )
        PetscCLEAR(self.obj); self.lgm = newlgm
        return self

    def getSize(self):
        cdef PetscInt n = 0
        CHKERR( ISLocalToGlobalMappingGetSize(self.lgm, &n) )
        return toInt(n)

    def getBlockSize(self):
        cdef PetscInt bs = 1
        CHKERR( ISLocalToGlobalMappingGetBlockSize(self.lgm, &bs) )
        return toInt(bs)

    def getIndices(self):
        cdef PetscInt size = 0
        cdef const PetscInt *indices = NULL
        CHKERR( ISLocalToGlobalMappingGetSize(
                self.lgm, &size) )
        CHKERR( ISLocalToGlobalMappingGetIndices(
                self.lgm, &indices) )
        cdef object oindices = None
        try:
            oindices = array_i(size, indices)
        finally:
            CHKERR( ISLocalToGlobalMappingRestoreIndices(
                    self.lgm, &indices) )
        return oindices

    def getBlockIndices(self):
        cdef PetscInt size = 0, bs = 1
        cdef const PetscInt *indices = NULL
        CHKERR( ISLocalToGlobalMappingGetSize(
                self.lgm, &size) )
        CHKERR( ISLocalToGlobalMappingGetBlockSize(
                self.lgm, &bs) )
        CHKERR( ISLocalToGlobalMappingGetBlockIndices(
                self.lgm, &indices) )
        cdef object oindices = None
        try:
            oindices = array_i(size//bs, indices)
        finally:
            CHKERR( ISLocalToGlobalMappingRestoreBlockIndices(
                    self.lgm, &indices) )
        return oindices

    def getInfo(self):
        cdef PetscInt i, nproc = 0, *procs = NULL,
        cdef PetscInt *numprocs = NULL, **indices = NULL
        cdef object neighs = { }
        CHKERR( ISLocalToGlobalMappingGetInfo(
                self.lgm, &nproc, &procs, &numprocs, &indices) )
        try:
            for i from 0 <= i < nproc:
                neighs[toInt(procs[i])] = array_i(numprocs[i], indices[i])
        finally:
            ISLocalToGlobalMappingRestoreInfo(
                self.lgm, &nproc, &procs, &numprocs, &indices)
        return neighs

    def getBlockInfo(self):
        cdef PetscInt i, nproc = 0, *procs = NULL,
        cdef PetscInt *numprocs = NULL, **indices = NULL
        cdef object neighs = { }
        CHKERR( ISLocalToGlobalMappingGetBlockInfo(
                self.lgm, &nproc, &procs, &numprocs, &indices) )
        try:
            for i from 0 <= i < nproc:
                neighs[toInt(procs[i])] = array_i(numprocs[i], indices[i])
        finally:
            ISLocalToGlobalMappingRestoreBlockInfo(
                self.lgm, &nproc, &procs, &numprocs, &indices)
        return neighs

    #

    def apply(self, indices, result=None):
        cdef PetscInt niidx = 0, *iidx = NULL
        cdef PetscInt noidx = 0, *oidx = NULL
        indices = iarray_i(indices, &niidx, &iidx)
        if result is None: result = empty_i(niidx)
        result  = oarray_i(result,  &noidx, &oidx)
        assert niidx == noidx, "incompatible array sizes"
        CHKERR( ISLocalToGlobalMappingApply(
            self.lgm, niidx, iidx, oidx) )
        return result

    def applyBlock(self, indices, result=None):
        cdef PetscInt niidx = 0, *iidx = NULL
        cdef PetscInt noidx = 0, *oidx = NULL
        indices = iarray_i(indices, &niidx, &iidx)
        if result is None: result = empty_i(niidx)
        result  = oarray_i(result,  &noidx, &oidx)
        assert niidx == noidx, "incompatible array sizes"
        CHKERR( ISLocalToGlobalMappingApplyBlock(
            self.lgm, niidx, iidx, oidx) )
        return result

    def applyIS(self, IS iset):
        cdef IS result = IS()
        CHKERR( ISLocalToGlobalMappingApplyIS(
            self.lgm, iset.iset, &result.iset) )
        return result

    def applyInverse(self, indices, mode=None):
        cdef PetscGLMapMode cmode = PETSC_IS_GTOLM_MASK
        if mode is not None: cmode = mode
        cdef PetscInt n = 0, *idx = NULL
        indices = iarray_i(indices, &n, &idx)
        cdef PetscInt nout = n, *idxout = NULL
        if cmode != PETSC_IS_GTOLM_MASK:
            CHKERR( ISGlobalToLocalMappingApply(
                    self.lgm, cmode, n, idx, &nout, NULL) )
        result = oarray_i(empty_i(nout), &nout, &idxout)
        CHKERR( ISGlobalToLocalMappingApply(
                self.lgm, cmode, n, idx, &nout, idxout) )
        return result

    def applyBlockInverse(self, indices, mode=None):
        cdef PetscGLMapMode cmode = PETSC_IS_GTOLM_MASK
        if mode is not None: cmode = mode
        cdef PetscInt n = 0, *idx = NULL
        indices = iarray_i(indices, &n, &idx)
        cdef PetscInt nout = n, *idxout = NULL
        if cmode != PETSC_IS_GTOLM_MASK:
            CHKERR( ISGlobalToLocalMappingApply(
                    self.lgm, cmode, n, idx, &nout, NULL) )
        result = oarray_i(empty_i(nout), &nout, &idxout)
        CHKERR( ISGlobalToLocalMappingApplyBlock(
                self.lgm, cmode, n, idx, &nout, idxout) )
        return result
    #

    property size:
        def __get__(self):
            return self.getSize()

    property block_size:
        def __get__(self):
            return self.getBlockSize()

    property indices:
        def __get__(self):
            return self.getIndices()

    property block_indices:
        def __get__(self):
            return self.getBlockIndices()

    property info:
        def __get__(self):
            return self.getInfo()

    property block_info:
        def __get__(self):
            return self.getBlockInfo()

# --------------------------------------------------------------------

del ISType
del GLMapMode
del LGMapType
# --------------------------------------------------------------------
