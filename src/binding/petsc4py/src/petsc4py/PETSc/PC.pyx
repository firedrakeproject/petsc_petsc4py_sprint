# --------------------------------------------------------------------

class PCType(object):
    """The preconditioner method."""
    NONE               = S_(PCNONE)
    JACOBI             = S_(PCJACOBI)
    SOR                = S_(PCSOR)
    LU                 = S_(PCLU)
    QR                 = S_(PCQR)
    SHELL              = S_(PCSHELL)
    BJACOBI            = S_(PCBJACOBI)
    VPBJACOBI          = S_(PCVPBJACOBI)
    MG                 = S_(PCMG)
    EISENSTAT          = S_(PCEISENSTAT)
    ILU                = S_(PCILU)
    ICC                = S_(PCICC)
    ASM                = S_(PCASM)
    GASM               = S_(PCGASM)
    KSP                = S_(PCKSP)
    COMPOSITE          = S_(PCCOMPOSITE)
    REDUNDANT          = S_(PCREDUNDANT)
    SPAI               = S_(PCSPAI)
    NN                 = S_(PCNN)
    CHOLESKY           = S_(PCCHOLESKY)
    PBJACOBI           = S_(PCPBJACOBI)
    MAT                = S_(PCMAT)
    HYPRE              = S_(PCHYPRE)
    PARMS              = S_(PCPARMS)
    FIELDSPLIT         = S_(PCFIELDSPLIT)
    TFS                = S_(PCTFS)
    ML                 = S_(PCML)
    GALERKIN           = S_(PCGALERKIN)
    EXOTIC             = S_(PCEXOTIC)
    CP                 = S_(PCCP)
    BFBT               = S_(PCBFBT)
    LSC                = S_(PCLSC)
    PYTHON             = S_(PCPYTHON)
    PFMG               = S_(PCPFMG)
    SYSPFMG            = S_(PCSYSPFMG)
    REDISTRIBUTE       = S_(PCREDISTRIBUTE)
    SVD                = S_(PCSVD)
    GAMG               = S_(PCGAMG)
    CHOWILUVIENNACL    = S_(PCCHOWILUVIENNACL)
    ROWSCALINGVIENNACL = S_(PCROWSCALINGVIENNACL)
    SAVIENNACL         = S_(PCSAVIENNACL)
    BDDC               = S_(PCBDDC)
    KACZMARZ           = S_(PCKACZMARZ)
    TELESCOPE          = S_(PCTELESCOPE)
    PATCH              = S_(PCPATCH)
    LMVM               = S_(PCLMVM)
    HMG                = S_(PCHMG)
    DEFLATION          = S_(PCDEFLATION)
    HPDDM              = S_(PCHPDDM)
    H2OPUS             = S_(PCH2OPUS)

class PCSide(object):
    """The manner in which the preconditioner is applied."""
    # native
    LEFT      = PC_LEFT
    RIGHT     = PC_RIGHT
    SYMMETRIC = PC_SYMMETRIC
    # aliases
    L = LEFT
    R = RIGHT
    S = SYMMETRIC

class PCASMType(object):
    """The *ASM* subtype."""
    NONE        = PC_ASM_NONE
    BASIC       = PC_ASM_BASIC
    RESTRICT    = PC_ASM_RESTRICT
    INTERPOLATE = PC_ASM_INTERPOLATE

class PCGASMType(object):
    """The *GASM* subtype."""
    NONE        = PC_GASM_NONE
    BASIC       = PC_GASM_BASIC
    RESTRICT    = PC_GASM_RESTRICT
    INTERPOLATE = PC_GASM_INTERPOLATE

class PCMGType(object):
    """The *MG* subtype."""
    MULTIPLICATIVE = PC_MG_MULTIPLICATIVE
    ADDITIVE       = PC_MG_ADDITIVE
    FULL           = PC_MG_FULL
    KASKADE        = PC_MG_KASKADE

class PCMGCycleType(object):
    """The *MG* cycle type."""
    V = PC_MG_CYCLE_V
    W = PC_MG_CYCLE_W

class PCGAMGType(object):
    """The *GAMG* subtype."""
    AGG       = S_(PCGAMGAGG)
    GEO       = S_(PCGAMGGEO)
    CLASSICAL = S_(PCGAMGCLASSICAL)

class PCCompositeType(object):
    """The composite type."""
    ADDITIVE                 = PC_COMPOSITE_ADDITIVE
    MULTIPLICATIVE           = PC_COMPOSITE_MULTIPLICATIVE
    SYMMETRIC_MULTIPLICATIVE = PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE
    SPECIAL                  = PC_COMPOSITE_SPECIAL
    SCHUR                    = PC_COMPOSITE_SCHUR

class PCFieldSplitSchurPreType(object):
    """The field split Schur subtype."""
    SELF                     = PC_FIELDSPLIT_SCHUR_PRE_SELF
    SELFP                    = PC_FIELDSPLIT_SCHUR_PRE_SELFP
    A11                      = PC_FIELDSPLIT_SCHUR_PRE_A11
    USER                     = PC_FIELDSPLIT_SCHUR_PRE_USER
    FULL                     = PC_FIELDSPLIT_SCHUR_PRE_FULL

class PCFieldSplitSchurFactType(object):
    """The field split Schur factorization type."""
    DIAG                     = PC_FIELDSPLIT_SCHUR_FACT_DIAG
    LOWER                    = PC_FIELDSPLIT_SCHUR_FACT_LOWER
    UPPER                    = PC_FIELDSPLIT_SCHUR_FACT_UPPER
    FULL                     = PC_FIELDSPLIT_SCHUR_FACT_FULL

class PCPatchConstructType(object):
    """The patch construction type."""
    STAR                     = PC_PATCH_STAR
    VANKA                    = PC_PATCH_VANKA
    PARDECOMP                = PC_PATCH_PARDECOMP
    USER                     = PC_PATCH_USER
    PYTHON                   = PC_PATCH_PYTHON

class PCHPDDMCoarseCorrectionType(object):
    """The *HPDDM* coarse correction type."""
    DEFLATED                 = PC_HPDDM_COARSE_CORRECTION_DEFLATED
    ADDITIVE                 = PC_HPDDM_COARSE_CORRECTION_ADDITIVE
    BALANCED                 = PC_HPDDM_COARSE_CORRECTION_BALANCED

class PCDeflationSpaceType(object):
    """The deflation space subtype."""
    HAAR                     = PC_DEFLATION_SPACE_HAAR
    DB2                      = PC_DEFLATION_SPACE_DB2
    DB4                      = PC_DEFLATION_SPACE_DB4
    DB8                      = PC_DEFLATION_SPACE_DB8
    DB16                     = PC_DEFLATION_SPACE_DB16
    BIORTH22                 = PC_DEFLATION_SPACE_BIORTH22
    MEYER                    = PC_DEFLATION_SPACE_MEYER
    AGGERGATION              = PC_DEFLATION_SPACE_AGGREGATION
    USER                     = PC_DEFLATION_SPACE_USER

class PCFailedReason(object):
    """The reason the preconditioner has failed."""
    SETUP_ERROR              = PC_SETUP_ERROR
    NOERROR                  = PC_NOERROR
    FACTOR_STRUCT_ZEROPIVOT  = PC_FACTOR_STRUCT_ZEROPIVOT
    FACTOR_NUMERIC_ZEROPIVOT = PC_FACTOR_NUMERIC_ZEROPIVOT
    FACTOR_OUTMEMORY         = PC_FACTOR_OUTMEMORY
    FACTOR_OTHER             = PC_FACTOR_OTHER
    SUBPC_ERROR              = PC_SUBPC_ERROR

# --------------------------------------------------------------------

cdef class PC(Object):
    """Preconditioners.
    
    `PC` is described in the `PETSc manual <petsc:sec_ksppc>`.
    Calling the `PC` with a vector as an argument will `apply` the
    preconditioner as shown in the example below.

    Examples
    --------
    >>> from petsc4py import PETSc
    >>> v = PETSc.Vec().createWithArray([1,2])
    >>> m = PETSc.Mat().createDense(2,array=[[1,0],[0,1]])
    >>> pc = PETSc.PC().create()
    >>> pc.setOperators(m)
    >>> u = pc(v) # Vec u is created internally, can also be passed as second argument

    See Also
    --------
    petsc.PC
    """

    Type = PCType
    Side = PCSide

    ASMType                   = PCASMType
    GASMType                  = PCGASMType
    MGType                    = PCMGType
    MGCycleType               = PCMGCycleType
    GAMGType                  = PCGAMGType
    CompositeType             = PCCompositeType
    SchurFactType             = PCFieldSplitSchurFactType
    SchurPreType              = PCFieldSplitSchurPreType
    PatchConstructType        = PCPatchConstructType
    HPDDMCoarseCorrectionType = PCHPDDMCoarseCorrectionType
    DeflationSpaceType        = PCDeflationSpaceType
    FailedReason              = PCFailedReason

    # --- xxx ---

    def __cinit__(self):
        self.obj = <PetscObject*> &self.pc
        self.pc = NULL

    def __call__(self, x, y=None):
        if y is None: # XXX do this better
            y = self.getOperators()[0].createVecLeft()
        self.apply(x, y)
        return y

    # --- xxx ---

    def view(self, Viewer viewer=None) -> None:
        """View the `PC` object.

        Collective.

        Parameters
        ----------
        viewer
            The visualization context.

        See Also
        --------
        petsc.PCView
        """        
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PCView(self.pc, vwr) )

    def destroy(self) -> Self:
        """Destroy the `PC` that was created with `create`.

        Collective.

        See Also
        --------
        petsc.PCDestroy
        """
        CHKERR( PCDestroy(&self.pc) )
        self.pc = NULL
        return self

    def create(self, comm=None) -> Self:
        """Create an empty `PC`. 
        
        Collective. The default preconditioner for sparse matrices is `ILU` or
        `ICC` with 0 fill on one process and block Jacobi (`BJACOBI`) with `ILU`
        or `ICC` in parallel. For dense matrices it is always `None`.

        Parameters
        ----------
        comm
            The communicator or None for `Sys.getDefaultComm`.

        See Also
        --------
        destroy, petsc.PCCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscPC newpc = NULL
        CHKERR( PCCreate(ccomm, &newpc) )
        PetscCLEAR(self.obj); self.pc = newpc
        return self

    def setType(self, pc_type: Type | str) -> None:
        """Set the preconditioner type.

        Collective.

        Parameters
        ----------
        pc_type
            The preconditioner type.

        See Also
        --------
        petsc_options, getType, petsc.TSSetType
        """
        cdef PetscPCType cval = NULL
        pc_type = str2bytes(pc_type, &cval)
        CHKERR( PCSetType(self.pc, cval) )

    def getType(self) -> str:
        """Return the preconditioner type.
        
        Not collective.

        See Also
        --------
        setType, petsc.PCGetType
        """
        cdef PetscPCType cval = NULL
        CHKERR( PCGetType(self.pc, &cval) )
        return bytes2str(cval)

    def setOptionsPrefix(self, prefix: str) -> None:
        """Set the prefix used for all the `PC` options.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix to prepend to all option names.

        See Also
        --------
        petsc_options, petsc.PCSetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( PCSetOptionsPrefix(self.pc, cval) )

    def getOptionsPrefix(self) -> str:
        """Return the prefix used for all the `PC` options.

        Not collective.

        See Also
        --------
        petsc_options, petsc.PCGetOptionsPrefix
        """
        cdef const char *cval = NULL
        CHKERR( PCGetOptionsPrefix(self.pc, &cval) )
        return bytes2str(cval)

    def appendOptionsPrefix(self, prefix: str) -> None:
        """Append to the prefix used for all the `PC` options.

        Logically Collective.

        Parameters
        ----------
        prefix
            The prefix to append to the current prefix.

        See Also
        --------
        petsc_options, petsc.PCAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( PCAppendOptionsPrefix(self.pc, cval) )

    def setFromOptions(self) -> None:
        """Set various `PC` parameters from user options.

        Collective.

        See Also
        --------
        petsc_options, petsc.PCSetFromOptions
        """
        CHKERR( PCSetFromOptions(self.pc) )

    def setOperators(self, Mat A=None, Mat P=None) -> None:
        """Sets the matrices associated with the linear system.

        Logically collective. Passing `None` for ``A`` or ``P`` removes the
        matrix that is currently used. PETSc does not reset the matrix entries
        of either ``A`` or ``P`` to zero after a linear solve; the user is
        completely responsible for matrix assembly. See `Mat.zeroEntries` to
        zero all elements of a matrix.

        Parameters
        ----------
        A
            the matrix which defines the linear system
        P
            the matrix to be used in constructing the preconditioner, usually
            the same as ``A``
        
        See Also
        --------
        petsc.PCSetOperators
        """
        cdef PetscMat amat=NULL
        if A is not None: amat = A.mat
        cdef PetscMat pmat=amat
        if P is not None: pmat = P.mat
        CHKERR( PCSetOperators(self.pc, amat, pmat) )

    def getOperators(self) -> tuple[Mat,Mat]:
        """Return the matrices associated with a linear system.

        Not collective.

        See Also
        --------
        setOperators, petsc.PCGetOperators
        """
        cdef Mat A = Mat(), P = Mat()
        CHKERR( PCGetOperators(self.pc, &A.mat, &P.mat) )
        PetscINCREF(A.obj)
        PetscINCREF(P.obj)
        return (A, P)

    def setUseAmat(self, flag: bool) -> None:
        """Set to indicate to apply `PC` to ``A`` and not ``P``. 

        Logically collective. Sets a flag to indicate that when the
        preconditioner needs to apply (part of) the operator during the
        preconditioning process, it applies to ``A`` provided to
        `TS.setRHSJacobian`, `TS.setIJacobian`, `SNES.setJacobian`,
        `KSP.setOperators` or `PC.setOperators` not the ``P``.

        Parameter
        ---------
        flag
            Set True to use ``A`` and False to use ``P``.

        See Also
        --------
        setOperators, petsc.PCSetUseAmat
        """
        cdef PetscBool cflag = PETSC_FALSE
        if flag:
            cflag = PETSC_TRUE
        CHKERR( PCSetUseAmat(self.pc, cflag) )

    def getUseAmat(self) -> bool:
        """Return the flag to indicate if `PC` is applied to ``A`` or ``P``.

        Logically collective.

        Returns
        -------
        flag : bool
            True if ``A`` is used and False if ``P``.

        See Also
        --------
        setUseAmat, petsc.PCGetUseAmat
        """
        cdef PetscBool cflag = PETSC_FALSE
        CHKERR( PCGetUseAmat(self.pc, &cflag) )
        return toBool(cflag)

    def setReusePreconditioner(self, flag: bool) -> None:
        """Set to indicate the preconditioner is to be reused.

        Logically collective. Normally if the ``A`` matrix inside a `PC`
        changes, the `PC` automatically updates itself using information from
        the changed matrix. Enable this option prevents this.

        Parameters
        ----------
        flag
            Set to `True` to use the reuse the current preconditioner and
            `False` to recompute on changes to the matrix.

        See Also
        --------
        setOperators, petsc.PCSetReusePreconditioner
        """
        cdef PetscBool cflag = PETSC_FALSE
        if flag:
            cflag = PETSC_TRUE
        CHKERR( PCSetReusePreconditioner(self.pc, cflag) )

    def setFailedReason(self, reason: FailedReason | str) ->  None:
        """Set the reason the `PC` terminated.

        Logically collective.

        Parameters
        ----------
        reason
            the reason the `PC` terminated

        See Also
        --------
        petsc.PCSetFailedReason
        """
        cdef PetscPCFailedReason val = reason
        CHKERR( PCSetFailedReason(self.pc, val) )

    def getFailedReason(self) -> FailedReason:
        """Return the reason the `PC` terminated.

        Logically collective. This is the maximum over reason over all ranks in
        the `PC` communicator.

        See Also
        --------
        petsc.PCGetFailedReason
        """
        cdef PetscPCFailedReason reason = PC_NOERROR
        CHKERR( PCGetFailedReason(self.pc, &reason) )
        return reason

    def getFailedReasonRank(self) -> FailedReason:
        """Return the reason the `PC` terminated on this rank.

        Not collective. Different ranks may have different reasons.

        See Also
        --------
        getFailedReason, petsc.PCGetFailedReasonRank
        """
        cdef PetscPCFailedReason reason = PC_NOERROR
        CHKERR( PCGetFailedReasonRank(self.pc, &reason) )
        return reason

    def setUp(self) -> None:
        """Set up the internal data structures for the `PC`.

        Collective.

        See Also
        --------
        petsc.PCSetUp
        """
        CHKERR( PCSetUp(self.pc) )

    def reset(self) -> None:
        """Reset the `PC`, removing any allocated vectors and matrices.

        Collective.

        See Also
        --------
        petsc.PCReset
        """        
        CHKERR( PCReset(self.pc) )

    def setUpOnBlocks(self) -> None:
        """Set up the `PC` for each block.

        Collective. For nested preconditioners such as `BJACOBI`, `setUp` is not
        called on each sub-`KSP` when `setUp` is called on the outer `PC`. This
        routine ensures it is called.

        See Also
        --------
        setUp, petsc.PCSetUpOnBlocks
        """
        CHKERR( PCSetUpOnBlocks(self.pc) )

    def apply(self, Vec x, Vec y) -> None:
        """Apply the `PC` to a vector.

        Collective.

        Parameters
        ----------
        x
            The input vector.
        y
            The output vector, cannot be the same as ``x``.

        See also
        --------
        petsc.PCApply
        """
        CHKERR( PCApply(self.pc, x.vec, y.vec) )

    def matApply(self, Mat x, Mat y) -> None:
        """Apply the `PC` to many vectors stored as `Mat.Type.DENSE`.

        Collective.

        Parameters
        ----------
        x
            The input matrix.
        y
            The output matrix, cannot be the same as ``x``.

        See also
        --------
        petsc.PCMatApply, petsc.PCApply      
        """
        CHKERR( PCMatApply(self.pc, x.mat, y.mat) )

    def applyTranspose(self, Vec x, Vec y) -> None:
        """Apply the transpose of the `PC` to a vector.

        Collective. For complex numbers this applies the non-Hermitian
        transpose.

        Parameters
        ----------
        x
            The input vector.
        y
            The output vector, cannot be the same as ``x``.

        See also
        --------
        petsc.PCApply
        """
        CHKERR( PCApplyTranspose(self.pc, x.vec, y.vec) )

    def applySymmetricLeft(self, Vec x, Vec y) -> None:
        """Apply the left part of a symmetric `PC` to a vector.

        Collective.

        Parameters
        ----------
        x
            The input vector.
        y
            The output vector, cannot be the same as ``x``.

        See also
        --------
        petsc.PCApplySymmetricLeft
        """
        CHKERR( PCApplySymmetricLeft(self.pc, x.vec, y.vec) )

    def applySymmetricRight(self, Vec x, Vec y) -> None:
        """Apply the right part of a symmetric `PC` to a vector.

        Collective.

        Parameters
        ----------
        x
            The input vector.
        y
            The output vector, cannot be the same as ``x``.

        See also
        --------
        petsc.PCApplySymmetricRight
        """
        CHKERR( PCApplySymmetricRight(self.pc, x.vec, y.vec) )

    # --- discretization space ---

    def getDM(self) -> DM:
        """Return the `DM` associated with the `PC`.

        Not collective. 

        See Also
        --------
        petsc.PCGetDM
        """
        cdef PetscDM newdm = NULL
        CHKERR( PCGetDM(self.pc, &newdm) )
        cdef DM dm = subtype_DM(newdm)()
        dm.dm = newdm
        PetscINCREF(dm.obj)
        return dm

    def setDM(self, DM dm) -> None:
        """Set the `DM` that may be used by some preconditioners.

        Logically collective.

        Parameters
        ----------
        dm
            The `DM` object.

        See Also
        --------
        petsc.PCSetDM
        """        
        CHKERR( PCSetDM(self.pc, dm.dm) )

    def setCoordinates(self, coordinates: Sequence[Sequence[float]]) -> None:
        """Set the coordinates for the nodes on the local process.

        Collective.

        Parameters
        ----------
        coordinates
            The two dimensional coordinate array.

        See Also
        --------
        petsc.PCSetCoordinates
        """
        cdef ndarray xyz = iarray(coordinates, NPY_PETSC_REAL)
        if PyArray_ISFORTRAN(xyz): xyz = PyArray_Copy(xyz)
        if PyArray_NDIM(xyz) != 2: raise ValueError(
            ("coordinates must have two dimensions: "
             "coordinates.ndim=%d") % (PyArray_NDIM(xyz)) )
        cdef PetscInt nvtx = <PetscInt> PyArray_DIM(xyz, 0)
        cdef PetscInt ndim = <PetscInt> PyArray_DIM(xyz, 1)
        cdef PetscReal *coords = <PetscReal*> PyArray_DATA(xyz)
        CHKERR( PCSetCoordinates(self.pc, ndim, nvtx, coords) )

    # --- Python ---

    def createPython(self, context: Any = None, comm: Comm | None = None) -> Self:
        """Create a preconditioner of Python type.

        Collective.

        Parameters
        ----------
        context
            An instance of the Python class implementing the required methods.
        comm
            The communicator associated with the object. Defaults to
            `Sys.getDefaultComm`.

        See Also
        --------
        petsc_python_pc, setType, setPythonContext, PC.Type.PYTHON
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscPC newpc = NULL
        CHKERR( PCCreate(ccomm, &newpc) )
        PetscCLEAR(self.obj); self.pc = newpc
        CHKERR( PCSetType(self.pc, PCPYTHON) )
        CHKERR( PCPythonSetContext(self.pc, <void*>context) )
        return self

    def setPythonContext(self, context: Any) -> None:
        """Set the instance of the Python class implementing the required Python methods.

        Not collective.

        See Also
        --------
        petsc_python_pc, getPythonContext
        """
        CHKERR( PCPythonSetContext(self.pc, <void*>context) )

    def getPythonContext(self) -> Any:
        """Return the instance of the Python class implementing the required Python methods.

        Not collective.

        See Also
        --------
        petsc_python_pc, setPythonContext
        """
        cdef void *context = NULL
        CHKERR( PCPythonGetContext(self.pc, &context) )
        if context == NULL: return None
        else: return <object> context

    def setPythonType(self, py_type: str) -> None:
        """Set the fully qualified Python name of the class to be used.

        Collective.

        See Also
        --------
        petsc_python_pc, setPythonContext, getPythonType, petsc.PCPythonSetType

        """
        cdef const char *cval = NULL
        py_type = str2bytes(py_type, &cval)
        CHKERR( PCPythonSetType(self.pc, cval) )

    def getPythonType(self) -> str:
        """Return the fully qualified Python name of the class used by the `PC`.

        Not collective.

        See Also
        --------
        petsc_python_pc, setPythonContext, setPythonType, petsc.PCPythonGetType
        """
        cdef const char *cval = NULL
        CHKERR( PCPythonGetType(self.pc, &cval) )
        return bytes2str(cval)

    # --- ASM ---

    def setASMType(self, asmtype):
        cdef PetscPCASMType cval = asmtype
        CHKERR( PCASMSetType(self.pc, cval) )

    def setASMOverlap(self, overlap):
        cdef PetscInt ival = asInt(overlap)
        CHKERR( PCASMSetOverlap(self.pc, ival) )

    def setASMLocalSubdomains(self, nsd, is_=None, is_local=None):
        cdef PetscInt n = asInt(nsd)
        cdef PetscInt i = 0
        cdef PetscIS *isets = NULL
        cdef PetscIS *isets_local = NULL
        if is_ is not None:
            assert len(is_) == nsd
            CHKERR( PetscMalloc(<size_t>n*sizeof(PetscIS), &isets) )
            for i in range(n):
                isets[i] = (<IS?>is_[i]).iset
        if is_local is not None:
            assert len(is_local) == nsd
            CHKERR( PetscMalloc(<size_t>n*sizeof(PetscIS), &isets_local) )
            for i in range(n):
                isets_local[i] = (<IS?>is_local[i]).iset
        CHKERR( PCASMSetLocalSubdomains(self.pc, n, isets, isets_local) )
        CHKERR( PetscFree(isets) )
        CHKERR( PetscFree(isets_local) )

    def setASMTotalSubdomains(self, nsd, is_=None, is_local=None):
        cdef PetscInt n = asInt(nsd)
        cdef PetscInt i = 0
        cdef PetscIS *isets = NULL
        cdef PetscIS *isets_local = NULL
        if is_ is not None:
            assert len(is_) == nsd
            CHKERR( PetscMalloc(<size_t>n*sizeof(PetscIS), &isets) )
            for i in range(n):
                isets[i] = (<IS?>is_[i]).iset
        if is_local is not None:
            assert len(is_local) == nsd
            CHKERR( PetscMalloc(<size_t>n*sizeof(PetscIS), &isets_local) )
            for i in range(n):
                isets_local[i] = (<IS?>is_local[i]).iset
        CHKERR( PCASMSetTotalSubdomains(self.pc, n, isets, isets_local) )
        CHKERR( PetscFree(isets) )
        CHKERR( PetscFree(isets_local) )

    def getASMSubKSP(self):
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        CHKERR( PCASMGetSubKSP(self.pc, &n, NULL, &p) )
        return [ref_KSP(p[i]) for i from 0 <= i <n]

    def setASMSortIndices(self, dosort):
        cdef PetscBool cdosort = asBool(dosort)
        CHKERR( PCASMSetSortIndices(self.pc, cdosort) )

    # --- GASM ---

    def setGASMType(self, gasmtype):
        cdef PetscPCGASMType cval = gasmtype
        CHKERR( PCGASMSetType(self.pc, cval) )

    def setGASMOverlap(self, overlap):
        cdef PetscInt ival = asInt(overlap)
        CHKERR( PCGASMSetOverlap(self.pc, ival) )

    # --- GAMG ---

    def setGAMGType(self, gamgtype):
        cdef PetscPCGAMGType cval = NULL
        gamgtype = str2bytes(gamgtype, &cval)
        CHKERR( PCGAMGSetType(self.pc, cval) )

    def setGAMGLevels(self, levels):
        cdef PetscInt ival = asInt(levels)
        CHKERR( PCGAMGSetNlevels(self.pc, ival) )

    def setGAMGSmooths(self, smooths):
        cdef PetscInt ival = asInt(smooths)
        CHKERR( PCGAMGSetNSmooths(self.pc, ival) )

    # --- Hypre ---

    def getHYPREType(self):
        cdef PetscPCHYPREType cval = NULL
        CHKERR( PCHYPREGetType(self.pc, &cval) )
        return bytes2str(cval)

    def setHYPREType(self, hypretype):
        cdef PetscPCHYPREType cval = NULL
        hypretype = str2bytes(hypretype, &cval)
        CHKERR( PCHYPRESetType(self.pc, cval) )

    def setHYPREDiscreteCurl(self, Mat mat):
        CHKERR( PCHYPRESetDiscreteCurl(self.pc, mat.mat) )

    def setHYPREDiscreteGradient(self, Mat mat):
        CHKERR( PCHYPRESetDiscreteGradient(self.pc, mat.mat) )

    def setHYPRESetAlphaPoissonMatrix(self, Mat mat):
        CHKERR( PCHYPRESetAlphaPoissonMatrix(self.pc, mat.mat) )

    def setHYPRESetBetaPoissonMatrix(self, Mat mat=None):
        cdef PetscMat pmat = NULL
        if mat is not None: pmat = mat.mat
        CHKERR( PCHYPRESetBetaPoissonMatrix(self.pc, pmat) )
    
    def setHYPRESetInterpolations(self, dim, Mat RT_Pi_Full=None, RT_Pi=None,
                                  Mat ND_Pi_Full=None, ND_Pi=None):
        cdef PetscMat RT_full_mat = NULL
        if RT_Pi_Full is not None: RT_full_mat = RT_Pi_Full.mat
        cdef PetscMat ND_full_mat = NULL
        if ND_Pi_Full is not None: ND_full_mat = ND_Pi_Full.mat
        cdef PetscInt idim = asInt(dim)
        cdef PetscMat *RT_Pi_mat = NULL
        if RT_Pi is not None:
            PetscMalloc(<size_t>dim*sizeof(PetscMat), &RT_Pi_mat)
            assert len(RT_Pi) == dim
            for i in range(dim):
                RT_Pi_mat[i] = (<Mat?>RT_Pi[i]).mat
        cdef PetscMat *ND_Pi_mat = NULL
        if ND_Pi is not None:
            PetscMalloc(<size_t>dim*sizeof(PetscMat), &ND_Pi_mat)
            assert len(ND_Pi) == dim
            for i in range(dim):
                ND_Pi_mat[dim] = (<Mat?>ND_Pi[i]).mat
        CHKERR (PCHYPRESetInterpolations(self.pc, idim, RT_full_mat, RT_Pi_mat,
                                         ND_full_mat, ND_Pi_mat))
        CHKERR (PetscFree(RT_Pi_mat))
        CHKERR (PetscFree(ND_Pi_mat))
       
       

    def setHYPRESetEdgeConstantVectors(self, Vec ozz, Vec zoz, Vec zzo=None):
        cdef PetscVec zzo_vec = NULL
        if zzo is not None: zzo_vec = zzo.vec
        CHKERR( PCHYPRESetEdgeConstantVectors(self.pc, ozz.vec, zoz.vec,
                                              zzo_vec) )

    def setHYPREAMSSetInteriorNodes(self, Vec interior):
        CHKERR(PCHYPREAMSSetInteriorNodes(self.pc, interior.vec))

    # --- Factor ---

    def setFactorSolverType(self, solver):
        cdef PetscMatSolverType cval = NULL
        solver = str2bytes(solver, &cval)
        CHKERR( PCFactorSetMatSolverType(self.pc, cval) )

    def getFactorSolverType(self):
        cdef PetscMatSolverType cval = NULL
        CHKERR( PCFactorGetMatSolverType(self.pc, &cval) )
        return bytes2str(cval)

    def setFactorSetUpSolverType(self):
        CHKERR( PCFactorSetUpMatSolverType(self.pc) )

    def setFactorOrdering(self, ord_type=None, nzdiag=None, reuse=None):
        cdef PetscMatOrderingType cval = NULL
        if ord_type is not None:
            ord_type = str2bytes(ord_type, &cval)
            CHKERR( PCFactorSetMatOrderingType(self.pc, cval) )
        cdef PetscReal rval = 0
        if nzdiag is not None:
            rval = asReal(nzdiag)
            CHKERR( PCFactorReorderForNonzeroDiagonal(self.pc, rval) )
        cdef PetscBool bval = PETSC_FALSE
        if reuse is not None:
            bval = PETSC_TRUE if reuse else PETSC_FALSE
            CHKERR( PCFactorSetReuseOrdering(self.pc, bval) )

    def setFactorPivot(self, zeropivot=None, inblocks=None):
        cdef PetscReal rval = 0
        if zeropivot is not None:
            rval = asReal(zeropivot)
            CHKERR( PCFactorSetZeroPivot(self.pc, rval) )
        cdef PetscBool bval = PETSC_FALSE
        if inblocks is not None:
            bval = PETSC_TRUE if inblocks else PETSC_FALSE
            CHKERR( PCFactorSetPivotInBlocks(self.pc, bval) )

    def setFactorShift(self, shift_type=None, amount=None):
        cdef PetscMatFactorShiftType cval = MAT_SHIFT_NONE
        if shift_type is not None:
            cval = matfactorshifttype(shift_type)
            CHKERR( PCFactorSetShiftType(self.pc, cval) )
        cdef PetscReal rval = 0
        if amount is not None:
            rval = asReal(amount)
            CHKERR( PCFactorSetShiftAmount(self.pc, rval) )

    def setFactorLevels(self, levels):
        cdef PetscInt ival = asInt(levels)
        CHKERR( PCFactorSetLevels(self.pc, ival) )

    def getFactorMatrix(self):
        cdef Mat mat = Mat()
        CHKERR( PCFactorGetMatrix(self.pc, &mat.mat) )
        PetscINCREF(mat.obj)
        return mat

    # --- FieldSplit ---

    def setFieldSplitType(self, ctype):
        cdef PetscPCCompositeType cval = ctype
        CHKERR( PCFieldSplitSetType(self.pc, cval) )

    def setFieldSplitIS(self, *fields):
        cdef object name = None
        cdef IS field = None
        cdef const char *cname = NULL
        for name, field in fields:
            name = str2bytes(name, &cname)
            CHKERR( PCFieldSplitSetIS(self.pc, cname, field.iset) )

    def setFieldSplitFields(self, bsize, *fields):
        cdef PetscInt bs = asInt(bsize)
        CHKERR( PCFieldSplitSetBlockSize(self.pc, bs) )
        cdef object name = None
        cdef object field = None
        cdef const char *cname = NULL
        cdef PetscInt nfields = 0, *ifields = NULL
        for name, field in fields:
            name = str2bytes(name, &cname)
            field = iarray_i(field, &nfields, &ifields)
            CHKERR( PCFieldSplitSetFields(self.pc, cname,
                                          nfields, ifields, ifields) )

    def getFieldSplitSubKSP(self):
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        cdef object subksp = None
        try:
            CHKERR( PCFieldSplitGetSubKSP(self.pc, &n, &p) )
            subksp = [ref_KSP(p[i]) for i from 0 <= i <n]
        finally:
            CHKERR( PetscFree(p) )
        return subksp

    def getFieldSplitSchurGetSubKSP(self):
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        cdef object subksp = None
        try:
            CHKERR( PCFieldSplitSchurGetSubKSP(self.pc, &n, &p) )
            subksp = [ref_KSP(p[i]) for i from 0 <= i <n]
        finally:
            CHKERR( PetscFree(p) )
        return subksp

    def setFieldSplitSchurFactType(self, ctype):
        cdef PetscPCFieldSplitSchurFactType cval = ctype
        CHKERR( PCFieldSplitSetSchurFactType(self.pc, cval) )

    def setFieldSplitSchurPreType(self, ptype, Mat pre=None):
        cdef PetscPCFieldSplitSchurPreType pval = ptype
        cdef PetscMat pmat = NULL
        if pre is not None: pmat = pre.mat
        CHKERR( PCFieldSplitSetSchurPre(self.pc, pval, pmat) )

    # --- COMPOSITE ---

    def setCompositeType(self, ctype):
        cdef PetscPCCompositeType cval = ctype
        CHKERR( PCCompositeSetType(self.pc, cval) )

    def getCompositePC(self, n):
        cdef PC pc = PC()
        cdef cn = asInt(n)
        CHKERR( PCCompositeGetPC(self.pc, cn, &pc.pc) )
        PetscINCREF(pc.obj)
        return pc

    def addCompositePCType(self, pc_type):
        cdef PetscPCType cval = NULL
        pc_type = str2bytes(pc_type, &cval)
        CHKERR( PCCompositeAddPCType(self.pc, cval) )

    # --- KSP ---

    def getKSP(self):
        cdef KSP ksp = KSP()
        CHKERR( PCKSPGetKSP(self.pc, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    # --- MG ---

    def getMGType(self):
        cdef PetscPCMGType cval = PC_MG_ADDITIVE
        CHKERR( PCMGGetType(self.pc, &cval) )
        return cval

    def setMGType(self, mgtype):
        cdef PetscPCMGType cval = mgtype
        CHKERR( PCMGSetType(self.pc, cval) )

    def getMGLevels(self):
        cdef PetscInt levels = 0
        CHKERR( PCMGGetLevels(self.pc, &levels) )
        return toInt(levels)

    def setMGLevels(self, levels):
        cdef PetscInt clevels = asInt(levels)
        CHKERR( PCMGSetLevels(self.pc, clevels, NULL) )

    def getMGCoarseSolve(self):
        cdef KSP ksp = KSP()
        CHKERR( PCMGGetCoarseSolve(self.pc, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    def setMGInterpolation(self, level, Mat mat):
        cdef PetscInt clevel = asInt(level)
        CHKERR( PCMGSetInterpolation(self.pc, clevel, mat.mat) )

    def getMGInterpolation(self, level):
        cdef PetscInt clevel = asInt(level)
        cdef Mat interpolation = Mat()
        CHKERR( PCMGGetInterpolation(self.pc, clevel, &interpolation.mat) )
        PetscINCREF(interpolation.obj)
        return interpolation

    def setMGRestriction(self, level, Mat mat):
        cdef PetscInt clevel = asInt(level)
        CHKERR( PCMGSetRestriction(self.pc, clevel, mat.mat) )

    def getMGRestriction(self, level):
        cdef PetscInt clevel = asInt(level)
        cdef Mat restriction = Mat()
        CHKERR( PCMGGetRestriction(self.pc, clevel, &restriction.mat) )
        PetscINCREF(restriction.obj)
        return restriction

    def setMGRScale(self, level, Vec rscale):
        cdef PetscInt clevel = asInt(level)
        CHKERR( PCMGSetRScale(self.pc, clevel, rscale.vec) )

    def getMGRScale(self, level):
        cdef PetscInt clevel = asInt(level)
        cdef Vec rscale = Vec()
        CHKERR( PCMGGetRScale(self.pc, clevel, &rscale.vec) )
        PetscINCREF(rscale.obj)
        return rscale

    def getMGSmoother(self, level):
        cdef PetscInt clevel = asInt(level)
        cdef KSP ksp = KSP()
        CHKERR( PCMGGetSmoother(self.pc, clevel, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    def getMGSmootherDown(self, level):
        cdef PetscInt clevel = asInt(level)
        cdef KSP ksp = KSP()
        CHKERR( PCMGGetSmootherDown(self.pc, clevel, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    def getMGSmootherUp(self, level):
        cdef PetscInt clevel = asInt(level)
        cdef KSP ksp = KSP()
        CHKERR( PCMGGetSmootherUp(self.pc, clevel, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    def setMGCycleType(self, cycle_type):
        cdef PetscPCMGCycleType ctype = cycle_type
        CHKERR( PCMGSetCycleType(self.pc, ctype) )

    def setMGCycleTypeOnLevel(self, level, cycle_type):
        cdef PetscInt clevel = asInt(level)
        cdef PetscPCMGCycleType ctype = cycle_type
        CHKERR( PCMGSetCycleTypeOnLevel(self.pc, clevel, ctype) )

    def setMGRhs(self, level, Vec rhs):
        cdef PetscInt clevel = asInt(level)
        CHKERR( PCMGSetRhs(self.pc, clevel, rhs.vec) )

    def setMGX(self, level, Vec x):
        cdef PetscInt clevel = asInt(level)
        CHKERR( PCMGSetX(self.pc, clevel, x.vec) )

    def setMGR(self, level, Vec r):
        cdef PetscInt clevel = asInt(level)
        CHKERR( PCMGSetR(self.pc, clevel, r.vec) )

    # --- BDDC ---

    def setBDDCDivergenceMat(self, Mat div, trans=False, IS l2l=None):
        cdef PetscBool ptrans = trans
        cdef PetscIS pl2l = NULL
        if l2l is not None: pl2l = l2l.iset
        CHKERR( PCBDDCSetDivergenceMat(self.pc, div.mat, ptrans, pl2l) )

    def setBDDCDiscreteGradient(self, Mat G, order=1, field=1, gord=True, conforming=True):
        cdef PetscInt porder = asInt(order)
        cdef PetscInt pfield = asInt(field)
        cdef PetscBool pgord = gord
        cdef PetscBool pconforming = conforming
        CHKERR( PCBDDCSetDiscreteGradient(self.pc, G.mat, porder, pfield, pgord, pconforming) )

    def setBDDCChangeOfBasisMat(self, Mat T, interior=False):
        cdef PetscBool pinterior = interior
        CHKERR( PCBDDCSetChangeOfBasisMat(self.pc, T.mat, pinterior) )

    def setBDDCPrimalVerticesIS(self, IS primv):
        CHKERR( PCBDDCSetPrimalVerticesIS(self.pc, primv.iset) )

    def setBDDCPrimalVerticesLocalIS(self, IS primv):
        CHKERR( PCBDDCSetPrimalVerticesLocalIS(self.pc, primv.iset) )

    def setBDDCCoarseningRatio(self, cratio):
        cdef PetscInt pcratio = asInt(cratio)
        CHKERR( PCBDDCSetCoarseningRatio(self.pc, pcratio) )

    def setBDDCLevels(self, levels):
        cdef PetscInt plevels = asInt(levels)
        CHKERR( PCBDDCSetLevels(self.pc, plevels) )

    def setBDDCDirichletBoundaries(self, IS bndr):
        CHKERR( PCBDDCSetDirichletBoundaries(self.pc, bndr.iset) )

    def setBDDCDirichletBoundariesLocal(self, IS bndr):
        CHKERR( PCBDDCSetDirichletBoundariesLocal(self.pc, bndr.iset) )

    def setBDDCNeumannBoundaries(self, IS bndr):
        CHKERR( PCBDDCSetNeumannBoundaries(self.pc, bndr.iset) )

    def setBDDCNeumannBoundariesLocal(self, IS bndr):
        CHKERR( PCBDDCSetNeumannBoundariesLocal(self.pc, bndr.iset) )

    def setBDDCDofsSplitting(self, isfields):
        isfields = [isfields] if isinstance(isfields, IS) else list(isfields)
        cdef Py_ssize_t i, n = len(isfields)
        cdef PetscIS  *cisfields = NULL
        cdef object tmp
        tmp = oarray_p(empty_p(n), NULL, <void**>&cisfields)
        for i from 0 <= i < n: cisfields[i] = (<IS?>isfields[i]).iset
        CHKERR( PCBDDCSetDofsSplitting(self.pc, <PetscInt>n, cisfields) )


    def setBDDCDofsSplittingLocal(self, isfields):
        isfields = [isfields] if isinstance(isfields, IS) else list(isfields)
        cdef Py_ssize_t i, n = len(isfields)
        cdef PetscIS  *cisfields = NULL
        cdef object tmp
        tmp = oarray_p(empty_p(n), NULL, <void**>&cisfields)
        for i from 0 <= i < n: cisfields[i] = (<IS?>isfields[i]).iset
        CHKERR( PCBDDCSetDofsSplittingLocal(self.pc, <PetscInt>n, cisfields) )

    # --- Patch ---
    def setPatchCellNumbering(self, Section sec not None):
        CHKERR( PCPatchSetCellNumbering(self.pc, sec.sec) )

    def setPatchDiscretisationInfo(self, dms, bs,
                                   cellNodeMaps,
                                   subspaceOffsets,
                                   ghostBcNodes,
                                   globalBcNodes):
        cdef PetscInt numSubSpaces = 0
        cdef PetscInt numGhostBcs = 0, numGlobalBcs = 0
        cdef PetscInt *nodesPerCell = NULL
        cdef const PetscInt **ccellNodeMaps = NULL
        cdef PetscDM *cdms = NULL
        cdef PetscInt *cbs = NULL
        cdef PetscInt *csubspaceOffsets = NULL
        cdef PetscInt *cghostBcNodes = NULL
        cdef PetscInt *cglobalBcNodes = NULL
        cdef PetscInt i = 0

        bs = iarray_i(bs, &numSubSpaces, &cbs)
        ghostBcNodes = iarray_i(ghostBcNodes, &numGhostBcs, &cghostBcNodes)
        globalBcNodes = iarray_i(globalBcNodes, &numGlobalBcs, &cglobalBcNodes)
        subspaceOffsets = iarray_i(subspaceOffsets, NULL, &csubspaceOffsets)

        CHKERR( PetscMalloc(<size_t>numSubSpaces*sizeof(PetscInt), &nodesPerCell) )
        CHKERR( PetscMalloc(<size_t>numSubSpaces*sizeof(PetscDM), &cdms) )
        CHKERR( PetscMalloc(<size_t>numSubSpaces*sizeof(PetscInt*), &ccellNodeMaps) )
        for i in range(numSubSpaces):
            cdms[i] = (<DM?>dms[i]).dm
            _, nodes = asarray(cellNodeMaps[i]).shape
            cellNodeMaps[i] = iarray_i(cellNodeMaps[i], NULL, <PetscInt**>&(ccellNodeMaps[i]))
            nodesPerCell[i] = asInt(nodes)

        # TODO: refactor on the PETSc side to take ISes?
        CHKERR( PCPatchSetDiscretisationInfo(self.pc, numSubSpaces,
                                             cdms, cbs, nodesPerCell,
                                             ccellNodeMaps, csubspaceOffsets,
                                             numGhostBcs, cghostBcNodes,
                                             numGlobalBcs, cglobalBcNodes) )
        CHKERR( PetscFree(nodesPerCell) )
        CHKERR( PetscFree(cdms) )
        CHKERR( PetscFree(ccellNodeMaps) )

    def setPatchComputeOperator(self, operator, args=None, kargs=None):
        if args is  None: args  = ()
        if kargs is None: kargs = {}
        context = (operator, args, kargs)
        self.set_attr("__patch_compute_operator__", context)
        CHKERR( PCPatchSetComputeOperator(self.pc, PCPatch_ComputeOperator, <void*>context) )

    def setPatchComputeOperatorInteriorFacets(self, operator, args=None, kargs=None):
        if args is  None: args  = ()
        if kargs is None: kargs = {}
        context = (operator, args, kargs)
        self.set_attr("__patch_compute_operator_interior_facets__", context)
        CHKERR( PCPatchSetComputeOperatorInteriorFacets(self.pc, PCPatch_ComputeOperatorInteriorFacets, <void*>context) )

    def setPatchComputeFunction(self, function, args=None, kargs=None):
        if args is  None: args  = ()
        if kargs is None: kargs = {}
        context = (function, args, kargs)
        self.set_attr("__patch_compute_function__", context)
        CHKERR( PCPatchSetComputeFunction(self.pc, PCPatch_ComputeFunction, <void*>context) )

    def setPatchComputeFunctionInteriorFacets(self, function, args=None, kargs=None):
        if args is  None: args  = ()
        if kargs is None: kargs = {}
        context = (function, args, kargs)
        self.set_attr("__patch_compute_function_interior_facets__", context)
        CHKERR( PCPatchSetComputeFunction(self.pc, PCPatch_ComputeFunctionInteriorFacets, <void*>context) )

    def setPatchConstructType(self, typ, operator=None, args=None, kargs=None):
        if args is  None: args  = ()
        if kargs is None: kargs = {}

        if typ in {PC.PatchConstructType.PYTHON, PC.PatchConstructType.USER} and operator is None:
            raise ValueError("Must provide operator for USER or PYTHON type")
        if operator is not None:
            context = (operator, args, kargs)
        else:
            context = None
        self.set_attr("__patch_construction_operator__", context)
        CHKERR( PCPatchSetConstructType(self.pc, typ, PCPatch_UserConstructOperator, <void*>context) )

    # --- HPDDM ---

    def setHPDDMAuxiliaryMat(self, IS uis, Mat uaux):
        CHKERR( PCHPDDMSetAuxiliaryMat(self.pc, uis.iset, uaux.mat, NULL, <void*>NULL) )

    def setHPDDMRHSMat(self, Mat B):
        CHKERR( PCHPDDMSetRHSMat(self.pc, B.mat) )

    def setHPDDMHasNeumannMat(self, has):
        cdef PetscBool phas = has
        CHKERR( PCHPDDMHasNeumannMat(self.pc, phas) )

    def setHPDDMCoarseCorrectionType(self, correction_type):
        cdef PetscPCHPDDMCoarseCorrectionType ctype = correction_type
        CHKERR( PCHPDDMSetCoarseCorrectionType(self.pc, ctype) )

    def getHPDDMCoarseCorrectionType(self):
        cdef PetscPCHPDDMCoarseCorrectionType cval = PC_HPDDM_COARSE_CORRECTION_DEFLATED
        CHKERR( PCHPDDMGetCoarseCorrectionType(self.pc, &cval) )
        return cval

    def getHPDDMSTShareSubKSP(self):
        cdef PetscBool cval = PETSC_FALSE
        CHKERR( PCHPDDMGetSTShareSubKSP(self.pc, &cval) )
        return toBool(cval)

    def setHPDDMDeflationMat(self, IS uis, Mat U):
        CHKERR( PCHPDDMSetDeflationMat(self.pc, uis.iset, U.mat) )

    # --- SPAI ---

    def setSPAIEpsilon(self, val):
        cdef PetscReal cval = asReal(val)
        CHKERR( PCSPAISetEpsilon(self.pc, cval) )

    def setSPAINBSteps(self, nbsteps):
        cdef PetscInt cval = asInt(nbsteps)
        CHKERR( PCSPAISetNBSteps(self.pc, cval) )

    def setSPAIMax(self, maxval):
        cdef PetscInt cval = asInt(maxval)
        CHKERR( PCSPAISetMax(self.pc, cval) )

    def setSPAIMaxNew(self, maxval):
        cdef PetscInt cval = asInt(maxval)
        CHKERR( PCSPAISetMaxNew(self.pc, cval) )

    def setSPAIBlockSize(self, n):
        cdef PetscInt cval = asInt(n)
        CHKERR( PCSPAISetBlockSize(self.pc, cval) )

    def setSPAICacheSize(self, size):
        cdef PetscInt cval = asInt(size)
        CHKERR( PCSPAISetCacheSize(self.pc, cval) )

    def setSPAIVerbose(self, level):
        cdef PetscInt cval = asInt(level)
        CHKERR( PCSPAISetVerbose(self.pc, cval) )

    def setSPAISp(self, sym):
        cdef PetscInt cval = asInt(sym)
        CHKERR( PCSPAISetSp(self.pc, cval) )

    # --- DEFLATION ---

    def setDeflationInitOnly(self, flg):
        cdef PetscBool cflg = asBool(flg)
        CHKERR( PCDeflationSetInitOnly(self.pc, cflg) )

    def setDeflationLevels(self, levels):
        cdef PetscInt clevels = asInt(levels)
        CHKERR( PCDeflationSetLevels(self.pc, clevels) )

    def setDeflationReductionFactor(self, red):
        cdef PetscInt cred = asInt(red)
        CHKERR( PCDeflationSetReductionFactor(self.pc, cred) )

    def setDeflationCorrectionFactor(self, fact):
        cdef PetscScalar cfact = asScalar(fact)
        CHKERR( PCDeflationSetCorrectionFactor(self.pc, fact) )

    def setDeflationSpaceToCompute(self, space_type, size):
        cdef PetscInt csize = asInt(size)
        cdef PetscPCDeflationSpaceType ctype = space_type
        CHKERR( PCDeflationSetSpaceToCompute(self.pc, space_type, csize) )

    def setDeflationSpace(self, Mat W, transpose):
        cdef PetscBool ctranspose = asBool(transpose)
        CHKERR( PCDeflationSetSpace(self.pc, W.mat, ctranspose) )

    def setDeflationProjectionNullSpaceMat(self, Mat mat):
        CHKERR( PCDeflationSetProjectionNullSpaceMat(self.pc, mat.mat) )

    def setDeflationCoarseMat(self, Mat mat):
        CHKERR( PCDeflationSetCoarseMat(self.pc, mat.mat) )

    def getDeflationCoarseKSP(self):
        cdef KSP ksp = KSP()
        CHKERR( PCDeflationGetCoarseKSP(self.pc, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    def getDeflationPC(self):
        cdef PC apc = PC()
        CHKERR( PCDeflationGetPC(self.pc, &apc.pc) )
        PetscINCREF(apc.obj)
        return apc

# --------------------------------------------------------------------

del PCType
del PCSide
del PCASMType
del PCGASMType
del PCMGType
del PCMGCycleType
del PCGAMGType
del PCCompositeType
del PCFieldSplitSchurPreType
del PCFieldSplitSchurFactType
del PCPatchConstructType
del PCHPDDMCoarseCorrectionType
del PCDeflationSpaceType
del PCFailedReason

# --------------------------------------------------------------------
