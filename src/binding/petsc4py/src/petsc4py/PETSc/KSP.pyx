# --------------------------------------------------------------------

class KSPType(object):
    """

    See Also
    --------
    petsc.KSP
    resolves to `https://petsc.org/release/docs/manualpages/KSP/KSPType/`__

    """
    RICHARDSON = S_(KSPRICHARDSON)
    CHEBYSHEV  = S_(KSPCHEBYSHEV)
    CG         = S_(KSPCG)
    GROPPCG    = S_(KSPGROPPCG)
    PIPECG     = S_(KSPPIPECG)
    PIPECGRR   = S_(KSPPIPECGRR)
    PIPELCG    = S_(KSPPIPELCG)
    PIPEPRCG   = S_(KSPPIPEPRCG)
    PIPECG2    = S_(KSPPIPECG2)
    CGNE       = S_(KSPCGNE)
    NASH       = S_(KSPNASH)
    STCG       = S_(KSPSTCG)
    GLTR       = S_(KSPGLTR)
    FCG        = S_(KSPFCG)
    PIPEFCG    = S_(KSPPIPEFCG)
    GMRES      = S_(KSPGMRES)
    PIPEFGMRES = S_(KSPPIPEFGMRES)
    FGMRES     = S_(KSPFGMRES)
    LGMRES     = S_(KSPLGMRES)
    DGMRES     = S_(KSPDGMRES)
    PGMRES     = S_(KSPPGMRES)
    TCQMR      = S_(KSPTCQMR)
    BCGS       = S_(KSPBCGS)
    IBCGS      = S_(KSPIBCGS)
    QMRCGS     = S_(KSPQMRCGS)
    FBCGS      = S_(KSPFBCGS)
    FBCGSR     = S_(KSPFBCGSR)
    BCGSL      = S_(KSPBCGSL)
    PIPEBCGS   = S_(KSPPIPEBCGS)
    CGS        = S_(KSPCGS)
    TFQMR      = S_(KSPTFQMR)
    CR         = S_(KSPCR)
    PIPECR     = S_(KSPPIPECR)
    LSQR       = S_(KSPLSQR)
    PREONLY    = S_(KSPPREONLY)
    NONE       = S_(KSPNONE)
    QCG        = S_(KSPQCG)
    BICG       = S_(KSPBICG)
    MINRES     = S_(KSPMINRES)
    SYMMLQ     = S_(KSPSYMMLQ)
    LCD        = S_(KSPLCD)
    PYTHON     = S_(KSPPYTHON)
    GCR        = S_(KSPGCR)
    PIPEGCR    = S_(KSPPIPEGCR)
    TSIRM      = S_(KSPTSIRM)
    CGLS       = S_(KSPCGLS)
    FETIDP     = S_(KSPFETIDP)
    HPDDM      = S_(KSPHPDDM)

class KSPNormType(object):
    # native
    NORM_DEFAULT          = KSP_NORM_DEFAULT
    NORM_NONE             = KSP_NORM_NONE
    NORM_PRECONDITIONED   = KSP_NORM_PRECONDITIONED
    NORM_UNPRECONDITIONED = KSP_NORM_UNPRECONDITIONED
    NORM_NATURAL          = KSP_NORM_NATURAL
    # aliases
    DEFAULT          = NORM_DEFAULT
    NONE = NO        = NORM_NONE
    PRECONDITIONED   = NORM_PRECONDITIONED
    UNPRECONDITIONED = NORM_UNPRECONDITIONED
    NATURAL          = NORM_NATURAL

class KSPConvergedReason(object):
    #iterating
    CONVERGED_ITERATING       = KSP_CONVERGED_ITERATING
    ITERATING                 = KSP_CONVERGED_ITERATING
    # converged
    CONVERGED_RTOL_NORMAL     = KSP_CONVERGED_RTOL_NORMAL
    CONVERGED_ATOL_NORMAL     = KSP_CONVERGED_ATOL_NORMAL
    CONVERGED_RTOL            = KSP_CONVERGED_RTOL
    CONVERGED_ATOL            = KSP_CONVERGED_ATOL
    CONVERGED_ITS             = KSP_CONVERGED_ITS
    CONVERGED_NEG_CURVE       = KSP_CONVERGED_NEG_CURVE
    CONVERGED_STEP_LENGTH     = KSP_CONVERGED_STEP_LENGTH
    CONVERGED_HAPPY_BREAKDOWN = KSP_CONVERGED_HAPPY_BREAKDOWN
    # diverged
    DIVERGED_NULL             = KSP_DIVERGED_NULL
    DIVERGED_MAX_IT           = KSP_DIVERGED_MAX_IT
    DIVERGED_DTOL             = KSP_DIVERGED_DTOL
    DIVERGED_BREAKDOWN        = KSP_DIVERGED_BREAKDOWN
    DIVERGED_BREAKDOWN_BICG   = KSP_DIVERGED_BREAKDOWN_BICG
    DIVERGED_NONSYMMETRIC     = KSP_DIVERGED_NONSYMMETRIC
    DIVERGED_INDEFINITE_PC    = KSP_DIVERGED_INDEFINITE_PC
    DIVERGED_NANORINF         = KSP_DIVERGED_NANORINF
    DIVERGED_INDEFINITE_MAT   = KSP_DIVERGED_INDEFINITE_MAT
    DIVERGED_PCSETUP_FAILED   = KSP_DIVERGED_PC_FAILED

# --------------------------------------------------------------------

cdef class KSP(Object):
    """

    Further info in

    `KSP: Linear System Solvers <https://petsc.org/release/docs/manual/ksp/>`__
    """

    Type            = KSPType
    NormType        = KSPNormType
    ConvergedReason = KSPConvergedReason

    # --- xxx ---

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ksp
        self.ksp = NULL

    def __call__(self, b, x=None):
        if x is None: # XXX do this better
            x = self.getOperators()[0].createVecLeft()
        self.solve(b, x)
        return x

    # --- xxx ---

    def view(self, Viewer viewer=None):
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( KSPView(self.ksp, vwr) )

    def destroy(self):
        CHKERR( KSPDestroy(&self.ksp) )
        return self

    def create(self, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscKSP newksp = NULL
        CHKERR( KSPCreate(ccomm, &newksp) )
        PetscCLEAR(self.obj); self.ksp = newksp
        return self

    def setType(self, ksp_type: KSP.Type | str) -> None:
        """Build the `KSP` data structure for a particular `KSP.Type`.

        Logically collective.

        Parameters
        ----------
        ksp_type
            KSP Type object

        Notes
        -----
        See `petsc.KSPType` for available methods (for instance, `KSP.CG` or
        `KSP.GMRES`).

        Normally, it is best to use the `KSP.SetFromOptions` command
        and then set the KSP type from the options database rather than
        by using this routine. Using the options database provides the
        user with maximum flexibility in evaluating the many different
        Krylov methods. The `KSP.SetType` routine is provided for those
        situations where it is necessary to set the iterative solver
        independently of the command line or options database. This
        might be the case, for example, when the choice of iterative
        solver changes during the execution of the program, and the
        user's application is taking responsibility for choosing the
        appropriate method. In other words, this routine is not for
        beginners.

        See also
        --------
        petsc.KSPSetType

        """
        cdef PetscKSPType cval = NULL
        ksp_type = str2bytes(ksp_type, &cval)
        CHKERR( KSPSetType(self.ksp, cval) )

    def getType(self) -> str:
        """Get the KSP type as a string from the `KSP` object.

        Not collective.

        See also
        --------
        petsc.KSPGetType

        """
        cdef PetscKSPType cval = NULL
        CHKERR( KSPGetType(self.ksp, &cval) )
        return bytes2str(cval)

    def setOptionsPrefix(self, prefix: str) -> None:
        """Set the prefix used for all `KSP` options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The options prefix.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name. The first character of all runtime options is
        AUTOMATICALLY the hyphen. For example, to distinguish between
        the runtime options for two different `KSP` contexts, one could
        call
        ```
        KSPSetOptionsPrefix(ksp1, "sys1_")
        KSPSetOptionsPrefix(ksp2, "sys2_")
        ```

        This would enable use of different options for each system,
        such as
        ```
        -sys1_ksp_type gmres -sys1_ksp_rtol 1.e-3
        -sys2_ksp_type bcgs  -sys2_ksp_rtol 1.e-4
        ```

        See also
        --------
        petsc.KSPSetOptionsPrefix

        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( KSPSetOptionsPrefix(self.ksp, cval) )

    def getOptionsPrefix(self) -> str:
        """Get the prefix used for all `KSP` options in the database.

        Not collective.

        See also
        --------
        petsc.KSPGetOptionsPrefix

        """
        cdef const char *cval = NULL
        CHKERR( KSPGetOptionsPrefix(self.ksp, &cval) )
        return bytes2str(cval)

    def appendOptionsPrefix(self, prefix: str) -> None:
        """Append to prefix used for all `KSP` options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The options prefix to append.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name. The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        See also
        --------
        petsc.KSPAppendOptionsPrefix

        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( KSPAppendOptionsPrefix(self.ksp, cval) )

    def setFromOptions(self) -> None:
        """Sets `KSP` options from the options database.

        Collective.

        This routine must be called before `KSP.SetUp` if the user is
        to be allowed to set the Krylov type.

        See also
        --------
        petsc_options, petsc.KSPSetFromOptions

        """
        CHKERR( KSPSetFromOptions(self.ksp) )

    # --- application context ---

    def setAppCtx(self, appctx: Any) -> None:
        """Set the optional user-defined context for the linear solver.

        Logically collective.

        Parameters
        ----------
        appctx
            The user defined context

        Notes
        -----
        The user context is a way for users to attach any information
        to the `KSP` that they may need later when interacting with the
        KSP.
        Use `KSP.getAppCtx` to get access to the context at a later
        time.

        See also
        --------
        petsc.KSPSetApplicationContext

        """
        self.set_attr('__appctx__', appctx)

    def getAppCtx(self) -> Any:
        """Get the user-defined context for the linear solver.

        Not collective

        See also
        --------
        petsc.KSPGetApplicationContext

        """
        return self.get_attr('__appctx__')

    # --- discretization space ---

    def getDM(self) -> DM:
        """Get the `DM` that may be used by some preconditioners.

        Not collective.

        See also
        --------
        KSP, DM, petsc.KSPGetDM

        """
        cdef PetscDM newdm = NULL
        CHKERR( KSPGetDM(self.ksp, &newdm) )
        cdef DM dm = subtype_DM(newdm)()
        dm.dm = newdm
        PetscINCREF(dm.obj)
        return dm

    def setDM(self, DM dm) -> None:
        """Set the `DM` that may be used by some preconditioners

        Logically collective.

        Parameters
        ----------
        dm
            the dm, cannot be ``None``

        Notes
        -----
        If this is used then the `KSP` will attempt to use the `DM` to
        create the matrix and use the routine set with
        `DM.setKSPComputeOperators`. Use ``KSP.setDMActive(False)``
        to instead use the matrix you have provided with
        `KSP.setOperators`.

        A `DM` can only be used for solving one problem at a time
        because information about the problem is stored on the DM, even
        when not using interfaces like `DM.setKSPComputeOperators`. Use
        `DM.clone` to get a distinct `DM` when solving different
        problems using the same function space.

        See also
        --------
        KSP, DM, DM.setKSPComputeOperators, KSP.setOperators, DM.clone,
        petsc.KSPSetDM

        """
        CHKERR( KSPSetDM(self.ksp, dm.dm) )

    def setDMActive(self, bint flag) -> None:
        """`DM` should be used to generate system matrix & RHS vector.

        Logically collective

        Parameters
        ----------
        flag
            Boolean whether to use the `DM` (or not)

        Notes
        -----
        By default `KSP.setDM` sets the `DM` as active, call
        ``KSP.setDMActive(False)`` after ``KSP.setDM(dm)`` to not
        have the `KSP` object use the `DM` to generate the matrices.

        See also
        --------
        KSP, DM, KSP.setDM, petsc.KSPSetDMActive

        """
        cdef PetscBool cflag = PETSC_FALSE
        if flag: cflag = PETSC_TRUE
        CHKERR( KSPSetDMActive(self.ksp, cflag) )

    # --- operators and preconditioner ---

    def setComputeRHS(
        self,
        rhs: KSPComputeRHSFunction,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None
    ) -> None:
        """Set routine to compute the right hand side of the linear system.

        Logically collective.

        Parameters
        ----------
        rhs
            Function which computes the right hand side.
        args
            Positional arguments for callback function ``rhs``.
        kargs
            Keyword arguments for callback function ``rhs``.

        Notes
        -----
        The routine you provide will be called EACH you call `KSP.solve`
        to prepare the new right hand side for that solve.

        See also
        --------
        KSP, KSP.solve, petsc.KSPSetComputeRHS

        """
        if args  is None: args  = ()
        if kargs is None: kargs = {}
        context = (rhs, args, kargs)
        self.set_attr('__rhs__', context)
        CHKERR( KSPSetComputeRHS(self.ksp, KSP_ComputeRHS, <void*>context) )

    def setComputeOperators(
        self,
        operators: KSPComputeOperatorsFunction,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None
    ) -> None:
        """Set routine to compute the linear operators.

        Logically collective.

        Parameters
        ----------
        operators
            Function which computes the operators.
        args
            Positional arguments for callback function ``operators``.
        kargs
            Keyword arguments for callback function ``operators``.

        Notes
        -----
        The user provided function `operators` will be called
        automatically at the very next call to `KSP.solve`. It will NOT
        be called at future `KSP.solve` calls unless either
        `KSP.setComputeOperators` or `KSP.setOperators` is called
        before that `KSP.solve` is called. This allows the same system
        to be solved several times with different right hand side
        functions, but is a confusing API since one might expect it to
        be called for each `KSP.solve`.

        To reuse the same preconditioner for the next `KSP.solve` and
        not compute a new one based on the most recently computed
        matrix call `KSP.setReusePreconditioner`.

        See also
        --------
        KSP, KSP.solve, KSP.setOperators, KSP.setReusePreconditioner,
        petsc.KSPSetComputeOperators

        """
        if args  is None: args  = ()
        if kargs is None: kargs = {}
        context = (operators, args, kargs)
        self.set_attr('__operators__', context)
        CHKERR( KSPSetComputeOperators(self.ksp, KSP_ComputeOps, <void*>context) )

    def setOperators(self, Mat A=None, Mat P=None) -> None:
        """Set matrix associated with the linear system.

        Sets the matrix associated with the linear system and a
        (possibly) different one from which the preconditioner will be
        built.

        Collective.

        Parameters
        ----------
        A
            Matrix that defines the linear system.
        P
            Matrix to be used in constructing the preconditioner,
            usually the same as ``A``.

        Notes
        -----
        If you know the operator ``A`` has a null space you can use
        `Mat.setNullSpace` and `Mat.setTransposeNullSpace` to supply the
        null space to ``A`` and the `KSP` solvers will automatically use
        that null space as needed during the solution process.

        All future calls to `KSP.setOperators` must use the same size
        matrices!

        Passing ``None`` for ``A`` or ``P`` removes the matrix that is
        currently used.

        See also
        --------
        KSP, KSP.solve, KSP.setComputeOperators,
        petsc.KSPSetOperators

        """
        cdef PetscMat amat=NULL
        if A is not None: amat = A.mat
        cdef PetscMat pmat=amat
        if P is not None: pmat = P.mat
        CHKERR( KSPSetOperators(self.ksp, amat, pmat) )

    def getOperators(self) -> tuple[Mat, Mat]:
        """Get the matrix associated with the linear system.

        Gets the matrix associated with the linear system and a
        (possibly) different one used to construct the preconditioner.

        Collective.

        Returns
        -------
        A: Mat
            Matrix that defines the linear system.
        P: Mat
            Matrix to be used in constructing the preconditioner,
            usually the same as ``A``.

        See also
        --------
        KSP, KSP.solve, KSP.setOperators, petsc.KSPGetOperators

        """
        cdef Mat A = Mat(), P = Mat()
        CHKERR( KSPGetOperators(self.ksp, &A.mat, &P.mat) )
        PetscINCREF(A.obj)
        PetscINCREF(P.obj)
        return (A, P)

    def setPC(self, PC pc) -> None:
        """Set the preconditioner.

        Sets the preconditioner to be used to calculate the application
        of the preconditioner on a vector.

        Collective.

        Parameters
        ----------
        pc
            The preconditioner object

        See also
        --------
        KSP, KSP.getPC, petsc.KSPSetPC

        """
        CHKERR( KSPSetPC(self.ksp, pc.pc) )

    def getPC(self) -> PC:
        """Return the preconditioner.

        Not collective.

        See also
        --------
        KSP, KSP.setPC, petsc.KSPGetPC

        """
        cdef PC pc = PC()
        CHKERR( KSPGetPC(self.ksp, &pc.pc) )
        PetscINCREF(pc.obj)
        return pc

    # --- tolerances and convergence ---

    def setTolerances(
        self,
        float rtol=None,
        float atol=None,
        float divtol=None,
        int max_it=None
    ) -> None:
        """Set various tolerances used by the KSP convergence testers.

        Sets the relative, absolute, divergence, and maximum iteration
        tolerances used by the default KSP convergence testers.

        Logically collective.

        Parameters
        ----------
        rtol
            The relative convergence tolerance, relative decrease in
            the (possibly preconditioned) residual norm.
        atol
            The absolute convergence tolerance absolute size of the
            (possibly preconditioned) residual norm.
        dtol
            The divergence tolerance, amount (possibly preconditioned)
            residual norm can increase before `KSP.convergedDefault`
            concludes that the method is diverging.
        max_it
            Maximum number of iterations to use.

        Notes
        -----
        Use `PETSC_DEFAULT` to retain the default value of any of the
        tolerances.

        See also
        --------
        petsc_options, KSP.getTolerances, KSP.convergedDefault,
        KSP.setConvergenceTest, petsc.KSPsetTolerances


        """
        cdef PetscReal crtol, catol, cdivtol
        crtol = catol = cdivtol = PETSC_DEFAULT;
        if rtol   is not None: crtol   = asReal(rtol)
        if atol   is not None: catol   = asReal(atol)
        if divtol is not None: cdivtol = asReal(divtol)
        cdef PetscInt cmaxits = PETSC_DEFAULT
        if max_it is not None: cmaxits = asInt(max_it)
        CHKERR( KSPSetTolerances(self.ksp, crtol, catol, cdivtol, cmaxits) )

    def getTolerances(self) -> tuple[float, float, float, int]:
        """Get various tolerances used by the KSP convergence tests.

        Gets the relative, absolute, divergence, and maximum iteration
        tolerances used by the default KSP convergence tests.

        Not collective.

        Returns
        -------
        rtol: float
            The relative convergence tolerance
        atol: float
            The absolute convergence tolerance
        dtol: float
            The divergence tolerance
        maxits: int
            Maximum number of iterations

        See also
        --------
        KSP.setTolerances, petsc.KSPGetTolerances

        """
        cdef PetscReal crtol=0, catol=0, cdivtol=0
        cdef PetscInt cmaxits=0
        CHKERR( KSPGetTolerances(self.ksp, &crtol, &catol, &cdivtol, &cmaxits) )
        return (toReal(crtol), toReal(catol), toReal(cdivtol), toInt(cmaxits))

    def setConvergenceTest(
        self,
        converged: KSPConvergenceTestFunction,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None
    ) -> None:
        """Set the function to be used to determine convergence.

        Logically collective.

        Parameters
        ----------
        converged
            Callback which computes the convergence.
        args
            Positional arguments for callback function ``converged``.
        kargs
            Keyword arguments for callback function ``converged``.

        Notes
        -----
        Must be called after the KSP type has been set so put this
        after a call to `KSP.setType`, or `KSP.setFromOptions`.

        The default convergence test, `KSP.convergedDefault`, aborts if
        the residual grows to more than 10000 times the initial
        residual.

        The default is a combination of relative and absolute
        tolerances. The residual value that is tested may be an
        approximation; routines that need exact values should compute
        them.

        In the default PETSc convergence test, the precise values of
        reason are macros such as ``KSP_CONVERGED_RTOL``, which are
        defined in ``petscksp.h``.

        See also
        --------
        KSP.convergedDefault, KSP.getConvergenceContext, KSP.setTolerances,
        KSP.getConvergenceTest, KSP.getAndClearConvergenceTest,
        petsc.KSPSetConvergenceTest

        """
        cdef PetscKSPNormType normtype = KSP_NORM_NONE
        cdef void* cctx = NULL
        if converged is not None:
            CHKERR( KSPSetConvergenceTest(
                    self.ksp, KSP_Converged, NULL, NULL) )
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__converged__', (converged, args, kargs))
        else:
            CHKERR( KSPGetNormType(self.ksp, &normtype) )
            if normtype != KSP_NORM_NONE:
                CHKERR( KSPConvergedDefaultCreate(&cctx) )
                CHKERR( KSPSetConvergenceTest(
                        self.ksp, KSPConvergedDefault,
                        cctx, KSPConvergedDefaultDestroy) )
            else:
                CHKERR( KSPSetConvergenceTest(
                        self.ksp, KSPConvergedSkip,
                        NULL, NULL) )
            self.set_attr('__converged__', None)

    def getConvergenceTest(self) -> KSPConvergenceTestFunction:
        """Gets the function to be used to determine convergence.

        Logically collective.

        See also
        --------
        KSP.convergedDefault, KSP.getConvergenceContext, KSP.setTolerances,
        KSP.getConvergenceTest, KSP.getAndClearConvergenceTest,
        petsc.KSPGetConvergenceTest
        """
        return self.get_attr('__converged__')

    def callConvergenceTest(self, int its, float rnorm):
        """Call the convergence test callback.

        Parameters
        ----------
        its
            Number of iterations.
        rnorm
            The residual norm.

        Notes
        -----
        This functionality is implemented in petsc4py.

        """
        cdef PetscInt  ival = asInt(its)
        cdef PetscReal rval = asReal(rnorm)
        cdef PetscKSPConvergedReason reason = KSP_CONVERGED_ITERATING
        CHKERR( KSPConvergenceTestCall(self.ksp, ival, rval, &reason) )
        return reason

    def setConvergenceHistory(
        self,
        length: int | None = None,
        reset: bool = False
    ) -> None:
        """Set the array used to hold the residual history.

        If set, this array will contain the residual norms computed at
        each iteration of the solver.

        Not collective.

        Parameters
        ----------
        length
            Length of array to store history in.
        reset
            ``True`` indicates the history counter is reset to zero for
            each new linear solve.

        Notes
        -----
        If ``length`` is not provided or ``None`` then a default array
        of length 10000 is allocated.

        If the array is not long enough then once the iterations is
        longer than the array length `KSP.solve` stops recording the
        history.

        See also
        --------
        KSP.getResidualHistory, petsc.KSPSetResidualHistory

        """
        cdef PetscReal *data = NULL
        cdef PetscInt   size = 10000
        cdef PetscBool flag = PETSC_FALSE
        if   length is True:     pass
        elif length is not None: size = asInt(length)
        if size < 0: size = 10000
        if reset: flag = PETSC_TRUE
        cdef object hist = oarray_r(empty_r(size), NULL, &data)
        self.set_attr('__history__', hist)
        CHKERR( KSPSetResidualHistory(self.ksp, data, size, flag) )

    def getConvergenceHistory(self) -> ndarray:
        """Get array containing the residual history.

        Not collective.

        See also
        --------
        KSP.setResidualHistory, petsc.KSPGetResidualHistory

        """
        cdef const PetscReal *data = NULL
        cdef PetscInt   size = 0
        CHKERR( KSPGetResidualHistory(self.ksp, &data, &size) )
        return array_r(size, data)

    def logConvergenceHistory(self, float rnorm):
        """Add residual to convergence history.

        Parameters
        ----------
        rnorm
            Residual norm to be added to convergence history.

        See also
        --------
        petsc.KSPLogResidualHistory

        """
        cdef PetscReal rval = asReal(rnorm)
        CHKERR( KSPLogResidualHistory(self.ksp, rval) )

    # --- monitoring ---

    def setMonitor(self, monitor, args=None, kargs=None):
        if monitor is None: return
        cdef object monitorlist = self.get_attr('__monitor__')
        if monitorlist is None:
            monitorlist = []
            self.set_attr('__monitor__', monitorlist)
            CHKERR( KSPMonitorSet(self.ksp, KSP_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self):
        return self.get_attr('__monitor__')

    def monitorCancel(self):
        CHKERR( KSPMonitorCancel(self.ksp) )
        self.set_attr('__monitor__', None)

    cancelMonitor = monitorCancel

    def monitor(self, its, rnorm):
        cdef PetscInt  ival = asInt(its)
        cdef PetscReal rval = asReal(rnorm)
        CHKERR( KSPMonitor(self.ksp, ival, rval) )

    # --- customization ---

    def setPCSide(self, side):
        CHKERR( KSPSetPCSide(self.ksp, side) )

    def getPCSide(self):
        cdef PetscPCSide side = PC_LEFT
        CHKERR( KSPGetPCSide(self.ksp, &side) )
        return side

    def setNormType(self, normtype):
        CHKERR( KSPSetNormType(self.ksp, normtype) )

    def getNormType(self):
        cdef PetscKSPNormType normtype = KSP_NORM_NONE
        CHKERR( KSPGetNormType(self.ksp, &normtype) )
        return normtype

    def setComputeEigenvalues(self, bint flag):
        cdef PetscBool compute = PETSC_FALSE
        if flag: compute = PETSC_TRUE
        CHKERR( KSPSetComputeEigenvalues(self.ksp, compute) )

    def getComputeEigenvalues(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( KSPGetComputeEigenvalues(self.ksp, &flag) )
        return toBool(flag)

    def setComputeSingularValues(self, bint flag):
        cdef PetscBool compute = PETSC_FALSE
        if flag: compute = PETSC_TRUE
        CHKERR( KSPSetComputeSingularValues(self.ksp, compute) )

    def getComputeSingularValues(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( KSPGetComputeSingularValues(self.ksp, &flag) )
        return toBool(flag)

    # --- initial guess ---

    def setInitialGuessNonzero(self, bint flag):
        cdef PetscBool guess_nonzero = PETSC_FALSE
        if flag: guess_nonzero = PETSC_TRUE
        CHKERR( KSPSetInitialGuessNonzero(self.ksp, guess_nonzero) )

    def getInitialGuessNonzero(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( KSPGetInitialGuessNonzero(self.ksp, &flag) )
        return toBool(flag)

    def setInitialGuessKnoll(self, bint flag):
        cdef PetscBool guess_knoll = PETSC_FALSE
        if flag: guess_knoll = PETSC_TRUE
        CHKERR( KSPSetInitialGuessKnoll(self.ksp, guess_knoll) )

    def getInitialGuessKnoll(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( KSPGetInitialGuessKnoll(self.ksp, &flag) )
        return toBool(flag)

    def setUseFischerGuess(self, model, size):
        cdef PetscInt ival1 = asInt(model)
        cdef PetscInt ival2 = asInt(size)
        CHKERR( KSPSetUseFischerGuess(self.ksp, ival1, ival2) )

    # --- solving ---

    def setUp(self):
        CHKERR( KSPSetUp(self.ksp) )

    def reset(self):
        CHKERR( KSPReset(self.ksp) )

    def setUpOnBlocks(self):
        CHKERR( KSPSetUpOnBlocks(self.ksp) )

    def solve(self, Vec b or None, Vec x or None):
        cdef PetscVec b_vec = NULL
        cdef PetscVec x_vec = NULL
        if b is not None: b_vec = b.vec
        if x is not None: x_vec = x.vec
        CHKERR( KSPSolve(self.ksp, b_vec, x_vec) )

    def solveTranspose(self, Vec b, Vec x):
        CHKERR( KSPSolveTranspose(self.ksp, b.vec, x.vec) )

    def matSolve(self, Mat B, Mat X):
        CHKERR( KSPMatSolve(self.ksp, B.mat, X.mat) )

    def matSolveTranspose(self, Mat B, Mat X):
        CHKERR( KSPMatSolveTranspose(self.ksp, B.mat, X.mat) )

    def setIterationNumber(self, its):
        cdef PetscInt ival = asInt(its)
        CHKERR( KSPSetIterationNumber(self.ksp, ival) )

    def getIterationNumber(self):
        cdef PetscInt ival = 0
        CHKERR( KSPGetIterationNumber(self.ksp, &ival) )
        return toInt(ival)

    def setResidualNorm(self, rnorm):
        cdef PetscReal rval = asReal(rnorm)
        CHKERR( KSPSetResidualNorm(self.ksp, rval) )

    def getResidualNorm(self):
        cdef PetscReal rval = 0
        CHKERR( KSPGetResidualNorm(self.ksp, &rval) )
        return toReal(rval)

    def setConvergedReason(self, reason):
        cdef PetscKSPConvergedReason val = reason
        CHKERR( KSPSetConvergedReason(self.ksp, val) )

    def getConvergedReason(self):
        cdef PetscKSPConvergedReason reason = KSP_CONVERGED_ITERATING
        CHKERR( KSPGetConvergedReason(self.ksp, &reason) )
        return reason

    def setErrorIfNotConverged(self, bint flag):
        cdef PetscBool ernc = PETSC_FALSE
        if flag: ernc = PETSC_TRUE
        CHKERR( KSPSetErrorIfNotConverged(self.ksp, ernc) )

    def getErrorIfNotConverged(self):
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( KSPGetErrorIfNotConverged(self.ksp, &flag) )
        return toBool(flag)

    def getRhs(self):
        cdef Vec vec = Vec()
        CHKERR( KSPGetRhs(self.ksp, &vec.vec) )
        PetscINCREF(vec.obj)
        return vec

    def getSolution(self):
        cdef Vec vec = Vec()
        CHKERR( KSPGetSolution(self.ksp, &vec.vec) )
        PetscINCREF(vec.obj)
        return vec

    def getWorkVecs(self, right=None, left=None):
        cdef bint R = right is not None
        cdef bint L = left  is not None
        cdef PetscInt i=0, nr=0, nl=0
        cdef PetscVec *vr=NULL, *vl=NULL
        if R: nr = asInt(right)
        if L: nl = asInt(left)
        cdef object vecsr = [] if R else None
        cdef object vecsl = [] if L else None
        CHKERR( KSPCreateVecs(self.ksp, nr, &vr, nl, &vr) )
        try:
            for i from 0 <= i < nr:
                vecsr.append(ref_Vec(vr[i]))
            for i from 0 <= i < nl:
                vecsl.append(ref_Vec(vl[i]))
        finally:
            if nr > 0 and vr != NULL:
                VecDestroyVecs(nr, &vr) # XXX errors?
            if nl > 0 and vl !=NULL:
                VecDestroyVecs(nl, &vl) # XXX errors?
        #
        if R and L: return (vecsr, vecsl)
        elif R:     return vecsr
        elif L:     return vecsl
        else:       return None

    def buildSolution(self, Vec x=None):
        if x is None: x = Vec()
        if x.vec == NULL:
            CHKERR( KSPGetSolution(self.ksp, &x.vec) )
            CHKERR( VecDuplicate(x.vec, &x.vec) )
        CHKERR( KSPBuildSolution(self.ksp, x.vec, NULL) )
        return x

    def buildResidual(self, Vec r=None):
        if r is None: r = Vec()
        if r.vec == NULL:
            CHKERR( KSPGetRhs(self.ksp, &r.vec) )
            CHKERR( VecDuplicate(r.vec, &r.vec) )
        CHKERR( KSPBuildResidual(self.ksp , NULL, r.vec, &r.vec) )
        return r

    def computeEigenvalues(self):
        cdef PetscInt its = 0
        cdef PetscInt neig = 0
        cdef PetscReal *rdata = NULL
        cdef PetscReal *idata = NULL
        CHKERR( KSPGetIterationNumber(self.ksp, &its) )
        cdef ndarray r = oarray_r(empty_r(its), NULL, &rdata)
        cdef ndarray i = oarray_r(empty_r(its), NULL, &idata)
        CHKERR( KSPComputeEigenvalues(self.ksp, its, rdata, idata, &neig) )
        eigen = empty_c(neig)
        eigen.real = r[:neig]
        eigen.imag = i[:neig]
        return eigen

    def computeExtremeSingularValues(self):
        cdef PetscReal smax = 0
        cdef PetscReal smin = 0
        CHKERR( KSPComputeExtremeSingularValues(self.ksp, &smax, &smin) )
        return smax, smin

    # --- GMRES ---

    def setGMRESRestart(self, restart):
        cdef PetscInt ival = asInt(restart)
        CHKERR( KSPGMRESSetRestart(self.ksp, ival) )

    # --- Python ---

    def createPython(self, context=None, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscKSP newksp = NULL
        CHKERR( KSPCreate(ccomm, &newksp) )
        PetscCLEAR(self.obj); self.ksp = newksp
        CHKERR( KSPSetType(self.ksp, KSPPYTHON) )
        CHKERR( KSPPythonSetContext(self.ksp, <void*>context) )
        return self

    def setPythonContext(self, context):
        CHKERR( KSPPythonSetContext(self.ksp, <void*>context) )

    def getPythonContext(self):
        cdef void *context = NULL
        CHKERR( KSPPythonGetContext(self.ksp, &context) )
        if context == NULL: return None
        else: return <object> context

    def setPythonType(self, py_type):
        cdef const char *cval = NULL
        py_type = str2bytes(py_type, &cval)
        CHKERR( KSPPythonSetType(self.ksp, cval) )

    def getPythonType(self):
        cdef const char *cval = NULL
        CHKERR( KSPPythonGetType(self.ksp, &cval) )
        return bytes2str(cval)

    # --- application context ---

    property appctx:
        def __get__(self):
            return self.getAppCtx()
        def __set__(self, value):
            self.setAppCtx(value)

    # --- discretization space ---

    property dm:
        def __get__(self):
            return self.getDM()
        def __set__(self, value):
            self.setDM(value)

    # --- vectors ---

    property vec_sol:
        def __get__(self):
            return self.getSolution()

    property vec_rhs:
        def __get__(self):
            return self.getRhs()

    # --- operators ---

    property mat_op:
        def __get__(self):
            return self.getOperators()[0]

    property mat_pc:
        def __get__(self):
            return self.getOperators()[1]

    # --- initial guess ---

    property guess_nonzero:
        def __get__(self):
            return self.getInitialGuessNonzero()
        def __set__(self, value):
            self.setInitialGuessNonzero(value)

    property guess_knoll:
        def __get__(self):
            return self.getInitialGuessKnoll()
        def __set__(self, value):
            self.setInitialGuessKnoll(value)

    # --- preconditioner ---

    property pc:
        def __get__(self):
            return self.getPC()

    property pc_side:
        def __get__(self):
            return self.getPCSide()
        def __set__(self, value):
            self.setPCSide(value)

    property norm_type:
        def __get__(self):
            return self.getNormType()
        def __set__(self, value):
            self.setNormType(value)

    # --- tolerances ---

    property rtol:
        def __get__(self):
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(rtol=value)

    property atol:
        def __get__(self):
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(atol=value)

    property divtol:
        def __get__(self):
            return self.getTolerances()[2]
        def __set__(self, value):
            self.setTolerances(divtol=value)

    property max_it:
        def __get__(self):
            return self.getTolerances()[3]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    # --- iteration ---

    property its:
        def __get__(self):
            return self.getIterationNumber()
        def __set__(self, value):
            self.setIterationNumber(value)

    property norm:
        def __get__(self):
            return self.getResidualNorm()
        def __set__(self, value):
            self.setResidualNorm(value)

    property history:
        def __get__(self):
            return self.getConvergenceHistory()

    # --- convergence ---

    property reason:
        def __get__(self):
            return self.getConvergedReason()
        def __set__(self, value):
            self.setConvergedReason(value)

    property iterating:
        def __get__(self):
            return self.reason == 0

    property converged:
        def __get__(self):
            return self.reason > 0

    property diverged:
        def __get__(self):
            return self.reason < 0

# --------------------------------------------------------------------

del KSPType
del KSPNormType
del KSPConvergedReason

# --------------------------------------------------------------------
