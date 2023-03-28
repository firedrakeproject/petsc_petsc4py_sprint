# -----------------------------------------------------------------------------

class TSType(object):
    # native
    EULER           = S_(TSEULER)
    BEULER          = S_(TSBEULER)
    BASICSYMPLECTIC = S_(TSBASICSYMPLECTIC)
    PSEUDO          = S_(TSPSEUDO)
    CN              = S_(TSCN)
    SUNDIALS        = S_(TSSUNDIALS)
    RK              = S_(TSRK)
    PYTHON          = S_(TSPYTHON)
    THETA           = S_(TSTHETA)
    ALPHA           = S_(TSALPHA)
    ALPHA2          = S_(TSALPHA2)
    GLLE            = S_(TSGLLE)
    GLEE            = S_(TSGLEE)
    SSP             = S_(TSSSP)
    ARKIMEX         = S_(TSARKIMEX)
    ROSW            = S_(TSROSW)
    EIMEX           = S_(TSEIMEX)
    MIMEX           = S_(TSMIMEX)
    BDF             = S_(TSBDF)
    RADAU5          = S_(TSRADAU5)
    MPRK            = S_(TSMPRK)
    DISCGRAD        = S_(TSDISCGRAD)
    # aliases
    FE = EULER
    BE = BEULER
    TH = THETA
    CRANK_NICOLSON = CN
    RUNGE_KUTTA    = RK

class TSRKType(object):
    RK1FE = S_(TSRK1FE)
    RK2A  = S_(TSRK2A)
    RK2B  = S_(TSRK2B)
    RK4   = S_(TSRK4)
    RK3BS = S_(TSRK3BS)
    RK3   = S_(TSRK3)
    RK5F  = S_(TSRK5F)
    RK5DP = S_(TSRK5DP)
    RK5BS = S_(TSRK5BS)
    RK6VR = S_(TSRK6VR)
    RK7VR = S_(TSRK7VR)
    RK8VR = S_(TSRK8VR)

class TSARKIMEXType(object):
    ARKIMEX1BEE   = S_(TSARKIMEX1BEE)
    ARKIMEXA2     = S_(TSARKIMEXA2)
    ARKIMEXL2     = S_(TSARKIMEXL2)
    ARKIMEXARS122 = S_(TSARKIMEXARS122)
    ARKIMEX2C     = S_(TSARKIMEX2C)
    ARKIMEX2D     = S_(TSARKIMEX2D)
    ARKIMEX2E     = S_(TSARKIMEX2E)
    ARKIMEXPRSSP2 = S_(TSARKIMEXPRSSP2)
    ARKIMEX3      = S_(TSARKIMEX3)
    ARKIMEXBPR3   = S_(TSARKIMEXBPR3)
    ARKIMEXARS443 = S_(TSARKIMEXARS443)
    ARKIMEX4      = S_(TSARKIMEX4)
    ARKIMEX5      = S_(TSARKIMEX5)

class TSProblemType(object):
    LINEAR    = TS_LINEAR
    NONLINEAR = TS_NONLINEAR

class TSEquationType(object):
    UNSPECIFIED               = TS_EQ_UNSPECIFIED
    EXPLICIT                  = TS_EQ_EXPLICIT
    ODE_EXPLICIT              = TS_EQ_ODE_EXPLICIT
    DAE_SEMI_EXPLICIT_INDEX1  = TS_EQ_DAE_SEMI_EXPLICIT_INDEX1
    DAE_SEMI_EXPLICIT_INDEX2  = TS_EQ_DAE_SEMI_EXPLICIT_INDEX2
    DAE_SEMI_EXPLICIT_INDEX3  = TS_EQ_DAE_SEMI_EXPLICIT_INDEX3
    DAE_SEMI_EXPLICIT_INDEXHI = TS_EQ_DAE_SEMI_EXPLICIT_INDEXHI
    IMPLICIT                  = TS_EQ_IMPLICIT
    ODE_IMPLICIT              = TS_EQ_ODE_IMPLICIT
    DAE_IMPLICIT_INDEX1       = TS_EQ_DAE_IMPLICIT_INDEX1
    DAE_IMPLICIT_INDEX2       = TS_EQ_DAE_IMPLICIT_INDEX2
    DAE_IMPLICIT_INDEX3       = TS_EQ_DAE_IMPLICIT_INDEX3
    DAE_IMPLICIT_INDEXHI      = TS_EQ_DAE_IMPLICIT_INDEXHI

class TSExactFinalTime(object):
    UNSPECIFIED = TS_EXACTFINALTIME_UNSPECIFIED
    STEPOVER    = TS_EXACTFINALTIME_STEPOVER
    INTERPOLATE = TS_EXACTFINALTIME_INTERPOLATE
    MATCHSTEP   = TS_EXACTFINALTIME_MATCHSTEP

class TSConvergedReason(object):
    # iterating
    CONVERGED_ITERATING      = TS_CONVERGED_ITERATING
    ITERATING                = TS_CONVERGED_ITERATING
    # converged
    CONVERGED_TIME           = TS_CONVERGED_TIME
    CONVERGED_ITS            = TS_CONVERGED_ITS
    CONVERGED_USER           = TS_CONVERGED_USER
    CONVERGED_EVENT          = TS_CONVERGED_EVENT
    # diverged
    DIVERGED_NONLINEAR_SOLVE = TS_DIVERGED_NONLINEAR_SOLVE
    DIVERGED_STEP_REJECTED   = TS_DIVERGED_STEP_REJECTED

# -----------------------------------------------------------------------------

cdef class TS(Object):
    """An abstract PETSc object that manages all time-steppers (ODE integrators).

    See Also
    --------
    `TS`
    """
    Type = TSType
    RKType = TSRKType
    ARKIMEXType = TSARKIMEXType
    ProblemType = TSProblemType
    EquationType = TSEquationType
    ExactFinalTime = TSExactFinalTime
    ExactFinalTimeOption = TSExactFinalTime
    ConvergedReason = TSConvergedReason

    # --- xxx ---

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ts
        self.ts = NULL

    # --- xxx ---

    def view(self, Viewer viewer: Viewer=None) -> None:
        """Print the `TS` object.

        Collective.

        Parameters
        ----------
        viewer
            The visualization context.

        Notes
        -----
        ``-ts_view`` calls TSView at the end of TSStep

        See Also
        --------
        TSView
        """
        cdef PetscViewer cviewer = NULL
        if viewer is not None: cviewer = viewer.vwr
        CHKERR( TSView(self.ts, cviewer) )

    def load(self, Viewer viewer: Viewer) -> None:
        """Load a `TS` that has been stored in binary with `view`.

        Parameters
        ----------
        viewer
            The visualization context.

        See Also
        --------
        TSLoad
        """
        CHKERR( TSLoad(self.ts, viewer.vwr) )

    def destroy(self) -> Self:
        """Destroy the `TS` that was created with `create`.

        See Also
        --------
        TSDestroy
        """
        CHKERR( TSDestroy(&self.ts) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create an empty `TS`. 
        
        The problem type can then be set with `setProblemType` and the type of
        solver can then be set with `setType`.

        Parameters
        ----------
        comm
            The communicator or None for `Sys.getDefaultComm`.

        See Also
        --------
        TSCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscTS newts = NULL
        CHKERR( TSCreate(ccomm, &newts) )
        PetscCLEAR(self.obj); self.ts = newts
        return self

    def clone(self) -> TS:
        """Clone a `TS` object.

        Collective. Creates a shallow copy of the `TS`.

        Returns
        -------
        ts
            The cloned `TS` object.

        See Also
        --------
        TSClone
        """
        cdef TS ts = TS()
        CHKERR( TSClone(self.ts, &ts.ts) )
        return ts

    def setType(self, ts_type: TS.Type | str) -> None:
        """Set the method to be used as the `TS` solver.

        Parameters
        ----------
        ts_type : TS.Type
            The solver type.

        Notes
        -----
        ``-ts_type`` sets the method from the commandline

        See Also
        --------
        TSSetType
        """
        cdef PetscTSType cval = NULL
        ts_type = str2bytes(ts_type, &cval)
        CHKERR( TSSetType(self.ts, cval) )

    def setRKType(self, ts_type: TS.RKType | str) -> None:
        """Set the type of the `TSRK` scheme.

        Parameters
        ----------
        ts_type
            the type of `TSRK` scheme.

        Notes
        -----
            ``-ts_rk_type`` sets `TSRK` scheme type from the commandline.

        See Also
        --------
        TSRKSetType
        """
        cdef PetscTSRKType cval = NULL
        ts_type = str2bytes(ts_type, &cval)
        CHKERR( TSRKSetType(self.ts, cval) )

    def setARKIMEXType(self, ts_type: TS.ARKIMEXType | str) -> None:
        """Set the type of `TSARKIMEX` scheme.

        Parameters
        ----------
        ts_type
            The type of `TSARKIMEX` scheme.

        Notes
        -----
            ``-ts_arkimex_type`` sets TSARKIMEX scheme type from the commandline.

        See Also
        --------
        TSARKIMEXSetType
        """
        cdef PetscTSARKIMEXType cval = NULL
        ts_type = str2bytes(ts_type, &cval)
        CHKERR( TSARKIMEXSetType(self.ts, cval) )

    def setARKIMEXFullyImplicit(self, flag: bool) -> None:
        """Solve both parts of the equation implicitly.

        Parameters
        ----------
        flag : bool
            Set to True for fully implicit.

        See Also
        --------
        TSARKIMEXSetFullyImplicit
        """
        cdef PetscBool bval = asBool(flag)
        CHKERR( TSARKIMEXSetFullyImplicit(self.ts, bval) )

    def getType(self) -> str:
        """Return the `TS` type.
        
        See Also
        --------
        TSGetType
        """
        cdef PetscTSType cval = NULL
        CHKERR( TSGetType(self.ts, &cval) )
        return bytes2str(cval)

    def getRKType(self) -> str:
        """Return the `TSRK` scheme.
        
        See Also
        --------
        TSRKGetType
        """
        cdef PetscTSRKType cval = NULL
        CHKERR( TSRKGetType(self.ts, &cval) )
        return bytes2str(cval)

    def getARKIMEXType(self) -> str:
        """Return the `TSARKIMEX` scheme.
        
        See Also
        --------
        TSARKIMEXGetType
        """
        cdef PetscTSARKIMEXType cval = NULL
        CHKERR( TSARKIMEXGetType(self.ts, &cval) )
        return bytes2str(cval)

    def setProblemType(self, ptype: TS.ProblemType) -> None:
        """Set the type of problem to be solved.
        
        Parameters
        ----------
        ptype
            The type of problem of the forms.

        See Also
        --------
        TSSetProblemType
        """
        CHKERR( TSSetProblemType(self.ts, ptype) )

    def getProblemType(self) -> TS.ProblemType:
        """Return the type of problem to be solved.

        See Also
        --------
        TSGetProblemType
        """
        cdef PetscTSProblemType ptype = TS_NONLINEAR
        CHKERR( TSGetProblemType(self.ts, &ptype) )
        return ptype

    def setEquationType(self, eqtype: TS.EquationType) -> None:
        """Set the type of the equation that `TS` is solving.

        Not collective.

        Parameters
        ----------
        eqtype
            The type of equation.
        
        See Also
        --------
        TSSetEquationType
        """
        CHKERR( TSSetEquationType(self.ts, eqtype) )

    def getEquationType(self) -> TS.EquationType:
        """Get the type of the equation that `TS` is solving.

        See Also
        --------
        TSGetEquationType
        """
        cdef PetscTSEquationType eqtype = TS_EQ_UNSPECIFIED
        CHKERR( TSGetEquationType(self.ts, &eqtype) )
        return eqtype

    def setOptionsPrefix(self, prefix : str) -> None:
        """Set the prefix used for all the `TS` options.

        Parameters
        ----------
        prefix
            The prefix to prepend to all option names.

        Notes
        -----
        A hyphen must not be given at the beginning of the prefix name.

        See Also
        --------
        TSSetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( TSSetOptionsPrefix(self.ts, cval) )

    def getOptionsPrefix(self) -> str:
        """Return the prefix used for all the `TS` options.

        See Also
        --------
        TSGetOptionsPrefix
        """
        cdef const char *cval = NULL
        CHKERR( TSGetOptionsPrefix(self.ts, &cval) )
        return bytes2str(cval)

    def appendOptionsPrefix(self, prefix: str) -> None:
        """Append to the prefix used for all the `TS` options.

        Parameters
        ----------
        prefix
            The prefix to append to the current prefix.

        Notes
        -----
        A hyphen must not be given at the beginning of the prefix name.

        See Also
        --------
        TSAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( TSAppendOptionsPrefix(self.ts, cval) )

    def setFromOptions(self) -> None:
        """Set various `TS` parameters from user options.

        See Also
        --------
        TSSetFromOptions, ``petsc_options``
        """
        CHKERR( TSSetFromOptions(self.ts) )

    # --- application context ---

    def setAppCtx(self, appctx: Any) -> None:
        """Set the application context.

        Parameters
        ----------
        appctx
            The application context.
        """
        self.set_attr('__appctx__', appctx)

    def getAppCtx(self) -> Any:
        """Return the application context."""
        return self.get_attr('__appctx__')

    # --- user RHS Function/Jacobian routines ---

    def setRHSFunction(
        self, 
        function: TSRHSFunction, 
        Vec f=None,
        args : tuple[Any, ...] | None = None,
        kargs : dict[str, Any] | None = None) -> None:
        """Set the routine for evaluating the function ``G`` in ``U_t = G(t,u)``.

        Parameters
        ----------
        function
            The right-hand-side function.
        f
            The vector into which the right-hand-side is computed.
        args
            Additional posititional arguments for ``function``.
        kargs
            Additional keyword arguments for ``function``.

        See Also
        --------
        TSSetRHSFunction
        """
        cdef PetscVec fvec=NULL
        if f is not None: fvec = f.vec
        if function is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (function, args, kargs)
            self.set_attr('__rhsfunction__', context)
            CHKERR( TSSetRHSFunction(self.ts, fvec, TS_RHSFunction, <void*>context) )
        else:
            CHKERR( TSSetRHSFunction(self.ts, fvec, NULL, NULL) )

    def setRHSJacobian(
        self, 
        jacobian: TSRHSJacobian,
        Mat J=None,
        Mat P=None,
        args : tuple[Any, ...] | None = None,
        kargs : dict[str, Any] | None = None) -> None:
        """Set the function to compute the Jacobian of ``G`` in ``U_t = G(U,t)``.

        Logically collective.

        Parameters
        ----------
        jacobian
            The right-hand-side function.
        J
            The matrix into which the jacobian is computed.
        P
            The matrix into which the preconditioner is computed.
        args
            Additional posititional arguments for ``jacobian``.
        kargs
            Additional keyword arguments for ``jacobian``.

        See Also
        --------
        TSSetRHSJacobian
        """
        cdef PetscMat Jmat=NULL
        if J is not None: Jmat = J.mat
        cdef PetscMat Pmat=Jmat
        if P is not None: Pmat = P.mat
        if jacobian is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (jacobian, args, kargs)
            self.set_attr('__rhsjacobian__', context)
            CHKERR( TSSetRHSJacobian(self.ts, Jmat, Pmat, TS_RHSJacobian, <void*>context) )
        else:
            CHKERR( TSSetRHSJacobian(self.ts, Jmat, Pmat, NULL, NULL) )

    def computeRHSFunction(self, t: float, Vec x, Vec f) -> None:
        """Evaluate the right-hand-side function.

        Parameters
        ----------
        t
            The time at which to evaluate the RHS.
        x
            The state vector.
        f
            The Vec into which the RHS is computed.
        
        See Also
        --------
        TSComputeRHSFunction
        """
        cdef PetscReal time = asReal(t)
        CHKERR( TSComputeRHSFunction(self.ts, time, x.vec, f.vec) )

    def computeRHSFunctionLinear(self, t, Vec x, Vec f) -> None:
        """Evaluate the right-hand-side via the user-provided Jacobian.

        Parameters
        ----------
        t
            The time at which to evaluate the RHS.
        x
            The state vector.
        f
            The Vec into which the RHS is computed.
    
        See Also
        --------
        TSComputeRHSFunctionLinear
        """
        cdef PetscReal time = asReal(t)
        CHKERR( TSComputeRHSFunctionLinear(self.ts, time, x.vec, f.vec, NULL) )

    def computeRHSJacobian(self, t, Vec x, Mat J, Mat P=None) -> None:
        """Compute the Jacobian matrix that has been set with `setRHSJacobian`.

        Collective.

        Parameters
        ----------
        t
            The time at which to evaluate the Jacobian.
        x
            The state vector.
        J
            The matrix into which the Jacobian is computed.
        P
            The optional matrix to use for building a preconditioner matrix.

        See Also
        --------
        TSComputeRHSJacobian
        """
        cdef PetscReal time = asReal(t)
        cdef PetscMat jmat = J.mat, pmat = J.mat
        if P is not None: pmat = P.mat
        CHKERR( TSComputeRHSJacobian(self.ts, time, x.vec, jmat, pmat) )

    def computeRHSJacobianConstant(self, t, Vec x, Mat J, Mat P=None) -> None:
        """Reuse a Jacobian that is time-independent.

        Collective.

        Parameters
        ----------
        t
            The time at which to evaluate the Jacobian.
        x
            The state vector.
        J
            A pointer to the stored Jacobian.
        P
            An optional pointer to the preconditioner matrix.

        See Also
        --------
        TSComputeRHSJacobianConstant
        """
        cdef PetscReal time = asReal(t)
        cdef PetscMat jmat = J.mat, pmat = J.mat
        if P is not None: pmat = P.mat
        CHKERR( TSComputeRHSJacobianConstant(self.ts, time, x.vec, jmat, pmat, NULL) )

    def getRHSFunction(self) -> tuple[Vec, TSRHSFunction]:
        """Return the vector where the right hand side is stored and the
        function used to compute it.

        Not collective.

        See Also
        --------
        TSGetRHSFunction
        """
        cdef Vec f = Vec()
        CHKERR( TSGetRHSFunction(self.ts, &f.vec, NULL, NULL) )
        PetscINCREF(f.obj)
        cdef object function = self.get_attr('__rhsfunction__')
        return (f, function)

    def getRHSJacobian(self) -> tuple[Mat, Mat, TSRHSJacobian]:
        """Return the Jacobian and the function used to compute them.

        Not collective, but parallel objects are returned if `TS` is parallel.

        See Also
        --------
        TSGetRHSJacobian
        """
        cdef Mat J = Mat(), P = Mat()
        CHKERR( TSGetRHSJacobian(self.ts, &J.mat, &P.mat, NULL, NULL) )
        PetscINCREF(J.obj); PetscINCREF(P.obj)
        cdef object jacobian = self.get_attr('__rhsjacobian__')
        return (J, P, jacobian)

    # --- user Implicit Function/Jacobian routines ---

    def setIFunction(
        self, 
        function: TSIFunction, 
        Vec f=None, 
        args : tuple[Any, ...] | None = None,
        kargs : dict[str, Any] | None = None) -> None:
        """Set the function representing the DAE to be solved.

        Logically collective.

        Parameters
        ----------
        function
            The right-hand-side function.
        f
            The vector to store values or ``None`` to be created internally.
        args
            Additional posititional arguments for ``function``.
        kargs
            Additional keyword arguments for ``function``.

        See Also
        --------
        TSSetIFunction
        """
        cdef PetscVec fvec=NULL
        if f is not None: fvec = f.vec
        if function is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (function, args, kargs)
            self.set_attr('__ifunction__', context)
            CHKERR( TSSetIFunction(self.ts, fvec, TS_IFunction, <void*>context) )
        else:
            CHKERR( TSSetIFunction(self.ts, fvec, NULL, NULL) )

    def setIJacobian(
        self,
        jacobian: TSIJacobian,
        Mat J=None,
        Mat P=None, 
        args : tuple[Any, ...] | None = None,
        kargs : dict[str, Any] | None = None) -> None:
        """Set the function to compute the Jacobian.
        
        Set the function to compute the matrix ``dF/dU + a*dF/dU_t`` where
        ``F(t,U,U_t)`` is the function provided with `setIFunction`.

        Logically collective.

        Parameters
        ----------
        jacobian
            The function which computes the Jacobian.
        J
            The matrix into which the Jacobian is computed.
        P
            The optional matrix to use for building a preconditioner matrix.
        args
            Additional posititional arguments for ``jacobian``.
        kargs
            Additional keyword arguments for ``jacobian``.

        See Also
        --------
        TSSetIJacobian
        """
        cdef PetscMat Jmat=NULL
        if J is not None: Jmat = J.mat
        cdef PetscMat Pmat=Jmat
        if P is not None: Pmat = P.mat
        if jacobian is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (jacobian, args, kargs)
            self.set_attr('__ijacobian__', context)
            CHKERR( TSSetIJacobian(self.ts, Jmat, Pmat, TS_IJacobian, <void*>context) )
        else:
            CHKERR( TSSetIJacobian(self.ts, Jmat, Pmat, NULL, NULL) )

    def setIJacobianP(
        self,
        jacobian,
        Mat J=None, 
        args : tuple[Any, ...] | None = None,
        kargs : dict[str, Any] | None = None) -> None:
        """Set the function that computes the Jacobian of ``F`` with respect to
        the parameters ``P`` where ``F(Udot,U,t) = G(U,P,t)``, as well as the
        location to store the matrix.

        Logically collective.

        Parameters
        ----------
        jacobian
            The function which computes the Jacobian.
        J
            The matrix into which the Jacobian is computed.
        args
            Additional posititional arguments for ``jacobian``.
        kargs
            Additional keyword arguments for ``jacobian``.
        
        See Also
        --------
        TSSetIJacobianP
        """
        cdef PetscMat Jmat=NULL
        if J is not None: Jmat = J.mat
        if jacobian is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (jacobian, args, kargs)
            self.set_attr('__ijacobianp__', context)
            CHKERR( TSSetIJacobianP(self.ts, Jmat, TS_IJacobianP, <void*>context) )
        else:
            CHKERR( TSSetIJacobianP(self.ts, Jmat, NULL, NULL) )

    def computeIFunction(self,
                         t, Vec x, Vec xdot,
                         Vec f, imex: bool=False) -> None:
        """Evaluate the DAE residual written in implicit form.

        Collective.

        Parameters
        ----------
        t
            The current time.
        x
            The state vector.
        xdot
            The time derivative of the state vector.
        f
            The vector into which the residual is stored.
        imex
            A flag which indicates if the RHS should be kept separate.

        See Also
        --------
        TSComputeIFunction
        """
        cdef PetscReal rval = asReal(t)
        cdef PetscBool bval = imex
        CHKERR( TSComputeIFunction(self.ts, rval, x.vec, xdot.vec,
                                   f.vec, bval) )

    def computeIJacobian(self,
                         t, Vec x, Vec xdot, a: float,
                         Mat J, Mat P=None, imex=False) -> None:
        """Evaluate the Jacobian of the DAE.

        Collective. If ``F(t,U,Udot)=0`` is the DAE, the required Jacobian is
        ``dF/dU + shift*dF/dUdot``

        Parameters
        ----------
        t
            The current time.
        x
            The state vector.
        xdot
            The time derivative of the state vector.
        a
            The shift to apply
        J
            The matrix into which the Jacobian is computed.
        P
            The optional matrix to use for building a preconditioner matrix.
        imex
            A flag which indicates if the RHS should be kept separate.

        See Also
        --------
        TSComputeIJacobian
        """
        cdef PetscReal rval1 = asReal(t)
        cdef PetscReal rval2 = asReal(a)
        cdef PetscBool bval  = imex
        cdef PetscMat jmat = J.mat, pmat = J.mat
        if P is not None: pmat = P.mat
        CHKERR( TSComputeIJacobian(self.ts, rval1, x.vec, xdot.vec, rval2,
                                   jmat, pmat, bval) )

    def computeIJacobianP(self,
                         t, Vec x, Vec xdot, a,
                         Mat J, imex=False) -> None:
        """Evaluate the Jacobian with respect to parameters.

        Collective.

        Parameters
        ----------
        t
            The current time.
        x
            The state vector.
        xdot
            The time derivative of the state vector.
        a
            The shift to apply
        J
            The matrix into which the Jacobian is computed.
        imex
            A flag which indicates if the RHS should be kept separate.

        See Also
        --------
        TSComputeIJacobianP
        """
        cdef PetscReal rval1 = asReal(t)
        cdef PetscReal rval2 = asReal(a)
        cdef PetscBool bval  = asBool(imex)
        cdef PetscMat jmat = J.mat
        CHKERR( TSComputeIJacobianP(self.ts, rval1, x.vec, xdot.vec, rval2,
                                   jmat, bval) )

    def getIFunction(self) -> tuple[Vec, TSIFunction]:
        """Return the vector and function which computes the implicit residual.

        Not collective.

        See Also
        --------
        TSGetIFunction
        """
        cdef Vec f = Vec()
        CHKERR( TSGetIFunction(self.ts, &f.vec, NULL, NULL) )
        PetscINCREF(f.obj)
        cdef object function = self.get_attr('__ifunction__')
        return (f, function)

    def getIJacobian(self) -> tuple[Mat, Mat, TSIJacobian]:
        """Return the matrices and function which computes the implicit Jacobian.

        Not collective.

        See Also
        --------
        TSGetIJacobian
        """
        cdef Mat J = Mat(), P = Mat()
        CHKERR( TSGetIJacobian(self.ts, &J.mat, &P.mat, NULL, NULL) )
        PetscINCREF(J.obj); PetscINCREF(P.obj)
        cdef object jacobian = self.get_attr('__ijacobian__')
        return (J, P, jacobian)

    def setI2Function(
        self,
        function: TSI2Function,
        Vec f=None,
        args : tuple[Any, ...] | None = None,
        kargs : dict[str, Any] | None = None) -> None:
        """Set the function to compute the 2nd order DAE.

        Logically collective.

        Parameters
        ----------
        function
            The right-hand-side function.
        f
            The vector to store values or ``None`` to be created internally.
        args
            Additional posititional arguments for ``function``.
        kargs
            Additional keyword arguments for ``function``.


        See Also
        --------
        TSSetI2Function
        """
        cdef PetscVec fvec=NULL
        if f is not None: fvec = f.vec
        if function is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (function, args, kargs)
            self.set_attr('__i2function__', context)
            CHKERR( TSSetI2Function(self.ts, fvec, TS_I2Function, <void*>context) )
        else:
            CHKERR( TSSetI2Function(self.ts, fvec, NULL, NULL) )

    def setI2Jacobian(
        self,
        jacobian: TSI2Jacobian,
        Mat J=None, 
        Mat P=None, 
        args=None,
        kargs=None) -> None:
        """Set the function to compute the Jacobian of the 2nd order DAE.
        
        Logically collective.

        Parameters
        ----------
        jacobian
            The function which computes the Jacobian.
        J
            The matrix into which the Jacobian is computed.
        P
            The optional matrix to use for building a preconditioner matrix.
        args
            Additional posititional arguments for ``jacobian``.
        kargs
            Additional keyword arguments for ``jacobian``.

        See Also
        --------
        TSSetI2Jacobian
        """
        cdef PetscMat Jmat=NULL
        if J is not None: Jmat = J.mat
        cdef PetscMat Pmat=Jmat
        if P is not None: Pmat = P.mat
        if jacobian is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (jacobian, args, kargs)
            self.set_attr('__i2jacobian__', context)
            CHKERR( TSSetI2Jacobian(self.ts, Jmat, Pmat, TS_I2Jacobian, <void*>context) )
        else:
            CHKERR( TSSetI2Jacobian(self.ts, Jmat, Pmat, NULL, NULL) )

    def computeI2Function(self, t: float, Vec x, Vec xdot, Vec xdotdot, Vec f) -> None:
        """Evaluate the DAE residual in implicit form.

        Collective.

        Parameters
        ----------
        t
            The current time.
        x
            The state vector.
        xdot
            The time derivative of the state vector.
        xdotxdot
            The second time derivative of the state vector.
        f
            The vector into which the residual is stored.

        See Also
        --------
        TSComputeI2Function
        """
        cdef PetscReal rval = asReal(t)
        CHKERR( TSComputeI2Function(self.ts, rval, x.vec, xdot.vec, xdotdot.vec,
                                   f.vec) )

    def computeI2Jacobian(
        self,
        t: float,
        Vec x, 
        Vec xdot, 
        Vec xdotdot, 
        v: float, 
        a: float, 
        Mat J, 
        Mat P=None) -> None:
        """Evaluate the Jacobian of the DAE.

        If ``F(t,U,V,A)=0`` is the DAE, the required Jacobian is ``dF/dU + v dF/dV + a dF/dA``.

        Collective.

        Parameters
        ----------
        t
            The current time.
        x
            The state vector.
        xdot
            The time derivative of the state vector.
        xdotxdot
            The second time derivative of the state vector.
        v
            The shift to apply to the first derivative.
        a
            The shift to apply to the second derivative.
        J
            The matrix into which the Jacobian is computed.
        P
            The optional matrix to use for building a preconditioner matrix.

        See Also
        --------
        TSComputeI2Jacobian
        """
        cdef PetscReal rval1 = asReal(t)
        cdef PetscReal rval2 = asReal(v)
        cdef PetscReal rval3 = asReal(a)
        cdef PetscMat jmat = J.mat, pmat = J.mat
        if P is not None: pmat = P.mat
        CHKERR( TSComputeI2Jacobian(self.ts, rval1, x.vec, xdot.vec, xdotdot.vec, rval2, rval3,
                                   jmat, pmat) )

    def getI2Function(self) -> tuple[Vec, TSI2Function]:
        """Return the vector and function which computes the residual.

        Not collective.

        See Also
        --------
        TSGetI2Function
        """
        cdef Vec f = Vec()
        CHKERR( TSGetI2Function(self.ts, &f.vec, NULL, NULL) )
        PetscINCREF(f.obj)
        cdef object function = self.get_attr('__i2function__')
        return (f, function)

    def getI2Jacobian(self) -> tuple[Mat, Mat, TSI2Jacobian]:
        """Return the matrices and function which computes the Jacobian.

        Not collective.

        See Also
        --------
        TSGetI2Jacobian
        """
        cdef Mat J = Mat(), P = Mat()
        CHKERR( TSGetI2Jacobian(self.ts, &J.mat, &P.mat, NULL, NULL) )
        PetscINCREF(J.obj); PetscINCREF(P.obj)
        cdef object jacobian = self.get_attr('__i2jacobian__')
        return (J, P, jacobian)

    # --- solution vector ---

    def setSolution(self, Vec u) -> None:
        """Set the initial solution vector.

        Logically collective.

        Parameters
        ----------
        u
            The solution vector.

        See Also
        --------
        TSSetSolution
        """
        CHKERR( TSSetSolution(self.ts, u.vec) )

    def getSolution(self) -> Vec:
        """Return the solution at the present timestep.
        
        Not collective, but the vector is parallel if the `TS` is parallel. It
        is valid to call this routine inside the function that you are
        evaluating in order to move to the new timestep. This vector is not
        changed until the solution at the next timestep has been calculated.

        See Also
        --------
        TSGetSolution
        """
        cdef Vec u = Vec()
        CHKERR( TSGetSolution(self.ts, &u.vec) )
        PetscINCREF(u.obj)
        return u

    def setSolution2(self, Vec u, Vec v) -> None:
        """Set the initial solution and its time derivative.

        Logically collective.

        Parameters
        ----------
        u
            the solution vector
        v
            the time derivative vector

        See Also
        --------
        TS2SetSolution
        """
        CHKERR( TS2SetSolution(self.ts, u.vec, v.vec) )

    def getSolution2(self) -> tuple[Vec, Vec]:
        """Return the solution and time derivative at the present timestep.
        
        Not collective, but vectors are parallel if `TS` is parallel. It is
        valid to call this routine inside the function that you are evaluating
        in order to move to the new timestep. These vectors are not changed
        until the solution at the next timestep has been calculated.

        See Also
        --------
        TS2GetSolution
        """
        cdef Vec u = Vec()
        cdef Vec v = Vec()
        CHKERR( TS2GetSolution(self.ts, &u.vec, &v.vec) )
        PetscINCREF(u.obj)
        PetscINCREF(v.obj)
        return (u, v)

    # --- time span ---

    def setTimeSpan(self, tspan: Sequence[float]) -> None:
        """Set the time span. 
        
        Collective. The solution will be computed and stored for each time
        requested in the span. The times must be all increasing and correspond
        to the intermediate points for time integration.
        `TS_EXACTFINALTIME_MATCHSTEP` must be used to make the last time step in
        each sub-interval match the intermediate points specified. The
        intermediate solutions are saved in a vector array that can be accessed
        with `getTimeSpanSolutions`. 

        Parameters
        ----------
        tspan
            the sequence of time points
        
        Notes
        -----
        ``-ts_time_span <t0,...tf>`` sets the time span from the commandline
        
        See Also
        --------
        TSSetTimeSpan
        """
        cdef PetscInt  nt = 0
        cdef PetscReal *rtspan = NULL
        cdef object tmp = oarray_r(tspan, &nt, &rtspan)
        CHKERR( TSSetTimeSpan(self.ts, nt, rtspan) )

    def getTimeSpan(self) -> NDArray[float]:
        """
        See Also
        --------
        TSGetTimeSpan
        """
        cdef const PetscReal *rtspan = NULL
        cdef PetscInt   nt = 0
        CHKERR( TSGetTimeSpan(self.ts, &nt, &rtspan) )
        cdef object tspan = array_r(nt, rtspan)
        return tspan

    def getTimeSpanSolutions(self) -> list[Vec]:
        """
        See Also
        --------
        TSGetTimeSpanSolutions
        """
        cdef PetscInt nt = 0
        cdef PetscVec *sols = NULL
        CHKERR( TSGetTimeSpanSolutions(self.ts, &nt, &sols) )
        cdef object sollist = None
        if sols != NULL:
            sollist = [ref_Vec(sols[i]) for i from 0 <= i < nt]
        return sollist

    # --- inner solver ---

    def getSNES(self) -> SNES:
        """
        See Also
        --------
        TSGetSNES
        """
        cdef SNES snes = SNES()
        CHKERR( TSGetSNES(self.ts, &snes.snes) )
        PetscINCREF(snes.obj)
        return snes

    def getKSP(self) -> KSP:
        """
        See Also
        --------
        TSGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( TSGetKSP(self.ts, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    # --- discretization space ---

    def getDM(self) -> DM:
        """
        See Also
        --------
        TSGetDM
        """
        cdef PetscDM newdm = NULL
        CHKERR( TSGetDM(self.ts, &newdm) )
        cdef DM dm = subtype_DM(newdm)()
        dm.dm = newdm
        PetscINCREF(dm.obj)
        return dm

    def setDM(self, DM dm) -> None:
        """
        See Also
        --------
        TSSetDM
        """
        CHKERR( TSSetDM(self.ts, dm.dm) )

    # --- customization ---

    def setTime(self, t) -> None:
        """
        See Also
        --------
        TSSetTime
        """
        cdef PetscReal rval = asReal(t)
        CHKERR( TSSetTime(self.ts, rval) )

    def getTime(self) -> float:
        """
        See Also
        --------
        TSGetTime
        """
        cdef PetscReal rval = 0
        CHKERR( TSGetTime(self.ts, &rval) )
        return toReal(rval)

    def getPrevTime(self) -> float:
        """
        See Also
        --------
        TSGetPrevTime
        """
        cdef PetscReal rval = 0
        CHKERR( TSGetPrevTime(self.ts, &rval) )
        return toReal(rval)

    def getSolveTime(self) -> float:
        """
        See Also
        --------
        TSGetSolveTime
        """
        cdef PetscReal rval = 0
        CHKERR( TSGetSolveTime(self.ts, &rval) )
        return toReal(rval)

    def setTimeStep(self, time_step) -> None:
        """
        See Also
        --------
        TSSetTimeStep
        """
        cdef PetscReal rval = asReal(time_step)
        CHKERR( TSSetTimeStep(self.ts, rval) )

    def getTimeStep(self) -> float:
        """
        See Also
        --------
        TSGetTimeStep
        """
        cdef PetscReal tstep = 0
        CHKERR( TSGetTimeStep(self.ts, &tstep) )
        return toReal(tstep)

    def setStepNumber(self, step_number) -> None:
        """
        See Also
        --------
        TSSetStepNumber
        """
        cdef PetscInt ival = asInt(step_number)
        CHKERR( TSSetStepNumber(self.ts, ival) )

    def getStepNumber(self) -> int:
        """
        See Also
        --------
        TSGetStepNumber
        """
        cdef PetscInt ival = 0
        CHKERR( TSGetStepNumber(self.ts, &ival) )
        return toInt(ival)

    def setMaxTime(self, max_time) -> None:
        """
        See Also
        --------
        TSSetMaxTime
        """
        cdef PetscReal rval = asReal(max_time)
        CHKERR( TSSetMaxTime(self.ts, rval) )

    def getMaxTime(self) -> float:
        """
        See Also
        --------
        TSGetMaxTime
        """
        cdef PetscReal rval = 0
        CHKERR( TSGetMaxTime(self.ts, &rval) )
        return toReal(rval)

    def setMaxSteps(self, max_steps) -> None:
        """
        See Also
        --------
        TSSetMaxSteps
        """
        cdef PetscInt  ival = asInt(max_steps)
        CHKERR( TSSetMaxSteps(self.ts, ival) )

    def getMaxSteps(self) -> int:
        """
        See Also
        --------
        TSGetMaxSteps
        """
        cdef PetscInt ival = 0
        CHKERR( TSGetMaxSteps(self.ts, &ival) )
        return toInt(ival)

    def getSNESIterations(self) -> int:
        """
        See Also
        --------
        TSGetSNESIterations
        """
        cdef PetscInt n = 0
        CHKERR( TSGetSNESIterations(self.ts, &n) )
        return toInt(n)

    def getKSPIterations(self) -> int:
        """
        See Also
        --------
        TSGetKSPIterations
        """
        cdef PetscInt n = 0
        CHKERR( TSGetKSPIterations(self.ts, &n) )
        return toInt(n)

    def setMaxStepRejections(self, n) -> None:
        """
        See Also
        --------
        TSSetMaxStepRejections
        """
        cdef PetscInt rej = asInt(n)
        CHKERR( TSSetMaxStepRejections(self.ts, rej))

    #def getMaxStepRejections(self):
    #    cdef PetscInt n = 0
    #    CHKERR( TSGetMaxStepRejections(self.ts, &n))
    #    return toInt(n)

    def getStepRejections(self) -> int:
        """
        See Also
        --------
        TSGetStepRejections
        """
        cdef PetscInt n = 0
        CHKERR( TSGetStepRejections(self.ts, &n) )
        return toInt(n)

    def setMaxSNESFailures(self, n) -> None:
        """
        See Also
        --------
        TSSetMaxSNESFailures
        """
        cdef PetscInt fails = asInt(n)
        CHKERR( TSSetMaxSNESFailures(self.ts, fails))

    #def getMaxSNESFailures(self, n):
    #    cdef PetscInt n = 0
    #    CHKERR( TSGetMaxSNESFailures(self.ts, &n))
    #    return toInt(n)

    def getSNESFailures(self) -> int:
        """
        See Also
        --------
        TSGetSNESFailures
        """
        cdef PetscInt n = 0
        CHKERR( TSGetSNESFailures(self.ts, &n) )
        return toInt(n)

    def setErrorIfStepFails(self, flag=True) -> None:
        """
        See Also
        --------
        TSSetErrorIfStepFails
        """
        cdef PetscBool bval = flag
        CHKERR( TSSetErrorIfStepFails(self.ts, bval))

    def setTolerances(self, rtol=None, atol=None) -> None:
        """
        See Also
        --------
        TSSetTolerances
        """
        cdef PetscReal rrtol = PETSC_DEFAULT
        cdef PetscReal ratol = PETSC_DEFAULT
        cdef PetscVec  vrtol = NULL
        cdef PetscVec  vatol = NULL
        if rtol is None:
            pass
        elif isinstance(rtol, Vec):
            vrtol = (<Vec>rtol).vec
        else:
            rrtol = asReal(rtol)
        if atol is None:
            pass
        elif isinstance(atol, Vec):
            vatol = (<Vec>atol).vec
        else:
            ratol = asReal(atol)
        CHKERR( TSSetTolerances(self.ts, ratol, vatol, rrtol, vrtol) )

    def getTolerances(self) ->tuple[float,float]:
        """
        See Also
        --------
        TSGetTolerances
        """
        cdef PetscReal rrtol = PETSC_DEFAULT
        cdef PetscReal ratol = PETSC_DEFAULT
        cdef PetscVec  vrtol = NULL
        cdef PetscVec  vatol = NULL
        CHKERR( TSGetTolerances(self.ts, &ratol, &vatol, &rrtol, &vrtol) )
        cdef object rtol = None
        if vrtol != NULL:
            rtol = ref_Vec(vrtol)
        else:
            rtol = toReal(rrtol)
        cdef object atol = None
        if vatol != NULL:
            atol = ref_Vec(vatol)
        else:
            atol = toReal(ratol)
        return (rtol, atol)

    def setExactFinalTime(self, option) -> None:
        """
        See Also
        --------
        TSSetExactFinalTime
        """
        cdef PetscTSExactFinalTimeOption oval = option
        CHKERR( TSSetExactFinalTime(self.ts, oval) )

    def setConvergedReason(self, reason) -> None:
        """
        See Also
        --------
        TSSetConvergedReason
        """
        cdef PetscTSConvergedReason cval = reason
        CHKERR( TSSetConvergedReason(self.ts, cval) )

    def getConvergedReason(self) -> TSConvergedReason:
        """
        See Also
        --------
        TSGetConvergedReason
        """
        cdef PetscTSConvergedReason reason = TS_CONVERGED_ITERATING
        CHKERR( TSGetConvergedReason(self.ts, &reason) )
        return reason

    # --- monitoring ---

    def setMonitor(self, monitor, args=None, kargs=None) -> None:
        """
        See Also
        --------
        TSMonitorSet
        """
        if monitor is None: return
        cdef object monitorlist = self.get_attr('__monitor__')
        if monitorlist is None:
            monitorlist = []
            self.set_attr('__monitor__', monitorlist)
            CHKERR( TSMonitorSet(self.ts, TS_Monitor, NULL, NULL) )
        if args  is None: args  = ()
        if kargs is None: kargs = {}
        context = (monitor, args, kargs)
        monitorlist.append(context)

    def getMonitor(self) -> TSMonitorFunction:
        """
        """
        return self.get_attr('__monitor__')

    def monitorCancel(self) -> None:
        """
        See Also
        --------
        TSMonitorCancel
        """
        self.set_attr('__monitor__', None)
        CHKERR( TSMonitorCancel(self.ts) )

    cancelMonitor = monitorCancel

    def monitor(self, step, time, Vec u=None) -> None:
        """
        See Also
        --------
        TSMonitor
        """
        cdef PetscInt  ival = asInt(step)
        cdef PetscReal rval = asReal(time)
        cdef PetscVec  uvec = NULL
        if u is not None: uvec = u.vec
        if uvec == NULL:
            CHKERR( TSGetSolution(self.ts, &uvec) )
        CHKERR( TSMonitor(self.ts, ival, rval, uvec) )

    # --- event handling ---

    def setEventHandler(self, direction, terminate, eventhandler, postevent=None, args=None, kargs=None) -> None:
        """
        See Also
        --------
        TSSetEventHandler
        """
        cdef PetscInt  ndirs = 0
        cdef PetscInt *idirs = NULL
        direction = iarray_i(direction, &ndirs, &idirs)

        cdef PetscInt   nterm = 0
        cdef PetscBool *iterm = NULL
        terminate = iarray_b(terminate, &nterm, &iterm)
        assert nterm == ndirs

        cdef PetscInt nevents = ndirs
        if eventhandler is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            self.set_attr('__eventhandler__', (eventhandler, args, kargs))
            if postevent is not None:
                self.set_attr('__postevent__', (postevent, args, kargs))
                CHKERR( TSSetEventHandler(self.ts, nevents, idirs, iterm, TS_EventHandler, TS_PostEvent, <void*>NULL) )
            else:
                self.set_attr('__postevent__', None)
                CHKERR( TSSetEventHandler(self.ts, nevents, idirs, iterm, TS_EventHandler, NULL, <void*>NULL) )
        else:
            CHKERR( TSSetEventHandler(self.ts, nevents, idirs, iterm, NULL, NULL, <void*>NULL) )

    def setEventTolerances(self, tol=None, vtol=None) -> None:
        """
        See Also
        --------
        TSSetEventTolerances
        """
        cdef PetscInt  nevents = 0
        cdef PetscReal tolr = PETSC_DEFAULT
        cdef PetscInt  ntolr = 0
        cdef PetscReal *vtolr = NULL
        if tol is not None:
            tolr = asReal(tol)
        if vtol is not None:
            CHKERR( TSGetNumEvents(self.ts, &nevents) )
            vtol = iarray_r(vtol, &ntolr,  &vtolr)
            assert ntolr == nevents
        CHKERR( TSSetEventTolerances(self.ts, tolr, vtolr) )

    def getNumEvents(self) -> int:
        """
        See Also
        --------
        TSGetNumEvents
        """
        cdef PetscInt nevents = 0
        CHKERR( TSGetNumEvents(self.ts, &nevents) )
        return toInt(nevents)

    # --- solving ---

    def setPreStep(self, prestep, args=None, kargs=None) -> None:
        """
        See Also
        --------
        TSSetPreStep
        """
        if prestep is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (prestep, args, kargs)
            self.set_attr('__prestep__', context)
            CHKERR( TSSetPreStep(self.ts, TS_PreStep) )
        else:
            self.set_attr('__prestep__', None)
            CHKERR( TSSetPreStep(self.ts, NULL) )

    def getPreStep(self) -> TSPreStepFunction:
        """
        """
        return self.get_attr('__prestep__')

    def setPostStep(self, poststep, args=None, kargs=None) -> None:
        """
        See Also
        --------
        TSSetPostStep
        """
        if poststep is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (poststep, args, kargs)
            self.set_attr('__poststep__', context)
            CHKERR( TSSetPostStep(self.ts, TS_PostStep) )
        else:
            self.set_attr('__poststep__', None)
            CHKERR( TSSetPostStep(self.ts, NULL) )

    def getPostStep(self) -> TSPostStepFunction:
        """
        """
        return self.get_attr('__poststep__')

    def setUp(self) -> None:
        """
        See Also
        --------
        TSSetUp
        """
        CHKERR( TSSetUp(self.ts) )

    def reset(self) -> None:
        """
        See Also
        --------
        TSReset
        """
        CHKERR( TSReset(self.ts) )

    def step(self) -> None:
        """
        See Also
        --------
        TSStep
        """
        CHKERR( TSStep(self.ts) )

    def restartStep(self) -> None:
        """
        See Also
        --------
        TSRestartStep
        """
        CHKERR( TSRestartStep(self.ts) )

    def rollBack(self) -> None:
        """
        See Also
        --------
        TSRollBack
        """
        CHKERR( TSRollBack(self.ts) )

    def solve(self, Vec u) -> None:
        """
        See Also
        --------
        TSSolve
        """
        CHKERR( TSSolve(self.ts, u.vec) )

    def interpolate(self, t, Vec u) -> None:
        """
        See Also
        --------
        TSInterpolate
        """
        cdef PetscReal rval = asReal(t)
        CHKERR( TSInterpolate(self.ts, rval, u.vec) )

    def setStepLimits(self, hmin, hmax) -> None:
        """
        See Also
        --------
        TSAdaptSetStepLimits
        """
        cdef PetscTSAdapt tsadapt = NULL
        cdef PetscReal hminr = toReal(hmin)
        cdef PetscReal hmaxr = toReal(hmax)
        TSGetAdapt(self.ts, &tsadapt)
        CHKERR( TSAdaptSetStepLimits(tsadapt, hminr, hmaxr) )

    def getStepLimits(self) -> tuple[float,float]:
        """
        See Also
        --------
        TSAdaptGetStepLimits
        """
        cdef PetscTSAdapt tsadapt = NULL
        cdef PetscReal hminr = 0.
        cdef PetscReal hmaxr = 0.
        TSGetAdapt(self.ts, &tsadapt)
        CHKERR( TSAdaptGetStepLimits(tsadapt, &hminr, &hmaxr) )
        return (asReal(hminr), asReal(hmaxr))

    # --- Adjoint methods ---

    def setSaveTrajectory(self) -> None:
        """
        See Also
        --------
        TSSetSaveTrajectory
        """
        CHKERR(TSSetSaveTrajectory(self.ts))

    def removeTrajectory(self) -> None:
        """
        See Also
        --------
        TSRemoveTrajectory
        """
        CHKERR(TSRemoveTrajectory(self.ts))

    def getCostIntegral(self) -> Vec:
        """
        See Also
        --------
        TSGetCostIntegral
        """
        cdef Vec cost = Vec()
        CHKERR( TSGetCostIntegral(self.ts, &cost.vec) )
        PetscINCREF(cost.obj)
        return cost

    def setCostGradients(self, vl, vm=None) -> None:
        """
        See Also
        --------
        TSSetCostGradients
        """
        cdef PetscInt n = 0;
        cdef PetscVec *vecl = NULL
        cdef PetscVec *vecm = NULL
        cdef mem1 = None, mem2 = None
        if isinstance(vl, Vec): vl = [vl]
        if isinstance(vm, Vec): vm = [vm]
        if vl is not None:
            n = <PetscInt>len(vl)
        elif vm is not None:
            n = <PetscInt>len(vm)
        if vl is not None:
            assert len(vl) == <Py_ssize_t>n
            mem1 = oarray_p(empty_p(n), NULL, <void**>&vecl)
            for i from 0 <= i < n:
                vecl[i] = (<Vec?>vl[i]).vec
        if vm is not None:
            assert len(vm) == <Py_ssize_t>n
            mem2 = oarray_p(empty_p(n), NULL, <void**>&vecm)
            for i from 0 <= i < n:
                vecm[i] = (<Vec?>vm[i]).vec
        self.set_attr('__costgradients_memory', (mem1, mem2))
        CHKERR( TSSetCostGradients(self.ts, n, vecl, vecm) )

    def getCostGradients(self) -> tuple[list[Vec],list[Vec]]:
        """
        See Also
        --------
        TSGetCostGradients
        """
        cdef PetscInt i = 0, n = 0
        cdef PetscVec *vecl = NULL
        cdef PetscVec *vecm = NULL
        CHKERR( TSGetCostGradients(self.ts, &n, &vecl, &vecm) )
        cdef object vl = None, vm = None
        if vecl != NULL:
            vl = [ref_Vec(vecl[i]) for i from 0 <= i < n]
        if vecm != NULL:
            vm = [ref_Vec(vecm[i]) for i from 0 <= i < n]
        return (vl, vm)

    def setRHSJacobianP(self, jacobianp, Mat A=None, args=None, kargs=None) -> None:
        """
        See Also
        --------
        TSSetRHSJacobianP
        """
        cdef PetscMat Amat=NULL
        if A is not None: Amat = A.mat
        if jacobianp is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (jacobianp, args, kargs)
            self.set_attr('__rhsjacobianp__', context)
            CHKERR( TSSetRHSJacobianP(self.ts, Amat, TS_RHSJacobianP, <void*>context) )
        else:
            CHKERR( TSSetRHSJacobianP(self.ts, Amat, NULL, NULL) )

    def createQuadratureTS(self, forward=True) -> TS:
        """
        See Also
        --------
        TSCreateQuadratureTS
        """
        cdef TS qts = TS()
        cdef PetscBool fwd = forward
        CHKERR( TSCreateQuadratureTS(self.ts, fwd, &qts.ts) )
        PetscINCREF(qts.obj)
        return qts

    def getQuadratureTS(self) -> tuple[bool, TS]:
        """
        See Also
        --------
        TSGetQuadratureTS
        """
        cdef TS qts = TS()
        cdef PetscBool fwd = PETSC_FALSE
        CHKERR( TSGetQuadratureTS(self.ts, &fwd, &qts.ts) )
        PetscINCREF(qts.obj)
        return (toBool(fwd), qts)

    def setRHSJacobianP(self, rhsjacobianp, Mat A=None, args=None, kargs=None) -> None:
        """
        See Also
        --------
        TSSetRHSJacobianP
        """
        cdef PetscMat Amat=NULL
        if A is not None: Amat = A.mat
        if rhsjacobianp is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (rhsjacobianp, args, kargs)
            self.set_attr('__rhsjacobianp__', context)
            CHKERR( TSSetRHSJacobianP(self.ts, Amat, TS_RHSJacobianP, <void*>context) )
        else:
            CHKERR( TSSetRHSJacobianP(self.ts, Amat, NULL, NULL) )

    def computeRHSJacobianP(self, t, Vec x, Mat J) -> None:
        """
        See Also
        --------
        TSComputeRHSJacobianP
        """
        cdef PetscReal rval = asReal(t)
        CHKERR( TSComputeRHSJacobianP(self.ts, rval, x.vec, J.mat) )

    def adjointSetSteps(self, adjoint_steps) -> None:
        """
        See Also
        --------
        TSAdjointSetSteps
        """
        cdef PetscInt ival = asInt(adjoint_steps)
        CHKERR( TSAdjointSetSteps(self.ts, ival) )

    def adjointSetUp(self) -> None:
        """
        See Also
        --------
        TSAdjointSetUp
        """
        CHKERR(TSAdjointSetUp(self.ts))

    def adjointSolve(self) -> None:
        """
        See Also
        --------
        TSAdjointSolve
        """
        CHKERR( TSAdjointSolve(self.ts) )

    def adjointStep(self) -> None:
        """
        See Also
        --------
        TSAdjointStep
        """
        CHKERR(TSAdjointStep(self.ts))

    def adjointReset(self) -> None:
        """
        See Also
        --------
        TSAdjointReset
        """
        CHKERR(TSAdjointReset(self.ts))

    # --- Python ---

    def createPython(self, context=None, comm=None) -> Self:
        """
        See Also
        --------
        TSCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscTS newts = NULL
        CHKERR( TSCreate(ccomm, &newts) )
        PetscCLEAR(self.obj); self.ts = newts
        CHKERR( TSSetType(self.ts, TSPYTHON) )
        CHKERR( TSPythonSetContext(self.ts, <void*>context) )
        return self

    def setPythonContext(self, context) -> None:
        """
        See Also
        --------
        TSPythonSetContext
        """
        CHKERR( TSPythonSetContext(self.ts, <void*>context) )

    def getPythonContext(self) -> Any:
        """
        See Also
        --------
        TSPythonGetContext
        """
        cdef void *context = NULL
        CHKERR( TSPythonGetContext(self.ts, &context) )
        if context == NULL: return None
        else: return <object> context

    def setPythonType(self, py_type) -> None:
        """
        See Also
        --------
        TSPythonSetType
        """
        cdef const char *cval = NULL
        py_type = str2bytes(py_type, &cval)
        CHKERR( TSPythonSetType(self.ts, cval) )

    def getPythonType(self) -> str:
        """
        See Also
        --------
        TSPythonGetType
        """
        cdef const char *cval = NULL
        CHKERR( TSPythonGetType(self.ts, &cval) )
        return bytes2str(cval)

    # --- Theta ---

    def setTheta(self, theta) -> None:
        """
        See Also
        --------
        TSThetaSetTheta
        """
        cdef PetscReal rval = asReal(theta)
        CHKERR( TSThetaSetTheta(self.ts, rval) )

    def getTheta(self) -> float:
        """
        See Also
        --------
        TSThetaGetTheta
        """
        cdef PetscReal rval = 0
        CHKERR( TSThetaGetTheta(self.ts, &rval) )
        return toReal(rval)

    def setThetaEndpoint(self, flag=True) -> None:
        """
        See Also
        --------
        TSThetaSetEndpoint
        """
        cdef PetscBool bval = flag
        CHKERR( TSThetaSetEndpoint(self.ts, bval) )

    def getThetaEndpoint(self) -> bool:
        """
        See Also
        --------
        TSThetaGetEndpoint
        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( TSThetaGetEndpoint(self.ts, &flag) )
        return toBool(flag)

    # --- Alpha ---

    def setAlphaRadius(self, radius) -> None:
        """
        See Also
        --------
        TSAlphaSetRadius
        """
        cdef PetscReal rval = asReal(radius)
        CHKERR( TSAlphaSetRadius(self.ts, rval) )

    def setAlphaParams(self, alpha_m=None,alpha_f=None, gamma=None) -> None:
        """
        See Also
        --------
        TSAlphaSetParams
        """
        cdef PetscReal rval1 = 0, rval2 = 0, rval3 = 0
        try: CHKERR( TSAlphaGetParams(self.ts, &rval1, &rval2, &rval3) )
        except PetscError: pass
        if alpha_m is not None: rval1 = asReal(alpha_m)
        if alpha_f is not None: rval2 = asReal(alpha_f)
        if gamma   is not None: rval3 = asReal(gamma)
        CHKERR( TSAlphaSetParams(self.ts,  rval1,  rval2,  rval3) )

    def getAlphaParams(self) -> tuple[float, float, float]:
        """
        See Also
        --------
        TSAlphaGetParams
        """
        cdef PetscReal rval1 = 0, rval2 = 0, rval3 = 0
        CHKERR( TSAlphaGetParams(self.ts, &rval1, &rval2, &rval3) )
        return (toReal(rval1), toReal(rval2), toReal(rval3))

    # --- application context ---

    property appctx:
        """
        """
        def __get__(self) -> Any:
            return self.getAppCtx()
        def __set__(self, value) -> None:
            self.setAppCtx(value)

    # --- discretization space ---

    property dm:
        """
        """
        def __get__(self) -> DM:
            return self.getDM()
        def __set__(self, value) -> None:
            self.setDM(value)

    # --- xxx ---

    property problem_type:
        """
        """
        def __get__(self) -> TS.ProblemType:
            return self.getProblemType()
        def __set__(self, value) -> None:
            self.setProblemType(value)

    property equation_type:
        """
        """
        def __get__(self) -> TS.EquationType:
            return self.getEquationType()
        def __set__(self, value) -> None:
            self.setEquationType(value)

    property snes:
        """
        """
        def __get__(self) -> SNES:
            return self.getSNES()

    property ksp:
        def __get__(self) -> KSP:
            return self.getKSP()

    property vec_sol:
        def __get__(self) -> Vec:
            return self.getSolution()

    # --- xxx ---

    property time:
        def __get__(self) -> float:
            return self.getTime()
        def __set__(self, value) -> None:
            self.setTime(value)

    property time_step:
        def __get__(self) -> None:
            return self.getTimeStep()
        def __set__(self, value):
            self.setTimeStep(value)

    property step_number:
        def __get__(self) -> int:
            return self.getStepNumber()
        def __set__(self, value) -> None:
            self.setStepNumber(value)

    property max_time:
        def __get__(self) -> float:
            return self.getMaxTime()
        def __set__(self, value) -> None:
            self.setMaxTime(value)

    property max_steps:
        def __get__(self) -> int:
            return self.getMaxSteps()
        def __set__(self, value) -> None:
            self.setMaxSteps(value)

    # --- convergence ---

    property rtol:
        def __get__(self) -> float:
            return self.getTolerances()[0]
        def __set__(self, value) -> None:
            self.setTolerances(rtol=value)

    property atol:
        def __get__(self) -> float:
            return self.getTolerances()[1]
        def __set__(self, value) -> None:
            self.setTolerances(atol=value)

    property reason:
        def __get__(self) -> TSConvergedReason:
            return self.getConvergedReason()
        def __set__(self, value) -> None:
            self.setConvergedReason(value)

    property iterating:
        def __get__(self) -> bool:
            return self.reason == 0

    property converged:
        def __get__(self) -> bool:
            return self.reason > 0

    property diverged:
        def __get__(self) -> bool:
            return self.reason < 0

# -----------------------------------------------------------------------------

del TSType
del TSRKType
del TSARKIMEXType
del TSProblemType
del TSEquationType
del TSExactFinalTime
del TSConvergedReason

# -----------------------------------------------------------------------------
