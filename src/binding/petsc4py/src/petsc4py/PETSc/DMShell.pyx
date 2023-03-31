cdef class DMShell(DM):
"""TODO:"""

    def create(self, comm: Comm | None = None) -> Self:
        """TODO:

        Collective.

        Parameters
        ----------
        comm
            The MPI communicator.

        See also
        --------
        petsc.DMShellCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM newdm = NULL
        CHKERR( DMShellCreate(ccomm, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def setMatrix(self, Mat mat) -> None:
        """Set a template matrix.

        Collective.

        Parameters
        ----------
        mat
            The template matrix.

        See also
        --------
        petsc.DMShellSetMatrix

        """
        CHKERR( DMShellSetMatrix(self.dm, mat.mat) )

    def setGlobalVector(self, Vec gv) -> None:
        """Set a template global vector.

        Logically collective.

        Parameters
        ----------
        gv
            Template vector.

        See also
        --------
        petsc.DMShellSetGlobalVector, setLocalVector

        """
        CHKERR( DMShellSetGlobalVector(self.dm, gv.vec) )

    def setLocalVector(self, Vec lv) -> None:
        """Set a template local vector.

        Logically collective.

        Parameters
        ----------
        lv
            Template vector.

        See also
        --------
        petsc.DMShellSetLocalVector, setGlobalVector

        """
        CHKERR( DMShellSetLocalVector(self.dm, lv.vec) )

    def setCreateGlobalVector(
        self,
        create_gvec,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine to create a global vector.

        Logically collective.

        Parameters
        ----------
        create_gvec
            The creation routine.
        args
            Additional positional arguments for ``create_gvec``.
        kargs
            Additional keyword arguments for ``create_gvec``.

        See also
        --------
        petsc.DMShellSetCreateGlobalVector, setCreateLocalVector

        """
        if create_gvec is not None:
            if args  is None: args = ()
            if kargs is None: kargs = {}
            context = (create_gvec, args, kargs)
            self.set_attr('__create_global_vector__', context)
            CHKERR( DMShellSetCreateGlobalVector(self.dm, DMSHELL_CreateGlobalVector) )
        else:
            CHKERR( DMShellSetCreateGlobalVector(self.dm, NULL) )

    def setCreateLocalVector(
        self,
        create_lvec,
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:# TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:# TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:# TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        # TODO:
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine to create a local vector.

        Logically collective.

        Parameters
        ----------
        create_lvec
            The creation routine.
        args
            Additional positional arguments for ``create_lvec``.
        kargs
            Additional keyword arguments for ``create_lvec``.

        See also
        --------
        petsc.DMShellSetCreateLocalVector, setCreateGlobalVector

        """
        if create_lvec is not None:
            if args  is None: args = ()
            if kargs is None: kargs = {}
            context = (create_lvec, args, kargs)
            self.set_attr('__create_local_vector__', context)
            CHKERR( DMShellSetCreateLocalVector(self.dm, DMSHELL_CreateLocalVector) )
        else:
            CHKERR( DMShellSetCreateLocalVector(self.dm, NULL) )

    def setGlobalToLocal(
        self,
        begin,
        end,
        begin_args: tuple[Any, ...] | None = None,
        begin_kargs: dict[str, Any] | None = None,
        end_args: tuple[Any, ...] | None = None,
        end_kargs: dict[str, Any] | None = None,
    ):
        """Set the routines used to perform a global to local scatter.

        If these functions are not provided but DMShellSetGlobalToLocalVecScatter()
        is called then DMGlobalToLocalBeginDefaultShell() / DMGlobalToLocalEndDefaultShell()
        are used to to perform the transfers.

        Logically collective on the DM.

        Parameters
        ----------
        dm
            The `DMShell`.
        begin
            The routine that begins the global to local scatter.
        end
            The routine that ends the global to local scatter.
        begin_args
            Additional positional arguments for ``begin``.
        begin_kargs
            Additional keyword arguments for ``begin``.
        end_args
            Additional positional arguments for ``end``.
        end_kargs
            Additional keyword arguments for ``end``.

        See also
        --------
        petsc.DMShellSetGlobalToLocal

        """
        cdef PetscDMShellXToYFunction cbegin = NULL, cend = NULL
        if begin is not None:
            if begin_args  is None: begin_args = ()
            if begin_kargs is None: begin_kargs = {}
            context = (begin, begin_args, begin_kargs)
            self.set_attr('__g2l_begin__', context)
            cbegin = &DMSHELL_GlobalToLocalBegin
        if end is not None:
            if end_args  is None: end_args = ()
            if end_kargs is None: end_kargs = {}
            context = (end, end_args, end_kargs)
            self.set_attr('__g2l_end__', context)
            cend = &DMSHELL_GlobalToLocalEnd
        CHKERR( DMShellSetGlobalToLocal(self.dm, cbegin, cend) )

    def setGlobalToLocalVecScatter(self, Scatter gtol):
        """Sets a VecScatter context for global to local communication.

        Logically collective on the DM.

        Parameters
        ----------
        gtol
            The global to local VecScatter context.

        See also
        --------
        petsc.DMShellSetGlobalToLocalVecScatter

        """
        CHKERR( DMShellSetGlobalToLocalVecScatter(self.dm, gtol.sct) )

    def setLocalToGlobal(
        self,
        begin,
        end,
        begin_args: tuple[Any, ...] | None = None,
        begin_kargs: dict[str, Any] | None = None,
        end_args: tuple[Any, ...] | None = None,
        end_kargs: dict[str, Any] | None = None,
    ):
        """Set the routines used to perform a local to global scatter.

        If these functions are not provided but DMShellSetLocalToGlobalVecScatter()
        is called then DMLocalToGlobalBeginDefaultShell() / DMLocalToGlobalEndDefaultShell()
        are used to to perform the transfers.

        Logically collective on the DM.

        Parameters
        ----------
        begin
            The routine that begins the local to global scatter.
        end
            The routine that ends the local to global scatter.
        begin_args
            Additional positional arguments for ``begin``.
        begin_kargs
            Additional keyword arguments for ``begin``.
        end_args
            Additional positional arguments for ``end``.
        end_kargs
            Additional keyword arguments for ``end``.

        See also
        --------
        petsc.DMShellSetLocalToGlobal

        """
        cdef PetscDMShellXToYFunction cbegin = NULL, cend = NULL
        if begin is not None:
            if begin_args  is None: begin_args = ()
            if begin_kargs is None: begin_kargs = {}
            context = (begin, begin_args, begin_kargs)
            self.set_attr('__l2g_begin__', context)
            cbegin = &DMSHELL_LocalToGlobalBegin
        if end is not None:
            if end_args  is None: end_args = ()
            if end_kargs is None: end_kargs = {}
            context = (end, end_args, end_kargs)
            self.set_attr('__l2g_end__', context)
            cend = &DMSHELL_LocalToGlobalEnd
        CHKERR( DMShellSetLocalToGlobal(self.dm, cbegin, cend) )

    def setLocalToGlobalVecScatter(self, Scatter ltog):
        """Set a VecScatter context for local to global communication.

        Logically collective on the DM.

        Parameters
        ----------
        ltog
            The local to global VecScatter context.

        See also
        --------
        petsc.DMShellSetLocalToGlobalVecScatter

        """
        CHKERR( DMShellSetLocalToGlobalVecScatter(self.dm, ltog.sct) )

    def setLocalToLocal(
        self,
        begin,
        end,
        begin_args=None,
        begin_kargs=None,
        end_args=None,
        end_kargs=None,
        # args: tuple[Any, ...] | None = None,
        # kargs: dict[str, Any] | None = None,
    ):
        """Set the routines used to perform a local to local scatter.

        Logically collective on the DM.

        Note
        If these functions are not provided but `setLocalToLocalVecScatter` is
        called then DMLocalToLocalBeginDefaultShell / DMLocalToLocalEndDefaultShell
        are used to to perform the transfers.

        Parameters
        ----------
        begin
            The routine that begins the local to local scatter.
        end
            The routine that ends the local to local scatter.
        begin_args
            Additional positional arguments for ``begin``.
        begin_kargs
            Additional keyword arguments for ``begin``.
        end_args
            Additional positional arguments for ``end``.
        end_kargs
            Additional keyword arguments for ``end``.

        See also
        --------
        petsc.DMShellSetLocalToLocal

        """
        cdef PetscDMShellXToYFunction cbegin = NULL, cend = NULL
        cbegin = NULL
        cend = NULL
        if begin is not None:
            if begin_args  is None: begin_args = ()
            if begin_kargs is None: begin_kargs = {}
            context = (begin, begin_args, begin_kargs)
            self.set_attr('__l2l_begin__', context)
            cbegin = &DMSHELL_LocalToLocalBegin
        if end is not None:
            if end_args  is None: end_args = ()
            if end_kargs is None: end_kargs = {}
            context = (end, end_args, end_kargs)
            self.set_attr('__l2l_end__', context)
            cend = &DMSHELL_LocalToLocalEnd
        CHKERR( DMShellSetLocalToLocal(self.dm, cbegin, cend) )

    def setLocalToLocalVecScatter(self, Scatter ltol):
        """Set a VecScatter context for local to local communication.

        Logically collective.

        Parameters
        ----------
        ltol
            The local to local VecScatter context.

        See also
        --------
        petsc.DMShellSetLocalToLocalVecScatter

        """
        CHKERR( DMShellSetLocalToLocalVecScatter(self.dm, ltol.sct) )

    def setCreateMatrix(
        self,
        create_matrix,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ):
        """Set the routine to create a matrix.

        Logically collective.

        Parameters
        ----------
        create_matrix
            The function to create a matrix.
        args
            Additional positional arguments for ``create_matrix``.
        kargs
            Additional keyword arguments for ``create_matrix``.

        See also
        --------
        petsc.DMShellSetCreateMatrix

        """
        if create_matrix is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (create_matrix, args, kargs)
            self.set_attr('__create_matrix__', context)
            CHKERR( DMShellSetCreateMatrix(self.dm, DMSHELL_CreateMatrix) )
        else:
            CHKERR( DMShellSetCreateMatrix(self.dm, NULL) )

    def setCoarsen(
        self,
        coarsen,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ):
        """Set the routine used to coarsen the `DMShell`.

        Logically collective.

        Parameters
        ----------
        coarsen
            The routine that coarsens the DM.
        args
            Additional positional arguments for ``coarsen``.
        kargs
            Additional keyword arguments for ``coarsen``.

        See also
        --------
        petsc.DMShellSetCoarsen, setRefine

        """
        if coarsen is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (coarsen, args, kargs)
            self.set_attr('__coarsen__', context)
            CHKERR( DMShellSetCoarsen(self.dm, DMSHELL_Coarsen) )
        else:
            CHKERR( DMShellSetCoarsen(self.dm, NULL) )

    def setRefine(
        self,
        refine,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ):
        """Set the routine used to refine the `DMShell`.

        Logically collective.

        Parameters
        ----------
        refine
            The routine that refines the DM.
        args
            Additional positional arguments for ``refine``.
        kargs
            Additional keyword arguments for ``refine``.

        See also
        --------
        petsc.DMShellSetRefine, setCoarsen

        """
        if refine is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (refine, args, kargs)
            self.set_attr('__refine__', context)
            CHKERR( DMShellSetRefine(self.dm, DMSHELL_Refine) )
        else:
            CHKERR( DMShellSetRefine(self.dm, NULL) )

    def setCreateInterpolation(
        self,
        create_interpolation,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ):
        """Set the routine used to create the interpolation operator.

        Logically collective.

        Parameters
        ----------
        create_interpolation
            The routine to create the interpolation.
        args
            Additional positional arguments for ``create_interpolation``.
        kargs
            Additional keyword arguments for ``create_interpolation``.

        See also
        --------
        petsc.DMShellSetCreateInterpolation

        """
        if create_interpolation is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (create_interpolation, args, kargs)
            self.set_attr('__create_interpolation__', context)
            CHKERR( DMShellSetCreateInterpolation(self.dm, DMSHELL_CreateInterpolation) )
        else:
            CHKERR( DMShellSetCreateInterpolation(self.dm, NULL) )

    def setCreateInjection(
        self,
        create_injection,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ):
        """Set the routine used to create the injection operator.

        Logically collective.

        Parameters
        ----------
        create_injection
            The routine to create the injection.
        args
            Additional positional arguments for ``create_injection``.
        kargs
            Additional keyword arguments for ``create_injection``.

        See also
        --------
        petsc.DMShellSetCreateInjection

        """
        if create_injection is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (create_injection, args, kargs)
            self.set_attr('__create_injection__', context)
            CHKERR( DMShellSetCreateInjection(self.dm, DMSHELL_CreateInjection) )
        else:
            CHKERR( DMShellSetCreateInjection(self.dm, NULL) )

    def setCreateRestriction(
        self,
        create_restriction,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine used to create the restriction operator.

        Logically collective.

        Parameters
        ----------
        create_restriction TODO: was it striction in docs
            The routine to create the restriction
        args
            Additional positional arguments for ``create_restriction``.
        kargs
            Additional keyword arguments for ``create_restriction``.

        See also
        --------
        petsc.DMShellSetCreateRestriction

        """
        if create_restriction is not None:
            if args  is None: args  = ()
            if kargs is None: kargs = {}
            context = (create_restriction, args, kargs)
            self.set_attr('__create_restriction__', context)
            CHKERR( DMShellSetCreateRestriction(self.dm, DMSHELL_CreateRestriction) )
        else:
            CHKERR( DMShellSetCreateRestriction(self.dm, NULL) )

    def setCreateFieldDecomposition(
        self,
        decomp,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine used to create a decomposition of fields for the `DMShell`.

        Logically collective.

        Parameters
        ----------
        decomp
            The routine to create the decomposition
        args
            Additional positional arguments for ``decomp``.
        kargs
            Additional keyword arguments for ``decomp``.

        See also
        --------
        petsc.DMShellSetCreateFieldDecomposition

        """
        if decomp is not None:
            if args  is None: args = ()
            if kargs is None: kargs = {}
            context = (decomp, args, kargs)
            self.set_attr('__create_field_decomp__', context)
            CHKERR( DMShellSetCreateFieldDecomposition(self.dm, DMSHELL_CreateFieldDecomposition) )
        else:
            CHKERR( DMShellSetCreateFieldDecomposition(self.dm, NULL) )

    def setCreateDomainDecomposition(
        self,
        decomp,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine used to create a domain decomposition for the `DMShell`.

        Logically collective.

        Parameters
        ----------
        decomp
            The routine to create the decomposition.
        args
            Additional positional arguments for ``decomp``.
        kargs
            Additional keyword arguments for ``decomp``.

        See also
        --------
        petsc.DMShellSetCreateDomainDecomposition

        """
        if decomp is not None:
            if args  is None: args = ()
            if kargs is None: kargs = {}
            context = (decomp, args, kargs)
            self.set_attr('__create_domain_decomp__', context)
            CHKERR( DMShellSetCreateDomainDecomposition(self.dm, DMSHELL_CreateDomainDecomposition) )
        else:
            CHKERR( DMShellSetCreateDomainDecomposition(self.dm, NULL) )

    def setCreateDomainDecompositionScatters(
        self,
        scatter,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine used to create the scatter contexts for domain decomposition with a `DMShell`.

        Logically collective.

        Parameters
        ----------
        scatter
            The routine to create the scatters.
        args
            Additional positional arguments for ``scatter``.
        kargs
            Additional keyword arguments for ``scatter``.

        See also
        --------
        petsc.DMShellSetCreateDomainDecompositionScatters

        """
        if scatter is not None:
            if args  is None: args = ()
            if kargs is None: kargs = {}
            context = (scatter, args, kargs)
            self.set_attr('__create_domain_decomp_scatters__', context)
            CHKERR( DMShellSetCreateDomainDecompositionScatters(self.dm, DMSHELL_CreateDomainDecompositionScatters) )
        else:
            CHKERR( DMShellSetCreateDomainDecompositionScatters(self.dm, NULL) )

    def setCreateSubDM(
        self,
        create_subdm,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """Set the routine used to create a sub DM from the `DMShell`.

        Logically collective.

        Parameters
        ----------
        subdm
            The routine to create the decomposition.
        args
            Additional positional arguments for ``subdm``.
        kargs
            Additional keyword arguments for ``subdm``.

        See also
        --------
        petsc.DMShellSetCreateSubDM

        """
        if create_subdm is not None:
            if args  is None: args = ()
            if kargs is None: kargs = {}
            context = (create_subdm, args, kargs)
            self.set_attr('__create_subdm__', context)
            CHKERR( DMShellSetCreateSubDM(self.dm, DMSHELL_CreateSubDM) )
        else:
            CHKERR( DMShellSetCreateSubDM(self.dm, NULL) )
