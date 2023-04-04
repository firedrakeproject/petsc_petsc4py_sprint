# --------------------------------------------------------------------

class DMDAStencilType(object):
    STAR = DMDA_STENCIL_STAR
    BOX  = DMDA_STENCIL_BOX

class DMDAInterpolationType(object):
    Q0 = DMDA_INTERPOLATION_Q0
    Q1 = DMDA_INTERPOLATION_Q1

class DMDAElementType(object):
    P1 = DMDA_ELEMENT_P1
    Q1 = DMDA_ELEMENT_Q1

# --------------------------------------------------------------------

cdef class DMDA(DM):
    """A DM object that is used to manage data for a structured grid."""

    StencilType       = DMDAStencilType
    InterpolationType = DMDAInterpolationType
    ElementType       = DMDAElementType

    #

    def create(
        self,
        dim: int | None = None,
        dof: int | None = None,
        sizes: tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)] | None = None,
        proc_sizes: tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)] | None = None,
        boundary_type: tuple[()] | tuple[(DM.BoundaryType,)] | tuple[(DM.BoundaryType, DM.BoundaryType)] | tuple[(DM.BoundaryType, DM.BoundaryType, DM.BoundaryType)] | None = None,
        stencil_type: StencilType | None = None,
        stencil_width: int | None = None,
        bint setup: bool | None = True,
        ownership_ranges: tuple[Sequence[int]] | tuple[Sequence[int], Sequence[int]] | tuple[Sequence[int], Sequence[int], Sequence[int]] | None = None,
        comm: Comm | None = None,
    ) -> Self:
        """Create a ``DMDA`` object.

        This routine performs the following steps of the C API:
        - ``petsc.DMDACreate``
        - ``petsc.DMSetDimension``
        - ``petsc.DMDASetDof``
        - ``petsc.DMDASetSizes``
        - ``petsc.DMDASetNumProcs``
        - ``petsc.DMDASetOwnershipRanges``
        - ``petsc.DMDASetBoundaryType``
        - ``petsc.DMDASetStencilType``
        - ``petsc.DMDASetStencilWidth``
        - ``petsc.DMSetUp`` (optionally)

        Collective.

        Parameters
        ----------
        dim
            TODO.
        dof
            The number of degrees of freedom.
        sizes
            TODO.
        proc_sizes
            TODO.
        boundary_type
            TODO.
        stencil_type
            The ghost/halo stencil type.
        stencil_width
            The width of the ghost/halo region.
        setup
            TODO.
        ownership_ranges
            TODO.
        comm
            TODO.

        See Also
        --------
        petsc.DMDACreate, petsc.DMSetDimension, petsc.DMDASetDof,
        petsc.DMDASetSizes, petsc.DMDASetNumProcs,
        petsc.DMDASetOwnershipRanges, petsc.DMDASetBoundaryType,
        petsc.DMDASetStencilType, petsc.DMDASetStencilWidth, petsc.DMSetUp

        """
        #
        cdef object arg = None
        try: arg = tuple(dim)
        except TypeError: pass
        else: dim, sizes = None, arg
        #
        cdef PetscInt ndim = PETSC_DECIDE
        cdef PetscInt ndof = PETSC_DECIDE
        cdef PetscInt M = 1, m = PETSC_DECIDE, *lx = NULL
        cdef PetscInt N = 1, n = PETSC_DECIDE, *ly = NULL
        cdef PetscInt P = 1, p = PETSC_DECIDE, *lz = NULL
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        cdef PetscDMDAStencilType stype  = DMDA_STENCIL_BOX
        cdef PetscInt             swidth = PETSC_DECIDE
        # grid and proc sizes
        cdef object gsizes = sizes
        cdef object psizes = proc_sizes
        cdef PetscInt gdim = PETSC_DECIDE
        cdef PetscInt pdim = PETSC_DECIDE
        if sizes is not None:
            gdim = asDims(gsizes, &M, &N, &P)
        if psizes is not None:
            pdim = asDims(psizes, &m, &n, &p)
        if gdim>=0 and pdim>=0:
            assert gdim == pdim
        # dim and dof
        if dim is not None: ndim = asInt(dim)
        if dof is not None: ndof = asInt(dof)
        if ndim==PETSC_DECIDE: ndim = gdim
        if ndof==PETSC_DECIDE: ndof = 1
        # vertex distribution
        if ownership_ranges is not None:
            ownership_ranges = asOwnershipRanges(ownership_ranges,
                                                 ndim, &m, &n, &p,
                                                 &lx, &ly, &lz)
        # periodicity, stencil type & width
        if boundary_type is not None:
            asBoundary(boundary_type, &btx, &bty, &btz)
        if stencil_type is not None:
            stype = asStencil(stencil_type)
        if stencil_width is not None:
            swidth = asInt(stencil_width)
        if setup and swidth == PETSC_DECIDE: swidth = 0
        # create the DMDA object
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM newda = NULL
        CHKERR( DMDACreateND(ccomm, ndim, ndof,
                             M, N, P, m, n, p, lx, ly, lz,
                             btx, bty, btz, stype, swidth,
                             &newda) )
        if setup and ndim > 0: CHKERR( DMSetUp(newda) )
        PetscCLEAR(self.obj); self.dm = newda
        return self

    def duplicate(
        self,
        dof: int | None = None,
        boundary_type: tuple[()] | tuple[(DM.BoundaryType,)] | tuple[(DM.BoundaryType, DM.BoundaryType)] | tuple[(DM.BoundaryType, DM.BoundaryType, DM.BoundaryType)] | None = None,
        stencil_type: StencilType | None = None,
        stencil_width: int | None = None,
    ) -> DMDA:
        """Duplicate a DMDA.

        This routine retrieves the information from the DMDA and recreates it.
        Parameters ``dof``, ``boundary_type``, ``stencil_type``,
        ``stencil_width`` will be overwritten, if provided.

        Collective.

        Parameters
        ----------
        dof
            The number of degrees of freedom.
        boundary_type
            TODO.
        stencil_type
            The ghost/halo stencil type.
        stencil_width
            The width of the ghost/halo region.

        See Also
        --------
        petsc.DMDAGetInfo, create, petsc.DMSetUp

        """
        cdef PetscInt ndim = 0, ndof = 0
        cdef PetscInt M = 1, N = 1, P = 1
        cdef PetscInt m = 1, n = 1, p = 1
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        cdef PetscDMDAStencilType  stype  = DMDA_STENCIL_BOX
        cdef PetscInt              swidth = PETSC_DECIDE
        CHKERR( DMDAGetInfo(self.dm,
                          &ndim,
                          &M, &N, &P,
                          &m, &n, &p,
                          &ndof, &swidth,
                          &btx, &bty, &btz,
                          &stype) )
        cdef const PetscInt *lx = NULL, *ly = NULL, *lz = NULL
        CHKERR( DMDAGetOwnershipRanges(self.dm, &lx, &ly, &lz) )
        cdef MPI_Comm comm = MPI_COMM_NULL
        CHKERR( PetscObjectGetComm(<PetscObject>self.dm, &comm) )
        #
        if dof is not None:
            ndof = asInt(dof)
        if boundary_type is not None:
            asBoundary(boundary_type, &btx, &bty, &btz)
        if stencil_type  is not None:
            stype = asStencil(stencil_type)
        if stencil_width is not None:
            swidth = asInt(stencil_width)
        #
        cdef DMDA da = DMDA()
        CHKERR( DMDACreateND(comm, ndim, ndof,
                             M, N, P, m, n, p, lx, ly, lz,
                             btx, bty, btz, stype, swidth,
                             &da.dm) )
        CHKERR( DMSetUp(da.dm) )
        return da

    #

    def setDim(self, dim: int) -> None:
        """Set the topological dimension.

        Collective.

        Parameters
        ----------
        dim
            Topological dimension.

        See Also
        --------
        petsc.DMSetDimension

        """
        return self.setDimension(dim)

    def getDim(self) -> int:
        """Return the topological dimension.

        Not collective.

        See Also
        --------
        petsc.DMGetDimension

        """
        return self.getDimension()

    def setDof(self, dof: int) -> None:
        """Set the number of degrees of freedom per vertex

        Not collective.

        Parameters
        ----------
        dof
            Number of degrees of freedom.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc.DMDASetDof

        """
        cdef PetscInt ndof = asInt(dof)
        CHKERR( DMDASetDof(self.dm, ndof) )

    def getDof(self) -> int:
        """Return the number of degrees of freedom per node.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscInt dof = 0
        CHKERR( DMDAGetInfo(self.dm,
                            NULL,
                            NULL, NULL, NULL,
                            NULL, NULL, NULL,
                            &dof, NULL,
                            NULL, NULL, NULL,
                            NULL) )
        return toInt(dof)

    def setSizes(self, sizes: tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)]) -> None:
        """Set the number of grid points in the three dimensional directions.

        Developer Note
        Since the dimension may not yet have been set the code cannot error check for non-positive Y and Z number of grid points

        Logically collective.

        Parameters
        ----------
        M
            the global X size
        N
            the global Y size
        P
            the global Z size

        See Also
        --------
        petsc.DMDASetSizes

        """
        cdef tuple gsizes = tuple(sizes)
        cdef PetscInt gdim = PETSC_DECIDE
        cdef PetscInt M = 1
        cdef PetscInt N = 1
        cdef PetscInt P = 1
        gdim = asDims(gsizes, &M, &N, &P)
        cdef PetscInt dim = PETSC_DECIDE
        CHKERR( DMDAGetDim(self.dm, &dim) )
        if dim == PETSC_DECIDE:
            CHKERR( DMSetDimension(self.dm, gdim) )
        CHKERR( DMDASetSizes(self.dm, M, N, P) )

    def getSizes(self) -> tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)]:
        """Return the global dimension in first/second/third direction.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscInt dim = 0
        cdef PetscInt M = PETSC_DECIDE
        cdef PetscInt N = PETSC_DECIDE
        cdef PetscInt P = PETSC_DECIDE
        CHKERR( DMDAGetInfo(self.dm,
                            &dim,
                            &M, &N, &P,
                            NULL, NULL, NULL,
                            NULL, NULL,
                            NULL, NULL, NULL,
                            NULL) )
        return toDims(dim, M, N, P)

    def setProcSizes(self, proc_sizes: tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)]) -> None:
        """Set the number of processes in each dimension.

        Logically collective.

        Parameters
        ----------
        m
            the number of X procs (or PETSC_DECIDE)
        n
            the number of Y procs (or PETSC_DECIDE)
        p
            the number of Z procs (or PETSC_DECIDE)

        See Also
        --------
        petsc.DMDASetNumProcs

        """
        cdef tuple psizes = tuple(proc_sizes)
        cdef PetscInt pdim = PETSC_DECIDE
        cdef PetscInt m = PETSC_DECIDE
        cdef PetscInt n = PETSC_DECIDE
        cdef PetscInt p = PETSC_DECIDE
        pdim = asDims(psizes, &m, &n, &p)
        cdef PetscInt dim = PETSC_DECIDE
        CHKERR( DMDAGetDim(self.dm, &dim) )
        if dim == PETSC_DECIDE:
            CHKERR( DMSetDimension(self.dm, pdim) )
        CHKERR( DMDASetNumProcs(self.dm, m, n, p) )

    def getProcSizes(self) -> tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)]:
        """Return the number of processes in first/second/third dimensions.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscInt dim = 0
        cdef PetscInt m = PETSC_DECIDE
        cdef PetscInt n = PETSC_DECIDE
        cdef PetscInt p = PETSC_DECIDE
        CHKERR( DMDAGetInfo(self.dm,
                            &dim,
                            NULL, NULL, NULL,
                            &m, &n, &p,
                            NULL, NULL,
                            NULL, NULL, NULL,
                            NULL) )
        return toDims(dim, m, n, p)

    def setBoundaryType(
        self,
        boundary_type: tuple[()] | tuple[(DM.BoundaryType,)] | tuple[(DM.BoundaryType, DM.BoundaryType)] | tuple[(DM.BoundaryType, DM.BoundaryType, DM.BoundaryType)]
    ) -> None:
        """Set the type of ghost nodes on domain boundaries.

        Parameters
        ----------
        bx
            bx,by,bz is one of DM_BOUNDARY_NONE, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC

        Not collective.

        See Also
        --------
        petsc.DMDASetBoundaryType

        """
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        asBoundary(boundary_type, &btx, &bty, &btz)
        CHKERR( DMDASetBoundaryType(self.dm, btx, bty, btz) )

    def getBoundaryType(self) -> tuple[()] | tuple[(DM.BoundaryType,)] | tuple[(DM.BoundaryType, DM.BoundaryType)] | tuple[(DM.BoundaryType, DM.BoundaryType, DM.BoundaryType)]:
        """Return the type of ghost nodes at boundary in each dimensions.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscInt dim = 0
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        CHKERR( DMDAGetInfo(self.dm,
                          &dim,
                          NULL, NULL, NULL,
                          NULL, NULL, NULL,
                          NULL, NULL,
                          &btx, &bty, &btz,
                          NULL) )
        return toDims(dim, btx, bty, btz)

    def setStencilType(self, stencil_type: StencilType) -> None:
        """Set the type of the communication stencil

        Logically collective.

        Parameters
        ----------
        stype
            The stencil type, use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.



        See Also
        --------
        petsc.DMDASetStencilType

        """
        cdef PetscDMDAStencilType stype = asStencil(stencil_type)
        CHKERR( DMDASetStencilType(self.dm, stype) )

    def getStencilType(self) -> StencilType:
        """Return the stencil type.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscDMDAStencilType stype = DMDA_STENCIL_BOX
        CHKERR( DMDAGetInfo(self.dm,
                          NULL,
                          NULL, NULL, NULL,
                          NULL, NULL, NULL,
                          NULL, NULL,
                          NULL, NULL, NULL,
                          &stype) )
        return stype

    def setStencilWidth(self, stencil_width: int) -> None:
        """Set the width of the communication stencil.

        Logically collective.

        Parameters
        ----------
        width - The stencil width

        See Also
        --------
        petsc.DMDASetStencilWidth

        """
        cdef PetscInt swidth = asInt(stencil_width)
        CHKERR( DMDASetStencilWidth(self.dm, swidth) )

    def getStencilWidth(self) -> int:
        """Return the stencil width.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscInt swidth = 0
        CHKERR( DMDAGetInfo(self.dm,
                            NULL,
                            NULL, NULL, NULL,
                            NULL, NULL, NULL,
                            NULL, &swidth,
                            NULL, NULL, NULL,
                            NULL) )
        return toInt(swidth)

    def setStencil(
        self,
        stencil_type: StencilType,
        stencil_width: int,
    ) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc.DMDASetStencilType, petsc.DMDASetStencilWidth

        """
        cdef PetscDMDAStencilType stype = asStencil(stencil_type)
        cdef PetscInt swidth = asInt(stencil_width)
        CHKERR( DMDASetStencilType(self.dm, stype) )
        CHKERR( DMDASetStencilWidth(self.dm, swidth) )

    def getStencil(self) -> tuple[StencilType, int]:
        """Return the stencil type and width.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc.DMDAGetInfo

        """
        cdef PetscDMDAStencilType stype = DMDA_STENCIL_BOX
        cdef PetscInt swidth = 0
        CHKERR( DMDAGetInfo(self.dm,
                            NULL,
                            NULL, NULL, NULL,
                            NULL, NULL, NULL,
                            NULL, &swidth,
                            NULL, NULL, NULL,
                            &stype) )
        return (toStencil(stype), toInt(swidth))

    #

    # TODO: ??????
    def getRanges(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc.DMDAGetCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetCorners(self.dm,
                               &x, &y, &z,
                               &m, &n, &p) )
        return ((toInt(x), toInt(x+m)),
                (toInt(y), toInt(y+n)),
                (toInt(z), toInt(z+p)))[:<Py_ssize_t>dim]

    # TODO: ??????
    def getGhostRanges(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc.DMDAGetGhostCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetGhostCorners(self.dm,
                                    &x, &y, &z,
                                    &m, &n, &p) )
        return ((toInt(x), toInt(x+m)),
                (toInt(y), toInt(y+n)),
                (toInt(z), toInt(z+p)))[:<Py_ssize_t>dim]

    def getOwnershipRanges(self) -> tuple[Sequence[int]] | tuple[Sequence[int], Sequence[int]] | tuple[Sequence[int], Sequence[int], Sequence[int]]:
        """Return the ranges of indices in the x, y and z direction that are owned by each process.

        Output Parameters
        lx - ownership along x direction (optional)
        ly - ownership along y direction (optional)
        lz - ownership along z direction (optional)
        Note
        These correspond to the optional final arguments passed to DMDACreate(), DMDACreate2d(), DMDACreate3d()

        In C you should not free these arrays, nor change the values in them. They will only have valid values while the DMDA they came from still exists (has not been destroyed).

        These numbers are NOT multiplied by the number of dof per node.

        Not collective.

        See Also
        --------
        petsc.DMDAGetOwnershipRanges

        """
        cdef PetscInt dim=0, m=0, n=0, p=0
        cdef const PetscInt *lx = NULL, *ly = NULL, *lz = NULL
        CHKERR( DMDAGetInfo(self.dm,
                            &dim,
                            NULL, NULL, NULL,
                            &m, &n, &p,
                            NULL, NULL,
                            NULL, NULL, NULL,
                            NULL) )
        CHKERR( DMDAGetOwnershipRanges(self.dm, &lx, &ly, &lz) )
        return toOwnershipRanges(dim, m, n, p, lx, ly, lz)

    # TODO: ??????????
    def getCorners(self):
        """Return the global (x,y,z) indices of the lower left corner and size of the local region, excluding ghost points.

        Output Parameters
        x - the corner index for the first dimension
        y - the corner index for the second dimension (only used in 2D and 3D problems)
        z - the corner index for the third dimension (only used in 3D problems)
        m - the width in the first dimension
        n - the width in the second dimension (only used in 2D and 3D problems)
        p - the width in the third dimension (only used in 3D problems)
        Note
        The corner information is independent of the number of degrees of freedom per node set with the DMDACreateXX() routine. Thus the x, y, z, and m, n, p can be thought of as coordinates on a logical grid, where each grid point has (potentially) several degrees of freedom. Any of y, z, n, and p can be passed in as NULL if not needed.

        Not collective.

        See Also
        --------
        petsc.DMDAGetCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetCorners(self.dm,
                               &x, &y, &z,
                               &m, &n, &p) )
        return ((toInt(x), toInt(y), toInt(z))[:<Py_ssize_t>dim],
                (toInt(m), toInt(n), toInt(p))[:<Py_ssize_t>dim])

    # TODO: ???????
    def getGhostCorners(self):
        """Return the global (x,y,z) indices of the lower left corner and size of the local region, including ghost points.

        Output Parameters
        x - the corner index for the first dimension
        y - the corner index for the second dimension (only used in 2D and 3D problems)
        z - the corner index for the third dimension (only used in 3D problems)
        m - the width in the first dimension
        n - the width in the second dimension (only used in 2D and 3D problems)
        p - the width in the third dimension (only used in 3D problems)
        Note
        The corner information is independent of the number of degrees of freedom per node set with the DMDACreateXX() routine. Thus the x, y, z, and m, n, p can be thought of as coordinates on a logical grid, where each grid point has (potentially) several degrees of freedom. Any of y, z, n, and p can be passed in as NULL if not needed.

        Not collective.

        See Also
        --------
        petsc.DMDAGetGhostCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetGhostCorners(self.dm,
                                    &x, &y, &z,
                                    &m, &n, &p) )
        return ((toInt(x), toInt(y), toInt(z))[:<Py_ssize_t>dim],
                (toInt(m), toInt(n), toInt(p))[:<Py_ssize_t>dim])

    #

    def setFieldName(self, field: int, name: str) -> None:
        """Set the names of individual field components in multicomponent vectors associated with a DMDA.

        It must be called after having called DMSetUp().

        Logically collective; name must contain a common value.

        Parameters
        ----------
        nf
            field number for the DMDA (0, 1, ... dof-1), where dof indicates the number of degrees of freedom per node within the DMDA
        names
            the name of the field (component)

        See Also
        --------
        petsc.DMDASetFieldName

        """
        cdef PetscInt ival = asInt(field)
        cdef const char *cval = NULL
        name = str2bytes(name, &cval)
        CHKERR( DMDASetFieldName(self.dm, ival, cval) )

    def getFieldName(self, field: int) -> str:
        """Return the names of individual field components in multicomponent vectors associated with a DMDA.

        It must be called after having called DMSetUp().

        Not collective; name will contain a common value.

        Parameters
        ----------
        da
            the distributed array
        nf
            field number for the DMDA (0, 1, ... dof-1), where dof indicates the number of degrees of freedom per node within the DMDA

        See Also
        --------
        petsc.DMDAGetFieldName

        """
        cdef PetscInt ival = asInt(field)
        cdef const char *cval = NULL
        CHKERR( DMDAGetFieldName(self.dm, ival, &cval) )
        return bytes2str(cval)

    #

    # TODO: ?????
    def getVecArray(self, Vec vec) -> Any:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc._DMDA_Vec_array

        """
        return _DMDA_Vec_array(self, vec)

    #

    def setUniformCoordinates(
        self,
        xmin: float | None = 0,
        xmax: float | None = 1,
        ymin: float | None = 0,
        ymax: float | None = 1,
        zmin: float | None = 0,
        zmax: float | None = 1,
    ) -> None:
        """Set the DMDA coordinates to be a uniform grid.

        Collective.

        Parameters
        ----------
        xmin
            The minimum in the ``x`` direction.
        xmax
            The maximum in the ``x`` direction.
        ymin
            The minimum in the ``y`` direction (value ignored for 1 dimensional problems).
        ymax
            The maximum in the ``y`` direction (value ignored for 1 dimensional problems).
        zmin
            The minimum in the ``z`` direction (value ignored for 1 or 2 dimensional problems).
        zmax
            The maximum in the ``z`` direction (value ignored for 1 or 2 dimensional problems).

        See Also
        --------
        petsc.DMDASetUniformCoordinates

        """
        cdef PetscReal _xmin = asReal(xmin), _xmax = asReal(xmax)
        cdef PetscReal _ymin = asReal(ymin), _ymax = asReal(ymax)
        cdef PetscReal _zmin = asReal(zmin), _zmax = asReal(zmax)
        CHKERR( DMDASetUniformCoordinates(self.dm,
                                          _xmin, _xmax,
                                          _ymin, _ymax,
                                          _zmin, _zmax) )

    def setCoordinateName(self, index: int, name: str) -> None:
        """Set the name of the coordinate directions associated with a DMDA, for example "x" or "y"

        It must be called after having called DMSetUp().

        Logically collective; name must contain a common value.

        Parameters
        ----------
        index
            coordinate number for the DMDA (0, 1, ... dim-1),
        name
            the name of the coordinate

        See Also
        --------
        petsc.DMDASetCoordinateName

        """
        cdef PetscInt ival = asInt(index)
        cdef const char *cval = NULL
        name = str2bytes(name, &cval)
        CHKERR( DMDASetCoordinateName(self.dm, ival, cval) )

    def getCoordinateName(self, index: int) -> str:
        """Return the name of a coordinate direction associated with a DMDA.

        Not collective; name will contain a common value; No Fortran Support

        Parameters
        ----------
        nf
            number for the DMDA (0, 1, ... dim-1)

        See Also
        --------
        petsc.DMDAGetCoordinateName

        """
        cdef PetscInt ival = asInt(index)
        cdef const char *cval = NULL
        CHKERR( DMDAGetCoordinateName(self.dm, ival, &cval) )
        return bytes2str(cval)

    #

    def createNaturalVec(self) -> Vec:
        """Create a parallel PETSc vector that will hold vector values in the natural numbering, rather than in the PETSc parallel numbering associated with the DMDA.

        The output parameter, g, is a regular PETSc vector that should be destroyed with a call to VecDestroy() when usage is finished.

        The number of local entries in the vector on each process is the same as in a vector created with DMCreateGlobalVector().

        Output Parameter
        g - the distributed global vector

        Collective.

        See Also
        --------
        petsc.DMDACreateNaturalVector

        """
        cdef Vec vn = Vec()
        CHKERR( DMDACreateNaturalVector(self.dm, &vn.vec) )
        return vn

    def globalToNatural(
        self,
        Vec vg,
        Vec vn,
        addv: InsertMode | None = None,
    ) -> None:
        """Map values from the global vector to a global vector in the "natural" grid ordering. Must be followed by DMDAGlobalToNaturalEnd() to complete the exchange.

        Output Parameter
        l - the natural ordering values
        l - the natural ordering values
        l - the global values in the natural ordering

        Notes
        The global and natural vectors used here need not be the same as those obtained from DMCreateGlobalVector() and DMDACreateNaturalVector(), BUT they must have the same parallel data layout; they could, for example, be obtained with VecDuplicate() from the DMDA originating vectors.

        You must call DMDACreateNaturalVector() before using this routine

        Neighbor-wise collective.

        Parameters
        ----------
        g
            the global vector
        mode
            one of INSERT_VALUES or ADD_VALUES

        See Also
        --------
        petsc.DMDAGlobalToNaturalBegin, petsc.DMDAGlobalToNaturalEnd

        """
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMDAGlobalToNaturalBegin(self.dm, vg.vec, im, vn.vec) )
        CHKERR( DMDAGlobalToNaturalEnd  (self.dm, vg.vec, im, vn.vec) )

    def naturalToGlobal(
        self,
        Vec vn,
        Vec vg,
        addv: InsertMode | None = None,
    ) -> None:
        """Map values from a global vector from "natural" to grid ordering.

        Maps values from a global vector in the "natural" ordering to a global
        vector in the DMDA grid ordering.

        The global and natural vectors used here need not be the same as those
        obtained from DMCreateGlobalVector() and DMDACreateNaturalVector(), BUT
        they must have the same parallel data layout; they could, for example, be obtained with VecDuplicate() from the DMDA originating vectors.

        Neighbor-wise collective.

        Parameters
        ----------
        vn
            The global vector in a natural ordering.
        vg
            the global vector in a grid ordering.
        addv
            Insertion mode.

        See Also
        --------
        petsc.DMDANaturalToGlobalBegin, petsc.DMDANaturalToGlobalEnd

        """
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMDANaturalToGlobalBegin(self.dm, vn.vec, im, vg.vec) )
        CHKERR( DMDANaturalToGlobalEnd  (self.dm, vn.vec, im, vg.vec) )

    #

    def getAO(self) -> AO:
        """Return the application ordering context for a distributed array.

        In this case, the AO maps to the natural grid ordering that would be used for the DMDA if only 1 processor were employed (ordering most rapidly in the x-direction, then y, then z). Multiple degrees of freedom are numbered for each node (rather than 1 component for the whole grid, then the next component, etc.)

        Collective.

        See Also
        --------
        petsc.DMDAGetAO

        """
        cdef AO ao = AO()
        CHKERR( DMDAGetAO(self.dm, &ao.ao) )
        PetscINCREF(ao.obj)
        return ao

    def getScatter(self) -> tuple[Scatter, Scatter]:
        """Return the global-to-local, and local-to-local scatter contexts.

        Collective.

        See Also
        --------
        petsc.DMDAGetScatter

        """
        cdef Scatter l2g = Scatter()
        cdef Scatter g2l = Scatter()
        CHKERR( DMDAGetScatter(self.dm, &l2g.sct, &g2l.sct) )
        PetscINCREF(l2g.obj)
        PetscINCREF(g2l.obj)
        return (l2g, g2l)

    #

    def setRefinementFactor(
        self,
        refine_x: int | None = 2,
        refine_y: int | None = 2,
        refine_z: int | None = 2,
    ) -> None:
        """Set the ratios for the DMDA grid refinement.

        Logically collective.

        Parameters
        ----------
        refine_x
            Ratio of fine grid to coarse in x direction.
        refine_y
            Ratio of fine grid to coarse in y direction.
        refine_z
            Ratio of fine grid to coarse in z direction.

        See Also
        --------
        petsc.DMDASetRefinementFactor

        """
        cdef PetscInt refine[3]
        refine[0] = asInt(refine_x)
        refine[1] = asInt(refine_y)
        refine[2] = asInt(refine_z)
        CHKERR( DMDASetRefinementFactor(self.dm,
                                      refine[0],
                                      refine[1],
                                      refine[2]) )

    def getRefinementFactor(self):
        """Return the ratios that the DMDA grid is refined in each dimension.

        Not collective.

        See Also
        --------
        petsc.DMDAGetRefinementFactor

        """
        cdef PetscInt i, dim = 0, refine[3]
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetRefinementFactor(self.dm,
                                      &refine[0],
                                      &refine[1],
                                      &refine[2]) )
        return tuple([toInt(refine[i]) for 0 <= i < dim])

    def setInterpolationType(self, interp_type: InterpolationType) -> None:
        """Set the type of interpolation.

        You should call this on the coarser of the two DMDAs you pass to
        `DM.createInterpolation`.

        Logically collective.

        Parameters
        ----------
        ctype - DMDA_Q1 and DMDA_Q0 are currently the only supported forms.

        See Also
        --------
        petsc.DMDASetInterpolationType

        """
        cdef PetscDMDAInterpolationType ival = dainterpolationtype(interp_type)
        CHKERR( DMDASetInterpolationType(self.dm, ival) )

    # TODO: or just int64/long? why is it not a string?
    def getInterpolationType(self) -> int:
        """Return the type of interpolation.

        Not collective.

        See Also
        --------
        petsc.DMDAGetInterpolationType

        """
        cdef PetscDMDAInterpolationType ival = DMDA_INTERPOLATION_Q0
        CHKERR( DMDAGetInterpolationType(self.dm, &ival) )
        return <long>ival

    #

    def setElementType(self, elem_type: ElementType | str) -> None:
        """Set the element type to be returned by DMDAGetElements()

        Output Parameter
        etype - the element type, currently either DMDA_ELEMENT_P1 or DMDA_ELEMENT_Q1

        Not collective.

        See Also
        --------
        petsc.DMDASetElementType

        """
        cdef PetscDMDAElementType ival = daelementtype(elem_type)
        CHKERR( DMDASetElementType(self.dm, ival) )

    # TODO: or just int? why is it not a string?
    def getElementType(self):
        """Return the element type to be returned by DMDAGetElements()

        Output Parameter
        etype - the element type, currently either DMDA_ELEMENT_P1 or DMDA_ELEMENT_Q1

        Not collective.

        See Also
        --------
        petsc.DMDAGetElementType

        """
        cdef PetscDMDAElementType ival = DMDA_ELEMENT_Q1
        CHKERR( DMDAGetElementType(self.dm, &ival) )
        return <long>ival

    def getElements(self, elem_type: ElementType | None = None) -> ArrayInt:
        """Return an array containing the indices (in local coordinates) of all the local elements.


        Call DMDARestoreElements() once you have finished accessing the elements.

        Each process uniquely owns a subset of the elements. That is no element is owned by two or more processes.

        If on each process you integrate over its owned elements and use ADD_VALUES in Vec/MatSetValuesLocal() then you'll obtain the correct result.

        Output Parameters
        nel - number of local elements
        nen - number of element nodes
        e - the local indices of the elements’ vertices

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        petsc.DMDAGetElements

        """
        cdef PetscInt dim=0
        cdef PetscDMDAElementType etype
        cdef PetscInt nel=0, nen=0
        cdef const PetscInt *elems=NULL
        cdef object elements
        CHKERR( DMDAGetDim(self.dm, &dim) )
        if elem_type is not None:
            etype = daelementtype(elem_type)
            CHKERR( DMDASetElementType(self.dm, etype) )
        try:
            CHKERR( DMDAGetElements(self.dm, &nel, &nen, &elems) )
            elements = array_i(nel*nen, elems)
            elements.shape = (toInt(nel), toInt(nen))
        finally:
            CHKERR( DMDARestoreElements(self.dm, &nel, &nen, &elems) )
        return elements

    #

    property dim:
        """The dimension."""
        def __get__(self) -> int:
            return self.getDim()

    property dof:
        """The number of dof associated with each stratum of the grid."""
        def __get__(self) -> int:
            return self.getDof()

    property sizes:
        def __get__(self) -> tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)]:
            return self.getSizes()

    property proc_sizes:
        """The number of ranks in each direction in the global decomposition."""
        def __get__(self) -> tuple[()] | tuple[(int,)] | tuple[(int, int)] | tuple[(int, int, int)]:
            return self.getProcSizes()

    property boundary_type:
        """Boundary types in each direction."""
        def __get__(self) -> tuple[(str,)] | tuple[(str, str)] | tuple[(str, str, str)]:
            return self.getBoundaryType()

    property stencil:
        """Stencil type and width."""
        def __get__(self) -> tuple[StencilType, int]:
            return self.getStencil()

    property stencil_type:
        """Elementwise ghost/halo stencil type."""
        def __get__(self) -> str:
            return self.getStencilType()

    property stencil_width:
        """Elementwise stencil width."""
        def __get__(self) -> int:
            return self.getStencilWidth()

    # TODO: fix type once determined above
    property ranges:
        def __get__(self) -> None:
            return self.getRanges()

    # TODO: fix type once determined above
    property ghost_ranges:
        def __get__(self) -> None:
            return self.getGhostRanges()

    # TODO: fix type once determined above
    property corners:
        def __get__(self) -> None:
            return self.getCorners()

    # TODO: fix type once determined above
    property ghost_corners:
        def __get__(self) -> None:
            return self.getGhostCorners()

    # backward compatibility
    createNaturalVector = createNaturalVec


# backward compatibility alias
DA = DMDA

# --------------------------------------------------------------------

del DMDAStencilType
del DMDAInterpolationType
del DMDAElementType

# --------------------------------------------------------------------
