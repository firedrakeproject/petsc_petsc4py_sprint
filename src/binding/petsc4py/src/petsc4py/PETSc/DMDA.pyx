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

    # TODO: please verify
    def create(
        self,
        dim: int | None = None,
        dof: int | None = None,
        sizes: tuple[int, ...] | None = None,
        proc_sizes: tuple[int, ...] | None = None,
        boundary_type: tuple[DM.BoundaryType, ...] | tuple[int, ...] | tuple[str, ...] | tuple[bool, ...] | None = None,
        stencil_type: StencilType | None = None,
        stencil_width: int | None = None,
        bint setup: bool | None = True,
        ownership_ranges: tuple[Sequence[int], ...] | None = None,
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
            The number of dimensions.
        dof
            The number of degrees of freedom.
        sizes
            The number of elements in each dimension.
        proc_sizes
            The number of processes in x, y, z dimensions.
        boundary_type
            The boundary types.
        stencil_type
            The ghost/halo stencil type.
        stencil_width
            The width of the ghost/halo region.
        setup
            Whether to call the setup routine after creating the object.
        ownership_ranges
            Local x, y, z element counts, of length equal to ``proc_sizes``,
            summing to ``sizes``.
        comm
            The MPI communicator.

        See Also
        --------
        petsc.DMDACreate, petsc.DMSetDimension, petsc.DMDASetDof
        petsc.DMDASetSizes, petsc.DMDASetNumProcs
        petsc.DMDASetOwnershipRanges, petsc.DMDASetBoundaryType
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
        boundary_type: tuple[DM.BoundaryType, ...] | tuple[int, ...] | tuple[str, ...] | tuple[bool, ...] | None = None,
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
            Boundary types.
        stencil_type
            The ghost/halo stencil type.
        stencil_width
            The width of the ghost/halo region.

        See Also
        --------
        create, petsc.DMDAGetInfo, petsc.DMSetUp

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
        getDim, petsc.DMSetDimension

        """
        return self.setDimension(dim)

    def getDim(self) -> int:
        """Return the topological dimension.

        Not collective.

        See Also
        --------
        setDim, petsc.DMGetDimension

        """
        return self.getDimension()

    def setDof(self, dof: int) -> None:
        """Set the number of degrees of freedom per vertex.

        Not collective.

        Parameters
        ----------
        dof
            The number of degrees of freedom.

        Parameters
        ----------
        TODO
            TODO.

        See Also
        --------
        getDof, petsc.DMDASetDof

        """
        cdef PetscInt ndof = asInt(dof)
        CHKERR( DMDASetDof(self.dm, ndof) )

    def getDof(self) -> int:
        """Return the number of degrees of freedom per node.

        Not collective.

        See Also
        --------
        setDof, petsc.DMDAGetInfo

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

    def setSizes(self, sizes: tuple[int, ...]) -> None:
        """Set the number of grid points in each direction.

        Logically collective.

        Parameters
        ----------
        sizes
            The global (x), (x, y), or (x, y, z) size.

        See Also
        --------
        getSizes, petsc.DMDASetSizes

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

    def getSizes(self) -> tuple[int, ...]:
        """Return the global dimension in one/two/three directions.

        Not collective.

        See Also
        --------
        setSizes, petsc.DMDAGetInfo

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

    def setProcSizes(self, proc_sizes: tuple[int, ...]) -> None:
        """Set the number of processes in each dimension.

        Logically collective.

        Parameters
        ----------
        proc_sizes
            The number of processes in (x), (x, y), or (x, y, z) dimensions.

        See Also
        --------
        getProcSizes, petsc.DMDASetNumProcs

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

    def getProcSizes(self) -> tuple[int, ...]:
        """Return the number of processes in each dimensions.

        Not collective.

        See Also
        --------
        setProcSizes, petsc.DMDAGetInfo

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
        boundary_type: tuple[DM.BoundaryType, ...] | tuple[int, ...] | tuple[str, ...] | tuple[bool, ...]
    ) -> None:
        """Set the type of ghost nodes on domain boundaries.

        Not collective.

        Parameters
        ----------
        boundary_type
            The boundary type in (x), (x, y), or (x, y, z) dimensions.

        See Also
        --------
        getBoundaryType, petsc.DMDASetBoundaryType

        """
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        asBoundary(boundary_type, &btx, &bty, &btz)
        CHKERR( DMDASetBoundaryType(self.dm, btx, bty, btz) )

    # TODO: expected a string, but this seems to be an int
    def getBoundaryType(self) -> tuple[int, ...]:
        """Return the type of ghost nodes at boundary in each dimension.

        Not collective.

        See Also
        --------
        setBoundaryType

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
        """Set the type of the communication stencil.

        Logically collective.

        Parameters
        ----------
        stype
            The stencil type.

        See Also
        --------
        getStencilType, setStencil, petsc.DMDASetStencilType

        """
        cdef PetscDMDAStencilType stype = asStencil(stencil_type)
        CHKERR( DMDASetStencilType(self.dm, stype) )

    def getStencilType(self) -> StencilType:
        """Return the stencil type.

        Not collective.

        See Also
        --------
        setStencilType, petsc.DMDAGetInfo

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
        stencil_width
            The stencil width.

        See Also
        --------
        getStencilWidth, setStencil, petsc.DMDASetStencilWidth

        """
        cdef PetscInt swidth = asInt(stencil_width)
        CHKERR( DMDASetStencilWidth(self.dm, swidth) )

    def getStencilWidth(self) -> int:
        """Return the stencil width.

        Not collective.

        See Also
        --------
        setStencilWidth

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
        """Set the stencil type and width.

        Not collective.

        Parameters
        ----------
        stencil_type
            The stencil type.
        stencil_width
            The stencil width.

        See Also
        --------
        setStencilWidth, setStencilType, petsc.DMDASetStencilType
        petsc.DMDASetStencilWidth

        """
        cdef PetscDMDAStencilType stype = asStencil(stencil_type)
        cdef PetscInt swidth = asInt(stencil_width)
        CHKERR( DMDASetStencilType(self.dm, stype) )
        CHKERR( DMDASetStencilWidth(self.dm, swidth) )

    def getStencil(self) -> tuple[StencilType, int]:
        """Return the stencil type and width.

        Not collective.

        See Also
        --------
        getStencilType, getStencilWidth

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

    # TODO: should we clarify how is getRanges different than getGhostRanges and getCorners?
    def getRanges(self) -> tuple[tuple[int, int], ...]:
        """Return the ranges of the local region in each dimension.

        Excluding ghost nodes.

        Not collective.

        See Also
        --------
        getGhostRanges, petsc.DMDAGetCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetCorners(self.dm,
                               &x, &y, &z,
                               &m, &n, &p) )
        return ((toInt(x), toInt(x+m)),
                (toInt(y), toInt(y+n)),
                (toInt(z), toInt(z+p)))[:<Py_ssize_t>dim]

    def getGhostRanges(self) -> tuple[tuple[int, int], ...]:
        """Return the ranges of local region in each dim. including ghost nodes.

        Not collective.

        See Also
        --------
        getRanges, petsc.DMDAGetGhostCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetGhostCorners(self.dm,
                                    &x, &y, &z,
                                    &m, &n, &p) )
        return ((toInt(x), toInt(x+m)),
                (toInt(y), toInt(y+n)),
                (toInt(z), toInt(z+p)))[:<Py_ssize_t>dim]

    # TODO: this seems to be a sequence and not a tuple, right?
    def getOwnershipRanges(self) -> tuple[Sequence[int], ...]:
        """Return the ranges of indices in each direction owned by each process.

        These numbers are not multiplied by the number of DOFs per node.

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

    def getCorners(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        """Return lower left corner and size of local region in each dimension.

        Returns the global (x,y,z) indices of the lower left corner (first
        tuple) and size of the local region (second tuple).

        Excluding ghost points.

        The corner information is independent of the number of degrees of
        freedom per node. Thus the returned values can be thought of as
        coordinates on a logical grid, where each grid point has (potentially)
        several degrees of freedom.

        Not collective.

        See Also
        --------
        getGhostCorners, petsc.DMDAGetCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMDAGetDim(self.dm, &dim) )
        CHKERR( DMDAGetCorners(self.dm,
                               &x, &y, &z,
                               &m, &n, &p) )
        return ((toInt(x), toInt(y), toInt(z))[:<Py_ssize_t>dim],
                (toInt(m), toInt(n), toInt(p))[:<Py_ssize_t>dim])

    def getGhostCorners(self) -> tuple[()] | tuple[tuple[int], tuple[int]] | tuple[tuple[int, int], tuple[int, int]] | tuple[tuple[int, int, int], tuple[int, int, int]]:
        """Return lower left corner and size of local region in each dimension.

        Returns the global (x,y,z) indices of the lower left corner (first
        tuple) and size of the local region (second tuple).

        Excluding ghost points.

        The corner information is independent of the number of degrees of
        freedom per node. Thus the returned values can be thought of as
        coordinates on a logical grid, where each grid point has (potentially)
        several degrees of freedom.

        Not collective.

        See Also
        --------
        getCorners, petsc.DMDAGetGhostCorners

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
        """Set the name of individual field components.

        Logically collective; ``name`` must contain a common value.

        Parameters
        ----------
        field
            The field number for the DMDA (``0``, ``1``, ..., ``dof-1``),
            where ``dof`` indicates the number of degrees of freedom per node
            within the `DMDA`.
        name
            The name of the field (component).

        See Also
        --------
        getFieldName, petsc.DMDASetFieldName

        """
        cdef PetscInt ival = asInt(field)
        cdef const char *cval = NULL
        name = str2bytes(name, &cval)
        CHKERR( DMDASetFieldName(self.dm, ival, cval) )

    def getFieldName(self, field: int) -> str:
        """Return the name of an individual field component.

        Not collective; name will contain a common value.

        Parameters
        ----------
        field
            The field number for the DMDA (``0``, ``1``, ..., ``dof-1``),
            where ``dof`` indicates the number of degrees of freedom per node
            within the `DMDA`.

        See Also
        --------
        setFieldName, petsc.DMDAGetFieldName

        """
        cdef PetscInt ival = asInt(field)
        cdef const char *cval = NULL
        CHKERR( DMDAGetFieldName(self.dm, ival, &cval) )
        return bytes2str(cval)

    #

    def getVecArray(self, Vec vec) -> Any:
        """Get access to the vector.

        Use via `with` context manager (PEP 343).

        Not collective.

        Parameters
        ----------
        vec
            The vector to which access is being requested.

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
            The minimum in the ``y`` direction (value ignored for 1 dimensional
             problems).
        ymax
            The maximum in the ``y`` direction (value ignored for 1 dimensional
             problems).
        zmin
            The minimum in the ``z`` direction (value ignored for 1 or 2
            dimensional problems).
        zmax
            The maximum in the ``z`` direction (value ignored for 1 or 2
            dimensional problems).

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
        """Set the name of the coordinate direction.

        Logically collective; ``name`` must contain a common value.

        Parameters
        ----------
        index
            The coordinate number for the DMDA (``0``, ``1``, ..., ``dim-1``).
        name
            The name of the coordinate.

        See Also
        --------
        getCoordinateName, petsc.DMDASetCoordinateName

        """
        cdef PetscInt ival = asInt(index)
        cdef const char *cval = NULL
        name = str2bytes(name, &cval)
        CHKERR( DMDASetCoordinateName(self.dm, ival, cval) )

    def getCoordinateName(self, index: int) -> str:
        """Return the name of a coordinate direction.

        Not collective; the returned value will contain a common value.

        Parameters
        ----------
        index
            The coordinate number for the DMDA (``0``, ``1``, ..., ``dim-1``).

        See Also
        --------
        setCoordinateName, petsc.DMDAGetCoordinateName

        """
        cdef PetscInt ival = asInt(index)
        cdef const char *cval = NULL
        CHKERR( DMDAGetCoordinateName(self.dm, ival, &cval) )
        return bytes2str(cval)

    #

    def createNaturalVec(self) -> Vec:
        """Create a vector that will hold values in the natural numbering.

        The number of local entries in the vector on each process is the same
        as in a vector created with `DM.createGlobalVector`.

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
        """Map values to the "natural" grid ordering.

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
        vg
            The global vector in a grid ordering.
        vn
            The global vector in a "natural" ordering.
        addv
            The insertion mode.

        See Also
        --------
        naturalToGlobal, petsc.DMDAGlobalToNaturalBegin
        petsc.DMDAGlobalToNaturalEnd

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
        """Map values the to grid ordering.

        Neighbor-wise collective.

        Parameters
        ----------
        vn
            The global vector in a natural ordering.
        vg
            the global vector in a grid ordering.
        addv
            The insertion mode.

        See Also
        --------
        globalToNatural, petsc.DMDANaturalToGlobalBegin
        petsc.DMDANaturalToGlobalEnd

        """
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMDANaturalToGlobalBegin(self.dm, vn.vec, im, vg.vec) )
        CHKERR( DMDANaturalToGlobalEnd  (self.dm, vn.vec, im, vg.vec) )

    #

    def getAO(self) -> AO:
        """Return the application ordering context for a distributed array.

        In this case, the AO maps to the natural grid ordering that would be
        used for the `DMDA` if only 1 processor were employed (ordering most
        rapidly in the x-direction, then y, then z). Multiple degrees of
        freedom are numbered for each node (rather than 1 component for the
        whole grid, then the next component, etc.).

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
        getRefinementFactor, petsc.DMDASetRefinementFactor

        """
        cdef PetscInt refine[3]
        refine[0] = asInt(refine_x)
        refine[1] = asInt(refine_y)
        refine[2] = asInt(refine_z)
        CHKERR( DMDASetRefinementFactor(self.dm,
                                      refine[0],
                                      refine[1],
                                      refine[2]) )

    def getRefinementFactor(self) -> tuple[int, ...]:
        """Return the ratios that the DMDA grid is refined in each dimension.

        Not collective.

        See Also
        --------
        setRefinementFactor, petsc.DMDAGetRefinementFactor

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
        interp_type
            The interpolation type.

        See Also
        --------
        getInterpolationType, petsc.DMDASetInterpolationType

        """
        cdef PetscDMDAInterpolationType ival = dainterpolationtype(interp_type)
        CHKERR( DMDASetInterpolationType(self.dm, ival) )

    def getInterpolationType(self) -> int:
        """Return the type of interpolation.

        Not collective.

        See Also
        --------
        setInterpolationType, petsc.DMDAGetInterpolationType

        """
        cdef PetscDMDAInterpolationType ival = DMDA_INTERPOLATION_Q0
        CHKERR( DMDAGetInterpolationType(self.dm, &ival) )
        return <long>ival

    #

    def setElementType(self, elem_type: ElementType | str) -> None:
        """Set the element type to be returned by `getElements`.

        Not collective.

        See Also
        --------
        getElementType, petsc.DMDASetElementType

        """
        cdef PetscDMDAElementType ival = daelementtype(elem_type)
        CHKERR( DMDASetElementType(self.dm, ival) )

    def getElementType(self) -> int:
        """Return the element type to be returned by `getElements`.

        Not collective.

        See Also
        --------
        setElementType, petsc.DMDAGetElementType

        """
        cdef PetscDMDAElementType ival = DMDA_ELEMENT_Q1
        CHKERR( DMDAGetElementType(self.dm, &ival) )
        return <long>ival

    def getElements(self, elem_type: ElementType | None = None) -> ArrayInt:
        """Return an array containing the indices of all the local elements.

        The elements are in local coordinates.

        Each process uniquely owns a subset of the elements. That is, no
        element is owned by two or more processes.

        Not collective.

        Parameters
        ----------
        elem_type
            The element type.

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
        """The global dimension in each direction."""
        def __get__(self) -> tuple[int, ...]:
            return self.getSizes()

    property proc_sizes:
        """The number of ranks in each direction in the global decomposition."""
        def __get__(self) -> tuple[int, ...]:
            return self.getProcSizes()

    # TODO: see not above at getBoundaryType
    property boundary_type:
        """Boundary types in each direction."""
        def __get__(self) -> tuple[int, ...]:
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

    property ranges:
        """Ranges of the local region in each dimension."""
        def __get__(self) -> tuple[tuple[int, int], ...]:
            return self.getRanges()

    property ghost_ranges:
        """Ranges of local region in each dim. including ghost nodes."""
        def __get__(self) -> tuple[tuple[int, int], ...]:
            return self.getGhostRanges()

    property corners:
        """The lower left corner and size of local region in each dimension."""
        def __get__(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
            return self.getCorners()

    property ghost_corners:
        """The lower left corner and size of local region in each dimension."""
        def __get__(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
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
