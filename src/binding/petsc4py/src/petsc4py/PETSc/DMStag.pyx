# --------------------------------------------------------------------

class DMStagStencilType(object):
    STAR = DMSTAG_STENCIL_STAR
    BOX  = DMSTAG_STENCIL_BOX
    NONE = DMSTAG_STENCIL_NONE

class DMStagStencilLocation(object):
    NULLLOC          = DMSTAG_NULL_LOCATION
    BACK_DOWN_LEFT   = DMSTAG_BACK_DOWN_LEFT
    BACK_DOWN        = DMSTAG_BACK_DOWN
    BACK_DOWN_RIGHT  = DMSTAG_BACK_DOWN_RIGHT
    BACK_LEFT        = DMSTAG_BACK_LEFT
    BACK             = DMSTAG_BACK
    BACK_RIGHT       = DMSTAG_BACK_RIGHT
    BACK_UP_LEFT     = DMSTAG_BACK_UP_LEFT
    BACK_UP          = DMSTAG_BACK_UP
    BACK_UP_RIGHT    = DMSTAG_BACK_UP_RIGHT
    DOWN_LEFT        = DMSTAG_DOWN_LEFT
    DOWN             = DMSTAG_DOWN
    DOWN_RIGHT       = DMSTAG_DOWN_RIGHT
    LEFT             = DMSTAG_LEFT
    ELEMENT          = DMSTAG_ELEMENT
    RIGHT            = DMSTAG_RIGHT
    UP_LEFT          = DMSTAG_UP_LEFT
    UP               = DMSTAG_UP
    UP_RIGHT         = DMSTAG_UP_RIGHT
    FRONT_DOWN_LEFT  = DMSTAG_FRONT_DOWN_LEFT
    FRONT_DOWN       = DMSTAG_FRONT_DOWN
    FRONT_DOWN_RIGHT = DMSTAG_FRONT_DOWN_RIGHT
    FRONT_LEFT       = DMSTAG_FRONT_LEFT
    FRONT            = DMSTAG_FRONT
    FRONT_RIGHT      = DMSTAG_FRONT_RIGHT
    FRONT_UP_LEFT    = DMSTAG_FRONT_UP_LEFT
    FRONT_UP         = DMSTAG_FRONT_UP
    FRONT_UP_RIGHT   = DMSTAG_FRONT_UP_RIGHT

# --------------------------------------------------------------------

cdef class DMStag(DM):
    """TODO"""
    StencilType       = DMStagStencilType
    StencilLocation   = DMStagStencilLocation

    def create(
        self,
        dim,
        dofs=None,
        sizes=None,
        boundary_types: None # DM.BoundaryType | None = None, # sequence? 0,1,2,3
        stencil_type: StencilType | None = None,
        stencil_width: int | None = None,
        proc_sizes=None,
        ownership_ranges=None,
        comm: Comm | None = None,
        setUp: bool | None = False,
    ) -> Self:
        """TODO
        Create an object to manage data living on the elements and vertices of a parallelized regular 1D grid.
        Create an object to manage data living on the elements, faces, and vertices of a parallelized regular 2D grid.
        Create an object to manage data living on the elements, faces, edges, and vertices of a parallelized regular 3D grid.

        Collective.

        notes

        Parameters
        ----------
        bndx
            boundary type: DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, or DM_BOUNDARY_GHOSTED
        M
            global number of elements
        dof0
            number of degrees of freedom per vertex/0-cell
        dof1
            number of degrees of freedom per element/1-cell

        dof0
            number of degrees of freedom per vertex/0-cell
        dof1
            number of degrees of freedom per face/1-cell
        dof2
            number of degrees of freedom per element/2-cell

        dof0
            number of degrees of freedom per vertex/0-cell
        dof1
            number of degrees of freedom per edge/1-cell
        dof2
            number of degrees of freedom per face/2-cell
        dof3
            number of degrees of freedom per element/3-cell

        stencil_type
            Ghost/halo region type.
        stencil_width
            Width, in elements, of halo/ghost region.
        lx
            array of local sizes, of length equal to the comm size, summing to M

        See also
        --------
        petsc.DMStagCreate1d, petsc.DMStagCreate2d, petsc.DMStagCreate3d,
        petsc.DMSetUp

        """
        # TODO: do all see alsos render? they didn't before
        # ndim
        cdef PetscInt ndim = asInt(dim)

        # sizes
        cdef object gsizes = sizes
        cdef PetscInt nsizes=PETSC_DECIDE, M=1, N=1, P=1
        if sizes is not None:
            nsizes = asStagDims(gsizes, &M, &N, &P)
            assert(nsizes==ndim)

        # dofs
        cdef object cdofs = dofs
        cdef PetscInt ndofs=PETSC_DECIDE, dof0=1, dof1=0, dof2=0, dof3=0
        if dofs is not None:
            ndofs = asDofs(cdofs, &dof0, &dof1, &dof2, &dof3)
            assert(ndofs==ndim+1)

        # boundary types
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        asBoundary(boundary_types, &btx, &bty, &btz)

        # stencil
        cdef PetscInt swidth = 0
        if stencil_width is not None:
            swidth = asInt(stencil_width)
        cdef PetscDMStagStencilType stype = DMSTAG_STENCIL_NONE
        if stencil_type is not None:
            stype = asStagStencil(stencil_type)

        # comm
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)

        # proc sizes
        cdef object psizes = proc_sizes
        cdef PetscInt nprocs=PETSC_DECIDE, m=PETSC_DECIDE, n=PETSC_DECIDE, p=PETSC_DECIDE
        if proc_sizes is not None:
            nprocs = asStagDims(psizes, &m, &n, &p)
            assert(nprocs==ndim)

        # ownership ranges
        cdef PetscInt *lx = NULL, *ly = NULL, *lz = NULL
        if ownership_ranges is not None:
            nranges = asStagOwnershipRanges(ownership_ranges, ndim, &m, &n, &p, &lx, &ly, &lz)
            assert(nranges==ndim)

        # create
        cdef PetscDM newda = NULL
        if dim == 1:
            CHKERR( DMStagCreate1d(ccomm, btx, M, dof0, dof1, stype, swidth, lx, &newda) )
        if dim == 2:
            CHKERR( DMStagCreate2d(ccomm, btx, bty, M, N, m, n, dof0, dof1, dof2, stype, swidth, lx, ly, &newda) )
        if dim == 3:
            CHKERR( DMStagCreate3d(ccomm, btx, bty, btz, M, N, P, m, n, p, dof0, dof1, dof2, dof3, stype, swidth, lx, ly, lz, &newda) )
        PetscCLEAR(self.obj); self.dm = newda
        if setUp:
            CHKERR( DMSetUp(self.dm) )
        return self

    # Setters

    def setStencilWidth(self,swidth: int) -> None:
        """Set elementwise stencil width.

        Logically collective; stencilWidth must contain common value.

        The width value is not used when `StencilType.NONE` is specified.

        Parameters
        ----------
        swidth
            Stencil/halo/ghost width in elements.

        See also
        --------
        petsc.DMStagSetStencilWidth

        """
        cdef PetscInt sw = asInt(swidth)
        CHKERR( DMStagSetStencilWidth(self.dm, sw) )

    def setStencilType(self, stenciltype: StencilType) -> None:
        """Set elementwise ghost/halo stencil type.

        Logically collective; stenciltype must contain common value.

        Parameters
        ----------
        stenciltype
            The elementwise ghost stencil type.

        See also
        --------
        petsc.DMStagSetStencilType

        """
        cdef PetscDMStagStencilType stype = asStagStencil(stenciltype)
        CHKERR( DMStagSetStencilType(self.dm, stype) )

    def setBoundaryTypes(self, boundary_types: DM.BoundaryType) -> None: # TODO: should be seq?
        """Set DMSTAG boundary types.

        Logically collective; boundaryType0, boundaryType1, and boundaryType2
        must contain common values.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

        Parameters
        ----------
        boundaryTypeX
            Boundary type for x direction.
        boundaryTypeY
            Boundary type for y direction, not set for one dimensional problems.
        boundaryTypeZ
            Boundary type for z direction, not set for one and two dimensional problems.


        See also
        --------
        petsc.DMStagSetBoundaryTypes

        """
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        asBoundary(boundary_types, &btx, &bty, &btz)
        CHKERR( DMStagSetBoundaryTypes(self.dm, btx, bty, btz) )

    def setDof(self, dofs) -> None:
        """Set dof/stratum.

        Logically collective; dof0, dof1, dof2, and dof3 must contain common values.

        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids.

        Parameters
        ----------
        dm
            The DMSTAG object.
        dof0
            The number of points per 0-cell (vertex/node).
        dof1
            The number of points per 1-cell (element in 1D, edge in 2D and 3D).
        dof2
            The number of points per 2-cell (element in 2D, face in 3D).
        dof3
            The number of points per 3-cell (element in 3D).

        See also
        --------
        petsc.DMStagSetDOF

        """
        cdef tuple gdofs = tuple(dofs)
        cdef PetscInt gdim=PETSC_DECIDE, dof0=1, dof1=0, dof2=0, dof3=0
        gdim = asDofs(gdofs, &dof0, &dof1, &dof2, &dof3)
        CHKERR( DMStagSetDOF(self.dm, dof0, dof1, dof2, dof3) )

    def setGlobalSizes(self, sizes) -> None:
        """Set global element counts in each direction.

        Logically collective; N0, N1, and N2 must contain common values.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

        Parameters
        ----------
        dm
            The DMSTAG object.
        N0
            Global elementwise size in the x direction.
        N1
            Global elementwise size in the y direction.
        N2
            Global elementwise size in the z direction.

        See also
        --------
        petsc.DMStagSetGlobalSizes

        """
        cdef tuple gsizes = tuple(sizes)
        cdef PetscInt gdim=PETSC_DECIDE, M=1, N=1, P=1
        gdim = asStagDims(gsizes, &M, &N, &P)
        CHKERR( DMStagSetGlobalSizes(self.dm, M, N, P) )

    def setProcSizes(self, sizes) -> None:
        """Set ranks in each direction in the global rank grid.

        Logically collective; nRanks0, nRanks1, and nRanks2 must contain common values.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

        Parameters
        ----------
        dm
            The DMSTAG object.
        nRanks0
            Number of ranks in the x direction.
        nRanks1
            Number of ranks in the y direction.
        nRanks2
            Number of ranks in the z direction.

        See also
        --------
        petsc.DMStagSetNumRanks

        """
        cdef tuple psizes = tuple(sizes)
        cdef PetscInt pdim=PETSC_DECIDE, m=PETSC_DECIDE, n=PETSC_DECIDE, p=PETSC_DECIDE
        pdim = asStagDims(psizes, &m, &n, &p)
        CHKERR( DMStagSetNumRanks(self.dm, m, n, p) )

    def setOwnershipRanges(self, ranges) -> None:
        """Set elements per rank in each direction.

        Logically collective; lx, ly, and lz must contain common values.

        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        Parameters
        ----------
        lx
            element counts for each rank in the x direction
        ly
            element counts for each rank in the y direction
        lz
            element counts for each rank in the z direction

        See also
        --------
        petsc.DMStagSetOwnershipRanges

        """
        cdef PetscInt dim=0, m=PETSC_DECIDE, n=PETSC_DECIDE, p=PETSC_DECIDE
        cdef PetscInt *lx = NULL, *ly = NULL, *lz = NULL
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetNumRanks(self.dm, &m, &n, &p) )
        ownership_ranges = asStagOwnershipRanges(ranges, dim, &m, &n, &p, &lx, &ly, &lz)
        CHKERR( DMStagSetOwnershipRanges(self.dm, lx, ly, lz) )

    # Getters

    def getDim(self):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        ``TODO``

        """
        return self.getDimension()

    def getEntriesPerElement(self) -> int:
        """Return number of entries per element in the local representation.

        Notes
        This is the natural block size for most local operations. In 1D it is equal to dof0++dof1,
        in 2D it is equal to dof0+2+2dof1++dof2, and in 3D it is equal to dof0+3+3dof1+3+3dof2++dof3.

        Output Parameter
        entriesPerElement
            Number of entries associated with each element in the local representation.

        Not collective.

        See also
        --------
        petsc.DMStagGetEntriesPerElement

        """
        cdef PetscInt epe=0
        CHKERR( DMStagGetEntriesPerElement(self.dm, &epe) )
        return toInt(epe)

    def getStencilWidth(self) -> int:
        """Return elementwise stencil width.

        Not collective.

        Output Parameter
        stencilWidth
            Stencil/halo/ghost width in elements.

        See also
        --------
        petsc.DMStagGetStencilWidth

        """
        cdef PetscInt swidth=0
        CHKERR( DMStagGetStencilWidth(self.dm, &swidth) )
        return toInt(swidth)

    def getDof(self):
        """Get number of DOF associated with each stratum of the grid.

        Not collective.

        Output Parameters
        dof0
            The number of points per 0-cell (vertex/node).
        dof1
            The number of points per 1-cell (element in 1D, edge in 2D and 3D).
        dof2
            The number of points per 2-cell (element in 2D, face in 3D).
        dof3
            The number of points per 3-cell (element in 3D).

        See also
        --------
        petsc.DMStagGetDOF, petsc.DMGetDimension

        """
        cdef PetscInt dim=0, dof0=0, dof1=0, dof2=0, dof3=0
        CHKERR( DMStagGetDOF(self.dm, &dof0, &dof1, &dof2, &dof3) )
        CHKERR( DMGetDimension(self.dm, &dim) )
        return toDofs(dim+1,dof0,dof1,dof2,dof3)

    def getCorners(self):
        """Return global element indices of the local region (excluding ghost points).

        Not collective.

        Notes
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        The number of extra partial elements is either 1 or 0. The value is 1
        on right, top, and front non-periodic domain (“physical”) boundaries,
        in the x, y, and z directions respectively, and otherwise 0.

        Output Parameters
        x
            starting element index in first direction
        y
            starting element index in second direction
        z
            starting element index in third direction
        m
            element width in first direction
        n
            element width in second direction
        p
            element width in third direction
        nExtrax
            number of extra partial elements in first direction
        nExtray
            number of extra partial elements in second direction
        nExtraz
            number of extra partial elements in third direction

        See also
        --------
        petsc.DMStagGetCorners, petsc.DMGetDimension

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0, nExtrax=0, nExtray=0, nExtraz=0
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetCorners(self.dm, &x, &y, &z, &m, &n, &p, &nExtrax, &nExtray, &nExtraz) )
        return (asInt(x), asInt(y), asInt(z))[:<Py_ssize_t>dim], (asInt(m), asInt(n), asInt(p))[:<Py_ssize_t>dim], (asInt(nExtrax), asInt(nExtray), asInt(nExtraz))[:<Py_ssize_t>dim]

    def getGhostCorners(self):
        """Return global element indices of the local region, including ghost points.

        Not collective.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        Output Parameters
        x
            the starting element index in the first direction
        y
            the starting element index in the second direction
        z
            the starting element index in the third direction
        m
            the element width in the first direction
        n
            the element width in the second direction
        p
            the element width in the third direction

        See also
        --------
        petsc.DMStagGetGhostCorners

        """
        cdef PetscInt dim=0, x=0, y=0, z=0, m=0, n=0, p=0
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetGhostCorners(self.dm, &x, &y, &z, &m, &n, &p) )
        return (asInt(x), asInt(y), asInt(z))[:<Py_ssize_t>dim], (asInt(m), asInt(n), asInt(p))[:<Py_ssize_t>dim]

    def getLocalSizes(self) -> tuple[] | tuple[int] | tuple[int, int] | tuple[int, int, int]:
        """Return local elementwise sizes.

        Not collective.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        Output Parameters
        m
            local element counts (excluding ghosts) in the x direction
        n
            local element counts (excluding ghosts) in the y direction
        p
            local element counts (excluding ghosts) in the z direction

        See also
        --------
        petsc.DMStagGetLocalSizes

        """
        cdef PetscInt dim=0, m=PETSC_DECIDE, n=PETSC_DECIDE, p=PETSC_DECIDE
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetLocalSizes(self.dm, &m, &n, &p) )
        return toStagDims(dim, m, n, p)

    def getGlobalSizes(self) -> tuple[] | tuple[int] | tuple[int, int] | tuple[int, int, int]:
        """Return global element counts.

        Not collective.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        Output Parameters
        M
            global element counts in the x direction
        N
            global element counts in the y direction
        P
            global element counts in the z direction


        See also
        --------
        petsc.DMStagGetGlobalSizes

        """
        cdef PetscInt dim=0, m=PETSC_DECIDE, n=PETSC_DECIDE, p=PETSC_DECIDE
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetGlobalSizes(self.dm, &m, &n, &p) )
        return toStagDims(dim, m, n, p)

    def getProcSizes(self) -> tuple[] | tuple[int] | tuple[int, int] | tuple[int, int, int]:
        """Return number of ranks in each direction in the global grid decomposition.

        Not collective.

        nRanks0
            number of ranks in the x direction in the grid decomposition
        nRanks1
            number of ranks in the y direction in the grid decomposition
        nRanks2
            number of ranks in the z direction in the grid decomposition

        See also
        --------
        petsc.DMStagGetNumRanks

        """
        cdef PetscInt dim=0, m=PETSC_DECIDE, n=PETSC_DECIDE, p=PETSC_DECIDE
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetNumRanks(self.dm, &m, &n, &p) )
        return toStagDims(dim, m, n, p)

    def getStencilType(self):
        """Return elementwise ghost/halo stencil type.

        Not collective.

        Output Parameter
        stencilType
            The elementwise ghost stencil type: DMSTAG_STENCIL_BOX,
            DMSTAG_STENCIL_STAR, or DMSTAG_STENCIL_NONE.

        See also
        --------
        petsc.DMStagGetStencilType

        """
        cdef PetscDMStagStencilType stype = DMSTAG_STENCIL_BOX
        CHKERR( DMStagGetStencilType(self.dm, &stype) )
        return toStagStencil(stype)

    def getOwnershipRanges(self):
        """Return elements per rank in each direction.

        Not collective.

        Notes
        These correspond to the optional final arguments passed to
        DMStagCreate1d(), DMStagCreate2d(), and DMStagCreate3d().

        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        In C you should not free these arrays, nor change the values in them.
        They will only have valid values while the DMSTAG they came from still
        exists (has not been destroyed).

        lx
            ownership along x direction (optional)
        ly
            ownership along y direction (optional)
        lz
            ownership along z direction (optional)

        See also
        --------
        petsc.DMStagGetOwnershipRanges

        """
        cdef PetscInt dim=0, m=0, n=0, p=0
        cdef const PetscInt *lx = NULL, *ly = NULL, *lz = NULL
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetNumRanks(self.dm, &m, &n, &p) )
        CHKERR( DMStagGetOwnershipRanges(self.dm, &lx, &ly, &lz) )
        return toStagOwnershipRanges(dim, m, n, p, lx, ly, lz)

    def getBoundaryTypes(self):
        """Return boundary types.

        Not collective.

        Output Parameters
        boundaryTypeX
            boundary type for x direction
        boundaryTypeY
            boundary type for y direction, not set for one dimensional problems
        boundaryTypeZ
            boundary type for z direction, not set for one and two dimensional problems

        See also
        --------
        petsc.DMStagGetBoundaryTypes

        """
        cdef PetscInt dim=0
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetBoundaryTypes(self.dm, &btx, &bty, &btz) )
        return toStagBoundaryTypes(dim, btx, bty, btz)

    def getIsFirstRank(self) -> tuple[] | tuple[int] | tuple[int, int] | tuple[int, int, int]:
        """Return boolean value for whether this rank is first in each direction in the rank grid.

        Not collective.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        Output Parameters
        isFirstRank0
            whether this rank is first in the x direction
        isFirstRank1
            whether this rank is first in the y direction
        isFirstRank2
            whether this rank is first in the z direction

        See also
        --------
        petsc.DMStagGetIsFirstRank

        """
        cdef PetscBool rank0=PETSC_FALSE, rank1=PETSC_FALSE, rank2=PETSC_FALSE
        cdef PetscInt dim=0
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetIsFirstRank(self.dm, &rank0, &rank1, &rank2) )
        return toStagDims(dim, rank0, rank1, rank2)

    def getIsLastRank(self) -> tuple[] | tuple[int] | tuple[int, int] | tuple[int, int, int]:
        """Return boolean value for whether this rank is last in each direction in the rank grid.

        Not collective.

        Note
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids. These arguments may be set to NULL in this case.

        Output Parameters
        isFirstRank0
            whether this rank is last in the x direction
        isFirstRank1
            whether this rank is last in the y direction
        isFirstRank2
            whether this rank is last in the z direction

        See also
        --------
        petsc.DMStagGetIsLastRank

        """
        cdef PetscBool rank0=PETSC_FALSE, rank1=PETSC_FALSE, rank2=PETSC_FALSE
        cdef PetscInt dim=0
        CHKERR( DMGetDimension(self.dm, &dim) )
        CHKERR( DMStagGetIsLastRank(self.dm, &rank0, &rank1, &rank2) )
        return toStagDims(dim, rank0, rank1, rank2)

    # Coordinate-related functions

    def setUniformCoordinatesExplicit(
        self,
        xmin: float | None = 0,
        xmax: float | None = 1,
        ymin: float | None = 0,
        ymax: float | None = 1,
        zmin: float | None = 0,
        zmax: float | None = 1,
    ) -> None:
        """Set `DMStag` coordinates to be a uniform grid, storing all values.

        Collective. TODO:?

        Notes
        DMSTAG supports 2 different types of coordinate DM: either another
        DMSTAG, or a DMPRODUCT. If the grid is orthogonal, using DMPRODUCT
        should be more efficient.

        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids.

        See the manual page for DMStagSetUniformCoordinates() for information
        on how coordinates for dummy cells outside the physical domain boundary
        are populated.

        Parameters
        ----------
        dm
            the DMSTAG object
        xmin
            minimum global coordinate value in the x direction
        xmax
            maximum global coordinate values in the x direction
        ymin
            minimum global coordinate value in the y direction
        ymax
            maximum global coordinate value in the y direction
        zmin
            minimum global coordinate value in the z direction
        zmax
            maximum global coordinate value in the z direction

        See also
        --------
        petsc.DMStagSetUniformCoordinatesExplicit

        """
        cdef PetscReal _xmin = asReal(xmin), _xmax = asReal(xmax)
        cdef PetscReal _ymin = asReal(ymin), _ymax = asReal(ymax)
        cdef PetscReal _zmin = asReal(zmin), _zmax = asReal(zmax)
        CHKERR( DMStagSetUniformCoordinatesExplicit(self.dm, _xmin, _xmax, _ymin, _ymax, _zmin, _zmax) )

    def setUniformCoordinatesProduct(
        self,
        xmin: float | None = 0,
        xmax: float | None = 1,
        ymin: float | None = 0,
        ymax: float | None = 1,
        zmin: float | None = 0,
        zmax: float | None = 1,
    ) -> None:
        """Create uniform coordinates, as a product of 1D arrays.

        Set the coordinate DM to be a DMPRODUCT of 1D DMSTAG objects, each of
        which have a coordinate DM (also a 1d DMSTAG) holding uniform
        coordinates.

        Collective.

        Notes
        -----
        Arguments corresponding to higher dimensions are ignored for 1D and 2D
        grids.

        The per-dimension 1-dimensional DMSTAG objects that comprise the
        product always have active 0-cells (vertices, element boundaries) and
        1-cells (element centers).

        See the manual page for DMStagSetUniformCoordinates() for information
        on how coordinates for dummy cells outside the physical domain boundary
        are populated.

        Parameters
        ----------
        xmin
            minimum global coordinate value in the x direction
        xmax
            maximum global coordinate value in the x direction
        ymin
            minimum global coordinate value in the y direction
        ymax
            maximum global coordinate value in the y direction
        zmin
            minimum global coordinate value in the z direction
        zmax
            maximum global coordinate value in the z direction

        See also
        --------
        petsc.DMStagSetUniformCoordinatesProduct

        """
        cdef PetscReal _xmin = asReal(xmin), _xmax = asReal(xmax)
        cdef PetscReal _ymin = asReal(ymin), _ymax = asReal(ymax)
        cdef PetscReal _zmin = asReal(zmin), _zmax = asReal(zmax)
        CHKERR( DMStagSetUniformCoordinatesProduct(self.dm, _xmin, _xmax, _ymin, _ymax, _zmin, _zmax) )

    def setUniformCoordinates(
        self,
        xmin: float | None = 0,
        xmax: float | None = 1,
        ymin: float | None = 0,
        ymax: float | None = 1,
        zmin: float | None = 0,
        zmax: float | None = 1,
    ) -> None:
        """Set DMSTAG coordinates to be a uniform grid.

        Collective.

        Notes
        DMSTAG supports 2 different types of coordinate DM: DMSTAG and
        DMPRODUCT. Arguments corresponding to higher dimensions are ignored for
        1D and 2D grids.

        Local coordinates are populated (using DMSetCoordinatesLocal()),
        linearly extrapolated to ghost cells, including those outside the
        physical domain. This is also done in case of periodic boundaries,
        meaning that the same global point may have different coordinates in
        different local representations, which are equivalent assuming a
        periodicity implied by the arguments to this function, i.e. two points
        are equivalent if their difference is a multiple of
        ((xmax−xmin)) in the x direction,
        ((ymax−−ymin)) in the y direction,
        and((zmax−−zmin)) in the z direction.

        Parameters
        ----------
        xmin
            minimum global coordinate value in the x direction
        xmax
            maximum global coordinate values in the x direction
        ymin
            minimum global coordinate value in the y direction
        ymax
            maximum global coordinate value in the y direction
        zmin
            minimum global coordinate value in the z direction
        zmax
            maximum global coordinate value in the z direction

        See also
        --------
        petsc.DMStagSetUniformCoordinates

        """
        cdef PetscReal _xmin = asReal(xmin), _xmax = asReal(xmax)
        cdef PetscReal _ymin = asReal(ymin), _ymax = asReal(ymax)
        cdef PetscReal _zmin = asReal(zmin), _zmax = asReal(zmax)
        CHKERR( DMStagSetUniformCoordinates(self.dm, _xmin, _xmax, _ymin, _ymax, _zmin, _zmax) )

    def setCoordinateDMType(self, dmtype: DM.Type) -> None:
        """Set DM type to store coordinates.

        Logically collective; dmtype must contain common value.

        Parameters
        ----------
        dmtype
            `DM.Type` for coordinates, either DMSTAG or DMPRODUCT.

        See also
        --------
        petsc.DMStagSetCoordinateDMType

        """
        cdef PetscDMType cval = NULL
        dmtype = str2bytes(dmtype, &cval)
        CHKERR( DMStagSetCoordinateDMType(self.dm, cval) )

    # Location slot related functions

    def getLocationSlot(self, loc, c: int) -> int:
        """Return index to use in accessing raw local arrays.

        Not collective.

        Notes
        Provides an appropriate index to use with DMStagVecGetArray() and
        friends. This is required so that the user doesn't need to know about
        the ordering of dof associated with each local element.

        Parameters
        ----------
        loc
            Location relative to an element.
        c
            Component.

        Output Parameter
        slot
            index to use

        See also
        --------
        petsc.DMStagGetLocationSlot

        """
        cdef PetscInt slot=0
        cdef PetscInt comp=asInt(c)
        cdef PetscDMStagStencilLocation sloc = asStagStencilLocation(loc)
        CHKERR( DMStagGetLocationSlot(self.dm, sloc, comp, &slot) )
        return toInt(slot)

    def getProductCoordinateLocationSlot(self, loc) -> None:
        """Return slot for use with local product coordinate arrays.

        Not collective.

        Notes
        High-level helper function to get slot indices for 1D coordinate DMs,
        for use with DMStagGetProductCoordinateArrays() and related functions.

        For loc, one should use DMSTAG_LEFT, DMSTAG_ELEMENT, or DMSTAG_RIGHT
        for “previous”, “center” and “next” locations, respectively, in each
        dimension. One can equivalently use DMSTAG_DOWN or DMSTAG_BACK in place
        of DMSTAG_LEFT, and DMSTAG_UP or DMSTACK_FRONT in place of DMSTAG_RIGHT.

        This function checks that the coordinates are actually set up so that
        using the slots from any of the 1D coordinate sub-DMs are valid for all
        the 1D coordinate sub-DMs.

        Parameters
        ----------
        dm
            The DMSTAG object.
        loc
            The grid location.
        Output Parameter
        slot
            the index to use in local arrays


        See also
        --------
        petsc.DMStagGetProductCoordinateLocationSlot

        """
        cdef PetscInt slot=0
        cdef PetscDMStagStencilLocation sloc = asStagStencilLocation(loc)
        CHKERR( DMStagGetProductCoordinateLocationSlot(self.dm, sloc, &slot) )
        return toInt(slot)

    def getLocationDof(self, loc) -> int:
        """Return number of DOF associated with a given point in a DMSTAG grid.

        Not collective.

        loc
            grid point (see DMStagStencilLocation)
        Output Parameter
        dof
            the number of DOF (components) living at loc in dm

        See also
        --------
        petsc.DMStagGetLocationDOF

        """
        cdef PetscInt dof=0
        cdef PetscDMStagStencilLocation sloc = asStagStencilLocation(loc)
        CHKERR( DMStagGetLocationDOF(self.dm, sloc, &dof) )
        return toInt(dof)

    # Random other functions

    def migrateVec(self, Vec vec, DM dmTo, Vec vecTo) -> None:
        """Transfer a vector associated with a DMSTAG to a vector associated with a compatible DMSTAG.

        Collective.

        Parameters
        ----------
        vec
            The source vector, compatible with dm.
        dmTo
            The compatible destination DMSTAG object.
        vecTo
            The destination vector, compatible with dmTo.
        Notes
        Extra dof are ignored, and unfilled dof are zeroed. Currently only
        implemented to migrate global vectors to global vectors. For the
        definition of compatibility of DMs, see DMGetCompatibility().

        See also
        --------
        petsc.DMStagMigrateVec

        """
        CHKERR( DMStagMigrateVec(self.dm, vec.vec, dmTo.dm, vecTo.vec ) )

    def createCompatibleDMStag(self, dofs) -> DM:
        """Create a compatible DMSTAG with different dof/stratum.

        Collective.

        Parameters
        ----------
        dof0
            Number of dof on the first stratum in the new DMSTAG.
        dof1
            Number of dof on the second stratum in the new DMSTAG.
        dof2
            Number of dof on the third stratum in the new DMSTAG.
        dof3
            Number of dof on the fourth stratum in the new DMSTAG.

        Notes
        DOF supplied for strata too big for the dimension are ignored; these
        may be set to 0. For example, for a 2-dimensional DMSTAG, dof2 sets the
        number of dof per element, and dof3 is unused. For a 3-dimensional
        DMSTAG, dof3 sets the number of DOF per element.

        In contrast to DMDACreateCompatibleDMDA(), coordinates are not reused.

        See also
        --------
        petsc.DMStagCreateCompatibleDMStag

        """
        cdef tuple gdofs = tuple(dofs)
        cdef PetscInt gdim=PETSC_DECIDE, dof0=1, dof1=0, dof2=0, dof3=0
        gdim = asDofs(gdofs, &dof0, &dof1, &dof2, &dof3)
        cdef PetscDM newda = NULL
        CHKERR( DMStagCreateCompatibleDMStag(self.dm, dof0, dof1, dof2, dof3, &newda) )
        cdef DM newdm = type(self)()
        PetscCLEAR(newdm.obj); newdm.dm = newda
        return newdm

    def VecSplitToDMDA(self, Vec vec, loc, c: int):
        """Create a DMDA and Vec from a subgrid of a DMSTAG and its Vec.

        Collective.

        Notes
        If a c value of -k is provided, the first k DOF for that position are extracted, padding with zero values if needed. If a non-negative value is provided, a single DOF is extracted.

        The caller is responsible for destroying the created DMDA and Vec.

        Parameters
        ----------
        vec
            Vec object associated with dm.
        loc
            Which subgrid to extract (see DMStagStencilLocation).
        c
            Which component to extract (see note below).
        Output Parameters
        pda
            the DMDA
        pdavec
            the new Vec

        See also
        --------
        petsc.DMStagVecSplitToDMDA

        """
        cdef PetscInt pc = asInt(c)
        cdef PetscDMStagStencilLocation sloc = asStagStencilLocation(loc)
        cdef PetscDM pda = NULL
        cdef PetscVec pdavec = NULL
        CHKERR( DMStagVecSplitToDMDA(self.dm, vec.vec, sloc, pc, &pda, &pdavec) )
        cdef DM da = DMDA()
        PetscCLEAR(da.obj); da.dm = pda
        cdef Vec davec = Vec()
        PetscCLEAR(davec.obj); davec.vec = pdavec
        return (da,davec)

    def getVecArray(self, Vec vec):
        """**Not implemented.**"""
        raise NotImplementedError('getVecArray for DMStag not yet implemented in petsc4py')

    def get1dCoordinatecArrays(self):
        """**Not implemented.**"""
        raise NotImplementedError('get1dCoordinatecArrays for DMStag not yet implemented in petsc4py')

    property dim:
        # TODO: docstring
        def __get__(self):
            return self.getDim()

    property dofs:
        # TODO: docstring
        def __get__(self):
            return self.getDof()

    property entries_per_element:
        # TODO: docstring
        def __get__(self):
            return self.getEntriesPerElement()

    property global_sizes:
        # TODO: docstring
        def __get__(self):
            return self.getGlobalSizes()

    property local_sizes:
        # TODO: docstring
        def __get__(self):
            return self.getLocalSizes()

    property proc_sizes:
        # TODO: docstring
        def __get__(self):
            return self.getProcSizes()

    property boundary_types:
        # TODO: docstring
        def __get__(self):
            return self.getBoundaryTypes()

    property stencil_type:
        # TODO: docstring
        def __get__(self):
            return self.getStencilType()

    property stencil_width:
        # TODO: docstring
        def __get__(self):
            return self.getStencilWidth()

    property corners:
        # TODO: docstring
        def __get__(self):
            return self.getCorners()

    property ghost_corners:
        # TODO: docstring
        def __get__(self):
            return self.getGhostCorners()


# --------------------------------------------------------------------

del DMStagStencilType
del DMStagStencilLocation

# --------------------------------------------------------------------
