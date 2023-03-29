# --------------------------------------------------------------------

class DMPlexReorderDefaultFlag(object):
    NOTSET = DMPLEX_REORDER_DEFAULT_NOTSET
    FALSE  = DMPLEX_REORDER_DEFAULT_FALSE
    TRUE   = DMPLEX_REORDER_DEFAULT_TRUE

# --------------------------------------------------------------------

cdef class DMPlex(DM):

    ReorderDefaultFlag = DMPlexReorderDefaultFlag

    #

    def create(self, comm: Comm | None = None) -> Self:
        """Create a `DMPlex` object, which encapsulates an unstructured mesh, or CW complex, which can be expressed using a Hasse Diagram.

        Collective.

        Parameters
        ----------
        comm
            The communicator for the `DMPlex` object.

        See Also
        --------
        `DM`, `DMPlex`, `DMType`, `DM.create`, `DM.setType`, petsc.DMPlexCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM newdm = NULL
        CHKERR( DMPlexCreate(ccomm, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createFromCellList(self, dim: int, cells: Sequence[int], coords: Sequence[float], interpolate: bool = True, comm: Comm | None = None) -> Self:
        """Create `DMPlex` from a list of vertices for each cell (common mesh generator output), but only process 0 takes in the input.

        Collective.

        Parameters
        ----------
        dim
            The topological dimension of the mesh.
        cells
            An array of number of cells times number of vertices on each cell, only on process 0.
        coords
            An array of number of vertices times spatial dimension, only on process 0.
        interpolate
            Flag indicating that intermediate mesh entities (faces, edges) should be created automatically.
        comm
            The communicator.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.buildFromCellList`, `DMPlex.buildCoordinatesFromCellList`, `DMPlex.createFromDAG`, `DMPlex.create`, petsc.DMPlexCreateFromCellListPetsc

        """
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        cdef PetscInt  cdim = asInt(dim)
        cdef PetscInt  numCells = 0
        cdef PetscInt  numCorners = 0
        cdef PetscInt  *cellVertices = NULL
        cdef PetscInt  numVertices = 0
        cdef PetscInt  spaceDim= 0
        cdef PetscReal *vertexCoords = NULL
        cdef int npy_flags = NPY_ARRAY_ALIGNED|NPY_ARRAY_NOTSWAPPED|NPY_ARRAY_CARRAY
        cells  = PyArray_FROM_OTF(cells,  NPY_PETSC_INT,  npy_flags)
        coords = PyArray_FROM_OTF(coords, NPY_PETSC_REAL, npy_flags)
        if PyArray_NDIM(cells) != 2: raise ValueError(
                ("cell indices must have two dimensions: "
                 "cells.ndim=%d") % (PyArray_NDIM(cells)) )
        if PyArray_NDIM(coords) != 2: raise ValueError(
                ("coords vertices must have two dimensions: "
                 "coords.ndim=%d") % (PyArray_NDIM(coords)) )
        numCells     = <PetscInt>   PyArray_DIM(cells,  0)
        numCorners   = <PetscInt>   PyArray_DIM(cells,  1)
        numVertices  = <PetscInt>   PyArray_DIM(coords, 0)
        spaceDim     = <PetscInt>   PyArray_DIM(coords, 1)
        cellVertices = <PetscInt*>  PyArray_DATA(cells)
        vertexCoords = <PetscReal*> PyArray_DATA(coords)
        CHKERR( DMPlexCreateFromCellListPetsc(ccomm, cdim, numCells, numVertices,
                                              numCorners, interp, cellVertices,
                                              spaceDim, vertexCoords, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createBoxMesh(self, faces: Sequence[int], lower: Sequence[int] | None = (0,0,0), upper: Sequence[int] | None = (1,1,1),
                      simplex: bool | None = True, periodic: Sequence | str | int | bool | None = False, interpolate: bool | None = True, comm: Comm | None = None) -> Self:
        """Create a mesh on the tensor product of unit intervals (box) using simplices or tensor cells (hexahedra).

        Collective.

        Parameters
        ----------
        faces
            Number of faces per dimension, or ``None`` for (1,) in 1D and (2, 2) in 2D and (1, 1, 1) in 3D.
        lower
            The lower left corner, defaults to (0, 0, 0).
        upper
            The upper right corner, defaults to (1, 1, 1).
        simplex
            ``True`` for simplices, ``False`` for tensor cells, defaults to ``True``.
        periodic
            The boundary type for the X,Y,Z direction, or ``None`` for `DM_BOUNDARY_NONE`, defaults to ``False``.
        interpolate
            Flag to create intermediate mesh pieces (edges, faces), defaults to ``True``.
        comm
            The communicator for the `DMPlex` object.

        See Also
        --------
        `DM`, `DMPlex`, `DM.setFromOptions`, `DMPlex.createFromFile`, `DMPlex.createHexCylinderMesh`, `DM.setType`, `DM.create`, petsc.DMPlexCreateBoxMesh

        """
        cdef Py_ssize_t i = 0
        cdef PetscInt dim = 0, *cfaces = NULL
        faces = iarray_i(faces, &dim, &cfaces)
        assert dim >= 1 and dim <= 3
        cdef PetscReal clower[3]
        clower[0] = clower[1] = clower[2] = 0
        for i from 0 <= i < dim: clower[i] = lower[i]
        cdef PetscReal cupper[3]
        cupper[0] = cupper[1] = cupper[2] = 1
        for i from 0 <= i < dim: cupper[i] = upper[i]
        cdef PetscDMBoundaryType btype[3];
        asBoundary(periodic, &btype[0], &btype[1], &btype[2])
        cdef PetscBool csimplex = simplex
        cdef PetscBool cinterp = interpolate
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM   newdm = NULL
        CHKERR( DMPlexCreateBoxMesh(ccomm, dim, csimplex, cfaces,
                                    clower, cupper, btype, cinterp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createBoxSurfaceMesh(self, faces, lower=(0,0,0), upper=(1,1,1),
                             interpolate=True, comm: Comm | None = None):
        """DMPlexCreateBoxSurfaceMesh - Creates a mesh on the surface of the tensor product of unit intervals (box) using tensor cells (hexahedra).

        Collective.

        Parameters
        ----------
        comm
            The communicator for the `DM` object
        dim
            The spatial dimension of the box, so the resulting mesh is has dimension `dim`-1
        faces
            Number of faces per dimension, or ``None`` for (1,) in 1D and (2, 2) in 2D and (1, 1, 1) in 3D
        lower
            The lower left corner, or ``None`` for (0, 0, 0)
        upper
            The upper right corner, or ``None`` for (1, 1, 1)
        interpolate
            Flag to create intermediate mesh pieces (edges, faces)

        Returns
        -------
        dm
            The `DM` object

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMSetFromOptions`, `DMPlexCreateBoxMesh`, `DMPlexCreateFromFile`, `DMSetType`, `DMCreate`

        See Also
        --------
        petsc.DMPlexCreateBoxSurfaceMesh

        """
        cdef Py_ssize_t i = 0
        cdef PetscInt dim = 0, *cfaces = NULL
        faces = iarray_i(faces, &dim, &cfaces)
        assert dim >= 1 and dim <= 3
        cdef PetscReal clower[3]
        clower[0] = clower[1] = clower[2] = 0
        for i from 0 <= i < dim: clower[i] = lower[i]
        cdef PetscReal cupper[3]
        cupper[0] = cupper[1] = cupper[2] = 1
        for i from 0 <= i < dim: cupper[i] = upper[i]
        cdef PetscBool cinterp = interpolate
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM   newdm = NULL
        CHKERR( DMPlexCreateBoxSurfaceMesh(ccomm, dim, cfaces, clower, cupper, cinterp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createFromFile(self, filename, plexname="unnamed", interpolate=True, comm: Comm | None = None):
        """DMPlexCreateFromFile - This takes a filename and produces a `DM`

        Collective.

        Parameters
        ----------
        comm
            The communicator
        filename
            A file name
        plexname
            The object name of the resulting `DM`, also used for intra-datafile lookup by some formats
        interpolate
            Flag to create intermediate mesh pieces (edges, faces)

        Returns
        -------
        dm
            The `DM`

        Options Database Key:
        . -dm_plex_create_from_hdf5_xdmf - use the `PETSC_VIEWER_HDF5_XDMF` format for reading HDF5

        Use -dm_plex_create_ prefix to pass options to the internal PetscViewer, e.g.
        $ -dm_plex_create_viewer_hdf5_collective

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexCreateFromDAG`, `DMPlexCreateFromCellListPetsc`, `DMPlex.create`, `PetscObjectSetName`, `DMView`, `DMLoad`

        See Also
        --------
        petsc.DMPlexCreateFromFile

        """
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        cdef const char *cfile = NULL
        cdef const char *pname = NULL
        filename = str2bytes(filename, &cfile)
        plexname = str2bytes(plexname, &pname)
        CHKERR( DMPlexCreateFromFile(ccomm, cfile, pname, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createCGNS(self, cgid, interpolate=True, comm: Comm | None = None):
        """DMPlexCreateCGNS - Create a `DMPlex` mesh from a CGNS file.

        Collective.

        Parameters
        ----------
        comm
            The MPI communicator
        filename
            The name of the CGNS file
        interpolate
            Create faces and edges in the mesh

        Returns
        -------
        dm
            The `DM` object representing the mesh

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.create`, `DMPlexCreateCGNS`, `DMPlexCreateExodus`

        See Also
        --------
        petsc.DMPlexCreateCGNS

        """
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        cdef PetscInt  ccgid = asInt(cgid)
        CHKERR( DMPlexCreateCGNS(ccomm, ccgid, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createCGNSFromFile(self, filename, interpolate=True, comm: Comm | None = None):
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        cdef const char *cfile = NULL
        filename = str2bytes(filename, &cfile)
        CHKERR( DMPlexCreateCGNSFromFile(ccomm, cfile, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createExodusFromFile(self, filename, interpolate=True, comm: Comm | None = None):
        """DMPlexCreateExodusFromFile - Create a `DMPlex` mesh from an ExodusII file.

        Collective.

        Parameters
        ----------
        comm
            The MPI communicator
        filename
            The name of the ExodusII file
        interpolate
            Create faces and edges in the mesh

        Returns
        -------
        dm
            The `DM` object representing the mesh

        .seealso: [](chapter_unstructured), `DM`, `PETSCVIEWEREXODUSII`, `DMPlex`, `DMCreate`, `DMPlexCreateExodus`

        See Also
        --------
        petsc.DMPlexCreateExodusFromFile

        """
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        cdef const char *cfile = NULL
        filename = str2bytes(filename, &cfile)
        CHKERR( DMPlexCreateExodusFromFile(ccomm, cfile, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createExodus(self, exoid, interpolate=True, comm: Comm | None = None):
        """DMPlexCreateExodus - Create a `DMPlex` mesh from an ExodusII file ID.

        Collective.

        Parameters
        ----------
        comm
            The MPI communicator
        exoid
            The ExodusII id associated with a exodus file and obtained using ex_open
        interpolate
            Create faces and edges in the mesh

        Returns
        -------
        dm
            The `DM` object representing the mesh

        .seealso: [](chapter_unstructured), `DM`, `PETSCVIEWEREXODUSII`, `DMPlex`, `DMPlex`, `DMCreate`

        See Also
        --------
        petsc.DMPlexCreateExodus

        """
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        cdef PetscInt  cexoid = asInt(exoid)
        CHKERR( DMPlexCreateExodus(ccomm, cexoid, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createGmsh(self, Viewer viewer, interpolate=True, comm: Comm | None = None):
        """DMPlexCreateGmsh - Create a `DMPlex` mesh from a Gmsh file viewer

        Collective.

        Parameters
        ----------
        comm
            The MPI communicator
        viewer
            The `PetscViewer` associated with a Gmsh file
        interpolate
            Create faces and edges in the mesh

        Returns
        -------
        dm
            The `DM` object representing the mesh

        Options Database Keys:
        + -dm_plex_gmsh_hybrid        - Force triangular prisms to use tensor order
        . -dm_plex_gmsh_periodic      - Read Gmsh periodic section and construct a periodic Plex
        . -dm_plex_gmsh_highorder     - Generate high-order coordinates
        . -dm_plex_gmsh_project       - Project high-order coordinates to a different space, use the prefix dm_plex_gmsh_project_ to define the space
        . -dm_plex_gmsh_use_regions   - Generate labels with region names
        . -dm_plex_gmsh_mark_vertices - Add vertices to generated labels
        . -dm_plex_gmsh_multiple_tags - Allow multiple tags for default labels
        - -dm_plex_gmsh_spacedim <d>  - Embedding space dimension, if different from topological dimension

        Note:
        The Gmsh file format is described in http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

        By default, the "Cell Sets", "Face Sets", and "Vertex Sets" labels are created, and only insert the first tag on a point. By using -dm_plex_gmsh_multiple_tags, all tags can be inserted. Instead, -dm_plex_gmsh_use_regions creates labels based on the region names from the PhysicalNames section, and all tags are used.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMCreate`

        See Also
        --------
        petsc.DMPlexCreateGmsh

        """
        cdef MPI_Comm  ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscBool interp = interpolate
        cdef PetscDM   newdm = NULL
        CHKERR( DMPlexCreateGmsh(ccomm, viewer.vwr, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def createCohesiveSubmesh(self, hasLagrange, value):
        """DMPlexCreateCohesiveSubmesh - Extract from a mesh with cohesive cells the hypersurface defined by one face of the cells. Optionally, a label can be given to restrict the cells.

        Parameters
        ----------
        dm
            The original mesh
        hasLagrange
            The mesh has Lagrange unknowns in the cohesive cells
        label
            A label name, or ``None``
        value
            A label value

        Returns
        -------
        subdm
            The surface mesh

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getSubpointMap`, `DMPlexCreateSubmesh`

        See Also
        --------
        petsc.DMPlexCreateCohesiveSubmesh

        """
        cdef PetscBool flag = hasLagrange
        cdef PetscInt cvalue = asInt(value)
        cdef DM subdm = DMPlex()
        CHKERR( DMPlexCreateCohesiveSubmesh(self.dm, flag, NULL, cvalue, &subdm.dm) )
        return subdm

    def getChart(self) -> tuple[int, int]:
        """Return the interval for all mesh points [``pStart``, ``pEnd``).

        Not collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.setChart`, petsc.DMPlexGetChart

        """
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        return toInt(pStart), toInt(pEnd)

    def setChart(self, pStart: int, pEnd: int) -> None:
        """Set the interval for all mesh points [``pStart``, ``pEnd``).

        Not collective.

        Parameters
        ----------
        pStart
            The first mesh point.
        pEnd
            The upper bound for mesh points.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getChart`, petsc.DMPlexSetChart

        """
        cdef PetscInt cStart = asInt(pStart)
        cdef PetscInt cEnd   = asInt(pEnd)
        CHKERR( DMPlexSetChart(self.dm, cStart, cEnd) )

    def getConeSize(self, p: int) -> int:
        """Return the number of in-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.setConeSize`, `DMPlex.setChart`, petsc.DMPlexGetConeSize

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt csize = 0
        CHKERR( DMPlexGetConeSize(self.dm, cp, &csize) )
        return toInt(csize)

    def setConeSize(self, p: int, size: int) -> None:
        """Set the number of in-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.
        size
            The cone size for point ``p``.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getConeSize`, `DMPlex.setChart`, petsc.DMPlexSetConeSize

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt csize = asInt(size)
        CHKERR( DMPlexSetConeSize(self.dm, cp, csize) )

    def getCone(self, p: int) -> ndarray[int]:
        """Return the points on the in-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getConeSize`, `DMPlex.setCone`, `DMPlex.getConeTuple`, `DMPlex.setChart`, `DMPlex.restoreCone`, petsc.DMPlexGetCone

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt        ncone = 0
        cdef const PetscInt *icone = NULL
        CHKERR( DMPlexGetConeSize(self.dm, cp, &ncone) )
        CHKERR( DMPlexGetCone(self.dm, cp, &icone) )
        return array_i(ncone, icone)

    def setCone(self, p: int, cone: Sequence[int], orientation: Sequence[int] | None = None) -> None:
        """Set the points on the in-edges for this point in the DAG; that is these are the points that cover the specific point.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.
        cone
            An array of points which are on the in-edges for point ``p``.
        orientation
            An array of orientations, defaults to ``None``.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getCone`, `DMPlex.setChart`, `DMPlex.setConeSize`, `DM.setUp`, `DMPlex.setSupport`, `DMPlex.setSupportSize`, petsc.DMPlexSetCone

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        #
        cdef PetscInt  ncone = 0
        cdef PetscInt *icone = NULL
        cone = iarray_i(cone, &ncone, &icone)
        CHKERR( DMPlexSetConeSize(self.dm, cp, ncone) )
        CHKERR( DMPlexSetCone(self.dm, cp, icone) )
        #
        cdef PetscInt  norie = 0
        cdef PetscInt *iorie = NULL
        if orientation is not None:
            orientation = iarray_i(orientation, &norie, &iorie)
            assert norie == ncone
            CHKERR( DMPlexSetConeOrientation(self.dm, cp, iorie) )

    def insertCone(self, p: int, conePos: int, conePoint: int) -> None:
        """DMPlexInsertCone - Insert a point into the in-edges for the point p in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.
        conePos
            The local index in the cone where the point should be put.
        conePoint
            The mesh point to insert.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getCone`, `DMPlex.setChart`, `DMPlex.setConeSize`, `DM.setUp`, petsc.DMPlexInsertCone

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt cconePos = asInt(conePos)
        cdef PetscInt cconePoint = asInt(conePoint)
        CHKERR( DMPlexInsertCone(self.dm,cp,cconePos,cconePoint) )

    def insertConeOrientation(self, p: int, conePos: int, coneOrientation: int) -> None:
        """Insert a point orientation for the in-edge for the point p in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`
        conePos
            The local index in the cone where the point should be put.
        coneOrientation
            The point orientation to insert.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getCone`, `DMPlex.setChart`, `DMPlex.setConeSize`, `DM.setUp`, petsc.DMPlexInsertConeOrientation

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt cconePos = asInt(conePos)
        cdef PetscInt cconeOrientation = asInt(coneOrientation)
        CHKERR( DMPlexInsertConeOrientation(self.dm, cp, cconePos, cconeOrientation) )

    def getConeOrientation(self, p: int) -> ndarray[int]:
        """Return the orientations on the in-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.

        See Also
        --------
        `DM`, `DMPlex`, `DMPolytopeType.composeOrientation`, `DMPolytopeType.composeOrientationInv`, `DMPlex.create`, `DMPlex.getCone`, `DMPlex.setCone`, `DMPlex.setChart`, petsc.DMPlexGetConeOrientation

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt        norie = 0
        cdef const PetscInt *iorie = NULL
        CHKERR( DMPlexGetConeSize(self.dm, cp, &norie) )
        CHKERR( DMPlexGetConeOrientation(self.dm, cp, &iorie) )
        return array_i(norie, iorie)

    def setConeOrientation(self, p: int, orientation: Sequence[int]) -> None:
        """Set the orientations on the in-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.
        orientation
            An array of orientations.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getConeOrientation`, `DMPlex.setCone`, `DMPlex.setChart`, `DMPlex.setConeSize`, `DM.setUp`, petsc.DMPlexSetConeOrientation

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt ncone = 0
        CHKERR( DMPlexGetConeSize(self.dm, cp, &ncone) )
        cdef PetscInt  norie = 0
        cdef PetscInt *iorie = NULL
        orientation = iarray_i(orientation, &norie, &iorie)
        assert norie == ncone
        CHKERR( DMPlexSetConeOrientation(self.dm, cp, iorie) )

    def setCellType(self, p: int, ctype: DM.DMPolytopeType) -> None:
        """Set the polytope type of a given cell.

        Not collective.

        Parameters
        ----------
        p
            The cell.
        ctype
            The polytope type of the cell.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getCellTypeLabel`, `DMPlex.getDepthLabel`, `DMPlex.getDepth`, `DMPlex.computeCellTypes`, `DM.createLabel`, petsc.DMPlexSetCellType

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscDMPolytopeType val = ctype
        CHKERR( DMPlexSetCellType(self.dm, cp, val) )

    def getCellType(self, p: int) -> int:
        """Get the polytope type of a given cell.

        Not collective.

        Parameters
        ----------
        p
            The cell.

        See Also
        --------
        `DM`, `DMPlex`, `DMPolytopeType`, `DMPlex.getCellTypeLabel`, `DMPlex.getDepthLabel`, `DMPlex.getDepth`, petsc.DMPlexGetCellType

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscDMPolytopeType ctype = DM_POLYTOPE_UNKNOWN
        CHKERR( DMPlexGetCellType(self.dm, cp, &ctype) )
        return toInt(ctype)

    def getCellTypeLabel(self) -> DMLabel:
        """Get the `DMLabel` recording the polytope type of each cell.

        Not collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getCellType`, `DMPlex.getDepthLabel`, `DM.createLabel`, petsc.DMPlexGetCellTypeLabel

        """
        cdef DMLabel label = DMLabel()
        CHKERR( DMPlexGetCellTypeLabel(self.dm, &label.dmlabel) )
        PetscINCREF(label.obj)
        return label

    def getSupportSize(self, p: int) -> int:
        """Return the number of out-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.setConeSize`, `DMPlex.setChart`, `DMPlex.getConeSize`, petsc.DMPlexGetSupportSize

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt ssize = 0
        CHKERR( DMPlexGetSupportSize(self.dm, cp, &ssize) )
        return toInt(ssize)

    def setSupportSize(self, p: int, size: int) -> None:
        """Set the number of out-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.
        size
            The support size for point ``p``.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.getSupportSize`, `DMPlex.setChart`, petsc.DMPlexSetSupportSize

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt ssize = asInt(size)
        CHKERR( DMPlexSetSupportSize(self.dm, cp, ssize) )

    def getSupport(self, p: int) -> ndarray[int]:
        """Return the points on the out-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getSupportSize`, `DMPlex.setSupport`, `DMPlex.getCone`, `DMPlex.setChart`, petsc.DMPlexGetSupport

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt        nsupp = 0
        cdef const PetscInt *isupp = NULL
        CHKERR( DMPlexGetSupportSize(self.dm, cp, &nsupp) )
        CHKERR( DMPlexGetSupport(self.dm, cp, &isupp) )
        return array_i(nsupp, isupp)

    def setSupport(self, p: int, supp: Sequence[int]) -> None:
        """Set the points on the out-edges for this point in the DAG, that is the list of points that this point covers.

        Not collective.

        Parameters
        ----------
        p
            The point, which must lie in the chart set with `DMPlex.setChart`.
        supp
            An array of points which are on the out-edges for point ``p``.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.setCone`, `DMPlex.setConeSize`, `DMPlex.create`, `DMPlex.getSupport`, `DMPlex.setChart`, `DMPlex.setSupportSize`, `DM.setUp`, petsc.DMPlexSetSupport

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscInt  nsupp = 0
        cdef PetscInt *isupp = NULL
        supp = iarray_i(supp, &nsupp, &isupp)
        CHKERR( DMPlexSetSupportSize(self.dm, cp, nsupp) )
        CHKERR( DMPlexSetSupport(self.dm, cp, isupp) )

    def getMaxSizes(self) -> tuple[int, int]:
        """Return the maximum number of in-edges (cone) and out-edges (support) for any point in the DAG.

        Not collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.setConeSize`, `DMPlex.setChart`, petsc.DMPlexGetMaxSizes

        """
        cdef PetscInt maxConeSize = 0, maxSupportSize = 0
        CHKERR( DMPlexGetMaxSizes(self.dm, &maxConeSize, &maxSupportSize) )
        return toInt(maxConeSize), toInt(maxSupportSize)

    def symmetrize(self) -> None:
        """Create support (out-edge) information from cone (in-edge) information.

        Not collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.setChart`, `DMPlex.setConeSize`, `DMPlex.setCone`, petsc.DMPlexSymmetrize

        """
        CHKERR( DMPlexSymmetrize(self.dm) )

    def stratify(self) -> None:
        """The DAG for most topologies is a graded poset (https://en.wikipedia.org/wiki/Graded_poset), and
        can be illustrated by a Hasse Diagram (https://en.wikipedia.org/wiki/Hasse_diagram). The strata group all points of the
        same grade, and this function calculates the strata. This grade can be seen as the height (or depth) of the point in
        the DAG.

        Collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.create`, `DMPlex.symmetrize`, `DMPlex.computeCellTypes`, petsc.DMPlexStratify

        """
        CHKERR( DMPlexStratify(self.dm) )

    def orient(self) -> None:
        """Give a consistent orientation to the input mesh.

        See Also
        --------
        `DM`, `DMPlex`, `DM.create`, petsc.DMPlexOrient

        """
        CHKERR( DMPlexOrient(self.dm) )

    def getCellNumbering(self) -> IS:
        """Get a global cell numbering for all cells on this process.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getVertexNumbering`, petsc.DMPlexGetCellNumbering

        """
        cdef IS iset = IS()
        CHKERR( DMPlexGetCellNumbering(self.dm, &iset.iset) )
        PetscINCREF(iset.obj)
        return iset

    def getVertexNumbering(self) -> IS:
        """Get a global vertex numbering for all vertices on this process.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getCellNumbering`, petsc.DMPlexGetVertexNumbering

        """
        cdef IS iset = IS()
        CHKERR( DMPlexGetVertexNumbering(self.dm, &iset.iset) )
        PetscINCREF(iset.obj)
        return iset

    def createPointNumbering(self) -> IS:
        """Create a global numbering for all points.

        Collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getCellNumbering`, petsc.DMPlexCreatePointNumbering

        """
        cdef IS iset = IS()
        CHKERR( DMPlexCreatePointNumbering(self.dm, &iset.iset) )
        return iset

    def getDepth(self) -> int:
        """Get the depth of the DAG representing this mesh.

        Not collective.

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getDepthLabel`, `DMPlex.getDepthStratum`, `DMPlex.getPointDepth`, `DMPlex.symmetrize`, petsc.DMPlexGetDepth

        """
        cdef PetscInt depth = 0
        CHKERR( DMPlexGetDepth(self.dm,&depth) )
        return toInt(depth)

    def getDepthStratum(self, svalue: int) -> tuple[int, int]:
        """Get the bounds [``start``, ``end``) for all points at a certain depth.

        Not collective.

        Parameters
        ----------
        svalue
            The requested depth

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getHeightStratum`, `DMPlex.getDepth`, `DMPlex.getDepthLabel`, `DMPlex.getPointDepth`, `DMPlex.symmetrize`, `DMPlex.interpolate`, petsc.DMPlexGetDepthStratum

        """
        cdef PetscInt csvalue = asInt(svalue), sStart = 0, sEnd = 0
        CHKERR( DMPlexGetDepthStratum(self.dm, csvalue, &sStart, &sEnd) )
        return (toInt(sStart), toInt(sEnd))

    def getHeightStratum(self, svalue: int) -> tuple[int, int]:
        """Get the bounds [``start``, ``end``) for all points at a certain height.

        Not collective.

        Parameters
        ----------
        svalue
            The requested height

        See Also
        --------
        `DM`, `DMPlex`, `DMPlex.getDepthStratum`, `DMPlex.getDepth`, `DMPlex.getPointHeight`, petsc.DMPlexGetHeightStratum

        """
        cdef PetscInt csvalue = asInt(svalue), sStart = 0, sEnd = 0
        CHKERR( DMPlexGetHeightStratum(self.dm, csvalue, &sStart, &sEnd) )
        return (toInt(sStart), toInt(sEnd))

    def getMeet(self, points: Sequence[int]) -> ndarray[int]:
        """Get an array for the meet of the set of points.

        Not collective.

        Parameters
        ----------
        points
            The input points.

        See Also
        --------
        `DM`, `DMPLEX`, `DMPlex.restoreMeet`, `DMPlex.getJoin`, petsc.DMPlexGetMeet

        """
        cdef PetscInt  numPoints = 0
        cdef PetscInt *ipoints = NULL
        cdef PetscInt  numCoveringPoints = 0
        cdef const PetscInt *coveringPoints = NULL
        points = iarray_i(points, &numPoints, &ipoints)
        CHKERR( DMPlexGetMeet(self.dm, numPoints, ipoints, &numCoveringPoints, &coveringPoints) )
        try:
            return array_i(numCoveringPoints, coveringPoints)
        finally:
            CHKERR( DMPlexRestoreMeet(self.dm, numPoints, ipoints, &numCoveringPoints, &coveringPoints) )

    def getJoin(self, points: Sequence[int]) -> ndarray[int]:
        """Get an array for the join of the set of points.

        Not collective.

        Parameters
        ----------
        points
            The input points.

        See Also
        --------
        `DM`, `DMPLEX`, `DMPlex.restoreJoin`, `DMPlex.getMeet`, petsc.DMPlexGetJoin

        """
        cdef PetscInt  numPoints = 0
        cdef PetscInt *ipoints = NULL
        cdef PetscInt  numCoveringPoints = 0
        cdef const PetscInt *coveringPoints = NULL
        points = iarray_i(points, &numPoints, &ipoints)
        CHKERR( DMPlexGetJoin(self.dm, numPoints, ipoints, &numCoveringPoints, &coveringPoints) )
        try:
            return array_i(numCoveringPoints, coveringPoints)
        finally:
            CHKERR( DMPlexRestoreJoin(self.dm, numPoints, ipoints, &numCoveringPoints, &coveringPoints) )

    def getFullJoin(self, points: Sequence[int]) -> ndarray[int]:
        """Get an array for the join of the set of points.

        Not collective.

        Parameters
        ----------
        points
            The input points.

        See Also
        --------
        `DM`, `DMPLEX`, `DMPlex.getJoin`, `DMPlex.restoreJoin`, `DMPlex.getMeet`, petsc.DMPlexGetFullJoin

        """
        cdef PetscInt  numPoints = 0
        cdef PetscInt *ipoints = NULL
        cdef PetscInt  numCoveringPoints = 0
        cdef const PetscInt *coveringPoints = NULL
        points = iarray_i(points, &numPoints, &ipoints)
        CHKERR( DMPlexGetFullJoin(self.dm, numPoints, ipoints, &numCoveringPoints, &coveringPoints) )
        try:
            return array_i(numCoveringPoints, coveringPoints)
        finally:
            CHKERR( DMPlexRestoreJoin(self.dm, numPoints, ipoints, &numCoveringPoints, &coveringPoints) )

    def getTransitiveClosure(self, p: int, useCone: bool | None = True) -> tuple[ndarray[int], ndarray[int]]:
        """Return the points on the transitive closure of the in-edges or out-edges for this point in the DAG.

        Not collective.

        Parameters
        ----------
        p
            The mesh point.
        useCone
            ``True`` for the closure, otherwise return the star, defaults to ``True``.

        See Also
        --------
        `DM`, `DMPLEX`, `DMPlex.restoreTransitiveClosure`, `DMPlex.create`, `DMPlex.setCone`, `DMPlex.setChart`, `DMPlex.getCone`, petsc.DMPlexGetTransitiveClosure

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( DMPlexGetChart(self.dm, &pStart, &pEnd) )
        assert cp>=pStart and cp<pEnd
        cdef PetscBool cuseCone = useCone
        cdef PetscInt  numPoints = 0
        cdef PetscInt *points = NULL
        CHKERR( DMPlexGetTransitiveClosure(self.dm, cp, cuseCone, &numPoints, &points) )
        try:
            out = array_i(2*numPoints,points)
        finally:
            CHKERR( DMPlexRestoreTransitiveClosure(self.dm, cp, cuseCone, &numPoints, &points) )
        return out[::2],out[1::2]

    def vecGetClosure(self, Section sec, Vec vec, p):
        """DMPlexVecRestoreClosure - Restore the array of the values on the closure of 'point'

        Not collective

        Parameters
        ----------
        dm
            The `DM`
        section
            The section describing the layout in `v`, or ``None`` to use the default section
        v
            The local vector
        point
            The point in the `DM`
        csize
            The number of values in the closure, or ``None``
        values
            The array of values, which is a borrowed array and should not be freed

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexVecGetClosure`, `DMPlexVecSetClosure`, `DMPlexMatSetClosure`

        See Also
        --------
        petsc.DMPlexVecRestoreClosure

        """
        cdef PetscInt cp = asInt(p), csize = 0
        cdef PetscScalar *cvals = NULL
        CHKERR( DMPlexVecGetClosure(self.dm, sec.sec, vec.vec, cp, &csize, &cvals) )
        try:
            closure = array_s(csize, cvals)
        finally:
            CHKERR( DMPlexVecRestoreClosure(self.dm, sec.sec, vec.vec, cp, &csize, &cvals) )
        return closure

    def getVecClosure(self, Section sec or None, Vec vec, point):
        """DMPlexVecRestoreClosure - Restore the array of the values on the closure of 'point'

        Not collective

        Parameters
        ----------
        dm
            The `DM`
        section
            The section describing the layout in `v`, or ``None`` to use the default section
        v
            The local vector
        point
            The point in the `DM`
        csize
            The number of values in the closure, or ``None``
        values
            The array of values, which is a borrowed array and should not be freed

        Note:
        The array values are discarded and not copied back into `v`. In order to copy values back to `v`, use `DMPlexVecSetClosure`

        Fortran Note:
        The `csize` argument is not present in the Fortran binding since it is internal to the array.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexVecGetClosure`, `DMPlexVecSetClosure`, `DMPlexMatSetClosure`

        See Also
        --------
        petsc.DMPlexVecRestoreClosure

        """
        cdef PetscSection csec = sec.sec if sec is not None else NULL
        cdef PetscInt cp = asInt(point), csize = 0
        cdef PetscScalar *cvals = NULL
        CHKERR( DMPlexVecGetClosure(self.dm, csec, vec.vec, cp, &csize, &cvals) )
        try:
            closure = array_s(csize, cvals)
        finally:
            CHKERR( DMPlexVecRestoreClosure(self.dm, csec, vec.vec, cp, &csize, &cvals) )
        return closure

    def setVecClosure(self, Section sec or None, Vec vec, point, values, addv=None):
        """DMPlexVecSetClosure - Set an array of the values on the closure of `point`

        Not collective

        Parameters
        ----------
        dm
            The `DM`
        section
            The section describing the layout in `v`, or ``None`` to use the default section
        v
            The local vector
        point
            The point in the `DM`
        values
            The array of values
        mode
            The insert mode. One of `INSERT_ALL_VALUES`, `ADD_ALL_VALUES`, `INSERT_VALUES`, `ADD_VALUES`, `INSERT_BC_VALUES`, and `ADD_BC_VALUES`,
        where `INSERT_ALL_VALUES` and `ADD_ALL_VALUES` also overwrite boundary conditions.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexVecGetClosure`, `DMPlexMatSetClosure`

        See Also
        --------
        petsc.DMPlexVecSetClosure

        """
        cdef PetscSection csec = sec.sec if sec is not None else NULL
        cdef PetscInt cp = asInt(point)
        cdef PetscInt csize = 0
        cdef PetscScalar *cvals = NULL
        cdef object tmp = iarray_s(values, &csize, &cvals)
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMPlexVecSetClosure(self.dm, csec, vec.vec, cp, cvals, im) )

    def setMatClosure(self, Section sec or None, Section gsec or None,
                      Mat mat, point, values, addv=None):
        """DMPlexMatSetClosure - Set an array of the values on the closure of 'point'

        Not collective

        Parameters
        ----------
        dm
            The `DM`
        section
            The section describing the layout in `v`, or ``None`` to use the default section
        globalSection
            The section describing the layout in `v`, or ``None`` to use the default global section
        A
            The matrix
        point
            The point in the `DM`
        values
            The array of values
        mode
            The insert mode, where `INSERT_ALL_VALUES` and `ADD_ALL_VALUES` also overwrite boundary conditions

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexMatSetClosureGeneral`, `DMPlexVecGetClosure`, `DMPlexVecSetClosure`

        See Also
        --------
        petsc.DMPlexMatSetClosure

        """
        cdef PetscSection csec  =  sec.sec if  sec is not None else NULL
        cdef PetscSection cgsec = gsec.sec if gsec is not None else NULL
        cdef PetscInt cp = asInt(point)
        cdef PetscInt csize = 0
        cdef PetscScalar *cvals = NULL
        cdef object tmp = iarray_s(values, &csize, &cvals)
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMPlexMatSetClosure(self.dm, csec, cgsec, mat.mat, cp, cvals, im) )

    def generate(self, DMPlex boundary, name=None, interpolate=True):
        """DMPlexGenerate - Generates a mesh.

        Not collective.

        Parameters
        ----------
        boundary
            The `DMPlex` boundary object
        name
            The mesh generation package name
        interpolate
            Flag to create intermediate mesh elements

        Returns
        -------
        mesh
            The `DMPlex` object

        Options Database Keys:
        +  -dm_plex_generate <name> - package to generate mesh, for example, triangle, ctetgen or tetgen
        -  -dm_generator <name> - package to generate mesh, for example, triangle, ctetgen or tetgen

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.create`, `DMRefine`

        See Also
        --------
        petsc.DMPlexGenerate

        """
        cdef PetscBool interp = interpolate
        cdef const char *cname = NULL
        if name: name = str2bytes(name, &cname)
        cdef PetscDM   newdm = NULL
        CHKERR( DMPlexGenerate(boundary.dm, cname, interp, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def setTriangleOptions(self, opts):
        """DMPlexTriangleSetOptions - Set the options used for the Triangle mesh generator

        Not collective.

        Inputs Parameters:
        + dm - The `DMPlex` object
        - opts - The command line options

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexTetgenSetOptions`, `DMPlexGenerate`

        See Also
        --------
        petsc.DMPlexTriangleSetOptions

        """
        cdef const char *copts = NULL
        opts = str2bytes(opts, &copts)
        CHKERR( DMPlexTriangleSetOptions(self.dm, copts) )

    def setTetGenOptions(self, opts):
        """DMPlexTetgenSetOptions - Set the options used for the Tetgen mesh generator

        Not collective.

        Inputs Parameters:
        + dm - The `DMPlex` object
        - opts - The command line options

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexTriangleSetOptions`, `DMPlexGenerate`

        See Also
        --------
        petsc.DMPlexTetgenSetOptions

        """
        cdef const char *copts = NULL
        opts = str2bytes(opts, &copts)
        CHKERR( DMPlexTetgenSetOptions(self.dm, copts) )

    def markBoundaryFaces(self, label, value=None):
        """DMPlexMarkBoundaryFaces - Mark all faces on the boundary

        Not collective.

        Parameters
        ----------
        dm
            The original `DM`
        val
            The marker value, or `PETSC_DETERMINE` to use some value in the closure (or 1 if none are found)

        Returns
        -------
        label
            The `DMLabel` marking boundary faces with the given value

        Note:
        This function will use the point `PetscSF` from the input `DM` to exclude points on the partition boundary from being marked, unless the partition overlap is greater than zero. If you also wish to mark the partition boundary, you can use `DMSetPointSF` to temporarily set it to ``None``, and then reset it to the original object after the call.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMLabelCreate`, `DMCreateLabel`

        See Also
        --------
        petsc.DMPlexMarkBoundaryFaces

        """
        cdef PetscInt ival = PETSC_DETERMINE
        if value is not None: ival = asInt(value)
        if not self.hasLabel(label):
            self.createLabel(label)
        cdef const char *cval = NULL
        label = str2bytes(label, &cval)
        cdef PetscDMLabel clbl = NULL
        CHKERR( DMGetLabel(self.dm, cval, &clbl) )
        CHKERR( DMPlexMarkBoundaryFaces(self.dm, ival, clbl) )

    def labelComplete(self, DMLabel label):
        """DMPlexLabelComplete - Starting with a label marking points on a surface, we add the transitive closure to the surface

        Parameters
        ----------
        dm
            The `DM`
        label
            A `DMLabel` marking the surface points

        Returns
        -------
        label
            A `DMLabel` marking all surface points in the transitive closure

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexLabelCohesiveComplete`

        See Also
        --------
        petsc.DMPlexLabelComplete

        """
        CHKERR( DMPlexLabelComplete(self.dm, label.dmlabel) )

    def labelCohesiveComplete(self, DMLabel label, DMLabel bdlabel, bdvalue, flip, DMPlex subdm):
        """DMPlexLabelCohesiveComplete - Starting with a label marking points on an internal surface, we add all other mesh pieces
        to complete the surface

        Parameters
        ----------
        dm
            The `DM`
        label
            A `DMLabel` marking the surface
        blabel
            A `DMLabel` marking the vertices on the boundary which will not be duplicated, or ``None`` to find them automatically
        bvalue
            Value of `DMLabel` marking the vertices on the boundary
        flip
            Flag to flip the submesh normal and replace points on the other side
        subdm
            The `DM` associated with the label, or ``None``

        Returns
        -------
        label
            A `DMLabel` marking all surface points

        Note:
        The vertices in blabel are called "unsplit" in the terminology from hybrid cell creation.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexConstructCohesiveCells`, `DMPlexLabelComplete`

        See Also
        --------
        petsc.DMPlexLabelCohesiveComplete

        """
        cdef PetscBool flg = flip
        cdef PetscInt  val = asInt(bdvalue)
        CHKERR( DMPlexLabelCohesiveComplete(self.dm, label.dmlabel, bdlabel.dmlabel, val, flg, subdm.dm) )

    def setAdjacencyUseAnchors(self, useAnchors=True):
        """DMPlexSetAdjacencyUseAnchors - Define adjacency in the mesh using the point-to-point constraints.

        Parameters
        ----------
        dm
            The `DM` object
        useAnchors
            Flag to use the constraints. If PETSC_TRUE, then constrained points are omitted from DMPlexGetAdjacency(), and their anchor points appear in their place.

        .seealso: `DMPlex`, `DMGetAdjacency`, `DMSetAdjacency`, `DMPlexDistribute`, `DMPlexPreallocateOperator`, `DMPlex.setAnchors`

        See Also
        --------
        petsc.DMPlexSetAdjacencyUseAnchors

        """
        cdef PetscBool flag = useAnchors
        CHKERR( DMPlexSetAdjacencyUseAnchors(self.dm, flag) )

    def getAdjacencyUseAnchors(self):
        """DMPlexGetAdjacencyUseAnchors - Query whether adjacency in the mesh uses the point-to-point constraints.

        Parameters
        ----------
        dm
            The `DM` object

        Returns
        -------
        useAnchors
            Flag to use the closure. If PETSC_TRUE, then constrained points are omitted from DMPlexGetAdjacency(), and their anchor points appear in their place.

        .seealso: `DMPlex`, `DMPlex.setAdjacencyUseAnchors`, `DMSetAdjacency`, `DMGetAdjacency`, `DMPlexDistribute`, `DMPlexPreallocateOperator`, `DMPlex.setAnchors`

        See Also
        --------
        petsc.DMPlexGetAdjacencyUseAnchors

        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( DMPlexGetAdjacencyUseAnchors(self.dm, &flag) )
        return toBool(flag)

    def getAdjacency(self, p):
        """DMPlexGetAdjacency - Return all points adjacent to the given point

        Parameters
        ----------
        dm
            The `DM` object
        p
            The point

        Returns
        -------
        adjSize
            The maximum size of `adj` if it is non-``None``, or `PETSC_DETERMINE`;
        on output the number of adjacent points
        adj
            Either ``None`` so that the array is allocated, or an existing array with size `adjSize`;
        on output contains the adjacent points

        Notes:
        The user must `PetscFree` the `adj` array if it was not passed in.

        .seealso: `DMPlex`, `DMSetAdjacency`, `DMPlexDistribute`, `DMCreateMatrix`, `DMPlexPreallocateOperator`

        See Also
        --------
        petsc.DMPlexGetAdjacency

        """
        cdef PetscInt cp = asInt(p)
        cdef PetscInt nadj = PETSC_DETERMINE
        cdef PetscInt *iadj = NULL
        CHKERR( DMPlexGetAdjacency(self.dm, cp, &nadj, &iadj) )
        try:
            adjacency = array_i(nadj, iadj)
        finally:
            CHKERR( PetscFree(iadj) )
        return adjacency

    def setPartitioner(self, Partitioner part):
        """DMPlexSetPartitioner - Set the mesh partitioner

        logically Collective

        Parameters
        ----------
        dm
            The `DM`
        part
            The partitioner

        Note:
        Any existing `PetscPartitioner` will be destroyed.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `PetscPartitioner`,`DMPlexDistribute`, `DMPlex.getPartitioner`, `PetscPartitionerCreate`

        See Also
        --------
        petsc.DMPlexSetPartitioner

        """
        CHKERR( DMPlexSetPartitioner(self.dm, part.part) )

    def getPartitioner(self):
        """DMPlexGetPartitioner - Get the mesh partitioner

        Not collective.

        Parameters
        ----------
        dm
            The `DM`

        Returns
        -------
        part
            The `PetscPartitioner`

        Note:
        This gets a borrowed reference, so the user should not destroy this `PetscPartitioner`.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `PetscPartitioner`, `PetscSection`, `DMPlexDistribute`, `DMPlex.setPartitioner`, `PetscPartitionerDMPlexPartition`, `PetscPartitionerCreate`

        See Also
        --------
        petsc.DMPlexGetPartitioner

        """
        cdef Partitioner part = Partitioner()
        CHKERR( DMPlexGetPartitioner(self.dm, &part.part) )
        PetscINCREF(part.obj)
        return part

    def rebalanceSharedPoints(self, entityDepth=0, useInitialGuess=True, parallel=True):
        """DMPlexRebalanceSharedPoints - Redistribute points in the plex that are shared in order to achieve better balancing. This routine updates the `PointSF` of the `DM` inplace.

        Input parameters:
        + dm               - The `DMPlex` object.
        . entityDepth      - depth of the entity to balance (0 -> balance vertices).
        . useInitialGuess  - whether to use the current distribution as initial guess (only used by ParMETIS).
        - parallel         - whether to use ParMETIS and do the partition in parallel or whether to gather the graph onto a single process and use METIS.

        Output parameter:
        . success          - whether the graph partitioning was successful or not, optional. Unsuccessful simply means no change to the partitioning

        Options Database Keys:
        +  -dm_plex_rebalance_shared_points_parmetis - Use ParMetis instead of Metis for the partitioner
        .  -dm_plex_rebalance_shared_points_use_initial_guess - Use current partition to bootstrap ParMetis partition
        .  -dm_plex_rebalance_shared_points_use_mat_partitioning - Use the MatPartitioning object to perform the partition, the prefix for those operations is -dm_plex_rebalance_shared_points_
        -  -dm_plex_rebalance_shared_points_monitor - Monitor the shared points rebalance process

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexDistribute`, `DMPlexCreateOverlap`

        See Also
        --------
        petsc.DMPlexRebalanceSharedPoints

        """
        cdef PetscInt centityDepth = asInt(entityDepth)
        cdef PetscBool cuseInitialGuess = asBool(useInitialGuess)
        cdef PetscBool cparallel = asBool(parallel)
        cdef PetscBool csuccess = PETSC_FALSE
        CHKERR( DMPlexRebalanceSharedPoints(self.dm, centityDepth, cuseInitialGuess, cparallel, &csuccess) )
        return toBool(csuccess)

    def distribute(self, overlap=0):
        """DMPlexDistribute - Distributes the mesh and any associated sections.

        Collective.

        Parameters
        ----------
        dm
            The original `DMPlex` object
        overlap
            The overlap of partitions, 0 is the default

        Returns
        -------
        sf
            The `PetscSF` used for point distribution, or ``None`` if not needed
        dmParallel
            The distributed `DMPlex` object

        Note:
        If the mesh was not distributed, the output `dmParallel` will be ``None``.

        The user can control the definition of adjacency for the mesh using `DMSetAdjacency`. They should choose the combination appropriate for the function
        representation on the mesh.

        .seealso: `DMPlex`, `DM`, `DMPlex.create`, `DMSetAdjacency`, `DMPlex.getOverlap`

        See Also
        --------
        petsc.DMPlexDistribute

        """
        cdef PetscDM dmParallel = NULL
        cdef PetscInt coverlap = asInt(overlap)
        cdef SF sf = SF()
        CHKERR( DMPlexDistribute(self.dm, coverlap, &sf.sf, &dmParallel) )
        if dmParallel != NULL:
            PetscCLEAR(self.obj); self.dm = dmParallel
            return sf

    def distributeOverlap(self, overlap=0):
        """DMPlexDistributeOverlap - Add partition overlap to a distributed non-overlapping `DM`.

        Collective.

        Parameters
        ----------
        dm
            The non-overlapping distributed `DMPlex` object
        overlap
            The overlap of partitions (the same on all ranks)

        Returns
        -------
        sf
            The `PetscSF` used for point distribution
        dmOverlap
            The overlapping distributed `DMPlex` object, or ``None``

        Options Database Keys:
        + -dm_plex_overlap_labels <name1,name2,...> - List of overlap label names
        . -dm_plex_overlap_values <int1,int2,...>   - List of overlap label values
        . -dm_plex_overlap_exclude_label <name>     - Label used to exclude points from overlap
        - -dm_plex_overlap_exclude_value <int>      - Label value used to exclude points from overlap

        Notes:
        If the mesh was not distributed, the return value is ``None``.

        The user can control the definition of adjacency for the mesh using `DMSetAdjacency`. They should choose the combination appropriate for the function
        representation on the mesh.

        .seealso: `DMPlex`, `PetscSF`, `DM`, `DMPlex.create`, `DMSetAdjacency`, `DMPlexDistribute`, `DMPlexCreateOverlapLabel`, `DMPlex.getOverlap`

        See Also
        --------
        petsc.DMPlexDistributeOverlap

        """
        cdef PetscInt coverlap = asInt(overlap)
        cdef SF sf = SF()
        cdef PetscDM dmOverlap = NULL
        CHKERR( DMPlexDistributeOverlap(self.dm, coverlap,
                                        &sf.sf, &dmOverlap) )
        PetscCLEAR(self.obj); self.dm = dmOverlap
        return sf

    def isDistributed(self):
        """DMPlexIsDistributed - Find out whether this `DM` is distributed, i.e. more than one rank owns some points.

        Collective.

        Parameters
        ----------
        dm
            The `DM` object

        Returns
        -------
        distributed
            Flag whether the `DM` is distributed

        Notes:
        This currently finds out whether at least two ranks have any DAG points.
        This involves `MPI_Allreduce` with one integer.
        The result is currently not stashed so every call to this routine involves this global communication.

        .seealso: `DMPlex`, `DMPlex`, `DMPlexDistribute`, `DMPlex.getOverlap`, `DMPlexIsInterpolated`

        See Also
        --------
        petsc.DMPlexIsDistributed

        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( DMPlexIsDistributed(self.dm, &flag) )
        return toBool(flag)

    def isSimplex(self):
        """DMPlexIsSimplex - Is the first cell in this mesh a simplex?

        Parameters
        ----------
        dm
            The `DMPlex` object

        Returns
        -------
        simplex
            Flag checking for a simplex

        Note:
        This just gives the first range of cells found. If the mesh has several cell types, it will only give the first.
        If the mesh has no cells, this returns ``False``.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getSimplexOrBoxCells`, `DMPlex.getCellType`, `DMPlex.getHeightStratum`, `DMPolytopeTypeGetNumVertices`

        See Also
        --------
        petsc.DMPlexIsSimplex

        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( DMPlexIsSimplex(self.dm, &flag) )
        return toBool(flag)

    def distributeGetDefault(self):
        """DMPlexDistributeGetDefault - Get flag indicating whether the `DM` should be distributed by default

        Not collective.

        Parameters
        ----------
        dm
            The `DM`

        Returns
        -------
        dist
            Flag for distribution

        .seealso: `DMPlex`, `DM`, `DMPlexDistributeSetDefault`, `DMPlexDistribute`

        See Also
        --------
        petsc.DMPlexDistributeGetDefault

        """
        cdef PetscBool dist = PETSC_FALSE
        CHKERR( DMPlexDistributeGetDefault(self.dm, &dist) )
        return toBool(dist)

    def distributeSetDefault(self, flag):
        """DMPlexDistributeSetDefault - Set flag indicating whether the `DM` should be distributed by default

        Logically collective.

        Parameters
        ----------
        dm
            The `DM`
        dist
            Flag for distribution

        .seealso: `DMPlex`, `DMPlexDistributeGetDefault`, `DMPlexDistribute`

        See Also
        --------
        petsc.DMPlexDistributeSetDefault

        """
        cdef PetscBool dist = asBool(flag)
        CHKERR( DMPlexDistributeSetDefault(self.dm, dist) )
        return

    def distributionSetName(self, name):
        """DMPlexDistributionSetName - Set the name of the specific parallel distribution

        Parameters
        ----------
        dm
            The `DM`
        name
            The name of the specific parallel distribution

        Note:
        If distribution name is set when saving, `DMPlexTopologyView` saves the plex's
        parallel distribution (i.e., partition, ownership, and local ordering of points) under
        this name. Conversely, if distribution name is set when loading, `DMPlexTopologyLoad`
        loads the parallel distribution stored in file under this name.

        .seealso: `DMPlex`, `DMPlexDistributionGetName`, `DMPlexTopologyView`, `DMPlexTopologyLoad`

        See Also
        --------
        petsc.DMPlexDistributionSetName

        """
        cdef const char *cname = NULL
        if name is not None:
            name = str2bytes(name, &cname)
        CHKERR( DMPlexDistributionSetName(self.dm, cname) )

    def distributionGetName(self):
        """DMPlexDistributionGetName - Retrieve the name of the specific parallel distribution

        Parameters
        ----------
        dm
            The `DM`

        Returns
        -------
        name
            The name of the specific parallel distribution

        Note:
        If distribution name is set when saving, `DMPlexTopologyView` saves the plex's
        parallel distribution (i.e., partition, ownership, and local ordering of points) under
        this name. Conversely, if distribution name is set when loading, `DMPlexTopologyLoad`
        loads the parallel distribution stored in file under this name.

        .seealso: `DMPlex`, `DMPlexDistributionSetName`, `DMPlexTopologyView`, `DMPlexTopologyLoad`

        See Also
        --------
        petsc.DMPlexDistributionGetName

        """
        cdef const char *cname = NULL
        CHKERR( DMPlexDistributionGetName(self.dm, &cname) )
        return bytes2str(cname)

    def isSimplex(self):
        """DMPlexIsSimplex - Is the first cell in this mesh a simplex?

        Parameters
        ----------
        dm
            The `DMPlex` object

        Returns
        -------
        simplex
            Flag checking for a simplex

        Note:
        This just gives the first range of cells found. If the mesh has several cell types, it will only give the first.
        If the mesh has no cells, this returns ``False``.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getSimplexOrBoxCells`, `DMPlex.getCellType`, `DMPlex.getHeightStratum`, `DMPolytopeTypeGetNumVertices`

        See Also
        --------
        petsc.DMPlexIsSimplex

        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( DMPlexIsSimplex(self.dm, &flag) )
        return toBool(flag)

    def interpolate(self):
        """DMPlexInterpolate - Take in a cell-vertex mesh and return one with all intermediate faces, edges, etc.

        Collective.

        Parameters
        ----------
        dm
            The `DMPlex` object with only cells and vertices

        Returns
        -------
        dmInt
            The complete `DMPlex` object

        Note:
        Labels and coordinates are copied.

        Developer Note:
        It sets plex->interpolated = `DMPLEX_INTERPOLATED_FULL`.

        .seealso: `DMPlex`, `DMPlexUninterpolate`, `DMPlexCreateFromCellListPetsc`, `DMPlexCopyCoordinates`

        See Also
        --------
        petsc.DMPlexInterpolate

        """
        cdef PetscDM newdm = NULL
        CHKERR( DMPlexInterpolate(self.dm, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm

    def uninterpolate(self):
        """DMPlexUninterpolate - Take in a mesh with all intermediate faces, edges, etc. and return a cell-vertex mesh

        Collective.

        Parameters
        ----------
        dm
            The complete `DMPlex` object

        Returns
        -------
        dmUnint
            The `DMPlex` object with only cells and vertices

        .seealso: `DMPlex`, `DMPlexInterpolate`, `DMPlexCreateFromCellListPetsc`, `DMPlexCopyCoordinates`

        See Also
        --------
        petsc.DMPlexUninterpolate

        """
        cdef PetscDM newdm = NULL
        CHKERR( DMPlexUninterpolate(self.dm, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm

    def distributeField(self, SF sf, Section sec, Vec vec,
                        Section newsec=None, Vec newvec=None):
        """DMPlexDistributeField - Distribute field data to match a given `PetscSF`, usually the `PetscSF` from mesh distribution

        Collective.

        Parameters
        ----------
        dm
            The `DMPlex` object
        pointSF
            The `PetscSF` describing the communication pattern
        originalSection
            The `PetscSection` for existing data layout
        originalVec
            The existing data in a local vector

        Returns
        -------
        newSection
            The `PetscSF` describing the new data layout
        newVec
            The new data in a local vector

        .seealso: `DMPlex`, `DMPlexDistribute`, `DMPlexDistributeFieldIS`, `DMPlexDistributeData`

        See Also
        --------
        petsc.DMPlexDistributeField

        """
        cdef MPI_Comm ccomm = MPI_COMM_NULL
        if newsec is None: newsec = Section()
        if newvec is None: newvec = Vec()
        if newsec.sec == NULL:
            CHKERR( PetscObjectGetComm(<PetscObject>sec.sec, &ccomm) )
            CHKERR( PetscSectionCreate(ccomm, &newsec.sec) )
        if newvec.vec == NULL:
            CHKERR( PetscObjectGetComm(<PetscObject>vec.vec, &ccomm) )
            CHKERR( VecCreate(ccomm, &newvec.vec) )
        CHKERR( DMPlexDistributeField(self.dm, sf.sf,
                                      sec.sec, vec.vec,
                                      newsec.sec, newvec.vec))
        return (newsec, newvec)

    def getMinRadius(self):
        """DMPlexGetMinRadius - Returns the minimum distance from any cell centroid to a face

        Not collective.

        Parameters
        ----------
        dm
            the `DMPlex`

        Returns
        -------
        minradius
            the minimum cell radius

        .seealso: `DMPlex`, `DMGetCoordinates`

        See Also
        --------
        petsc.DMPlexGetMinRadius

        """
        cdef PetscReal cminradius = 0.
        CHKERR( DMPlexGetMinRadius(self.dm, &cminradius))
        return asReal(cminradius)

    def createCoarsePointIS(self):
        """DMPlexCreateCoarsePointIS - Creates an `IS` covering the coarse `DM` chart with the fine points as data

        Collective.

        Parameters
        ----------
        dm
            The coarse `DM`

        Returns
        -------
        fpointIS
            The `IS` of all the fine points which exist in the original coarse mesh

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `IS`, `DMRefine`, `DMPlex.setRefinementUniform`, `DMPlex.getSubpointIS`

        See Also
        --------
        petsc.DMPlexCreateCoarsePointIS

        """
        cdef IS fpoint = IS()
        CHKERR( DMPlexCreateCoarsePointIS(self.dm, &fpoint.iset) )
        return fpoint

    def createSection(self, numComp, numDof,
                      bcField=None, bcComps=None, bcPoints=None,
                      IS perm=None):
        """DMPlexCreateSection - Create a `PetscSection` based upon the dof layout specification provided.

        Not collective.

        Parameters
        ----------
        dm
            The `DMPlex` object
        label
            The label indicating the mesh support of each field, or NULL for the whole mesh
        numComp
            An array of size numFields that holds the number of components for each field
        numDof
            An array of size numFields*(dim+1) which holds the number of dof for each field on a mesh piece of dimension d
        numBC
            The number of boundary conditions
        bcField
            An array of size numBC giving the field number for each boundary condition
        bcComps
            [Optional] An array of size numBC giving an `IS` holding the field components to which each boundary condition applies
        bcPoints
            An array of size numBC giving an `IS` holding the `DMPlex` points to which each boundary condition applies
        perm
            Optional permutation of the chart, or NULL

        Returns
        -------
        section
            The `PetscSection` object

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.create`, `PetscSectionCreate`, `PetscSectionSetPermutation`

        See Also
        --------
        petsc.DMPlexCreateSection

        """
        # topological dimension
        cdef PetscInt dim = 0
        CHKERR( DMGetDimension(self.dm, &dim) )
        # components and DOFs
        cdef PetscInt ncomp = 0, ndof = 0
        cdef PetscInt *icomp = NULL, *idof = NULL
        numComp = iarray_i(numComp, &ncomp, &icomp)
        numDof  = iarray_i(numDof, &ndof, &idof)
        assert ndof == ncomp*(dim+1)
        # boundary conditions
        cdef PetscInt nbc = 0, i = 0
        cdef PetscInt *bcfield = NULL
        cdef PetscIS *bccomps  = NULL
        cdef PetscIS *bcpoints = NULL
        if bcField is not None:
            bcField = iarray_i(bcField, &nbc, &bcfield)
            if bcComps is not None:
                bcComps = list(bcComps)
                assert len(bcComps) == nbc
                tmp1 = oarray_p(empty_p(nbc), NULL, <void**>&bccomps)
                for i from 0 <= i < nbc:
                    bccomps[i] = (<IS?>bcComps[<Py_ssize_t>i]).iset
            if bcPoints is not None:
                bcPoints = list(bcPoints)
                assert len(bcPoints) == nbc
                tmp2 = oarray_p(empty_p(nbc), NULL, <void**>&bcpoints)
                for i from 0 <= i < nbc:
                    bcpoints[i] = (<IS?>bcPoints[<Py_ssize_t>i]).iset
            else:
                raise ValueError("bcPoints is a required argument")
        else:
            assert bcComps  is None
            assert bcPoints is None
        # optional chart permutations
        cdef PetscIS cperm = NULL
        if perm is not None: cperm = perm.iset
        # create section
        cdef Section sec = Section()
        CHKERR( DMPlexCreateSection(self.dm, NULL, icomp, idof,
                                    nbc, bcfield, bccomps, bcpoints,
                                    cperm, &sec.sec) )
        return sec

    def getPointLocal(self, point):
        """DMPlexGetPointLocal - get location of point data in local `Vec`

        Not collective.

        Parameters
        ----------
        dm
            `DM` defining the topological space
        point
            topological point

        Returns
        -------
        start
            start of point data
        end
            end of point data

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getPointLocalField`, `DMGetLocalSection`, `PetscSectionGetOffset`, `PetscSectionGetDof`, `DMPlexPointLocalRead`, `DMPlexPointLocalRead`, `DMPlexPointLocalRef`

        See Also
        --------
        petsc.DMPlexGetPointLocal

        """
        cdef PetscInt start = 0, end = 0
        cdef PetscInt cpoint = asInt(point)
        CHKERR( DMPlexGetPointLocal(self.dm, cpoint, &start, &end) )
        return toInt(start), toInt(end)

    def getPointLocalField(self, point, field):
        """DMPlexGetPointLocalField - get location of point field data in local Vec

        Not collective.

        Parameters
        ----------
        dm
            `DM` defining the topological space
        point
            topological point
        field
            the field number

        Returns
        -------
        start
            start of point data
        end
            end of point data

        Note:
        This is a half open interval [start, end)

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getPointLocal`, `DMGetLocalSection`, `PetscSectionGetOffset`, `PetscSectionGetDof`, `DMPlexPointLocalRead`, `DMPlexPointLocalRead`, `DMPlexPointLocalRef`

        See Also
        --------
        petsc.DMPlexGetPointLocalField

        """
        cdef PetscInt start = 0, end = 0
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        CHKERR( DMPlexGetPointLocalField(self.dm, cpoint, cfield, &start, &end) )
        return toInt(start), toInt(end)

    def getPointGlobal(self, point):
        """DMPlexGetPointGlobal - get location of point data in global Vec

        Not collective.

        Parameters
        ----------
        dm
            `DM` defining the topological space
        point
            topological point

        Returns
        -------
        start
            start of point data; returns -(globalStart+1) if point is not owned
        end
            end of point data; returns -(globalEnd+1) if point is not owned

        Note:
        This is a half open interval [start, end)

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getPointGlobalField`, `DMGetLocalSection`, `PetscSectionGetOffset`, `PetscSectionGetDof`, `DMPlexPointGlobalRead`, `DMPlex.getPointLocal`, `DMPlexPointGlobalRead`, `DMPlexPointGlobalRef`

        See Also
        --------
        petsc.DMPlexGetPointGlobal

        """
        cdef PetscInt start = 0, end = 0
        cdef PetscInt cpoint = asInt(point)
        CHKERR( DMPlexGetPointGlobal(self.dm, cpoint, &start, &end) )
        return toInt(start), toInt(end)

    def getPointGlobalField(self, point, field):
        """DMPlexGetPointGlobalField - get location of point field data in global `Vec`

        Not collective.

        Parameters
        ----------
        dm
            `DM` defining the topological space
        point
            topological point
        field
            the field number

        Returns
        -------
        start
            start of point data; returns -(globalStart+1) if point is not owned
        end
            end of point data; returns -(globalEnd+1) if point is not owned

        Note:
        This is a half open interval [start, end)

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlex.getPointGlobal`, `DMGetLocalSection`, `PetscSectionGetOffset`, `PetscSectionGetDof`, `DMPlexPointGlobalRead`, `DMPlex.getPointLocal`, `DMPlexPointGlobalRead`, `DMPlexPointGlobalRef`

        See Also
        --------
        petsc.DMPlexGetPointGlobalField

        """
        cdef PetscInt start = 0, end = 0
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        CHKERR( DMPlexGetPointGlobalField(self.dm, cpoint, cfield, &start, &end) )
        return toInt(start), toInt(end)

    def createClosureIndex(self, Section sec or None):
        """DMPlexCreateClosureIndex - Calculate an index for the given `PetscSection` for the closure operation on the `DM`

        Not collective.

        Parameters
        ----------
        dm
            The `DM`
        section
            The section describing the layout in the local vector, or NULL to use the default section

        Note:
        This should greatly improve the performance of the closure operations, at the cost of additional memory.

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `PetscSection`, `DMPlexVecGetClosure`, `DMPlexVecRestoreClosure`, `DMPlexVecSetClosure`, `DMPlexMatSetClosure`

        See Also
        --------
        petsc.DMPlexCreateClosureIndex

        """
        cdef PetscSection csec = sec.sec if sec is not None else NULL
        CHKERR( DMPlexCreateClosureIndex(self.dm, csec) )

    #

    def setRefinementUniform(self, refinementUniform=True):
        """DMPlexSetRefinementUniform - Set the flag for uniform refinement

        Parameters
        ----------
        dm
            The `DM`
        refinementUniform
            The flag for uniform refinement

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMRefine`, `DMPlex.getRefinementUniform`, `DMPlex.getRefinementLimit`, `DMPlex.setRefinementLimit`

        See Also
        --------
        petsc.DMPlexSetRefinementUniform

        """
        cdef PetscBool flag = refinementUniform
        CHKERR( DMPlexSetRefinementUniform(self.dm, flag) )

    def getRefinementUniform(self):
        """DMPlexGetRefinementUniform - Retrieve the flag for uniform refinement

        Parameters
        ----------
        dm
            The `DM`

        Returns
        -------
        refinementUniform
            The flag for uniform refinement

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMRefine`, `DMPlex.setRefinementUniform`, `DMPlex.getRefinementLimit`, `DMPlex.setRefinementLimit`

        See Also
        --------
        petsc.DMPlexGetRefinementUniform

        """
        cdef PetscBool flag = PETSC_FALSE
        CHKERR( DMPlexGetRefinementUniform(self.dm, &flag) )
        return toBool(flag)

    def setRefinementLimit(self, refinementLimit):
        """DMPlexSetRefinementLimit - Set the maximum cell volume for refinement

        Parameters
        ----------
        dm
            The `DM`
        refinementLimit
            The maximum cell volume in the refined mesh

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMRefine`, `DMPlex.getRefinementLimit`, `DMPlex.getRefinementUniform`, `DMPlex.setRefinementUniform`

        See Also
        --------
        petsc.DMPlexSetRefinementLimit

        """
        cdef PetscReal rval = asReal(refinementLimit)
        CHKERR( DMPlexSetRefinementLimit(self.dm, rval) )

    def getRefinementLimit(self):
        """DMPlexGetRefinementLimit - Retrieve the maximum cell volume for refinement

        Parameters
        ----------
        dm
            The `DM`

        Returns
        -------
        refinementLimit
            The maximum cell volume in the refined mesh

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMRefine`, `DMPlex.setRefinementLimit`, `DMPlex.getRefinementUniform`, `DMPlex.setRefinementUniform`

        See Also
        --------
        petsc.DMPlexGetRefinementLimit

        """
        cdef PetscReal rval = 0.0
        CHKERR( DMPlexGetRefinementLimit(self.dm, &rval) )
        return toReal(rval)

    def getOrdering(self, otype):
        """DMPlexGetOrdering - Calculate a reordering of the mesh

        Collective.

        Parameters
        ----------
        dm
            The DMPlex object
        otype
            type of reordering, see `MatOrderingType`
        label
            [Optional] Label used to segregate ordering into sets, or ``None``

        Returns
        -------
        perm
            The point permutation as an `IS`, `perm`[old point number] = new point number

        .seealso: `DMPlex`, `DMPlexPermute`, `MatOrderingType`, `MatGetOrdering`

        See Also
        --------
        petsc.DMPlexGetOrdering

        """
        cdef PetscMatOrderingType cval = NULL
        cdef PetscDMLabel label = NULL
        otype = str2bytes(otype, &cval)
        cdef IS perm = IS()
        CHKERR( DMPlexGetOrdering(self.dm, cval, label, &perm.iset) )
        return perm

    def permute(self, IS perm):
        """DMPlexPermute - Reorder the mesh according to the input permutation

        Collective.

        Parameters
        ----------
        dm
            The `DMPlex` object
        perm
            The point permutation, `perm`[old point number] = new point number

        Returns
        -------
        pdm
            The permuted `DM`

        .seealso: `DMPlex`, `MatPermute`

        See Also
        --------
        petsc.DMPlexPermute

        """
        cdef DMPlex dm = <DMPlex>type(self)()
        CHKERR( DMPlexPermute(self.dm, perm.iset, &dm.dm) )
        return dm

    def reorderGetDefault(self):
        """DMPlexReorderGetDefault - Get flag indicating whether the DM should be reordered by default

        Not collective.

        Parameters
        ----------
        dm
            The `DM`

        Returns
        -------
        reorder
            Flag for reordering

        .seealso: `DMPlexReorderSetDefault`

        See Also
        --------
        petsc.DMPlexReorderGetDefault

        """
        cdef PetscDMPlexReorderDefaultFlag reorder = DMPLEX_REORDER_DEFAULT_NOTSET
        CHKERR( DMPlexReorderGetDefault(self.dm, &reorder) )
        return reorder

    def reorderSetDefault(self, flag):
        """DMPlexReorderSetDefault - Set flag indicating whether the DM should be reordered by default

        Logically collective.

        Parameters
        ----------
        dm
            The `DM`
        reorder
            Flag for reordering

        .seealso: `DMPlexReorderGetDefault`

        See Also
        --------
        petsc.DMPlexReorderSetDefault

        """
        cdef PetscDMPlexReorderDefaultFlag reorder = flag
        CHKERR( DMPlexReorderSetDefault(self.dm, reorder) )
        return

    #

    def computeCellGeometryFVM(self, cell):
        """DMPlexComputeCellGeometryFVM - Compute the volume for a given cell

        Collective.

        Parameters
        ----------
        dm
            the `DMPlex`
        cell
            the cell

        Returns
        -------
        volume
            the cell volume
        centroid
            the cell centroid
        normal
            the cell normal, if appropriate

        .seealso: `DMPlex`, `DMGetCoordinateSection`, `DMGetCoordinates`

        See Also
        --------
        petsc.DMPlexComputeCellGeometryFVM

        """
        cdef PetscInt cdim = 0
        cdef PetscInt ccell = asInt(cell)
        CHKERR( DMGetCoordinateDim(self.dm, &cdim) )
        cdef PetscReal vol = 0, centroid[3], normal[3]
        CHKERR( DMPlexComputeCellGeometryFVM(self.dm, ccell, &vol, centroid, normal) )
        return (toReal(vol), array_r(cdim, centroid), array_r(cdim, normal))

    def constructGhostCells(self, labelName=None):
        """DMPlexConstructGhostCells - Construct ghost cells which connect to every boundary face

        Collective.

        Parameters
        ----------
        dm
            The original `DM`
        labelName
            The label specifying the boundary faces, or "Face Sets" if this is ``None``

        Returns
        -------
        numGhostCells
            The number of ghost cells added to the `DM`
        dmGhosted
            The new `DM`

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMCreate`

        See Also
        --------
        petsc.DMPlexConstructGhostCells

        """
        cdef const char *cname = NULL
        labelName = str2bytes(labelName, &cname)
        cdef PetscInt numGhostCells = 0
        cdef PetscDM dmGhosted = NULL
        CHKERR( DMPlexConstructGhostCells(self.dm, cname, &numGhostCells, &dmGhosted))
        PetscCLEAR(self.obj); self.dm = dmGhosted
        return toInt(numGhostCells)

    # Metric

    def metricSetFromOptions(self):
        CHKERR( DMPlexMetricSetFromOptions(self.dm) )

    def metricSetUniform(self, uniform):
        """DMPlexMetricSetUniform - Record whether a metric is uniform

        Input parameters:
        + dm      - The DM
        - uniform - Is the metric uniform?

        .seealso: `DMPlexMetricIsUniform`, `DMPlexMetricSetIsotropic`, `DMPlexMetricSetRestrictAnisotropyFirst`

        See Also
        --------
        petsc.DMPlexMetricSetUniform

        """
        cdef PetscBool bval = asBool(uniform)
        CHKERR( DMPlexMetricSetUniform(self.dm, bval) )

    def metricIsUniform(self):
        """DMPlexMetricIsUniform - Is a metric uniform?

        Input parameters:
        . dm      - The DM

        Output parameters:
        . uniform - Is the metric uniform?

        .seealso: `DMPlexMetricSetUniform`, `DMPlexMetricIsIsotropic`, `DMPlexMetricRestrictAnisotropyFirst`

        See Also
        --------
        petsc.DMPlexMetricIsUniform

        """
        cdef PetscBool uniform = PETSC_FALSE
        CHKERR( DMPlexMetricIsUniform(self.dm, &uniform) )
        return toBool(uniform)

    def metricSetIsotropic(self, isotropic):
        """DMPlexMetricSetIsotropic - Record whether a metric is isotropic

        Input parameters:
        + dm        - The DM
        - isotropic - Is the metric isotropic?

        .seealso: `DMPlexMetricIsIsotropic`, `DMPlexMetricSetUniform`, `DMPlexMetricSetRestrictAnisotropyFirst`

        See Also
        --------
        petsc.DMPlexMetricSetIsotropic

        """
        cdef PetscBool bval = asBool(isotropic)
        CHKERR( DMPlexMetricSetIsotropic(self.dm, bval) )

    def metricIsIsotropic(self):
        """DMPlexMetricIsIsotropic - Is a metric isotropic?

        Input parameters:
        . dm        - The DM

        Output parameters:
        . isotropic - Is the metric isotropic?

        .seealso: `DMPlexMetricSetIsotropic`, `DMPlexMetricIsUniform`, `DMPlexMetricRestrictAnisotropyFirst`

        See Also
        --------
        petsc.DMPlexMetricIsIsotropic

        """
        cdef PetscBool isotropic = PETSC_FALSE
        CHKERR( DMPlexMetricIsIsotropic(self.dm, &isotropic) )
        return toBool(isotropic)

    def metricSetRestrictAnisotropyFirst(self, restrictAnisotropyFirst):
        """DMPlexMetricSetRestrictAnisotropyFirst - Record whether anisotropy should be restricted before normalization

        Input parameters:
        + dm                      - The DM
        - restrictAnisotropyFirst - Should anisotropy be normalized first?

        .seealso: `DMPlexMetricSetIsotropic`, `DMPlexMetricRestrictAnisotropyFirst`

        See Also
        --------
        petsc.DMPlexMetricSetRestrictAnisotropyFirst

        """
        cdef PetscBool bval = asBool(restrictAnisotropyFirst)
        CHKERR( DMPlexMetricSetRestrictAnisotropyFirst(self.dm, bval) )

    def metricRestrictAnisotropyFirst(self):
        """DMPlexMetricRestrictAnisotropyFirst - Is anisotropy restricted before normalization or after?

        Input parameters:
        . dm                      - The DM

        Output parameters:
        . restrictAnisotropyFirst - Is anisotropy be normalized first?

        .seealso: `DMPlexMetricIsIsotropic`, `DMPlexMetricSetRestrictAnisotropyFirst`

        See Also
        --------
        petsc.DMPlexMetricRestrictAnisotropyFirst

        """
        cdef PetscBool restrictAnisotropyFirst = PETSC_FALSE
        CHKERR( DMPlexMetricRestrictAnisotropyFirst(self.dm, &restrictAnisotropyFirst) )
        return toBool(restrictAnisotropyFirst)

    def metricSetNoInsertion(self, noInsert):
        """DMPlexMetricSetNoInsertion - Should node insertion and deletion be turned off?

        Input parameters:
        + dm       - The DM
        - noInsert - Should node insertion and deletion be turned off?

        .seealso: `DMPlexMetricNoInsertion`, `DMPlexMetricSetNoSwapping`, `DMPlexMetricSetNoMovement`, `DMPlexMetricSetNoSurf`

        See Also
        --------
        petsc.DMPlexMetricSetNoInsertion

        """
        cdef PetscBool bval = asBool(noInsert)
        CHKERR( DMPlexMetricSetNoInsertion(self.dm, bval) )

    def metricNoInsertion(self):
        """DMPlexMetricNoInsertion - Are node insertion and deletion turned off?

        Input parameters:
        . dm       - The DM

        Output parameters:
        . noInsert - Are node insertion and deletion turned off?

        .seealso: `DMPlexMetricSetNoInsertion`, `DMPlexMetricNoSwapping`, `DMPlexMetricNoMovement`, `DMPlexMetricNoSurf`

        See Also
        --------
        petsc.DMPlexMetricNoInsertion

        """
        cdef PetscBool noInsert = PETSC_FALSE
        CHKERR( DMPlexMetricNoInsertion(self.dm, &noInsert) )
        return toBool(noInsert)

    def metricSetNoSwapping(self, noSwap):
        """DMPlexMetricSetNoSwapping - Should facet swapping be turned off?

        Input parameters:
        + dm     - The DM
        - noSwap - Should facet swapping be turned off?

        .seealso: `DMPlexMetricNoSwapping`, `DMPlexMetricSetNoInsertion`, `DMPlexMetricSetNoMovement`, `DMPlexMetricSetNoSurf`

        See Also
        --------
        petsc.DMPlexMetricSetNoSwapping

        """
        cdef PetscBool bval = asBool(noSwap)
        CHKERR( DMPlexMetricSetNoSwapping(self.dm, bval) )

    def metricNoSwapping(self):
        """DMPlexMetricNoSwapping - Is facet swapping turned off?

        Input parameters:
        . dm     - The DM

        Output parameters:
        . noSwap - Is facet swapping turned off?

        .seealso: `DMPlexMetricSetNoSwapping`, `DMPlexMetricNoInsertion`, `DMPlexMetricNoMovement`, `DMPlexMetricNoSurf`

        See Also
        --------
        petsc.DMPlexMetricNoSwapping

        """
        cdef PetscBool noSwap = PETSC_FALSE
        CHKERR( DMPlexMetricNoSwapping(self.dm, &noSwap) )
        return toBool(noSwap)

    def metricSetNoMovement(self, noMove):
        """DMPlexMetricSetNoMovement - Should node movement be turned off?

        Input parameters:
        + dm     - The DM
        - noMove - Should node movement be turned off?

        .seealso: `DMPlexMetricNoMovement`, `DMPlexMetricSetNoInsertion`, `DMPlexMetricSetNoSwapping`, `DMPlexMetricSetNoSurf`

        See Also
        --------
        petsc.DMPlexMetricSetNoMovement

        """
        cdef PetscBool bval = asBool(noMove)
        CHKERR( DMPlexMetricSetNoMovement(self.dm, bval) )

    def metricNoMovement(self):
        """DMPlexMetricNoMovement - Is node movement turned off?

        Input parameters:
        . dm     - The DM

        Output parameters:
        . noMove - Is node movement turned off?

        .seealso: `DMPlexMetricSetNoMovement`, `DMPlexMetricNoInsertion`, `DMPlexMetricNoSwapping`, `DMPlexMetricNoSurf`

        See Also
        --------
        petsc.DMPlexMetricNoMovement

        """
        cdef PetscBool noMove = PETSC_FALSE
        CHKERR( DMPlexMetricNoMovement(self.dm, &noMove) )
        return toBool(noMove)

    def metricSetNoSurf(self, noSurf):
        """DMPlexMetricSetNoSurf - Should surface modification be turned off?

        Input parameters:
        + dm     - The DM
        - noSurf - Should surface modification be turned off?

        .seealso: `DMPlexMetricNoSurf`, `DMPlexMetricSetNoMovement`, `DMPlexMetricSetNoInsertion`, `DMPlexMetricSetNoSwapping`

        See Also
        --------
        petsc.DMPlexMetricSetNoSurf

        """
        cdef PetscBool bval = asBool(noSurf)
        CHKERR( DMPlexMetricSetNoSurf(self.dm, bval) )

    def metricNoSurf(self):
        """DMPlexMetricNoSurf - Is surface modification turned off?

        Input parameters:
        . dm     - The DM

        Output parameters:
        . noSurf - Is surface modification turned off?

        .seealso: `DMPlexMetricSetNoSurf`, `DMPlexMetricNoMovement`, `DMPlexMetricNoInsertion`, `DMPlexMetricNoSwapping`

        See Also
        --------
        petsc.DMPlexMetricNoSurf

        """
        cdef PetscBool noSurf = PETSC_FALSE
        CHKERR( DMPlexMetricNoSurf(self.dm, &noSurf) )
        return toBool(noSurf)

    def metricSetVerbosity(self, verbosity):
        """DMPlexMetricSetVerbosity - Set the verbosity of the mesh adaptation package

        Input parameters:
        + dm        - The DM
        - verbosity - The verbosity, where -1 is silent and 10 is maximum

        .seealso: `DMPlexMetricGetVerbosity`, `DMPlexMetricSetNumIterations`

        See Also
        --------
        petsc.DMPlexMetricSetVerbosity

        """
        cdef PetscInt ival = asInt(verbosity)
        CHKERR( DMPlexMetricSetVerbosity(self.dm, ival) )

    def metricGetVerbosity(self):
        """DMPlexMetricGetVerbosity - Get the verbosity of the mesh adaptation package

        Input parameters:
        . dm        - The DM

        Output parameters:
        . verbosity - The verbosity, where -1 is silent and 10 is maximum

        .seealso: `DMPlexMetricSetVerbosity`, `DMPlexMetricGetNumIterations`

        See Also
        --------
        petsc.DMPlexMetricGetVerbosity

        """
        cdef PetscInt verbosity = 0
        CHKERR( DMPlexMetricGetVerbosity(self.dm, &verbosity) )
        return toInt(verbosity)

    def metricSetNumIterations(self, numIter):
        """DMPlexMetricSetNumIterations - Set the number of parallel adaptation iterations

        Input parameters:
        + dm      - The DM
        - numIter - the number of parallel adaptation iterations

        .seealso: `DMPlexMetricSetVerbosity`, `DMPlexMetricGetNumIterations`

        See Also
        --------
        petsc.DMPlexMetricSetNumIterations

        """
        cdef PetscInt ival = asInt(numIter)
        CHKERR( DMPlexMetricSetNumIterations(self.dm, ival) )

    def metricGetNumIterations(self):
        """DMPlexMetricGetNumIterations - Get the number of parallel adaptation iterations

        Input parameters:
        . dm      - The DM

        Output parameters:
        . numIter - the number of parallel adaptation iterations

        .seealso: `DMPlexMetricSetNumIterations`, `DMPlexMetricGetVerbosity`

        See Also
        --------
        petsc.DMPlexMetricGetNumIterations

        """
        cdef PetscInt numIter = 0
        CHKERR( DMPlexMetricGetNumIterations(self.dm, &numIter) )
        return toInt(numIter)

    def metricSetMinimumMagnitude(self, h_min):
        """DMPlexMetricSetMinimumMagnitude - Set the minimum tolerated metric magnitude

        Input parameters:
        + dm    - The DM
        - h_min - The minimum tolerated metric magnitude

        .seealso: `DMPlexMetricGetMinimumMagnitude`, `DMPlexMetricSetMaximumMagnitude`

        See Also
        --------
        petsc.DMPlexMetricSetMinimumMagnitude

        """
        cdef PetscReal rval = asReal(h_min)
        CHKERR( DMPlexMetricSetMinimumMagnitude(self.dm, rval) )

    def metricGetMinimumMagnitude(self):
        """DMPlexMetricGetMinimumMagnitude - Get the minimum tolerated metric magnitude

        Input parameters:
        . dm    - The DM

        Output parameters:
        . h_min - The minimum tolerated metric magnitude

        .seealso: `DMPlexMetricSetMinimumMagnitude`, `DMPlexMetricGetMaximumMagnitude`

        See Also
        --------
        petsc.DMPlexMetricGetMinimumMagnitude

        """
        cdef PetscReal h_min = 0
        CHKERR( DMPlexMetricGetMinimumMagnitude(self.dm, &h_min) )
        return toReal(h_min)

    def metricSetMaximumMagnitude(self, h_max):
        """DMPlexMetricSetMaximumMagnitude - Set the maximum tolerated metric magnitude

        Input parameters:
        + dm    - The DM
        - h_max - The maximum tolerated metric magnitude

        .seealso: `DMPlexMetricGetMaximumMagnitude`, `DMPlexMetricSetMinimumMagnitude`

        See Also
        --------
        petsc.DMPlexMetricSetMaximumMagnitude

        """
        cdef PetscReal rval = asReal(h_max)
        CHKERR( DMPlexMetricSetMaximumMagnitude(self.dm, rval) )

    def metricGetMaximumMagnitude(self):
        """DMPlexMetricGetMaximumMagnitude - Get the maximum tolerated metric magnitude

        Input parameters:
        . dm    - The DM

        Output parameters:
        . h_max - The maximum tolerated metric magnitude

        .seealso: `DMPlexMetricSetMaximumMagnitude`, `DMPlexMetricGetMinimumMagnitude`

        See Also
        --------
        petsc.DMPlexMetricGetMaximumMagnitude

        """
        cdef PetscReal h_max = 0
        CHKERR( DMPlexMetricGetMaximumMagnitude(self.dm, &h_max) )
        return toReal(h_max)

    def metricSetMaximumAnisotropy(self, a_max):
        """DMPlexMetricSetMaximumAnisotropy - Set the maximum tolerated metric anisotropy

        Input parameters:
        + dm    - The DM
        - a_max - The maximum tolerated metric anisotropy

        .seealso: `DMPlexMetricGetMaximumAnisotropy`, `DMPlexMetricSetMaximumMagnitude`

        See Also
        --------
        petsc.DMPlexMetricSetMaximumAnisotropy

        """
        cdef PetscReal rval = asReal(a_max)
        CHKERR( DMPlexMetricSetMaximumAnisotropy(self.dm, rval) )

    def metricGetMaximumAnisotropy(self):
        """DMPlexMetricGetMaximumAnisotropy - Get the maximum tolerated metric anisotropy

        Input parameters:
        . dm    - The DM

        Output parameters:
        . a_max - The maximum tolerated metric anisotropy

        .seealso: `DMPlexMetricSetMaximumAnisotropy`, `DMPlexMetricGetMaximumMagnitude`

        See Also
        --------
        petsc.DMPlexMetricGetMaximumAnisotropy

        """
        cdef PetscReal a_max = 0
        CHKERR( DMPlexMetricGetMaximumAnisotropy(self.dm, &a_max) )
        return toReal(a_max)

    def metricSetTargetComplexity(self, targetComplexity):
        """DMPlexMetricSetTargetComplexity - Set the target metric complexity

        Input parameters:
        + dm               - The DM
        - targetComplexity - The target metric complexity

        .seealso: `DMPlexMetricGetTargetComplexity`, `DMPlexMetricSetNormalizationOrder`

        See Also
        --------
        petsc.DMPlexMetricSetTargetComplexity

        """
        cdef PetscReal rval = asReal(targetComplexity)
        CHKERR( DMPlexMetricSetTargetComplexity(self.dm, rval) )

    def metricGetTargetComplexity(self):
        """DMPlexMetricGetTargetComplexity - Get the target metric complexity

        Input parameters:
        . dm               - The DM

        Output parameters:
        . targetComplexity - The target metric complexity

        .seealso: `DMPlexMetricSetTargetComplexity`, `DMPlexMetricGetNormalizationOrder`

        See Also
        --------
        petsc.DMPlexMetricGetTargetComplexity

        """
        cdef PetscReal targetComplexity = 0
        CHKERR( DMPlexMetricGetTargetComplexity(self.dm, &targetComplexity) )
        return toReal(targetComplexity)

    def metricSetNormalizationOrder(self, p):
        """DMPlexMetricSetNormalizationOrder - Set the order p for L-p normalization

        Input parameters:
        + dm - The DM
        - p  - The normalization order

        .seealso: `DMPlexMetricGetNormalizationOrder`, `DMPlexMetricSetTargetComplexity`

        See Also
        --------
        petsc.DMPlexMetricSetNormalizationOrder

        """
        cdef PetscReal rval = asReal(p)
        CHKERR( DMPlexMetricSetNormalizationOrder(self.dm, rval) )

    def metricGetNormalizationOrder(self):
        """DMPlexMetricGetNormalizationOrder - Get the order p for L-p normalization

        Input parameters:
        . dm - The DM

        Output parameters:
        . p - The normalization order

        .seealso: `DMPlexMetricSetNormalizationOrder`, `DMPlexMetricGetTargetComplexity`

        See Also
        --------
        petsc.DMPlexMetricGetNormalizationOrder

        """
        cdef PetscReal p = 0
        CHKERR( DMPlexMetricGetNormalizationOrder(self.dm, &p) )
        return toReal(p)

    def metricSetGradationFactor(self, beta):
        """DMPlexMetricSetGradationFactor - Set the metric gradation factor

        Input parameters:
        + dm   - The DM
        - beta - The metric gradation factor

        .seealso: `DMPlexMetricGetGradationFactor`, `DMPlexMetricSetHausdorffNumber`

        See Also
        --------
        petsc.DMPlexMetricSetGradationFactor

        """
        cdef PetscReal rval = asReal(beta)
        CHKERR( DMPlexMetricSetGradationFactor(self.dm, rval) )

    def metricGetGradationFactor(self):
        """DMPlexMetricGetGradationFactor - Get the metric gradation factor

        Input parameters:
        . dm   - The DM

        Output parameters:
        . beta - The metric gradation factor

        .seealso: `DMPlexMetricSetGradationFactor`, `DMPlexMetricGetHausdorffNumber`

        See Also
        --------
        petsc.DMPlexMetricGetGradationFactor

        """
        cdef PetscReal beta = 0
        CHKERR( DMPlexMetricGetGradationFactor(self.dm, &beta) )
        return toReal(beta)

    def metricSetHausdorffNumber(self, hausd):
        """DMPlexMetricSetHausdorffNumber - Set the metric Hausdorff number

        Input parameters:
        + dm    - The DM
        - hausd - The metric Hausdorff number

        .seealso: `DMPlexMetricSetGradationFactor`, `DMPlexMetricGetHausdorffNumber`

        See Also
        --------
        petsc.DMPlexMetricSetHausdorffNumber

        """
        cdef PetscReal rval = asReal(hausd)
        CHKERR( DMPlexMetricSetHausdorffNumber(self.dm, rval) )

    def metricGetHausdorffNumber(self):
        """DMPlexMetricGetHausdorffNumber - Get the metric Hausdorff number

        Input parameters:
        . dm    - The DM

        Output parameters:
        . hausd - The metric Hausdorff number

        .seealso: `DMPlexMetricGetGradationFactor`, `DMPlexMetricSetHausdorffNumber`

        See Also
        --------
        petsc.DMPlexMetricGetHausdorffNumber

        """
        cdef PetscReal hausd = 0
        CHKERR( DMPlexMetricGetHausdorffNumber(self.dm, &hausd) )
        return toReal(hausd)

    def metricCreate(self, field=0):
        """DMPlexMetricCreate - Create a Riemannian metric field

        Input parameters:
        + dm     - The DM
        - f      - The field number to use

        Output parameter:
        . metric - The metric

        Notes:

        It is assumed that the DM is comprised of simplices.

        Command line options for Riemannian metrics:

        + -dm_plex_metric_isotropic                 - Is the metric isotropic?
        . -dm_plex_metric_uniform                   - Is the metric uniform?
        . -dm_plex_metric_restrict_anisotropy_first - Should anisotropy be restricted before normalization?
        . -dm_plex_metric_h_min                     - Minimum tolerated metric magnitude
        . -dm_plex_metric_h_max                     - Maximum tolerated metric magnitude
        . -dm_plex_metric_a_max                     - Maximum tolerated anisotropy
        . -dm_plex_metric_p                         - L-p normalization order
        - -dm_plex_metric_target_complexity         - Target metric complexity

        Switching between remeshers can be achieved using

        . -dm_adaptor <pragmatic/mmg/parmmg>        - specify dm adaptor to use

        Further options that are only relevant to Mmg and ParMmg:

        + -dm_plex_metric_gradation_factor          - Maximum ratio by which edge lengths may grow during gradation
        . -dm_plex_metric_num_iterations            - Number of parallel mesh adaptation iterations for ParMmg
        . -dm_plex_metric_no_insert                 - Should node insertion/deletion be turned off?
        . -dm_plex_metric_no_swap                   - Should facet swapping be turned off?
        . -dm_plex_metric_no_move                   - Should node movement be turned off?
        - -dm_plex_metric_verbosity                 - Choose a verbosity level from -1 (silent) to 10 (maximum).

        .seealso: `DMPlexMetricCreateUniform`, `DMPlexMetricCreateIsotropic`

        See Also
        --------
        petsc.DMPlexMetricCreate

        """
        cdef PetscInt ival = asInt(field)
        cdef Vec metric = Vec()
        CHKERR( DMPlexMetricCreate(self.dm, ival, &metric.vec) )
        return metric

    def metricCreateUniform(self, alpha, field=0):
        """DMPlexMetricCreateUniform - Construct a uniform isotropic metric

        Input parameters:
        + dm     - The DM
        . f      - The field number to use
        - alpha  - Scaling parameter for the diagonal

        Output parameter:
        . metric - The uniform metric

        .seealso: `DMPlexMetricCreate`, `DMPlexMetricCreateIsotropic`

        See Also
        --------
        petsc.DMPlexMetricCreateUniform

        """
        cdef PetscInt  ival = asInt(field)
        cdef PetscReal rval = asReal(alpha)
        cdef Vec metric = Vec()
        CHKERR( DMPlexMetricCreateUniform(self.dm, ival, rval, &metric.vec) )
        return metric

    def metricCreateIsotropic(self, Vec indicator, field=0):
        """DMPlexMetricCreateIsotropic - Construct an isotropic metric from an error indicator

        Input parameters:
        + dm        - The DM
        . f         - The field number to use
        - indicator - The error indicator

        Output parameter:
        . metric    - The isotropic metric

        .seealso: `DMPlexMetricCreate`, `DMPlexMetricCreateUniform`

        See Also
        --------
        petsc.DMPlexMetricCreateIsotropic

        """
        cdef PetscInt  ival = asInt(field)
        cdef Vec metric = Vec()
        CHKERR( DMPlexMetricCreateIsotropic(self.dm, ival, indicator.vec, &metric.vec) )
        return metric

    def metricDeterminantCreate(self, field=0):
        """DMPlexMetricDeterminantCreate - Create the determinant field for a Riemannian metric

        Input parameters:
        + dm          - The DM of the metric field
        - f           - The field number to use

        Output parameter:
        + determinant - The determinant field
        - dmDet       - The corresponding DM

        .seealso: DMPlexMetricCreateUniform(), DMPlexMetricCreateIsotropic(), DMPlexMetricCreate()

        See Also
        --------
        petsc.DMPlexMetricDeterminantCreate

        """
        cdef PetscInt  ival = asInt(field)
        cdef Vec determinant = Vec()
        cdef DM dmDet = DM()
        CHKERR( DMPlexMetricDeterminantCreate(self.dm, ival, &determinant.vec, &dmDet.dm) )
        return (determinant, dmDet)

    def metricEnforceSPD(self, Vec metric, Vec ometric, Vec determinant, restrictSizes=False, restrictAnisotropy=False):
        """DMPlexMetricEnforceSPD - Enforce symmetric positive-definiteness of a metric

        Input parameters:
        + dm                 - The DM
        . metricIn           - The metric
        . restrictSizes      - Should maximum/minimum metric magnitudes be enforced?
        - restrictAnisotropy - Should maximum anisotropy be enforced?

        Output parameter:
        + metricOut          - The metric
        - determinant        - Its determinant

        Notes:

        Relevant command line options:

        + -dm_plex_metric_isotropic - Is the metric isotropic?
        . -dm_plex_metric_uniform   - Is the metric uniform?
        . -dm_plex_metric_h_min     - Minimum tolerated metric magnitude
        . -dm_plex_metric_h_max     - Maximum tolerated metric magnitude
        - -dm_plex_metric_a_max     - Maximum tolerated anisotropy

        .seealso: `DMPlexMetricNormalize`, `DMPlexMetricIntersection`

        See Also
        --------
        petsc.DMPlexMetricEnforceSPD

        """
        cdef PetscBool bval_rs = asBool(restrictSizes)
        cdef PetscBool bval_ra = asBool(restrictAnisotropy)
        cdef DM dmDet = DM()
        CHKERR( DMPlexMetricEnforceSPD(self.dm, metric.vec, bval_rs, bval_ra, ometric.vec, determinant.vec) )
        return (ometric, determinant)

    def metricNormalize(self, Vec metric, Vec ometric, Vec determinant, restrictSizes=True, restrictAnisotropy=True):
        """DMPlexMetricNormalize - Apply L-p normalization to a metric

        Input parameters:
        + dm                 - The DM
        . metricIn           - The unnormalized metric
        . restrictSizes      - Should maximum/minimum metric magnitudes be enforced?
        - restrictAnisotropy - Should maximum metric anisotropy be enforced?

        Output parameter:
        . metricOut          - The normalized metric

        Notes:

        Relevant command line options:

        + -dm_plex_metric_isotropic                 - Is the metric isotropic?
        . -dm_plex_metric_uniform                   - Is the metric uniform?
        . -dm_plex_metric_restrict_anisotropy_first - Should anisotropy be restricted before normalization?
        . -dm_plex_metric_h_min                     - Minimum tolerated metric magnitude
        . -dm_plex_metric_h_max                     - Maximum tolerated metric magnitude
        . -dm_plex_metric_a_max                     - Maximum tolerated anisotropy
        . -dm_plex_metric_p                         - L-p normalization order
        - -dm_plex_metric_target_complexity         - Target metric complexity

        .seealso: `DMPlexMetricEnforceSPD`, `DMPlexMetricIntersection`

        See Also
        --------
        petsc.DMPlexMetricNormalize

        """
        cdef PetscBool bval_rs = asBool(restrictSizes)
        cdef PetscBool bval_ra = asBool(restrictAnisotropy)
        CHKERR( DMPlexMetricNormalize(self.dm, metric.vec, bval_rs, bval_ra, ometric.vec, determinant.vec) )
        return (ometric, determinant)

    def metricAverage2(self, Vec metric1, Vec metric2, Vec metricAvg):
        """DMPlexMetricAverage2 - Compute the unweighted average of two metrics

        Parameters
        ----------
        dm
            The DM
        metric1
            The first metric to be averaged
        metric2
            The second metric to be averaged

        Returns
        -------
        metricAvg
            The averaged metric

        .seealso: `DMPlexMetricAverage`, `DMPlexMetricAverage3`

        See Also
        --------
        petsc.DMPlexMetricAverage2

        """
        CHKERR( DMPlexMetricAverage2(self.dm, metric1.vec, metric2.vec, metricAvg.vec) )
        return metricAvg

    def metricAverage3(self, Vec metric1, Vec metric2, Vec metric3, Vec metricAvg):
        """DMPlexMetricAverage3 - Compute the unweighted average of three metrics

        Parameters
        ----------
        dm
            The DM
        metric1
            The first metric to be averaged
        metric2
            The second metric to be averaged
        metric3
            The third metric to be averaged

        Returns
        -------
        metricAvg
            The averaged metric

        .seealso: `DMPlexMetricAverage`, `DMPlexMetricAverage2`

        See Also
        --------
        petsc.DMPlexMetricAverage3

        """
        CHKERR( DMPlexMetricAverage3(self.dm, metric1.vec, metric2.vec, metric3.vec, metricAvg.vec) )
        return metricAvg

    def metricIntersection2(self, Vec metric1, Vec metric2, Vec metricInt):
        """DMPlexMetricIntersection2 - Compute the intersection of two metrics

        Parameters
        ----------
        dm
            The DM
        metric1
            The first metric to be intersected
        metric2
            The second metric to be intersected

        Returns
        -------
        metricInt
            The intersected metric

        .seealso: `DMPlexMetricIntersection`, `DMPlexMetricIntersection3`

        See Also
        --------
        petsc.DMPlexMetricIntersection2

        """
        CHKERR( DMPlexMetricIntersection2(self.dm, metric1.vec, metric2.vec, metricInt.vec) )
        return metricInt

    def metricIntersection3(self, Vec metric1, Vec metric2, Vec metric3, Vec metricInt):
        """DMPlexMetricIntersection3 - Compute the intersection of three metrics

        Parameters
        ----------
        dm
            The DM
        metric1
            The first metric to be intersected
        metric2
            The second metric to be intersected
        metric3
            The third metric to be intersected

        Returns
        -------
        metricInt
            The intersected metric

        .seealso: `DMPlexMetricIntersection`, `DMPlexMetricIntersection2`

        See Also
        --------
        petsc.DMPlexMetricIntersection3

        """
        CHKERR( DMPlexMetricIntersection3(self.dm, metric1.vec, metric2.vec, metric3.vec, metricInt.vec) )
        return metricInt

    def computeGradientClementInterpolant(self, Vec locX, Vec locC):
        """DMPlexComputeGradientClementInterpolant - This function computes the L2 projection of the cellwise gradient of a function u onto P1

        Collective.

        Parameters
        ----------
        dm
            The `DM`
        locX
            The coefficient vector u_h

        Returns
        -------
        locC
            A `Vec` which holds the Clement interpolant of the gradient

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMProjectFunction`, `DMComputeL2Diff`, `DMPlexComputeL2FieldDiff`, `DMComputeL2GradientDiff`

        See Also
        --------
        petsc.DMPlexComputeGradientClementInterpolant

        """
        CHKERR( DMPlexComputeGradientClementInterpolant(self.dm, locX.vec, locC.vec) )
        return locC

    # View

    def topologyView(self, Viewer viewer):
        """DMPlexTopologyView - Saves a `DMPlex` topology into a file

        Collective.

        Parameters
        ----------
        dm
            The `DM` whose topology is to be saved
        viewer
            The `PetscViewer` to save it in

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMView`, `DMPlexCoordinatesView`, `DMPlexLabelsView`, `DMPlexTopologyLoad`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexTopologyView

        """
        CHKERR( DMPlexTopologyView(self.dm, viewer.vwr))

    def coordinatesView(self, Viewer viewer):
        """DMPlexCoordinatesView - Saves `DMPlex` coordinates into a file

        Collective.

        Parameters
        ----------
        dm
            The `DM` whose coordinates are to be saved
        viewer
            The `PetscViewer` for saving

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMView`, `DMPlexTopologyView`, `DMPlexLabelsView`, `DMPlexCoordinatesLoad`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexCoordinatesView

        """
        CHKERR( DMPlexCoordinatesView(self.dm, viewer.vwr))

    def labelsView(self, Viewer viewer):
        """DMPlexLabelsView - Saves `DMPlex` labels into a file

        Collective.

        Parameters
        ----------
        dm
            The `DM` whose labels are to be saved
        viewer
            The `PetscViewer` for saving

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMView`, `DMPlexTopologyView`, `DMPlexCoordinatesView`, `DMPlexLabelsLoad`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexLabelsView

        """
        CHKERR( DMPlexLabelsView(self.dm, viewer.vwr))

    def sectionView(self, Viewer viewer, DM sectiondm):
        """DMPlexSectionView - Saves a section associated with a `DMPlex`

        Collective.

        Parameters
        ----------
        dm
            The `DM` that contains the topology on which the section to be saved is defined
        viewer
            The `PetscViewer` for saving
        sectiondm
            The `DM` that contains the section to be saved

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMView`, `DMPlexTopologyView`, `DMPlexCoordinatesView`, `DMPlexLabelsView`, `DMPlexGlobalVectorView`, `DMPlexLocalVectorView`, `PetscSectionView`, `DMPlexSectionLoad`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexSectionView

        """
        CHKERR( DMPlexSectionView(self.dm, viewer.vwr, sectiondm.dm))

    def globalVectorView(self, Viewer viewer, DM sectiondm, Vec vec):
        """DMPlexGlobalVectorView - Saves a global vector

        Collective.

        Parameters
        ----------
        dm
            The `DM` that represents the topology
        viewer
            The `PetscViewer` to save data with
        sectiondm
            The `DM` that contains the global section on which vec is defined
        vec
            The global vector to be saved

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexTopologyView`, `DMPlexSectionView`, `DMPlexLocalVectorView`, `DMPlexGlobalVectorLoad`, `DMPlexLocalVectorLoad`

        See Also
        --------
        petsc.DMPlexGlobalVectorView

        """
        CHKERR( DMPlexGlobalVectorView(self.dm, viewer.vwr, sectiondm.dm, vec.vec))

    def localVectorView(self, Viewer viewer, DM sectiondm, Vec vec):
        """DMPlexLocalVectorView - Saves a local vector

        Collective.

        Parameters
        ----------
        dm
            The `DM` that represents the topology
        viewer
            The `PetscViewer` to save data with
        sectiondm
            The `DM` that contains the local section on which `vec` is defined; may be the same as `dm`
        vec
            The local vector to be saved

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexTopologyView`, `DMPlexSectionView`, `DMPlexGlobalVectorView`, `DMPlexGlobalVectorLoad`, `DMPlexLocalVectorLoad`

        See Also
        --------
        petsc.DMPlexLocalVectorView

        """
        CHKERR( DMPlexLocalVectorView(self.dm, viewer.vwr, sectiondm.dm, vec.vec))

    # Load

    def topologyLoad(self, Viewer viewer):
        """DMPlexTopologyLoad - Loads a topology into a `DMPlex`

        Collective.

        Parameters
        ----------
        dm
            The `DM` into which the topology is loaded
        viewer
            The `PetscViewer` for the saved topology

        Returns
        -------
        globalToLocalPointSF
            The `PetscSF` that pushes points in [0, N) to the associated points in the loaded `DMPlex`, where N is the global number of points; ``None`` if unneeded

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMLoad`, `DMPlexCoordinatesLoad`, `DMPlexLabelsLoad`, `DMView`, `PetscViewerHDF5Open`, `PetscViewerPushFormat`,
        `PetscViewer`, `PetscSF`

        See Also
        --------
        petsc.DMPlexTopologyLoad

        """
        cdef SF sf = SF()
        CHKERR( DMPlexTopologyLoad(self.dm, viewer.vwr, &sf.sf))
        return sf

    def coordinatesLoad(self, Viewer viewer, SF sfxc):
        """DMPlexCoordinatesLoad - Loads coordinates into a `DMPlex`

        Collective.

        Parameters
        ----------
        dm
            The `DM` into which the coordinates are loaded
        viewer
            The `PetscViewer` for the saved coordinates
        globalToLocalPointSF
            The `PetscSF` returned by `DMPlexTopologyLoad` when loading dm from viewer

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMLoad`, `DMPlexTopologyLoad`, `DMPlexLabelsLoad`, `DMView`, `PetscViewerHDF5Open`, `PetscViewerPushFormat`,
        `PetscSF`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexCoordinatesLoad

        """
        CHKERR( DMPlexCoordinatesLoad(self.dm, viewer.vwr, sfxc.sf))

    def labelsLoad(self, Viewer viewer, SF sfxc):
        """DMPlexLabelsLoad - Loads labels into a `DMPlex`

        Collective.

        Parameters
        ----------
        dm
            The `DM` into which the labels are loaded
        viewer
            The `PetscViewer` for the saved labels
        globalToLocalPointSF
            The `PetscSF` returned by `DMPlexTopologyLoad` when loading `dm` from viewer

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMLoad`, `DMPlexTopologyLoad`, `DMPlexCoordinatesLoad`, `DMView`, `PetscViewerHDF5Open`, `PetscViewerPushFormat`,
        `PetscSF`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexLabelsLoad

        """
        CHKERR( DMPlexLabelsLoad(self.dm, viewer.vwr, sfxc.sf))

    def sectionLoad(self, Viewer viewer, DM sectiondm, SF sfxc):
        """DMPlexSectionLoad - Loads section into a `DMPlex`

        Collective.

        Parameters
        ----------
        dm
            The `DM` that represents the topology
        viewer
            The `PetscViewer` that represents the on-disk section (sectionA)
        sectiondm
            The `DM` into which the on-disk section (sectionA) is migrated
        globalToLocalPointSF
            The `PetscSF` returned by `DMPlexTopologyLoad(`) when loading dm from viewer

        Returns
        -------
        globalDofSF
            The `PetscSF` that migrates any on-disk `Vec` data associated with sectionA into a global `Vec` associated with the `sectiondm`'s global section (``None`` if not needed)
        localDofSF
            The `PetscSF` that migrates any on-disk `Vec` data associated with sectionA into a local `Vec` associated with the `sectiondm`'s local section (``None`` if not needed)

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMLoad`, `DMPlexTopologyLoad`, `DMPlexCoordinatesLoad`, `DMPlexLabelsLoad`, `DMPlexGlobalVectorLoad`, `DMPlexLocalVectorLoad`, `PetscSectionLoad`, `DMPlexSectionView`, `PetscSF`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexSectionLoad

        """
        cdef SF gsf = SF()
        cdef SF lsf = SF()
        CHKERR( DMPlexSectionLoad(self.dm, viewer.vwr, sectiondm.dm, sfxc.sf, &gsf.sf, &lsf.sf))
        return gsf, lsf

    def globalVectorLoad(self, Viewer viewer, DM sectiondm, SF sf, Vec vec):
        """DMPlexGlobalVectorLoad - Loads on-disk vector data into a global vector

        Collective.

        Parameters
        ----------
        dm
            The `DM` that represents the topology
        viewer
            The `PetscViewer` that represents the on-disk vector data
        sectiondm
            The `DM` that contains the global section on which vec is defined
        sf
            The `PetscSF` that migrates the on-disk vector data into vec
        vec
            The global vector to set values of

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexTopologyLoad`, `DMPlexSectionLoad`, `DMPlexLocalVectorLoad`, `DMPlexGlobalVectorView`, `DMPlexLocalVectorView`,
        `PetscSF`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexGlobalVectorLoad

        """
        CHKERR( DMPlexGlobalVectorLoad(self.dm, viewer.vwr, sectiondm.dm, sf.sf, vec.vec))

    def localVectorLoad(self, Viewer viewer, DM sectiondm, SF sf, Vec vec):
        """
DMPlexLocalVectorLoad - Loads on-disk vector data into a local vector

        Collective.

        Parameters
        ----------
        dm
            The `DM` that represents the topology
        viewer
            The `PetscViewer` that represents the on-disk vector data
        sectiondm
            The `DM` that contains the local section on which vec is defined
        sf
            The `PetscSF` that migrates the on-disk vector data into vec
        vec
            The local vector to set values of

        .seealso: [](chapter_unstructured), `DM`, `DMPlex`, `DMPlexTopologyLoad`, `DMPlexSectionLoad`, `DMPlexGlobalVectorLoad`, `DMPlexGlobalVectorView`, `DMPlexLocalVectorView`,
        `PetscSF`, `PetscViewer`

        See Also
        --------
        petsc.DMPlexLocalVectorLoad

        """
        CHKERR( DMPlexLocalVectorLoad(self.dm, viewer.vwr, sectiondm.dm, sf.sf, vec.vec))

# --------------------------------------------------------------------

del DMPlexReorderDefaultFlag

# --------------------------------------------------------------------
