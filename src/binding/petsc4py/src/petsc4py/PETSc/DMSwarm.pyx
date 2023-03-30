# --------------------------------------------------------------------

class DMSwarmType(object):
    BASIC = DMSWARM_BASIC
    PIC = DMSWARM_PIC

class DMSwarmMigrateType(object):
    MIGRATE_BASIC = DMSWARM_MIGRATE_BASIC
    MIGRATE_DMCELLNSCATTER = DMSWARM_MIGRATE_DMCELLNSCATTER
    MIGRATE_DMCELLEXACT = DMSWARM_MIGRATE_DMCELLEXACT
    MIGRATE_USER = DMSWARM_MIGRATE_USER

class DMSwarmCollectType(object):
    COLLECT_BASIC = DMSWARM_COLLECT_BASIC
    COLLECT_DMDABOUNDINGBOX = DMSWARM_COLLECT_DMDABOUNDINGBOX
    COLLECT_GENERAL = DMSWARM_COLLECT_GENERAL
    COLLECT_USER = DMSWARM_COLLECT_USER

class DMSwarmPICLayoutType(object):
    LAYOUT_REGULAR = DMSWARMPIC_LAYOUT_REGULAR
    LAYOUT_GAUSS = DMSWARMPIC_LAYOUT_GAUSS
    LAYOUT_SUBDIVISION = DMSWARMPIC_LAYOUT_SUBDIVISION


cdef class DMSwarm(DM):
    """
    A DM object used to represent arrays of data (fields) of arbitrary type.
    """
    Type = DMSwarmType
    MigrateType = DMSwarmMigrateType
    CollectType = DMSwarmCollectType
    PICLayoutType = DMSwarmPICLayoutType

    def create(self, comm: Comm | None = None) -> Self:
        """Create an empty DM object and set its type to `DMSWARM`. TODO:

        DMs are the abstract objects in PETSc that mediate between meshes and
        discretizations and the algebraic solvers, time integrators, and
        optimization algorithms.

        Collective.

        Parameters
        ----------
        comm
            The communicator for the DM object.

        See also
        --------
        petsc.DMCreate, petsc.DMSetType

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM newdm = NULL
        CHKERR( DMCreate(ccomm, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        CHKERR( DMSetType(self.dm, DMSWARM) )
        return self

    def createGlobalVectorFromField(self, fieldname: str) -> Vec:
        """Create a `Vec` object sharing the array associated with a given field.

        The vector must be returned using a matching call to
        `destroyGlobalVectorFromField`.

        Collective.

        Parameters
        ----------
        fieldname
            The textual name given to a registered field.

        See also
        --------
        petsc.DMSwarmCreateGlobalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef Vec vg = Vec()
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmCreateGlobalVectorFromField(self.dm, cfieldname, &vg.vec) )
        return vg

    def destroyGlobalVectorFromField(self, fieldname: str) -> None:
        """Destroy the `Vec` object which share the array associated with a given field.

        Collective.

        Parameters
        ----------
        fieldname
            The textual name given to a registered field.

        See also
        --------
        petsc.DMSwarmDestroyGlobalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef PetscVec vec = NULL
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmDestroyGlobalVectorFromField(self.dm, cfieldname, &vec) )

    def createLocalVectorFromField(self, fieldname: str) -> Vec:
        """Create a `Vec` object sharing the array associated with a given field.

        The vector must be returned using a matching call to `destroyLocalVectorFromField`.

        Collective.

        Parameters
        ----------
        fieldname
            The textual name given to a registered field.

        See also
        --------
        petsc.DMSwarmCreateLocalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef Vec vl = Vec()
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmCreateLocalVectorFromField(self.dm, cfieldname, &vl.vec) )
        return vl

    def destroyLocalVectorFromField(self, fieldname: str) -> None:
        """Destroy the `Vec` object which share the array associated with a given field.

        Collective.

        Parameters
        ----------
        fieldname
            The textual name given to a registered field.

        See also
        --------
        petsc.DMSwarmDestroyLocalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef PetscVec vec
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmDestroyLocalVectorFromField(self.dm, cfieldname, &vec) )

    def initializeFieldRegister(self):
        """Initiate the registration of fields to a `DMSwarm`.

        After all fields have been registered, you must call `finalizeFieldRegister`.

        Collective.

        See also
        --------
        petsc.DMSwarmInitializeFieldRegister

        """
        CHKERR( DMSwarmInitializeFieldRegister(self.dm) )

    def finalizeFieldRegister(self) -> None:
        """Finalize the registration of fields to a `DMSwarm`.

        Collective.

        See also
        --------
        petsc.DMSwarmFinalizeFieldRegister

        """
        CHKERR( DMSwarmFinalizeFieldRegister(self.dm) )

    def setLocalSizes(self, nlocal: int, buffer: int) -> Self:
        """Set the length of all registered fields on the `DMSwarm`.

        Not collective.

        Parameters
        ----------
nlocal - the length of each registered field
buffer - the length of the buffer used to efficient dynamic re-sizing

        See also
        --------
        petsc.DMSwarmSetLocalSizes

        """
        cdef PetscInt cnlocal = asInt(nlocal)
        cdef PetscInt cbuffer = asInt(buffer)
        CHKERR( DMSwarmSetLocalSizes(self.dm, cnlocal, cbuffer) )
        return self

    def registerField(self, fieldname, blocksize, dtype=ScalarType): #TODO:
        """Register a field to a `DMSwarm` with a native PETSc data type.

        Collective.

        Parameters
        ----------
        fieldname
            The textual name to identify this field.
        blocksize
            The number of each data type.
        type
            A valid PETSc data type (PETSC_CHAR, PETSC_SHORT, PETSC_INT, PETSC_FLOAT, PETSC_REAL, PETSC_LONG).

        See also
        --------
        petsc.DMSwarmRegisterPetscDatatypeField

        """
        cdef const char *cfieldname = NULL
        cdef PetscInt cblocksize = asInt(blocksize)
        cdef PetscDataType ctype  = PETSC_DATATYPE_UNKNOWN
        if dtype == IntType:     ctype  = PETSC_INT
        if dtype == RealType:    ctype = PETSC_REAL
        if dtype == ScalarType:  ctype = PETSC_SCALAR
        if dtype == ComplexType: ctype = PETSC_COMPLEX
        assert ctype != PETSC_DATATYPE_UNKNOWN
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmRegisterPetscDatatypeField(self.dm, cfieldname, cblocksize, ctype) )

    # TODO: return type
    def getField(self, fieldname):
        """get/return access to the underlying array storing all entries associated with a registered field.

        The array must be returned using a matching call to `restoreField`.

        Not collective.

        Parameters
        ----------
        fieldname
            The textual name to identify this field.

        Returns
        -------
        TODO:
        blocksize
            the number of each data type
        type
            the data type
        data
            pointer to raw array


        See also
        --------
        petsc.DMSwarmGetField

        """
        cdef const char *cfieldname = NULL
        cdef PetscInt blocksize = 0
        cdef PetscDataType ctype = PETSC_DATATYPE_UNKNOWN
        cdef PetscReal *data = NULL
        cdef PetscInt nlocal = 0
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmGetField(self.dm, cfieldname, &blocksize, &ctype, <void**> &data) )
        CHKERR( DMSwarmGetLocalSize(self.dm, &nlocal) )
        cdef int typenum = -1
        if ctype == PETSC_INT:     typenum = NPY_PETSC_INT
        if ctype == PETSC_REAL:    typenum = NPY_PETSC_REAL
        if ctype == PETSC_SCALAR:  typenum = NPY_PETSC_SCALAR
        if ctype == PETSC_COMPLEX: typenum = NPY_PETSC_COMPLEX
        assert typenum != -1
        cdef npy_intp s = <npy_intp> nlocal * blocksize
        return <object> PyArray_SimpleNewFromData(1, &s, typenum, data)

    def restoreField(self, fieldname: str):
        """Restore access to the underlying array storing all entries associated with a registered field.

        Not collective.

        Parameters
        ----------
        fieldname
            TODO the textual name to identify this field
        Output parameters
        blocksize TODO the number of each data type
        type TODO the data type
        data TODO pointer to raw array

        See also
        --------
        petsc.DMSwarmRestoreField

        """
        cdef const char *cfieldname = NULL
        cdef PetscInt blocksize = 0
        cdef PetscDataType ctype = PETSC_DATATYPE_UNKNOWN
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmRestoreField(self.dm, cfieldname, &blocksize, &ctype, <void**> 0) )

    def vectorDefineField(self, fieldname):
        """Set the field from which to define a `Vec` object when `createLocalVector`, or `createGlobalVector` is called.

        Collective.

        Parameters
        ----------
        fieldname
            The textual name given to a registered field.

        See also
        --------
        petsc.DMSwarmVectorDefineField

        """
        cdef const char *cval = NULL
        fieldname = str2bytes(fieldname, &cval)
        CHKERR( DMSwarmVectorDefineField(self.dm, cval) )

    def addPoint(self):
        """Add space for one new point in the `DMSwarm`.

        Not collective.

        See also
        --------
        petsc.DMSwarmAddPoint

        """
        CHKERR( DMSwarmAddPoint(self.dm) )

    def addNPoints(self, npoints: int) -> None:
        """Add space for a number of new points in the `DMSwarm`.

        Not collective.

        Parameters
        ----------
        npoints
            The number of new points to add.


        See also
        --------
        petsc.DMSwarmAddNPoints

        """
        cdef PetscInt cnpoints = asInt(npoints)
        CHKERR( DMSwarmAddNPoints(self.dm, cnpoints) )

    def removePoint(self) -> None:
        """Remove the last point from the `DMSwarm`.

        Not collective.

        See also
        --------
        petsc.DMSwarmRemovePoint

        """
        CHKERR( DMSwarmRemovePoint(self.dm) )

    def removePointAtIndex(self, index: int) -> None:
        """Remove a specific point from the DMSWARM

        Not collective.

        Parameters
        ----------
        idx
            Index of point to remove

        See also
        --------
        petsc.DMSwarmRemovePointAtIndex

        """
        cdef PetscInt cindex = asInt(index)
        CHKERR( DMSwarmRemovePointAtIndex(self.dm, cindex) )

    def copyPoint(self, pi: int, pj: int) -> None:
        """Copy point pi to point pj in the `DMSwarm`.

        Not collective.

        Parameters
        ----------
        pi
            The index of the point to copy (source).
        pj
            The point index where the copy should be located (destination).

        See also
        --------
        petsc.DMSwarmCopyPoint

        """
        cdef PetscInt cpi = asInt(pi)
        cdef PetscInt cpj = asInt(pj)
        CHKERR( DMSwarmCopyPoint(self.dm, cpi, cpj) )

    def getLocalSize(self) -> int:
        """Return the local length of fields registered.

        Not collective.

        See also
        --------
        petsc.DMSwarmGetLocalSize

        """
        cdef PetscInt size = asInt(0)
        CHKERR( DMSwarmGetLocalSize(self.dm, &size) )
        return toInt(size)

    def getSize(self) -> int:
        """Return the total length of fields registered.

        Collective.

        See also
        --------
        petsc.DMSwarmGetSize

        """
        cdef PetscInt size = asInt(0)
        CHKERR( DMSwarmGetSize(self.dm, &size) )
        return toInt(size)

    def migrate(self, remove_sent_points=False) -> None: #TODO: remove_sent_points = bool | None = False ?
        """Relocate points defined in the `DMSwarm` to other MPI-ranks.

        Collective.

        Parameters
        ----------
        remove_sent_points
            Flag indicating if sent points should be removed from the current MPI-rank.

        See also
        --------
        petsc.DMSwarmMigrate

        """
        cdef PetscBool remove_pts = asBool(remove_sent_points)
        CHKERR( DMSwarmMigrate(self.dm, remove_pts) )

    def collectViewCreate(self) -> None:
        """Apply a collection method and gathers points in neighbour ranks into the `DMSwarm`.

        Collective.

        See also
        --------
        petsc.DMSwarmCollectViewCreate

        """
        CHKERR( DMSwarmCollectViewCreate(self.dm) )

    def collectViewDestroy(self) -> None:
        """Reset the `DMSwarm` to the size prior to calling `collectViewCreate`.

        Collective.

        See also
        --------
        petsc.DMSwarmCollectViewDestroy

        """
        CHKERR( DMSwarmCollectViewDestroy(self.dm) )

    def setCellDM(self, DM dm):
        """Attach a `DM` to a `DMSwarm`.

        Collective.

        Parameters
        ----------
        dm
            The `DM` to attach to the `DMSwarm`.

        See also
        --------
        petsc.DMSwarmSetCellDM

        """
        CHKERR( DMSwarmSetCellDM(self.dm, dm.dm) )

    def getCellDM(self) -> DM:
        """Return `DM` cell attached to `DMSwarm`.

        Collective.

        See also
        --------
        petsc.DMSwarmGetCellDM

        """
        cdef PetscDM newdm = NULL
        CHKERR( DMSwarmGetCellDM(self.dm, &newdm) )
        cdef DM dm = subtype_DM(newdm)()
        dm.dm = newdm
        PetscINCREF(dm.obj)
        return dm

    def setType(self, dmswarm_type) -> None:
        """Set particular flavor of `DMSwarm`.

        Collective.

        Parameters
        ----------
        stype
            TODO the DMSWARM type (e.g. DMSWARM_PIC)

        See also
        --------
        petsc.DMSwarmSetType

        """
        cdef PetscDMSwarmType cval = dmswarm_type
        CHKERR( DMSwarmSetType(self.dm, cval) )

    def setPointsUniformCoordinates(self, min, max, npoints, mode=None) -> Self:
        """Set point coordinates in a `DMSwarm` on a regular (ijk) grid.

        Collective.

        Parameters
        ----------
        min
            TODO minimum coordinate values in the x, y, z directions (array of length dim)
        max
            TODO maximum coordinate values in the x, y, z directions (array of length dim)
        npoints
            TODO number of points in each spatial direction (array of length dim)
        mode
            TODO indicates whether to append points to the swarm (ADD_VALUES), or over-ride existing points (INSERT_VALUES)

        See also
        --------
        petsc.DMSwarmSetPointsUniformCoordinates

        """
        cdef PetscInt dim = asInt(0)
        CHKERR( DMGetDimension(self.dm, &dim) )
        cdef PetscReal cmin[3]
        cmin[0] = cmin[1] = cmin[2] = asReal(0.)
        for i from 0 <= i < dim: cmin[i] = min[i]
        cdef PetscReal cmax[3]
        cmax[0] = cmax[1] = cmax[2] = asReal(0.)
        for i from 0 <= i < dim: cmax[i] = max[i]
        cdef PetscInt cnpoints[3]
        cnpoints[0] = cnpoints[1] = cnpoints[2] = asInt(0)
        for i from 0 <= i < dim: cnpoints[i] = npoints[i]
        cdef PetscInsertMode cmode = insertmode(mode)
        CHKERR( DMSwarmSetPointsUniformCoordinates(self.dm, cmin, cmax, cnpoints, cmode) )
        return self

    def setPointCoordinates(self, coordinates, redundant=False, mode=None) -> None:
        """Set point coordinates in a `DMSwarm` from a user defined list.

        Collective.

        Parameters
        ----------
        npoints
            TODO: the number of points to insert
        coor
            TODO: the coordinate values
        redundant
            TODO: if set to PETSC_TRUE, it is assumed that npoints and coor are only valid on rank 0 and should be broadcast to other ranks
        mode
            TODO: indicates whether to append points to the swarm (ADD_VALUES), or over-ride existing points (INSERT_VALUES)

        See also
        --------
        petsc.DMSwarmSetPointCoordinates

        """
        cdef ndarray xyz = iarray(coordinates, NPY_PETSC_REAL)
        if PyArray_ISFORTRAN(xyz): xyz = PyArray_Copy(xyz)
        if PyArray_NDIM(xyz) != 2: raise ValueError(
            ("coordinates must have two dimensions: "
             "coordinates.ndim=%d") % (PyArray_NDIM(xyz)) )
        cdef PetscInt cnpoints = <PetscInt> PyArray_DIM(xyz, 0)
        cdef PetscBool credundant = asBool(redundant)
        cdef PetscInsertMode cmode = insertmode(mode)
        cdef PetscReal *coords = <PetscReal*> PyArray_DATA(xyz)
        CHKERR( DMSwarmSetPointCoordinates(self.dm, cnpoints, coords, credundant, cmode) )

    # TODO: layoutType: is not none
    def insertPointUsingCellDM(self, layoutType: None, fill_param: int) -> None:
        """Insert point coordinates within each cell.

        Not collective.

        Parameters
        ----------
        layout_type
            Method used to fill each cell with the cell DM.
        fill_param
            Parameter controlling how many points per cell are added (the meaning of this parameter is dependent on the layout type).

        See also
        --------
        petsc.DMSwarmInsertPointsUsingCellDM

        """
        cdef PetscDMSwarmPICLayoutType clayoutType = layoutType
        cdef PetscInt cfill_param = asInt(fill_param)
        CHKERR( DMSwarmInsertPointsUsingCellDM(self.dm, clayoutType, cfill_param) )

    def setPointCoordinatesCellwise(self, coordinates) -> None:
        """Insert point coordinates (defined over the reference cell) within each cell.

        Not collective.

        Parameters
        ----------
        celldm
            TODO the cell DM
        npoints
            TODO the number of points to insert in each cell
        xi
            TODO the coordinates (defined in the local coordinate system for each cell) to insert

        See also
        --------
        petsc.DMSwarmSetPointCoordinatesCellwise

        """
        cdef ndarray xyz = iarray(coordinates, NPY_PETSC_REAL)
        if PyArray_ISFORTRAN(xyz): xyz = PyArray_Copy(xyz)
        if PyArray_NDIM(xyz) != 2: raise ValueError(
            ("coordinates must have two dimensions: "
             "coordinates.ndim=%d") % (PyArray_NDIM(xyz)) )
        cdef PetscInt cnpoints = <PetscInt> PyArray_DIM(xyz, 0)
        cdef PetscReal *coords = <PetscReal*> PyArray_DATA(xyz)
        CHKERR( DMSwarmSetPointCoordinatesCellwise(self.dm, cnpoints, coords) )

    def viewFieldsXDMF(self, filename: str, fieldnames: Sequence[str]) -> None:
        """Write a selection of `DMSwarm` fields to an XDMF3 file.

        Collective.

        Parameters
        ----------
        filename
            The file name of the XDMF file (must have the extension .xmf).
        fieldnames
            Array containing the textual name of fields to write.

        See also
        --------
        petsc.DMSwarmViewFieldsXDMF

        """
        cdef const char *cval = NULL
        cdef const char *cfilename = NULL
        filename = str2bytes(filename, &cfilename)
        cdef PetscInt cnfields = <PetscInt> len(fieldnames)
        cdef const char** cfieldnames = NULL
        cdef object tmp = oarray_p(empty_p(cnfields), NULL, <void**>&cfieldnames)
        fieldnames = list(fieldnames)
        for i from 0 <= i < cnfields:
            fieldnames[i] = str2bytes(fieldnames[i], &cval)
            cfieldnames[i] = cval
        CHKERR( DMSwarmViewFieldsXDMF(self.dm, cfilename, cnfields, cfieldnames ) )

    def viewXDMF(self, filename: str) -> None:
        """Write this `DMSwarm` fields to an XDMF3 file.

        Collective.

        Parameters
        ----------
        filename
            The file name of the XDMF file (must have the extension .xmf).

        See also
        --------
        petsc.DMSwarmViewXDMF

        """
        cdef const char *cval = NULL
        filename = str2bytes(filename, &cval)
        CHKERR( DMSwarmViewXDMF(self.dm, cval) )

    def sortGetAccess(self) -> None:
        """Setup up a `DMSwarm` point sort context for efficient traversal of points within a cell.

        You must call `sortRestoreAccess` when you no longer need access to the sort context.

        Not collective.

        See also
        --------
        petsc.DMSwarmSortGetAccess

        """
        CHKERR( DMSwarmSortGetAccess(self.dm) )

    def sortRestoreAccess(self) -> None:
        """Invalidate the `DMSwarm` point sorting context.

        Not collective.

        See also
        --------
        petsc.DMSwarmSortRestoreAccess

        """
        CHKERR( DMSwarmSortRestoreAccess(self.dm) )

    def sortGetPointsPerCell(self, e: int) -> list[int]:
        """Create an array of point indices for all points in a cell.

        Not collective.

        Parameters
        ----------
        e
            The index of the cell.

        See also
        --------
        petsc.DMSwarmSortGetPointsPerCell

        """
        cdef PetscInt ce = asInt(e)
        cdef PetscInt cnpoints = asInt(0)
        cdef PetscInt *cpidlist = NULL
        cdef list pidlist = []
        CHKERR( DMSwarmSortGetPointsPerCell(self.dm, ce, &cnpoints, &cpidlist) )
        npoints = asInt(cnpoints)
        for i from 0 <= i < npoints: pidlist.append(asInt(cpidlist[i]))
        return pidlist

    def sortGetNumberOfPointsPerCell(self, e: int) -> int:
        """Return the number of points in a cell.

        Not collective.

        Parameters
        ----------
        e
            The index of the cell.

        See also
        --------
        petsc.DMSwarmSortGetNumberOfPointsPerCell

        """
        cdef PetscInt ce = asInt(e)
        cdef PetscInt npoints = asInt(0)
        CHKERR( DMSwarmSortGetNumberOfPointsPerCell(self.dm, ce, &npoints) )
        return toInt(npoints)

    def sortGetIsValid(self) -> bool:
        """Return flag indicating whether the sort context is up-to-date.

        Returns the `isvalid` flag associated with a `DMSwarm` point sorting context.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSortGetIsValid

        """
        cdef PetscBool isValid = asBool(False)
        CHKERR( DMSwarmSortGetIsValid(self.dm, &isValid) )
        return toBool(isValid)

    def sortGetSizes(self) -> tuple(int, int):
        """Return the sizes associated with a `DMSwarm` point sorting context.

        Not collective.

        Returns
        ----------
        ncells : int
            Number of cells within the sort context.
        npoints : int
            Number of points used to create the sort context.

        See also
        --------
        petsc.DMSwarmSortGetSizes

        """
        cdef PetscInt ncells = asInt(0)
        cdef PetscInt npoints = asInt(0)
        CHKERR( DMSwarmSortGetSizes(self.dm, &ncells, &npoints) )
        return (toInt(ncells), toInt(npoints))

    def projectFields(self, fieldnames, reuse=False): #TODO:
        """Project a set of swarm fields onto the cell `DM`.

        Collective.

        Parameters
        ----------
        nfields - the number of swarm fields to project
        fieldnames - the textual names of the swarm fields to project
        fields - an array of ``Vec``s of length nfields
        reuse - flag indicating whether the array and contents of fields should be re-used or internally allocated

        See also
        --------
        petsc.DMSwarmProjectFields

        """
        cdef PetscBool creuse = asBool(reuse)
        cdef const char *cval = NULL
        cdef PetscInt cnfields = <PetscInt> len(fieldnames)
        cdef const char** cfieldnames = NULL
        cdef object tmp = oarray_p(empty_p(cnfields), NULL, <void**>&cfieldnames)
        cdef PetscVec *cfieldvecs
        fieldnames = list(fieldnames)
        for i from 0 <= i < cnfields:
            fieldnames[i] = str2bytes(fieldnames[i], &cval)
            cfieldnames[i] = cval
        CHKERR( DMSwarmProjectFields(self.dm, cnfields, cfieldnames, &cfieldvecs, creuse) )
        cdef list fieldvecs = []
        for i from 0 <= i < cnfields:
            newVec = Vec()
            newVec.vec = cfieldvecs[i]
            fieldvecs.append(newVec)
        return fieldvecs


del DMSwarmType
del DMSwarmMigrateType
del DMSwarmCollectType
del DMSwarmPICLayoutType
