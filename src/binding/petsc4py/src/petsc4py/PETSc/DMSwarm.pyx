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

    Type = DMSwarmType
    MigrateType = DMSwarmMigrateType
    CollectType = DMSwarmCollectType
    PICLayoutType = DMSwarmPICLayoutType

    def create(self, comm=None):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def createGlobalVectorFromField(self, fieldname):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmCreateGlobalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef Vec vg = Vec()
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmCreateGlobalVectorFromField(self.dm, cfieldname, &vg.vec) )
        return vg

    def destroyGlobalVectorFromField(self, fieldname):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmDestroyGlobalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef PetscVec vec = NULL
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmDestroyGlobalVectorFromField(self.dm, cfieldname, &vec) )

    def createLocalVectorFromField(self, fieldname):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmCreateLocalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef Vec vl = Vec()
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmCreateLocalVectorFromField(self.dm, cfieldname, &vl.vec) )
        return vl

    def destroyLocalVectorFromField(self, fieldname):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmDestroyLocalVectorFromField

        """
        cdef const char *cfieldname = NULL
        cdef PetscVec vec
        fieldname = str2bytes(fieldname, &cfieldname)
        CHKERR( DMSwarmDestroyLocalVectorFromField(self.dm, cfieldname, &vec) )

    def initializeFieldRegister(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmInitializeFieldRegister

        """
        CHKERR( DMSwarmInitializeFieldRegister(self.dm) )

    def finalizeFieldRegister(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmFinalizeFieldRegister

        """
        CHKERR( DMSwarmFinalizeFieldRegister(self.dm) )

    def setLocalSizes(self, nlocal, buffer):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSetLocalSizes

        """
        cdef PetscInt cnlocal = asInt(nlocal)
        cdef PetscInt cbuffer = asInt(buffer)
        CHKERR( DMSwarmSetLocalSizes(self.dm, cnlocal, cbuffer) )
        return self

    def registerField(self, fieldname, blocksize, dtype=ScalarType):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def getField(self, fieldname):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def restoreField(self, fieldname):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmVectorDefineField

        """
        cdef const char *cval = NULL
        fieldname = str2bytes(fieldname, &cval)
        CHKERR( DMSwarmVectorDefineField(self.dm, cval) )

    def addPoint(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmAddPoint

        """
        CHKERR( DMSwarmAddPoint(self.dm) )

    def addNPoints(self, npoints):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmAddNPoints

        """
        cdef PetscInt cnpoints = asInt(npoints)
        CHKERR( DMSwarmAddNPoints(self.dm, cnpoints) )

    def removePoint(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmRemovePoint

        """
        CHKERR( DMSwarmRemovePoint(self.dm) )

    def removePointAtIndex(self, index):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmRemovePointAtIndex

        """
        cdef PetscInt cindex = asInt(index)
        CHKERR( DMSwarmRemovePointAtIndex(self.dm, cindex) )

    def copyPoint(self, pi, pj):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmCopyPoint

        """
        cdef PetscInt cpi = asInt(pi)
        cdef PetscInt cpj = asInt(pj)
        CHKERR( DMSwarmCopyPoint(self.dm, cpi, cpj) )

    def getLocalSize(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmGetLocalSize

        """
        cdef PetscInt size = asInt(0)
        CHKERR( DMSwarmGetLocalSize(self.dm, &size) )
        return toInt(size)

    def getSize(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmGetSize

        """
        cdef PetscInt size = asInt(0)
        CHKERR( DMSwarmGetSize(self.dm, &size) )
        return toInt(size)

    def migrate(self, remove_sent_points=False):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmMigrate

        """
        cdef PetscBool remove_pts = asBool(remove_sent_points)
        CHKERR( DMSwarmMigrate(self.dm, remove_pts) )

    def collectViewCreate(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmCollectViewCreate

        """
        CHKERR( DMSwarmCollectViewCreate(self.dm) )

    def collectViewDestroy(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmCollectViewDestroy

        """
        CHKERR( DMSwarmCollectViewDestroy(self.dm) )

    def setCellDM(self, DM dm):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSetCellDM

        """
        CHKERR( DMSwarmSetCellDM(self.dm, dm.dm) )

    def getCellDM(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def setType(self, dmswarm_type):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSetType

        """
        cdef PetscDMSwarmType cval = dmswarm_type
        CHKERR( DMSwarmSetType(self.dm, cval) )

    def setPointsUniformCoordinates(self, min, max, npoints, mode=None):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def setPointCoordinates(self, coordinates, redundant=False, mode=None):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def insertPointUsingCellDM(self, layoutType, fill_param):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmInsertPointsUsingCellDM

        """
        cdef PetscDMSwarmPICLayoutType clayoutType = layoutType
        cdef PetscInt cfill_param = asInt(fill_param)
        CHKERR( DMSwarmInsertPointsUsingCellDM(self.dm, clayoutType, cfill_param) )

    def setPointCoordinatesCellwise(self, coordinates):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def viewFieldsXDMF(self, filename, fieldnames):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def viewXDMF(self, filename):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmViewXDMF

        """
        cdef const char *cval = NULL
        filename = str2bytes(filename, &cval)
        CHKERR( DMSwarmViewXDMF(self.dm, cval) )

    def sortGetAccess(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSortGetAccess

        """
        CHKERR( DMSwarmSortGetAccess(self.dm) )

    def sortRestoreAccess(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSortRestoreAccess

        """
        CHKERR( DMSwarmSortRestoreAccess(self.dm) )

    def sortGetPointsPerCell(self, e):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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

    def sortGetNumberOfPointsPerCell(self, e):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSortGetNumberOfPointsPerCell

        """
        cdef PetscInt ce = asInt(e)
        cdef PetscInt npoints = asInt(0)
        CHKERR( DMSwarmSortGetNumberOfPointsPerCell(self.dm, ce, &npoints) )
        return toInt(npoints)

    def sortGetIsValid(self):
        """TODO

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

    def sortGetSizes(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMSwarmSortGetSizes

        """
        cdef PetscInt ncells = asInt(0)
        cdef PetscInt npoints = asInt(0)
        CHKERR( DMSwarmSortGetSizes(self.dm, &ncells, &npoints) )
        return (toInt(ncells), toInt(npoints))

    def projectFields(self, fieldnames, reuse=False):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
