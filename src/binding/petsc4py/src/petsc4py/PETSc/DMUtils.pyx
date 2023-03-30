
cdef class DMInterpolation:

    cdef PetscDMInterpolation dminterp

    def __cinit__(self):
        self.dminterp = NULL

    def __dealloc__(self):
        self.destroy()

    def create(self, comm=None):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_SELF)
        cdef PetscDMInterpolation new = NULL
        CHKERR( DMInterpolationCreate(ccomm, &new) )
        self.dminterp = new

    def destroy(self):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationDestroy

        """
        CHKERR( DMInterpolationDestroy(&self.dminterp))

    def evaluate(self, DM dm, Vec x):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationEvaluate

        """
        cdef Vec v = Vec()
        CHKERR( DMInterpolationEvaluate(self.dminterp, dm.dm, x.vec, v.vec ) )
        return v

    def getCoordinates(self):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationGetCoordinates

        """
        cdef Vec coords = Vec()
        CHKERR( DMInterpolationGetCoordinates(self.dminterp, &coords.vec) )
        return coords

    def getDim(self):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationGetDim

        """
        cdef PetscInt cdim = 0
        CHKERR( DMInterpolationGetDim(self.dminterp, &cdim) )
        return toInt(cdim)

    def getDof(self):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationGetDof

        """
        cdef PetscInt cdof = 0
        CHKERR( DMInterpolationGetDof(self.dminterp, &cdof) )
        return toInt(cdof)

    def setDim(self, dim):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationSetDim

        """
        cdef PetscInt cdim = asInt(dim)
        CHKERR( DMInterpolationSetDim(self.dminterp, cdim) )

    def setDof(self, dof):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationSetDof

        """
        cdef PetscInt cdof = asInt(dof)
        CHKERR( DMInterpolationSetDof(self.dminterp, cdof) )

    def setUp(self, DM dm, redundantPoints=False, ignoreOutsideDomain=False):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationSetUp

        """
        cdef PetscBool credundantPoints = asBool(redundantPoints)
        cdef PetscBool cignoreOutsideDomain = asBool(ignoreOutsideDomain)
        CHKERR( DMInterpolationSetUp(self.dminterp, dm.dm, credundantPoints, cignoreOutsideDomain) )

    def getVector(self):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationGetVector

        """
        cdef Vec vec = Vec()
        CHKERR( DMInterpolationGetVector(self.dminterp, &vec.vec))
        return vec

    def restoreVector(self, Vec vec):
        """TODO.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMInterpolationRestoreVector

        """
        CHKERR( DMInterpolationRestoreVector(self.dminterp, &vec.vec) )
        return vec