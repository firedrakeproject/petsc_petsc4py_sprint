
cdef class DMInterpolation:

    cdef PetscDMInterpolation dminterp

    def __cinit__(self):
        self.dminterp = NULL

    def __dealloc__(self):
        self.destroy()

    def create(self, comm=None):
        """Creates a DMInterpolationInfo context.

        Collective.

        Parameters
        ----------
comm - the communicator


        See also
        --------
        petsc.DMInterpolationCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_SELF)
        cdef PetscDMInterpolation new = NULL
        CHKERR( DMInterpolationCreate(ccomm, &new) )
        self.dminterp = new

    def destroy(self):
        """Destroys a DMInterpolationInfo context.

        Collective.

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
        """Using the input from dm and x, calculates interpolated field values at the interpolation points.

ctx - The DMInterpolationInfo context
dm - The DM
x - The local vector containing the field to be interpolated
Output Parameter
v - The vector containing the interpolated values

A suitable v can be obtained using DMInterpolationGetVector().

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
        """Gets a Vec with the coordinates of each interpolation point

Output Parameter
coordinates - the coordinates of interpolation points
Note
The local vector entries correspond to interpolation points lying on this process, according to the associated DM. This is a borrowed vector that the user should not destroy.



        Collective.

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
        """Gets the spatial dimension for the interpolation context

.

        Not collective.

Input Parameter
ctx - the context
Output Parameter
dim - the spatial dimension

        See also
        --------
        petsc.DMInterpolationGetDim

        """
        cdef PetscInt cdim = 0
        CHKERR( DMInterpolationGetDim(self.dminterp, &cdim) )
        return toInt(cdim)

    def getDof(self):
        """Gets the number of fields interpolated at a point for the interpolation context

.

        Not collective.

Input Parameter
ctx - the context
Output Parameter
dof - the number of fields

        See also
        --------
        petsc.DMInterpolationGetDof

        """
        cdef PetscInt cdof = 0
        CHKERR( DMInterpolationGetDof(self.dminterp, &cdof) )
        return toInt(cdof)

    def setDim(self, dim):
        """Sets the spatial dimension for the interpolation context

.

        Not collective.

Input Parameters
ctx - the context
dim - the spatial dimension

        See also
        --------
        petsc.DMInterpolationSetDim

        """
        cdef PetscInt cdim = asInt(dim)
        CHKERR( DMInterpolationSetDim(self.dminterp, cdim) )

    def setDof(self, dof):
        """Sets the number of fields interpolated at a point for the interpolation context

.

        Not collective.

        Parameters
        ----------
ctx - the context
dof - the number of fields

        See also
        --------
        petsc.DMInterpolationSetDof

        """
        cdef PetscInt cdof = asInt(dof)
        CHKERR( DMInterpolationSetDof(self.dminterp, cdof) )

    def setUp(self, DM dm, redundantPoints=False, ignoreOutsideDomain=False):
        """Compute spatial indices for point location during interpolation

.

        Collective.

        Parameters
        ----------
ctx - the context
dm - the DM for the function space used for interpolation
redundantPoints - If PETSC_TRUE, all processes are passing in the same array of points. Otherwise, points need to be communicated among processes.
ignoreOutsideDomain - If PETSC_TRUE, ignore points outside the domain, otherwise return an error

        See also
        --------
        petsc.DMInterpolationSetUp

        """
        cdef PetscBool credundantPoints = asBool(redundantPoints)
        cdef PetscBool cignoreOutsideDomain = asBool(ignoreOutsideDomain)
        CHKERR( DMInterpolationSetUp(self.dminterp, dm.dm, credundantPoints, cignoreOutsideDomain) )

    def getVector(self):
        """Gets a Vec which can hold all the interpolated field values

.

        Collective.

Input Parameter
ctx - the context
Output Parameter
v - a vector capable of holding the interpolated field values
Note
This vector should be returned using DMInterpolationRestoreVector().

        See also
        --------
        petsc.DMInterpolationGetVector

        """
        cdef Vec vec = Vec()
        CHKERR( DMInterpolationGetVector(self.dminterp, &vec.vec))
        return vec

    def restoreVector(self, Vec vec):
        """Returns a Vec which can hold all the interpolated field values

.

        Collective.

Input Parameters
ctx - the context
v - a vector capable of holding the interpolated field values

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