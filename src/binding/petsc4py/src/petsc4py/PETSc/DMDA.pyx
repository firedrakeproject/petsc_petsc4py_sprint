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

    StencilType       = DMDAStencilType
    InterpolationType = DMDAInterpolationType
    ElementType       = DMDAElementType

    #

    def create(self, dim=None, dof=None,
               sizes=None, proc_sizes=None, boundary_type=None,
               stencil_type=None, stencil_width=None,
               bint setup=True, ownership_ranges=None, comm=None):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDACreateND

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

    def duplicate(self, dof=None, boundary_type=None,
                  stencil_type=None, stencil_width=None):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDACreateND, petsc.DMSetUp

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

    def setDim(self, dim):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.TODO

        """
        return self.setDimension(dim)

    def getDim(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.TODO

        """
        return self.getDimension()

    def setDof(self, dof):
        """Set the number of degrees of freedom per vertex

Not Collective

Input Parameters
da - The DMDA
dof - Number of degrees of freedom


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetDof

        """
        cdef PetscInt ndof = asInt(dof)
        CHKERR( DMDASetDof(self.dm, ndof) )

    def getDof(self):
        """
        ------------
        getINFO
        Gets information about a given distributed array.

Not Collective

Input Parameter
da - the distributed array
Output Parameters
dim - dimension of the distributed array (1, 2, or 3)
M - global dimension in first direction of the array
N - global dimension in second direction of the array
P - global dimension in third direction of the array
m - corresponding number of procs in first dimension
n - corresponding number of procs in second dimension
p - corresponding number of procs in third dimension
dof - number of degrees of freedom per node
s - stencil width
bx - type of ghost nodes at boundary in first dimension
by - type of ghost nodes at boundary in second dimension
bz - type of ghost nodes at boundary in third dimension
st - stencil type, either DMDA_STENCIL_STAR or DMDA_STENCIL_BOX
Note
Use NULL (NULL_INTEGER in Fortran) in place of any output parameter that is not of interest.



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setSizes(self, sizes):
        """Set the number of grid points in the three dimensional directions

Logically Collective

Input Parameters
da - the DMDA
M - the global X size
N - the global Y size
P - the global Z size
Developer Note
Since the dimension may not yet have been set the code cannot error check for non-positive Y and Z number of grid points



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getSizes(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setProcSizes(self, proc_sizes):
        """Set the number of processes in each dimension

Logically Collective

Input Parameters
da - the DMDA
m - the number of X procs (or PETSC_DECIDE)
n - the number of Y procs (or PETSC_DECIDE)
p - the number of Z procs (or PETSC_DECIDE)


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getProcSizes(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setBoundaryType(self, boundary_type):
        """Set the type of ghost nodes on domain boundaries.

Not Collective

Input Parameters
da - The DMDA
bx,by,bz - One of DM_BOUNDARY_NONE, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetBoundaryType

        """
        cdef PetscDMBoundaryType btx = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType bty = DM_BOUNDARY_NONE
        cdef PetscDMBoundaryType btz = DM_BOUNDARY_NONE
        asBoundary(boundary_type, &btx, &bty, &btz)
        CHKERR( DMDASetBoundaryType(self.dm, btx, bty, btz) )

    def getBoundaryType(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setStencilType(self, stencil_type):
        """Set the type of the communication stencil

Logically Collective

Input Parameters
da - The DMDA
stype - The stencil type, use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetStencilType

        """
        cdef PetscDMDAStencilType stype = asStencil(stencil_type)
        CHKERR( DMDASetStencilType(self.dm, stype) )

    def getStencilType(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setStencilWidth(self, stencil_width):
        """Set the width of the communication stencil

Logically Collective

Input Parameters
da - The DMDA
width - The stencil width


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetStencilWidth

        """
        cdef PetscInt swidth = asInt(stencil_width)
        CHKERR( DMDASetStencilWidth(self.dm, swidth) )

    def getStencilWidth(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setStencil(self, stencil_type, stencil_width):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetStencilType, petsc.DMDASetStencilWidth

        """
        cdef PetscDMDAStencilType stype = asStencil(stencil_type)
        cdef PetscInt swidth = asInt(stencil_width)
        CHKERR( DMDASetStencilType(self.dm, stype) )
        CHKERR( DMDASetStencilWidth(self.dm, swidth) )

    def getStencil(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getRanges(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getGhostRanges(self):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getOwnershipRanges(self):
        """Return the ranges of indices in the x, y and z direction that are owned by each process

Not Collective

Input Parameter
da - the DMDA object
Output Parameters
lx - ownership along x direction (optional)
ly - ownership along y direction (optional)
lz - ownership along z direction (optional)
Note
These correspond to the optional final arguments passed to DMDACreate(), DMDACreate2d(), DMDACreate3d()

In C you should not free these arrays, nor change the values in them. They will only have valid values while the DMDA they came from still exists (has not been destroyed).

These numbers are NOT multiplied by the number of dof per node.



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getCorners(self):
        """Return the global (x,y,z) indices of the lower left corner and size of the local region, excluding ghost points.

Not Collective

Input Parameter
da - the distributed array
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

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def getGhostCorners(self):
        """Return the global (x,y,z) indices of the lower left corner and size of the local region, including ghost points.

Not Collective

Input Parameter
da - the distributed array
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

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setFieldName(self, field, name):
        """Set the names of individual field components in multicomponent vectors associated with a DMDA.

Logically Collective; name must contain a common value

Input Parameters
da - the distributed array
nf - field number for the DMDA (0, 1, … dof-1), where dof indicates the number of degrees of freedom per node within the DMDA
names - the name of the field (component)
Note
It must be called after having called DMSetUp().



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetFieldName

        """
        cdef PetscInt ival = asInt(field)
        cdef const char *cval = NULL
        name = str2bytes(name, &cval)
        CHKERR( DMDASetFieldName(self.dm, ival, cval) )

    def getFieldName(self, field):
        """Return the names of individual field components in multicomponent vectors associated with a DMDA.

Not Collective; name will contain a common value

Input Parameters
da - the distributed array
nf - field number for the DMDA (0, 1, … dof-1), where dof indicates the number of degrees of freedom per node within the DMDA
Output Parameter
names - the name of the field (component)
Note
It must be called after having called DMSetUp().



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDAGetFieldName

        """
        cdef PetscInt ival = asInt(field)
        cdef const char *cval = NULL
        CHKERR( DMDAGetFieldName(self.dm, ival, &cval) )
        return bytes2str(cval)

    #

    def getVecArray(self, Vec vec):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc._DMDA_Vec_array

        """
        return _DMDA_Vec_array(self, vec)

    #

    def setUniformCoordinates(self,
                              xmin=0, xmax=1,
                              ymin=0, ymax=1,
                              zmin=0, zmax=1):
        """Set a DMDA coordinates to be a uniform grid

Collective

Input Parameters
da - the distributed array object
xmin,xmax - extremes in the x direction
ymin,ymax - extremes in the y direction (value ignored for 1 dimensional problems)
zmin,zmax - extremes in the z direction (value ignored for 1 or 2 dimensional problems)


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setCoordinateName(self, index, name):
        """Set the name of the coordinate directions associated with a DMDA, for example “x” or “y”

Logically Collective; name must contain a common value; No Fortran Support

Input Parameters
dm - the DMDA
nf - coordinate number for the DMDA (0, 1, … dim-1),
name - the name of the coordinate
Note
Must be called after having called DMSetUp().



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetCoordinateName

        """
        cdef PetscInt ival = asInt(index)
        cdef const char *cval = NULL
        name = str2bytes(name, &cval)
        CHKERR( DMDASetCoordinateName(self.dm, ival, cval) )

    def getCoordinateName(self, index):
        """Return the name of a coordinate direction associated with a DMDA.

Not Collective; name will contain a common value; No Fortran Support

Input Parameters
dm - the DMDA
nf - number for the DMDA (0, 1, … dim-1)
Output Parameter
names - the name of the coordinate direction
Note
It must be called after having called DMSetUp().



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDAGetCoordinateName

        """
        cdef PetscInt ival = asInt(index)
        cdef const char *cval = NULL
        CHKERR( DMDAGetCoordinateName(self.dm, ival, &cval) )
        return bytes2str(cval)

    #

    def createNaturalVec(self):
        """Create a parallel PETSc vector that will hold vector values in the natural numbering, rather than in the PETSc parallel numbering associated with the DMDA.

Collective

Input Parameter
da - the distributed array
Output Parameter
g - the distributed global vector
Notes
The output parameter, g, is a regular PETSc vector that should be destroyed with a call to VecDestroy() when usage is finished.

The number of local entries in the vector on each process is the same as in a vector created with DMCreateGlobalVector().



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDACreateNaturalVector

        """
        cdef Vec vn = Vec()
        CHKERR( DMDACreateNaturalVector(self.dm, &vn.vec) )
        return vn

    def globalToNatural(self, Vec vg, Vec vn, addv=None):
        """Map values from the global vector to a global vector in the “natural” grid ordering. Must be followed by DMDAGlobalToNaturalEnd() to complete the exchange.

Neighbor-wise Collective

Input Parameters
da - the distributed array context
g - the global vector
mode - one of INSERT_VALUES or ADD_VALUES
Output Parameter
l - the natural ordering values
l - the natural ordering values
l - the global values in the natural ordering

Notes
The global and natural vectors used here need not be the same as those obtained from DMCreateGlobalVector() and DMDACreateNaturalVector(), BUT they must have the same parallel data layout; they could, for example, be obtained with VecDuplicate() from the DMDA originating vectors.

You must call DMDACreateNaturalVector() before using this routine



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDAGlobalToNaturalBegin, petsc.DMDAGlobalToNaturalEnd

        """
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMDAGlobalToNaturalBegin(self.dm, vg.vec, im, vn.vec) )
        CHKERR( DMDAGlobalToNaturalEnd  (self.dm, vg.vec, im, vn.vec) )

    def naturalToGlobal(self, Vec vn, Vec vg, addv=None):
        """Map values from a global vector in the “natural” ordering to a global vector in the PETSc DMDA grid ordering. Must be followed by DMDANaturalToGlobalEnd() to complete the exchange.

Neighbor-wise Collective

Input Parameters
da - the distributed array context
g - the global vector in a natural ordering
mode - one of INSERT_VALUES or ADD_VALUES
Output Parameter
l - the values in the DMDA ordering
Notes
The global and natural vectors used here need not be the same as those obtained from DMCreateGlobalVector() and DMDACreateNaturalVector(), BUT they must have the same parallel data layout; they could, for example, be obtained with VecDuplicate() from the DMDA originating vectors.



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDANaturalToGlobalBegin, petsc.DMDANaturalToGlobalEnd

        """
        cdef PetscInsertMode im = insertmode(addv)
        CHKERR( DMDANaturalToGlobalBegin(self.dm, vn.vec, im, vg.vec) )
        CHKERR( DMDANaturalToGlobalEnd  (self.dm, vn.vec, im, vg.vec) )

    #

    def getAO(self):
        """Return the application ordering context for a distributed array.

Collective

Input Parameter
da - the distributed array
Output Parameter
ao - the application ordering context for DMDA
Notes
In this case, the AO maps to the natural grid ordering that would be used for the DMDA if only 1 processor were employed (ordering most rapidly in the x-direction, then y, then z). Multiple degrees of freedom are numbered for each node (rather than 1 component for the whole grid, then the next component, etc.)

Do NOT call AODestroy() on the ao returned by this function.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDAGetAO

        """
        cdef AO ao = AO()
        CHKERR( DMDAGetAO(self.dm, &ao.ao) )
        PetscINCREF(ao.obj)
        return ao

    def getScatter(self):
        """Return the global-to-local, and local-to-local vector scatter contexts for a distributed array.

Collective

Input Parameter
da - the distributed array
Output Parameters
gtol - global-to-local scatter context (may be NULL)
ltol - local-to-local scatter context (may be NULL)
Note
The output contexts are valid only as long as the input da is valid. If you delete the da, the scatter contexts will become invalid.



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setRefinementFactor(self,
                            refine_x=2,
                            refine_y=2,
                            refine_z=2):
        """Set the ratios that the DMDA grid is refined

Logically Collective

Input Parameters
da - the DMDA object
refine_x - ratio of fine grid to coarse in x direction (2 by default)
refine_y - ratio of fine grid to coarse in y direction (2 by default)
refine_z - ratio of fine grid to coarse in z direction (2 by default)
Options Database Keys
-da_refine_x refine_x - refinement ratio in x direction
-da_refine_y rafine_y - refinement ratio in y direction
-da_refine_z refine_z - refinement ratio in z direction
-da_refine - refine the DMDA object n times when it is created.
Note
Pass PETSC_IGNORE to leave a value unchanged



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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
        """Return the ratios that the DMDA grid is refined

Not Collective

Input Parameter
da - the DMDA object
Output Parameters
refine_x - ratio of fine grid to coarse in x direction (2 by default)
refine_y - ratio of fine grid to coarse in y direction (2 by default)
refine_z - ratio of fine grid to coarse in z direction (2 by default)
Note
Pass NULL for values you do not need



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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

    def setInterpolationType(self, interp_type):
        """Set the type of interpolation that will be returned by DMCreateInterpolation()

Logically Collective

Input Parameters
da - initial distributed array
ctype - DMDA_Q1 and DMDA_Q0 are currently the only supported forms
Note
You should call this on the coarser of the two DMDA you pass to DMCreateInterpolation()



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetInterpolationType

        """
        cdef PetscDMDAInterpolationType ival = dainterpolationtype(interp_type)
        CHKERR( DMDASetInterpolationType(self.dm, ival) )

    def getInterpolationType(self):
        """Return the type of interpolation that will be used by DMCreateInterpolation()

Not Collective

Input Parameter
da - distributed array
Output Parameter
ctype - interpolation type (DMDA_Q1 and DMDA_Q0 are currently the only supported forms)


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDAGetInterpolationType

        """
        cdef PetscDMDAInterpolationType ival = DMDA_INTERPOLATION_Q0
        CHKERR( DMDAGetInterpolationType(self.dm, &ival) )
        return <long>ival

    #

    def setElementType(self, elem_type):
        """Set the element type to be returned by DMDAGetElements()

Not Collective

Input Parameter
da - the DMDA object
Output Parameter
etype - the element type, currently either DMDA_ELEMENT_P1 or DMDA_ELEMENT_Q1


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDASetElementType

        """
        cdef PetscDMDAElementType ival = daelementtype(elem_type)
        CHKERR( DMDASetElementType(self.dm, ival) )

    def getElementType(self):
        """Return the element type to be returned by DMDAGetElements()

Not Collective

Input Parameter
da - the DMDA object
Output Parameter
etype - the element type, currently either DMDA_ELEMENT_P1 or DMDA_ELEMENT_Q1


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMDAGetElementType

        """
        cdef PetscDMDAElementType ival = DMDA_ELEMENT_Q1
        CHKERR( DMDAGetElementType(self.dm, &ival) )
        return <long>ival

    def getElements(self, elem_type=None):
        """Return an array containing the indices (in local coordinates) of all the local elements

Not Collective; No Fortran Support

Input Parameter
dm - the DMDA object
Output Parameters
nel - number of local elements
nen - number of element nodes
e - the local indices of the elements’ vertices
Notes
Call DMDARestoreElements() once you have finished accessing the elements.

Each process uniquely owns a subset of the elements. That is no element is owned by two or more processes.

If on each process you integrate over its owned elements and use ADD_VALUES in Vec/MatSetValuesLocal() then you’ll obtain the correct result.



        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
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
        def __get__(self):
            return self.getDim()

    property dof:
        def __get__(self):
            return self.getDof()

    property sizes:
        def __get__(self):
            return self.getSizes()

    property proc_sizes:
        def __get__(self):
            return self.getProcSizes()

    property boundary_type:
        def __get__(self):
            return self.getBoundaryType()

    property stencil:
        def __get__(self):
            return self.getStencil()

    property stencil_type:
        def __get__(self):
            return self.getStencilType()

    property stencil_width:
        def __get__(self):
            return self.getStencilWidth()

    property ranges:
        def __get__(self):
            return self.getRanges()

    property ghost_ranges:
        def __get__(self):
            return self.getGhostRanges()

    property corners:
        def __get__(self):
            return self.getCorners()

    property ghost_corners:
        def __get__(self):
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
