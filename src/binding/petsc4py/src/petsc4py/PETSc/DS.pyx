# --------------------------------------------------------------------

class DSType(object):
    BASIC = S_(PETSCDSBASIC)

# --------------------------------------------------------------------

cdef class DS(Object):

    Type = DSType

    #

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ds
        self.ds  = NULL

    def view(self, Viewer viewer=None):
        """Views a PetscDS

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSView(PetscDS prob, PetscViewer v)

Collective

Input Parameters
prob - the PetscDS object to view
v - the viewer


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscDSView(self.ds, vwr) )

    def destroy(self):
        """Destroys a PetscDS object

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSDestroy(PetscDS *ds)

Collective

Input Parameter
prob - the PetscDS object to destroy


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSDestroy

        """
        CHKERR( PetscDSDestroy(&self.ds) )
        return self

    def create(self, comm=None):
        """Creates an empty PetscDS object. The type can then be set with PetscDSSetType().

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSCreate(MPI_Comm comm, PetscDS *ds)

Collective

Input Parameter
comm - The communicator for the PetscDS object
Output Parameter
ds - The PetscDS object


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDS newds = NULL
        CHKERR( PetscDSCreate(ccomm, &newds) )
        PetscCLEAR(self.obj); self.ds = newds
        return self

    def setType(self, ds_type):
        """Builds a particular PetscDS

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSSetType(PetscDS prob, PetscDSType name)

Collective; No Fortran Support

Input Parameters
prob - The PetscDS object
name - The PetscDSType
Options Database Key
-petscds_type - Sets the PetscDS type; use -help for a list of available types


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSSetType

        """
        cdef PetscDSType cval = NULL
        ds_type = str2bytes(ds_type, &cval)
        CHKERR( PetscDSSetType(self.ds, cval) )

    def getType(self):
        """Gets the PetscDSType name (as a string) from the PetscDS

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetType(PetscDS prob, PetscDSType *name)

Not Collective; No Fortran Support

Input Parameter
prob - The PetscDS
Output Parameter
name - The PetscDSType name


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetType

        """
        cdef PetscDSType cval = NULL
        CHKERR( PetscDSGetType(self.ds, &cval) )
        return bytes2str(cval)

    def setFromOptions(self):
        """sets parameters in a PetscDS from the options database

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSSetFromOptions(PetscDS prob)

Collective

Input Parameter
prob - the PetscDS object to set options for
Options Database Keys
-petscds_type - Set the PetscDS type
-petscds_view - View the PetscDS
-petscds_jac_pre - Turn formation of a separate Jacobian preconditioner on or off
-bc_ - Specify a list of label ids for a boundary condition
-bc__comp - Specify a list of field components to constrain for a boundary condition


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSSetFromOptions, petsc_options

        """
        CHKERR( PetscDSSetFromOptions(self.ds) )

    def setUp(self):
        """Construct data structures for the PetscDS

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSSetUp(PetscDS prob)

Collective

Input Parameter
prob - the PetscDS object to setup


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSSetUp

        """
        CHKERR( PetscDSSetUp(self.ds) )
        return self

    #

    def getSpatialDimension(self):
        """Returns the spatial dimension of the PetscDS, meaning the topological dimension of the discretizations

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetSpatialDimension(PetscDS prob, PetscInt *dim)

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
dim - The spatial dimension


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetSpatialDimension

        """
        cdef PetscInt dim = 0
        CHKERR( PetscDSGetSpatialDimension(self.ds, &dim) )
        return toInt(dim)

    def getCoordinateDimension(self):
        """Returns the coordinate dimension of the PetscDS, meaning the dimension of the space into which the discretiaztions are embedded

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetCoordinateDimension(PetscDS prob, PetscInt *dimEmbed)

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
dimEmbed - The coordinate dimension


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetCoordinateDimension

        """
        cdef PetscInt dim = 0
        CHKERR( PetscDSGetCoordinateDimension(self.ds, &dim) )
        return toInt(dim)

    def getNumFields(self):
        """Returns the number of fields in the PetscDS

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetNumFields(PetscDS prob, PetscInt *Nf)

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
Nf - The number of fields


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetNumFields

        """
        cdef PetscInt nf = 0
        CHKERR( PetscDSGetNumFields(self.ds, &nf) )
        return toInt(nf)

    def getFieldIndex(self, Object disc):
        """Returns the index of the given field

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetFieldIndex(PetscDS prob, PetscObject disc, PetscInt *f)

Not Collective

Input Parameters
prob - The PetscDS object
disc - The discretization object
Output Parameter
f - The field number


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetFieldIndex

        """
        cdef PetscInt field = 0
        CHKERR( PetscDSGetFieldIndex(self.ds, disc.obj[0], &field) )
        return toInt(field)

    def getTotalDimensions(self):
        """Returns the total size of the approximation space for this system

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetTotalDimension(PetscDS prob, PetscInt *dim)

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
dim - The total problem dimension


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetTotalDimension

        """
        cdef PetscInt tdim = 0
        CHKERR( PetscDSGetTotalDimension(self.ds, &tdim) )
        return toInt(tdim)

    def getTotalComponents(self):
        """Returns the total number of components in this system

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetTotalComponents(PetscDS prob, PetscInt *Nc)

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
dim - The total number of components


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetTotalComponents

        """
        cdef PetscInt tcmp = 0
        CHKERR( PetscDSGetTotalComponents(self.ds, &tcmp) )
        return toInt(tcmp)

    def getDimensions(self):
        """Returns the size of the approximation space for each field on an evaluation point

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetDimensions(PetscDS prob, PetscInt *dimensions[])

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
dimensions - The number of dimensions


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetDimensions

        """
        cdef PetscInt nf = 0, *dims = NULL
        CHKERR( PetscDSGetNumFields(self.ds, &nf) )
        CHKERR( PetscDSGetDimensions(self.ds, &dims) )
        return array_i(nf, dims)

    def getComponents(self):
        """Returns the number of components for each field on an evaluation point

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSGetComponents(PetscDS prob, PetscInt *components[])

Not Collective

Input Parameter
prob - The PetscDS object
Output Parameter
components - The number of components


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSGetComponents

        """
        cdef PetscInt nf = 0, *cmps = NULL
        CHKERR( PetscDSGetNumFields(self.ds, &nf) )
        CHKERR( PetscDSGetComponents(self.ds, &cmps) )
        return array_i(nf, cmps)

    def setDiscretisation(self, f, disc):
        """Sets the discretization object for the given field

Synopsis
#include "petscds.h"
PetscErrorCode PetscDSSetDiscretization(PetscDS prob, PetscInt f, PetscObject disc)

Not Collective

Input Parameters
prob - The PetscDS object
f - The field number
disc - The discretization object


        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDSSetDiscretization

        """
        cdef PetscInt cf = asInt(f)
        cdef FE fe = disc
        CHKERR( PetscDSSetDiscretization(self.ds, cf, <PetscObject> fe.fe) )



# --------------------------------------------------------------------

del DSType

# --------------------------------------------------------------------
