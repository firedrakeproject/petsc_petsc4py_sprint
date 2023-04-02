# --------------------------------------------------------------------

class AOType(object):
    BASIC          = S_(AOBASIC)
    ADVANCED       = S_(AOADVANCED)
    MAPPING        = S_(AOMAPPING)
    MEMORYSCALABLE = S_(AOMEMORYSCALABLE)

# --------------------------------------------------------------------

cdef class AO(Object):

    Type = AOType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ao
        self.ao = NULL

    def view(self, Viewer viewer=None):
        """Display an application ordering.

Collective

Input Parameters
ao - the application ordering context
viewer - viewer used for display
Options Database Key
-ao_view - calls AOView() at end of AOCreate()
Notes
The available visualization contexts include

PETSC_VIEWER_STDOUT_SELF - standard output (default)
PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where only the first processor opens the file. All other processors send their data to the first processor to print.
The user can open an alternative visualization context with PetscViewerASCIIOpen() - output to a specified file.



        See also
        --------
        petsc.AOView

        """
        cdef PetscViewer cviewer = NULL
        if viewer is not None: cviewer = viewer.vwr
        CHKERR( AOView(self.ao, cviewer) )

    def destroy(self):
        """Destroy an application ordering.

Collective

Input Parameter
ao - the application ordering context


        See also
        --------
        petsc.AODestroy

        """
        CHKERR( AODestroy(&self.ao) )
        return self

    def createBasic(self, app, petsc=None, comm=None):
        """Create a basic application ordering using two integer arrays.

Collective

Input Parameters
comm - MPI communicator that is to share AO
napp - size of integer arrays
myapp - integer array that defines an ordering
mypetsc - integer array that defines another ordering (may be NULL to indicate the natural ordering, that is 0,1,2,3,…)
Output Parameter
aoout - the new application ordering
Note
The arrays myapp and mypetsc must contain the all the integers 0 to napp-1 with no duplicates; that is there cannot be any “holes” in the indices. Use AOCreateMapping() or AOCreateMappingIS() if you wish to have “holes” in the indices.

Creates a basic application ordering using two IS index sets.

Collective

Input Parameters
isapp - index set that defines an ordering
ispetsc - index set that defines another ordering (may be NULL to use the natural ordering)
Output Parameter
aoout - the new application ordering
Note
The index sets isapp and ispetsc must contain the all the integers 0 to napp-1 (where napp is the length of the index sets) with no duplicates; that is there cannot be any “holes”



        See also
        --------
        petsc.AOCreateBasicIS, petsc.AOCreateBasic

        """
        cdef PetscIS isapp = NULL, ispetsc = NULL
        cdef PetscInt napp = 0, *idxapp = NULL,
        cdef PetscInt npetsc = 0, *idxpetsc = NULL
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscAO newao = NULL
        if isinstance(app, IS):
            isapp = (<IS>app).iset
            if petsc is not None:
                ispetsc = (<IS?>petsc).iset
            CHKERR( AOCreateBasicIS(isapp, ispetsc, &newao) )
        else:
            app = iarray_i(app, &napp, &idxapp)
            if petsc is not None:
                petsc = iarray_i(petsc, &npetsc, &idxpetsc)
                assert napp == npetsc, "incompatible array sizes"
            CHKERR( AOCreateBasic(ccomm, napp, idxapp, idxpetsc, &newao) )
        PetscCLEAR(self.obj); self.ao = newao
        return self

    def createMemoryScalable(self, app, petsc=None, comm=None):
        """Create a memory scalable application ordering using two integer arrays.

Collective

Input Parameters
comm - MPI communicator that is to share the AO
napp - size of integer arrays
myapp - integer array that defines an ordering
mypetsc - integer array that defines another ordering (may be NULL to indicate the natural ordering, that is 0,1,2,3,…)
Output Parameter
aoout - the new application ordering
Note
The arrays myapp and mypetsc must contain the all the integers 0 to napp-1 with no duplicates; that is there cannot be any “holes” in the indices. Use AOCreateMapping() or AOCreateMappingIS() if you wish to have “holes” in the indices. Comparing with AOCreateBasic(), this routine trades memory with message communication.

Creates a memory scalable application ordering using two index sets.

Collective

Input Parameters
isapp - index set that defines an ordering
ispetsc - index set that defines another ordering (may be NULL to use the natural ordering)
Output Parameter
aoout - the new application ordering
Notes
The index sets isapp and ispetsc must contain the all the integers 0 to napp-1 (where napp is the length of the index sets) with no duplicates; that is there cannot be any “holes”.

Comparing with AOCreateBasicIS(), this routine trades memory with message communication.



        See also
        --------
        petsc.AOCreateMemoryScalableIS, petsc.AOCreateMemoryScalable

        """
        cdef PetscIS isapp = NULL, ispetsc = NULL
        cdef PetscInt napp = 0, *idxapp = NULL,
        cdef PetscInt npetsc = 0, *idxpetsc = NULL
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscAO newao = NULL
        if isinstance(app, IS):
            isapp = (<IS>app).iset
            if petsc is not None:
                ispetsc = (<IS?>petsc).iset
            CHKERR( AOCreateMemoryScalableIS(isapp, ispetsc, &newao) )
        else:
            app = iarray_i(app, &napp, &idxapp)
            if petsc is not None:
                petsc = iarray_i(petsc, &npetsc, &idxpetsc)
                assert napp == npetsc, "incompatible array sizes"
            CHKERR( AOCreateMemoryScalable(ccomm, napp, idxapp, idxpetsc, &newao) )
        PetscCLEAR(self.obj); self.ao = newao
        return self

    def createMapping(self, app, petsc=None, comm=None):
        """Create an application mapping using two integer arrays.

Input Parameters
comm - MPI communicator that is to share the AO
napp - size of integer arrays
myapp - integer array that defines an ordering
mypetsc - integer array that defines another ordering (may be NULL to indicate the identity ordering)
Output Parameter
aoout - the new application mapping
Options Database Key
-ao_view - call AOView() at the conclusion of AOCreateMapping()
Note
The arrays myapp and mypetsc need NOT contain the all the integers 0 to napp-1, that is there CAN be “holes” in the indices. Use AOCreateBasic() or AOCreateBasicIS() if they do not have holes for better performance.


Creates an application mapping using two index sets.

Input Parameters
comm - MPI communicator that is to share AO
isapp - index set that defines an ordering
ispetsc - index set that defines another ordering, maybe NULL for identity IS
Output Parameter
aoout - the new application ordering
Options Database Key
-ao_view - call AOView() at the conclusion of AOCreateMappingIS()
Note
The index sets isapp and ispetsc need NOT contain the all the integers 0 to N-1, that is there CAN be “holes” in the indices. Use AOCreateBasic() or AOCreateBasicIS() if they do not have holes for better performance.



        See also
        --------
        petsc.AOCreateMappingIS, petsc.AOCreateMapping

        """
        cdef PetscIS isapp = NULL, ispetsc = NULL
        cdef PetscInt napp = 0, *idxapp = NULL,
        cdef PetscInt npetsc = 0, *idxpetsc = NULL
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscAO newao = NULL
        if isinstance(app, IS):
            isapp = (<IS>app).iset
            if petsc is not None:
                ispetsc = (<IS?>petsc).iset
            CHKERR( AOCreateMappingIS(isapp, ispetsc, &newao) )
        else:
            app = iarray_i(app, &napp, &idxapp)
            if petsc is not None:
                petsc = iarray_i(petsc, &npetsc, &idxpetsc)
                assert napp == npetsc, "incompatible array sizes"
            CHKERR( AOCreateMapping(ccomm, napp, idxapp, idxpetsc, &newao) )
        PetscCLEAR(self.obj); self.ao = newao
        return self

    def getType(self):
        """Return the AO type name (as a string) from the AO.

Not Collective

Input Parameter
ao - The vector
Output Parameter
type - The AO type name


        See also
        --------
        petsc.AOGetType

        """
        cdef PetscAOType cval = NULL
        CHKERR( AOGetType(self.ao, &cval) )
        return bytes2str(cval)

    def app2petsc(self, indices):
        """Map a set of integers in the application-defined ordering to the PETSc ordering.

Collective

Input Parameters
ao - the application ordering context
n - the number of integers
ia - the integers; these are replaced with their mapped value
Output Parameter
ia - the mapped integers
Notes
Any integers in ia[] that are negative are left unchanged. This allows one to convert, for example, neighbor lists that use negative entries to indicate nonexistent neighbors due to boundary conditions, etc.

Integers that are out of range are mapped to -1


Maps an index set in the application-defined ordering to the PETSc ordering.

Collective

Input Parameters
ao - the application ordering context
is - the index set; this is replaced with its mapped values
Output Parameter
is - the mapped index set
Notes
The index set cannot be of type stride or block

Any integers in is that are negative are left unchanged. This allows one to convert, for example, neighbor lists that use negative entries to indicate nonexistent neighbors due to boundary conditions, etc.




        See also
        --------
        petsc.AOApplicationToPetscIS, petsc.AOApplicationToPetsc

        """
        cdef PetscIS iset = NULL
        cdef PetscInt nidx = 0, *idx = NULL
        if isinstance(indices, IS):
            iset = (<IS>indices).iset
            CHKERR( AOApplicationToPetscIS(self.ao, iset) )
        else:
            indices = oarray_i(indices, &nidx, &idx)
            CHKERR( AOApplicationToPetsc(self.ao, nidx, idx) )
        return indices

    def petsc2app(self, indices):
        """Map a set of integers in the PETSc ordering to the application-defined ordering.

Collective

Input Parameters
ao - the application ordering context
n - the number of integers
ia - the integers; these are replaced with their mapped value
Output Parameter
ia - the mapped integers
Note
Any integers in ia[] that are negative are left unchanged. This allows one to convert, for example, neighbor lists that use negative entries to indicate nonexistent neighbors due to boundary conditions, etc.

Integers that are out of range are mapped to -1

Maps an index set in the PETSc ordering to the application-defined ordering.

Collective

Input Parameters
ao - the application ordering context
is - the index set; this is replaced with its mapped values
Output Parameter
is - the mapped index set
Notes
The index set cannot be of type stride or block

Any integers in is that are negative are left unchanged. This allows one to convert, for example, neighbor lists that use negative entries to indicate nonexistent neighbors due to boundary conditions etc.



        See also
        --------
        petsc.AOPetscToApplicationIS, petsc.AOPetscToApplication

        """
        cdef PetscIS iset = NULL
        cdef PetscInt nidx = 0, *idx = NULL
        if isinstance(indices, IS):
            iset = (<IS>indices).iset
            CHKERR( AOPetscToApplicationIS(self.ao, iset) )
        else:
            indices = oarray_i(indices, &nidx, &idx)
            CHKERR( AOPetscToApplication(self.ao, nidx, idx) )
        return indices

# --------------------------------------------------------------------

del AOType

# --------------------------------------------------------------------
