# --------------------------------------------------------------------

class PartitionerType(object):
    PARMETIS        = S_(PETSCPARTITIONERPARMETIS)
    PTSCOTCH        = S_(PETSCPARTITIONERPTSCOTCH)
    CHACO           = S_(PETSCPARTITIONERCHACO)
    SIMPLE          = S_(PETSCPARTITIONERSIMPLE)
    SHELL           = S_(PETSCPARTITIONERSHELL)
    GATHER          = S_(PETSCPARTITIONERGATHER)
    MATPARTITIONING = S_(PETSCPARTITIONERMATPARTITIONING)

# --------------------------------------------------------------------

cdef class Partitioner(Object):
    """TODO."""

    Type = PartitionerType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.part
        self.part = NULL

    def view(self, Viewer viewer=None):
        """Views a PetscPartitioner

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerView(PetscPartitioner part, PetscViewer v)

Collective

Input Parameters
part - the PetscPartitioner object to view
v - the viewer
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerView

        """

        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscPartitionerView(self.part, vwr) )

    def destroy(self):
        """Destroys a PetscPartitioner object

.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerDestroy

        """
        CHKERR( PetscPartitionerDestroy(&self.part) )
        return self

    def create(self, comm=None):
        """Creates an empty PetscPartitioner object. The type can then be set with PetscPartitionerSetType().

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerCreate(MPI_Comm comm, PetscPartitioner *part)

Collective

Input Parameter
comm - The communicator for the PetscPartitioner object
Output Parameter
part - The PetscPartitioner object
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscPartitioner newpart = NULL
        CHKERR( PetscPartitionerCreate(ccomm, &newpart) )
        PetscCLEAR(self.obj); self.part = newpart
        return self

    def setType(self, part_type):
        """Builds a particular PetscPartitioner

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerSetType(PetscPartitioner part, PetscPartitionerType name)

Collective

Input Parameters
part - The PetscPartitioner object
name - The kind of partitioner
Options Database Key
-petscpartitioner_type - Sets the PetscPartitioner type
Note
 PETSCPARTITIONERCHACO    - The Chaco partitioner (--download-chaco)
 PETSCPARTITIONERPARMETIS - The ParMetis partitioner (--download-parmetis)
 PETSCPARTITIONERSHELL    - A shell partitioner implemented by the user
 PETSCPARTITIONERSIMPLE   - A simple partitioner that divides cells into equal, contiguous chunks
 PETSCPARTITIONERGATHER   - Gathers all cells onto process 0

.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerSetType

        """
        cdef PetscPartitionerType cval = NULL
        part_type = str2bytes(part_type, &cval)
        CHKERR( PetscPartitionerSetType(self.part, cval) )

    def getType(self):
        """Gets the PetscPartitioner type name (as a string) from the object.

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerGetType(PetscPartitioner part, PetscPartitionerType *name)

Not Collective

Input Parameter
part - The PetscPartitioner
Output Parameter
name - The PetscPartitioner type name
.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerGetType

        """
        cdef PetscPartitionerType cval = NULL
        CHKERR( PetscPartitionerGetType(self.part, &cval) )
        return bytes2str(cval)

    def setFromOptions(self):
        """sets parameters in a PetscPartitioner from the options database

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerSetFromOptions(PetscPartitioner part)

Collective

Input Parameter
part - the PetscPartitioner object to set options for
Options Database Keys
-petscpartitioner_type - Sets the PetscPartitioner type; use -help for a list of available types
-petscpartitioner_use_vertex_weights - Uses weights associated with the graph vertices
-petscpartitioner_view_graph - View the graph each time PetscPartitionerPartition is called. Viewer can be customized, see PetscOptionsGetViewer()
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerSetFromOptions

        """
        CHKERR( PetscPartitionerSetFromOptions(self.part) )

    def setUp(self):
        """Construct data structures for the PetscPartitioner

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerSetUp(PetscPartitioner part)

Collective

Input Parameter
part - the PetscPartitioner object to setup
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerSetUp

        """
        CHKERR( PetscPartitionerSetUp(self.part) )

    def reset(self):
        """Resets data structures for the PetscPartitioner

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerReset(PetscPartitioner part)

Collective

Input Parameter
part - the PetscPartitioner object to reset
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerReset

        """
        CHKERR( PetscPartitionerReset(self.part) )

    def setShellPartition(self, numProcs, sizes=None, points=None):
        """Set an artificial partition for a mesh

Synopsis
#include "petscpartitioner.h"
PetscErrorCode PetscPartitionerShellSetPartition(PetscPartitioner part, PetscInt size, const PetscInt sizes[], const PetscInt points[])

Collective

Input Parameters
part - The PetscPartitioner
size - The number of partitions
sizes - array of length size (or NULL) providing the number of points in each partition
points - array of length sum(sizes) (may be NULL iff sizes is NULL), a permutation of the points that groups those assigned to each partition in order (i.e., partition 0 first, partition 1 next, etc.)
Note
It is safe to free the sizes and points arrays after use in this routine.

.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscPartitionerShellSetPartition

        """
        cdef PetscInt cnumProcs = asInt(numProcs)
        cdef PetscInt *csizes = NULL
        cdef PetscInt *cpoints = NULL
        cdef PetscInt nsize = 0
        if sizes is not None:
            sizes = iarray_i(sizes, &nsize, &csizes)
            if nsize != cnumProcs:
                raise ValueError("sizes array should have %d entries (has %d)" %
                                 numProcs, toInt(nsize))
            if points is None:
                raise ValueError("Must provide both sizes and points arrays")
        if points is not None:
            points = iarray_i(points, NULL, &cpoints)
        CHKERR( PetscPartitionerShellSetPartition(self.part, cnumProcs,
                                                  csizes, cpoints) )

# --------------------------------------------------------------------

del PartitionerType

# --------------------------------------------------------------------
