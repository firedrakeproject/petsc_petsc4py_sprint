# --------------------------------------------------------------------

class MatPartitioningType(object):
    PARTITIONINGCURRENT  = S_(MATPARTITIONINGCURRENT)
    PARTITIONINGAVERAGE  = S_(MATPARTITIONINGAVERAGE)
    PARTITIONINGSQUARE   = S_(MATPARTITIONINGSQUARE)
    PARTITIONINGPARMETIS = S_(MATPARTITIONINGPARMETIS)
    PARTITIONINGCHACO    = S_(MATPARTITIONINGCHACO)
    PARTITIONINGPARTY    = S_(MATPARTITIONINGPARTY)
    PARTITIONINGPTSCOTCH = S_(MATPARTITIONINGPTSCOTCH)
    PARTITIONINGHIERARCH = S_(MATPARTITIONINGHIERARCH)

# --------------------------------------------------------------------

cdef class MatPartitioning(Object):
    """TODO."""

    Type = MatPartitioningType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.part
        self.part = NULL

    def __call__(self):
        return self.getValue()

    def view(self, Viewer viewer=None):
        """Prints the partitioning data structure.

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningView(MatPartitioning part, PetscViewer viewer)

Collective

Input Parameters
part - the partitioning context
viewer - optional visualization context
Note
The available visualization contexts include

PETSC_VIEWER_STDOUT_SELF - standard output (default)
PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where only the first processor opens the file. All other processors send their data to the first processor to print.
The user can open alternative visualization contexts with

PetscViewerASCIIOpen() - output to a specified file
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningView

        """
        assert self.obj != NULL
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( MatPartitioningView(self.part, vwr) )

    def destroy(self):
        """Destroys the partitioning context.

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningDestroy(MatPartitioning *part)

Collective

Input Parameter
part - the partitioning context
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningDestroy

        """
        CHKERR( MatPartitioningDestroy(&self.part) )
        return self

    def create(self, comm=None):
        """Creates a partitioning context.

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningCreate(MPI_Comm comm, MatPartitioning *newp)

Collective

Input Parameter
comm - MPI communicator
Output Parameter
newp - location to put the context
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        CHKERR( MatPartitioningCreate(ccomm, &self.part) )
        return self

    def setType(self, matpartitioning_type):
        """Sets the type of partitioner to use

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningSetType(MatPartitioning part, MatPartitioningType type)

Collective

Input Parameters
part - the partitioning context.
type - a known method
Options Database Key
-mat_partitioning_type - (for instance, parmetis), use -help for a list of available methods or see MatPartitioningType
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningSetType

        """
        cdef PetscMatPartitioningType cval = NULL
        matpartitioning_type = str2bytes(matpartitioning_type, &cval)
        CHKERR( MatPartitioningSetType(self.part, cval) )

    def getType(self):
        """Gets the Partitioning method type and name (as a string) from the partitioning context.

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningGetType(MatPartitioning partitioning, MatPartitioningType *type)

Not Collective

Input Parameter
partitioning - the partitioning context
Output Parameter
type - partitioner type
.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningGetType

        """
        cdef PetscMatPartitioningType cval = NULL
        CHKERR( MatPartitioningGetType(self.part, &cval) )
        return bytes2str(cval)

    def setFromOptions(self):
        """Sets various partitioning options from the options database for the partitioning object

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningSetFromOptions(MatPartitioning part)

Collective

Input Parameter
part - the partitioning context.
Options Database Keys
-mat_partitioning_type - (for instance, parmetis), use -help for a list of available methods
-mat_partitioning_nparts - number of subgraphs
Note
If the partitioner has not been set by the user it uses one of the installed partitioner such as ParMetis. If there are no installed partitioners it does no repartioning.

.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningSetFromOptions

        """
        CHKERR( MatPartitioningSetFromOptions(self.part) )

    def setAdjacency(self, Mat adj):
        """Sets the adjacency graph (matrix) of the thing to be partitioned.

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningSetAdjacency(MatPartitioning part, Mat adj)

Collective

Input Parameters
part - the partitioning context
adj - the adjacency matrix, this can be any MatType but the natural representation is MATMPIADJ
.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningSetAdjacency

        """
        CHKERR( MatPartitioningSetAdjacency(self.part, adj.mat) )

    def apply(self, IS partitioning):
        """Gets a partitioning for the graph represented by a sparse matrix.

Synopsis
#include "petscmat.h"
PetscErrorCode MatPartitioningApply(MatPartitioning matp, IS *partitioning)

Collective

Input Parameter
matp - the matrix partitioning object
Output Parameter
partitioning - the partitioning. For each local node this tells the processor number that that node is assigned to.
Options Database Keys
-mat_partitioning_type - set the partitioning package or algorithm to use
-mat_partitioning_view - display information about the partitioning object
The user can define additional partitionings; see MatPartitioningRegister().

.

         collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.MatPartitioningApply

        """
        CHKERR( MatPartitioningApply(self.part, &partitioning.iset) )

# --------------------------------------------------------------------

del MatPartitioningType

# --------------------------------------------------------------------
