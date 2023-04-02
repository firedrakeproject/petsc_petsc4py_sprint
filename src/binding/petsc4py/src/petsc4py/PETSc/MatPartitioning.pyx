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

    def view(self, Viewer viewer=None) -> None:
        """Print the partitioning data structure.

        The available visualization contexts include
        PETSC_VIEWER_STDOUT_SELF - standard output (default)
        PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where only the
        first processor opens the file. All other processors send their data to
        the first processor to print.

        The user can open alternative visualization contexts with
        PetscViewerASCIIOpen() - output to a specified file

        Collective.

        Parameters
        ----------
        viewer - optional visualization context

        See also
        --------
        petsc.MatPartitioningView

        """
        assert self.obj != NULL
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( MatPartitioningView(self.part, vwr) )

    def destroy(self) -> Self:
        """Destroy the partitioning context.

        Collective.

        See also
        --------
        petsc.MatPartitioningDestroy

        """
        CHKERR( MatPartitioningDestroy(&self.part) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create a partitioning context.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator.

        See also
        --------
        petsc.MatPartitioningCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        CHKERR( MatPartitioningCreate(ccomm, &self.part) )
        return self

    def setType(self, matpartitioning_type: Type) -> None:
        """Set the type of partitioner to use

        Options Database Key
        -mat_partitioning_type - (for instance, parmetis), use -help for a list of available methods or see MatPartitioningType

        Collective.

        Parameters
        ----------
        type
            A known method.

        See also
        --------
        petsc.MatPartitioningSetType

        """
        cdef PetscMatPartitioningType cval = NULL
        matpartitioning_type = str2bytes(matpartitioning_type, &cval)
        CHKERR( MatPartitioningSetType(self.part, cval) )

    def getType(self) -> Type:
        """Return the Partitioning method type and name (as a string) from the partitioning context.

        Not collective.

        See also
        --------
        petsc.MatPartitioningGetType

        """
        cdef PetscMatPartitioningType cval = NULL
        CHKERR( MatPartitioningGetType(self.part, &cval) )
        return bytes2str(cval)

    def setFromOptions(self):
        """Set various partitioning options from the options database for the partitioning object

        If the partitioner has not been set by the user it uses one of the
        installed partitioner such as ParMetis. If there are no installed
        partitioners it does no repartioning.

        Collective.

        Options Database Keys
        -mat_partitioning_type - (for instance, parmetis), use -help for a list
         of available methods
        -mat_partitioning_nparts - number of subgraphs

        See also
        --------
        petsc.MatPartitioningSetFromOptions, petsc_options

        """
        CHKERR( MatPartitioningSetFromOptions(self.part) )

    def setAdjacency(self, Mat adj) -> None:
        """Set the adjacency graph (matrix) of the thing to be partitioned.

        Collective.

        Parameters
        ----------
        adj
            The adjacency matrix, this can be any MatType but the natural
            representation is MATMPIADJ.

        See also
        --------
        petsc.MatPartitioningSetAdjacency

        """
        CHKERR( MatPartitioningSetAdjacency(self.part, adj.mat) )

    def apply(self, IS partitioning) -> None:
        """Return a partitioning for the graph represented by a sparse matrix.

        ---------------------------------
        Output Parameter
        ---------------------------------
        partitioning - the partitioning. For each local node this tells the
        processor number that that node is assigned to.

        Options Database Keys
        -mat_partitioning_type - set the partitioning package or algorithm to use
        -mat_partitioning_view - display information about the partitioning object
        The user can define additional partitionings; see MatPartitioningRegister().

        Collective.

        Parameters
        ----------
        matp
            The matrix partitioning object.

        See also
        --------
        petsc.MatPartitioningApply

        """
        CHKERR( MatPartitioningApply(self.part, &partitioning.iset) )

# --------------------------------------------------------------------

del MatPartitioningType

# --------------------------------------------------------------------
