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
    """A graph partitioner."""

    Type = PartitionerType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.part
        self.part = NULL

    def view(self, Viewer viewer=None) -> None:
        """View the partitioner.

        Collective.

        Parameters
        ----------
        viewer
            A `Viewer` to display the graph.

        See also
        --------
        petsc.PetscPartitionerView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscPartitionerView(self.part, vwr) )

    def destroy(self) -> Self:
        """Destroy the partitioner object.

        Collective.

        See also
        --------
        petsc.PetscPartitionerDestroy

        """
        CHKERR( PetscPartitionerDestroy(&self.part) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create an empty partitioner object.

        The type can be set with `setType`.

        Collective.

        Parameters
        ----------
        comm
            The MPI communicator.

        See also
        --------
        petsc.PetscPartitionerCreate, setType

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscPartitioner newpart = NULL
        CHKERR( PetscPartitionerCreate(ccomm, &newpart) )
        PetscCLEAR(self.obj); self.part = newpart
        return self

    def setType(self, part_type: Type | str) -> None:
        """Build a particular type of the partitioner.

        Collective.

        Parameters
        ----------
        part_type
            The kind of partitioner.

        See also
        --------
        petsc.PetscPartitionerSetType, getType

        """
        cdef PetscPartitionerType cval = NULL
        part_type = str2bytes(part_type, &cval)
        CHKERR( PetscPartitionerSetType(self.part, cval) )

    def getType(self) -> Type:
        """Return the partitioner type.

        Not collective.

        See also
        --------
        petsc.PetscPartitionerGetType, setType

        """
        cdef PetscPartitionerType cval = NULL
        CHKERR( PetscPartitionerGetType(self.part, &cval) )
        return bytes2str(cval)

    def setFromOptions(self) -> None:
        """Set parameters in the partitioner from the options database.

        Collective.

        See also
        --------
        petsc.PetscPartitionerSetFromOptions, petsc_options

        """
        CHKERR( PetscPartitionerSetFromOptions(self.part) )

    def setUp(self) -> None:
        """Construct data structures for the partitioner.

        Collective.

        See also
        --------
        petsc.PetscPartitionerSetUp

        """
        CHKERR( PetscPartitionerSetUp(self.part) )

    def reset(self) -> None:
        """Reset data structures of the partitioner.

        Collective.

        See also
        --------
        petsc.PetscPartitionerReset

        """
        CHKERR( PetscPartitionerReset(self.part) )

    def setShellPartition(
        self,
        numProcs: int,
        sizes: Sequence[int] | None = None,
        points: Sequence[int] | None = None,
    ) -> None:
        """Set a custom partition for a mesh.

        Collective.

        Parameters
        ----------
        sizes
            The number of points in each partition.
        points
            A permutation of the points that groups those assigned to each
            partition in order (i.e., partition 0 first, partition 1 next,
            etc.).

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
