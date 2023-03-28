# --------------------------------------------------------------------

class SFType(object):
    BASIC      = S_(PETSCSFBASIC)
    NEIGHBOR   = S_(PETSCSFNEIGHBOR)
    ALLGATHERV = S_(PETSCSFALLGATHERV)
    ALLGATHER  = S_(PETSCSFALLGATHER)
    GATHERV    = S_(PETSCSFGATHERV)
    GATHER     = S_(PETSCSFGATHER)
    ALLTOALL   = S_(PETSCSFALLTOALL)
    WINDOW     = S_(PETSCSFWINDOW)

# --------------------------------------------------------------------

cdef class SF(Object):
    """SF object for setting up and managing the communication of certain
    entries of arrays and Vec between MPI ranks.
    """

    Type = SFType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.sf
        self.sf  = NULL

    def __dealloc__(self):
        CHKERR( PetscSFDestroy(&self.sf) )
        self.sf = NULL

    def view(self, Viewer viewer=None) -> None:
        """View a star forest.

        Collective.

        Parameters
        ----------
        viewer
            A `Viewer` to display the graph.

        See also
        --------
        petsc.PetscSFView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscSFView(self.sf, vwr) )

    def destroy(self) -> Self:
        """Destroy star forest.

        Collective.

        See also
        --------
        petsc.PetscSFDestroy

        """
        CHKERR( PetscSFDestroy(&self.sf) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create a star forest communication context.

        Collective.

        Parameters
        ----------
        comm
            The communicator on which the star forest will operate.

        See also
        --------
        petsc.PetscSFCreate

        """

        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscSF newsf = NULL
        CHKERR( PetscSFCreate(ccomm, &newsf) )
        PetscCLEAR(self.obj); self.sf = newsf
        return self

    def setType(self, sf_type: SF.Type | str) -> None:
        """Set the communication implementation.

        Collective.

        Parameters
        ----------
        sf_type
            The star forest type.

        See also
        --------
        petsc.PetscSFSetType

        """
        cdef PetscSFType cval = NULL
        sf_type = str2bytes(sf_type, &cval)
        CHKERR( PetscSFSetType(self.sf, cval) )

    def getType(self) -> str:
        """Return the type name of the communication implementation.

        Collective.

        See also
        --------
        petsc.PetscSFGetType

        """
        cdef PetscSFType cval = NULL
        CHKERR( PetscSFGetType(self.sf, &cval) )
        return bytes2str(cval)

    def setFromOptions(self) -> None:
        """Set options using the options database.

        Logically collective.

        See also
        --------
        petsc.PetscSFSetFromOptions, petsc_options

        """
        CHKERR( PetscSFSetFromOptions(self.sf) )

    def setUp(self) -> None:
        """Set up communication structures.

        Collective.

        See also
        --------
        petsc.PetscSFSetUp

        """
        CHKERR( PetscSFSetUp(self.sf) )

    def reset(self) -> None:
        """Reset a star forest so that different sizes or neighbors can be used.

        Collective.

        See also
        --------
        petsc.PetscSFReset

        """
        CHKERR( PetscSFReset(self.sf) )

    #

    def getGraph(self) -> tuple[int, ndarray, ndarray]:
        """Get graph.
        *nleaves* can be determined from the size of local.

        Not collective.

        Returns
        -------
        int
            Number of root vertices on the current process (these are possible
            targets for other process to attach leaves).
        ndarray
            Locations of leaves in leafdata buffers.
        ndarray
            Remote locations of root vertices for each leaf on the current
            process.

        See also
        --------
        petsc.PetscSFGetGraph

        """
        cdef PetscInt nroots = 0, nleaves = 0
        cdef const PetscInt *ilocal = NULL
        cdef const PetscSFNode *iremote = NULL
        CHKERR( PetscSFGetGraph(self.sf, &nroots, &nleaves, &ilocal, &iremote) )
        if ilocal == NULL:
            local = arange(0, nleaves, 1)
        else:
            local = array_i(nleaves, ilocal)
        remote = array_i(nleaves*2, <const PetscInt*>iremote)
        remote = remote.reshape(nleaves, 2)
        return toInt(nroots), local, remote

    def setGraph(self, nroots: int, local: Sequence[int], remote: Sequence[int]) -> None:
        """Set graph.

        The *nleaves* argument is determined from the size of local and/or
        remote.

        Collective.

        Parameters
        -------
        nroots
            Number of root vertices on the current process (these are possible
            targets for other process to attach leaves).
        local
            Locations of leaves in leafdata buffers, pass `None` for contiguous
            storage.
        remote
            Remote locations of root vertices for each leaf on the current
            process. Should be ``2*nleaves`` long as (rank, index) pairs.

        See also
        --------
        petsc.PetscSFSetGraph

        """

        cdef PetscInt cnroots = asInt(nroots)
        cdef PetscInt nleaves = 0
        cdef PetscInt nremote = 0
        cdef PetscInt *ilocal = NULL
        cdef PetscSFNode* iremote = NULL
        remote = iarray_i(remote, &nremote, <PetscInt**>&iremote)
        if local is not None:
            local = iarray_i(local, &nleaves, &ilocal)
            assert 2*nleaves == nremote
        else:
            assert nremote % 2 == 0
            nleaves = nremote // 2
        CHKERR( PetscSFSetGraph(self.sf, cnroots, nleaves, ilocal, PETSC_COPY_VALUES, iremote, PETSC_COPY_VALUES) )

    def setRankOrder(self, flag: bool) -> None:
        """Sort multi-points for gathers and scatters by rank order.

        Logically collective.

        Parameters
        ----------
        flag
            `True` to sort, `False` to skip sorting.

        See also
        --------
        petsc.PetscSFSetRankOrder

        """
        cdef PetscBool bval = asBool(flag)
        CHKERR( PetscSFSetRankOrder(self.sf, bval) )

    def getMulti(self) -> SF:
        """Get the inner SF implementing gathers and scatters.

        Collective.

        See also
        --------
        petsc.PetscSFGetMultiSF

        """
        cdef SF sf = SF()
        CHKERR( PetscSFGetMultiSF(self.sf, &sf.sf) )
        PetscINCREF(sf.obj)
        return sf

    def createInverse(self) -> SF:
        """Create the inverse map given a PetscSF in which all vertices have
        degree ``1``.

        Collective.

        See also
        --------
        petsc.PetscSFCreateInverseSF

        """
        cdef SF sf = SF()
        CHKERR( PetscSFCreateInverseSF(self.sf, &sf.sf) )
        return sf

    def computeDegree(self) -> ndarray:
        """Compute the degree for each root vertex.

        Collective.

        See also
        --------
        petsc.PetscSFComputeDegreeBegin, petsc.PetscSFComputeDegreeEnd

        """
        cdef const PetscInt *cdegree = NULL
        cdef PetscInt nroots
        CHKERR( PetscSFComputeDegreeBegin(self.sf, &cdegree) )
        CHKERR( PetscSFComputeDegreeEnd(self.sf, &cdegree) )
        CHKERR( PetscSFGetGraph(self.sf, &nroots, NULL, NULL, NULL) )
        degree = array_i(nroots, cdegree)
        return degree

    def createEmbeddedRootSF(self, selected: Sequence[int]) -> SF:
        """Remove edges from all but the selected roots, does not remap indices.

        Collective.

        Parameters
        ----------
        selected
            Indices of the selected roots on this process.

        See also
        --------
        petsc.PetscSFCreateEmbeddedRootSF

        """
        cdef PetscInt nroots = asInt(len(selected))
        cdef PetscInt *cselected = NULL
        selected = iarray_i(selected, &nroots, &cselected)
        cdef SF sf = SF()
        CHKERR( PetscSFCreateEmbeddedRootSF(self.sf, nroots, cselected, &sf.sf) )
        return sf

    def createEmbeddedLeafSF(self, selected: Sequence[int]) -> SF:
        """Remove edges from all but the selected leaves.

        Does not remap indices.

        Collective.

        Parameters
        ----------
        selected
            Indices of the selected roots on this process.

        See also
        --------
        petsc.PetscSFCreateEmbeddedLeafSF

        """
        cdef PetscInt nleaves = asInt(len(selected))
        cdef PetscInt *cselected = NULL
        selected = iarray_i(selected, &nleaves, &cselected)
        cdef SF sf = SF()
        CHKERR( PetscSFCreateEmbeddedLeafSF(self.sf, nleaves, cselected, &sf.sf) )
        return sf

    def createSectionSF(self, Section rootSection, remoteOffsets: Sequence[int] | None, Section leafSection) -> SF:
        """Create an expanded `SF` of dofs

        Assumes the input `SF` relates points.

        Collective.

        Parameters
        ----------
        rootSection
            Data layout of remote points for outgoing data (this is usually
            the serial section).
        remoteOffsets
            Offsets for point data on remote processes (these are offsets from
            the root section), or `None`.
        leafSection
            Data layout of local points for incoming data (this is the
            distributed section).

        See also
        --------
        petsc.PetscSFCreateSectionSF

        """
        cdef SF sectionSF = SF()
        cdef PetscInt noffsets = 0
        cdef PetscInt *cremoteOffsets = NULL
        if remoteOffsets is not None:
            remoteOffsets = iarray_i(remoteOffsets, &noffsets, &cremoteOffsets)
        CHKERR( PetscSFCreateSectionSF(self.sf, rootSection.sec, cremoteOffsets,
                                       leafSection.sec, &sectionSF.sf) )
        return sectionSF

    def distributeSection(self, Section rootSection, Section leafSection=None) -> tuple[ndarray, Section]:
        """Create a new, reorganized `Section` reorganized.

        Moves from the root to the leaves of the `SF`.

        Collective.

        Parameters:
        -----------
        rootSection
            Section defined on root space.
        leafSection
            Section defined on the leaf space.

        See also
        --------
        petsc.PetscSFDistributeSection

        """
        cdef PetscInt lpStart
        cdef PetscInt lpEnd
        cdef PetscInt *cremoteOffsets = NULL
        cdef ndarray remoteOffsets
        cdef MPI_Comm ccomm = def_Comm(self.comm, PETSC_COMM_DEFAULT)
        if leafSection is None:
            leafSection = Section()
        if leafSection.sec == NULL:
            CHKERR( PetscSectionCreate(ccomm, &leafSection.sec) )
        CHKERR( PetscSFDistributeSection(self.sf, rootSection.sec,
                                         &cremoteOffsets, leafSection.sec) )
        CHKERR( PetscSectionGetChart(leafSection.sec, &lpStart, &lpEnd) )
        remoteOffsets = array_i(lpEnd-lpStart, cremoteOffsets)
        CHKERR( PetscFree(cremoteOffsets) )
        return (remoteOffsets, leafSection)

    def compose(self, SF sf) -> SF:
        """Compose a new `SF`.

        Puts the `sf` under `self` in a top (roots) down (leaves) view.

        Collective.

        Parameters
        ----------
        sf
            `SF` to put under `self`.

        See also
        --------
        petsc.PetscSFCompose

        """
        cdef SF csf = SF()
        CHKERR( PetscSFCompose(self.sf, sf.sf, &csf.sf))
        return csf

    def bcastBegin(self, unit: Datatype, ndarray rootdata, ndarray leafdata, op: Op) -> None:
        """Begin pointwise broadcast.

        Root values are reduced to leaf values. This call has to be concluded
        with a call to `bcastEnd`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        rootdata
            Buffer to broadcast.
        leafdata
            Buffer to be reduced with values from each leaf's respective root.
        op
            MPI reduction operation.

        See also
        --------
        petsc.PetscSFBcastBegin

        """

        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFBcastBegin(self.sf, dtype, <const void*>PyArray_DATA(rootdata),
                                  <void*>PyArray_DATA(leafdata), cop) )

    def bcastEnd(self, unit: Datatype, ndarray rootdata, ndarray leafdata, op: Op) -> None:
        """End a broadcast & reduce operation started with `bcastBegin`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        rootdata
            Buffer to broadcast.
        leafdata
            Buffer to be reduced with values from each leaf's respective root.
        op
            MPI reduction operation.

        See also
        --------
        petsc.PetscSFBcastEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFBcastEnd(self.sf, dtype, <const void*>PyArray_DATA(rootdata),
                                <void*>PyArray_DATA(leafdata), cop) )

    def reduceBegin(self, unit: Datatype, ndarray leafdata, ndarray rootdata, op: Op) -> None:
        """Begin reduction of leafdata into rootdata.

        This call has to be completed with call to `reduceEnd`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        leafdata
            Values to reduce.
        rootdata
            Result of reduction of values from all leaves of each root.
        op
            MPI reduction operation.

        See also
        --------
        petsc.PetscSFReduceBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFReduceBegin(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                   <void*>PyArray_DATA(rootdata), cop) )

    def reduceEnd(self, unit: Datatype, ndarray leafdata, ndarray rootdata, op: Op) -> None:
        """End a reduction operation started with `reduceBegin`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        leafdata
            Values to reduce.
        rootdata
            Result of reduction of values from all leaves of each root.
        op
            MPI reduction operation.

        See also
        --------
        petsc.PetscSFReduceEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFReduceEnd(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                 <void*>PyArray_DATA(rootdata), cop) )

    def scatterBegin(self, unit: Datatype, ndarray multirootdata, ndarray leafdata) -> None:
        """Begin pointwise scatter operation.

        Operation is from multi-roots to leaves.
        This call has to be completed with scatterEnd.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        multirootdata
            Root buffer to send to each leaf, one unit of data per leaf.
        leafdata
            Leaf data to be updated with personal data from each respective root.

        See also
        --------
        petsc.PetscSFScatterBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFScatterBegin(self.sf, dtype, <const void*>PyArray_DATA(multirootdata),
                                    <void*>PyArray_DATA(leafdata)) )

    def scatterEnd(self, unit: Datatype, ndarray multirootdata, ndarray leafdata) -> None:
        """End scatter operation that was started with `scatterBegin`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        multirootdata
            Root buffer to send to each leaf, one unit of data per leaf.
        leafdata
            Leaf data to be updated with personal data from each respective root.

        See also
        --------
        petsc.PetscSFScatterEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFScatterEnd(self.sf, dtype, <const void*>PyArray_DATA(multirootdata),
                                  <void*>PyArray_DATA(leafdata)) )

    def gatherBegin(self, unit: Datatype, ndarray leafdata, ndarray multirootdata) -> None:
        """Begin pointwise gather of all leaves into multi-roots.

        This call has to be completed with `gatherEnd`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        leafdata
            Leaf data to gather to roots.
        multirootdata
            Root buffer to gather into, amount of space per root is
            equal to its degree.

        See also
        --------
        petsc.PetscSFGatherBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFGatherBegin(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                   <void*>PyArray_DATA(multirootdata)) )

    def gatherEnd(self, unit: Datatype, ndarray leafdata, ndarray multirootdata) -> None:
        """End gather operation that was started with `gatherBegin`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        leafdata
            Leaf data to gather to roots.
        multirootdata
            Root buffer to gather into, amount of space per root is
            equal to its degree.

        See also
        --------
        petsc.PetscSFGatherEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFGatherEnd(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                 <void*>PyArray_DATA(multirootdata)) )

    def fetchAndOpBegin(self, unit: Datatype, rootdata: ndarray, leafdata: ndarray, leafupdate: ndarray, op: Op) -> None:
        """Begin fetch and update operation.

        This operation fetches values from root and updates atomically
        by applying operation using the leaf value.

        This call has to be completed with `fetchAndOpEnd`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        rootdata
            Root values to be updated, input state is seen by first process
            to perform an update.
        leafdata
            Leaf values to use in reduction.
        leafupdate
            State at each leaf's respective root immediately prior to my atomic
            update.
        op
            MPI reduction operation.

        See also
        --------
        petsc.PetscSFFetchAndOpBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFFetchAndOpBegin(self.sf, dtype, <void*>PyArray_DATA(rootdata),
                                       <const void*>PyArray_DATA(leafdata),
                                       <void*>PyArray_DATA(leafupdate), cop) )

    def fetchAndOpEnd(self, unit: Datatype, rootdata: ndarray, leafdata: ndarray, leafupdate: ndarray, op: Op) -> None:
        """End operation started in a matching call to `fetchAndOpBegin`.

        Collective.

        Parameters
        ----------
        unit
            MPI datatype.
        rootdata
            Root values to be updated, input state is seen by first process
            to perform an update.
        leafdata
            Leaf values to use in reduction.
        leafupdate
            State at each leaf's respective root immediately prior to my atomic
            update.
        op
            MPI reduction operation.

        See also
        --------
        petsc.PetscSFFetchAndOpEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFFetchAndOpEnd(self.sf, dtype, <void*>PyArray_DATA(rootdata),
                                     <const void*>PyArray_DATA(leafdata),
                                     <void*>PyArray_DATA(leafupdate), cop) )

# --------------------------------------------------------------------

del SFType

# --------------------------------------------------------------------

# TODO: remove
cdef Datatype
cdef Op