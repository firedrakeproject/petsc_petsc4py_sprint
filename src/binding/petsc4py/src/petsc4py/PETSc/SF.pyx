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

    Type = SFType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.sf
        self.sf  = NULL

    def __dealloc__(self):
        CHKERR( PetscSFDestroy(&self.sf) )
        self.sf = NULL

    def view(self, Viewer viewer=None) -> None:
        """View a star forest
        Collective

        Parameters
        ----------
        viewer
            Viewer to display graph, for example `PETSC_VIEWER_STDOUT_WORLD`

        See also
        --------
        PetscSFView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscSFView(self.sf, vwr) )

    def destroy(self) -> Self:
        """Destroy star forest
        Collective

        See also
        --------
        PetscSFDestroy

        """
        CHKERR( PetscSFDestroy(&self.sf) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create a star forest communication context
        Collective

        Parameters
        ----------
        comm
            The communicator on which the star forest will operate

        See also
        --------
        PetscSFCreate

        """

        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscSF newsf = NULL
        CHKERR( PetscSFCreate(ccomm, &newsf) )
        PetscCLEAR(self.obj); self.sf = newsf
        return self

    def setType(self, sf_type : str) -> None:
        """Set the communication implementation
        Collective
        https://petsc.org/release/docs/manualpages/PetscSF/PetscSFSetType/

        Parameters
        ----------
        sf_type
            a known method
            TODO: 
            PETSCSFWINDOW - MPI-2/3 one-sided
            PETSCSFBASIC - basic implementation using MPI-1 two-sided

        See also
        --------
        PetscSFSetType

        """
        cdef PetscSFType cval = NULL
        sf_type = str2bytes(sf_type, &cval)
        CHKERR( PetscSFSetType(self.sf, cval) )

    def getType(self) -> str:
        """Return the type name of the communication implementation
        Collective

        See also
        --------
        PetscSFGetType

        """
        cdef PetscSFType cval = NULL
        CHKERR( PetscSFGetType(self.sf, &cval) )
        return bytes2str(cval)

    def setFromOptions(self):
        """Set options using the options database
        Logically collective
        
        TODO:

        See also
        --------
        PetscSFSetFromOptions

        """
        CHKERR( PetscSFSetFromOptions(self.sf) )

    def setUp(self):
        """Set up communication structures
        Collective

        See also
        --------
        PetscSFSetUp

        """
        CHKERR( PetscSFSetUp(self.sf) )

    def reset(self):
        """Reset a star forest so that different sizes or neighbors can be used
        Collective

        See also
        --------
        PetscSFReset

        """
        CHKERR( PetscSFReset(self.sf) )

    #

    def getGraph(self):
        """Get graph.
        *nleaves* can be determined from the size of local
        Not collective

        Returns
        -------
        int
            Number of root vertices on the current process (these are possible targets for other process to attach leaves)
        ndarray
            Locations of leaves in leafdata buffers
        ndarray
            Remote locations of root vertices for each leaf on the current process

        See also
        --------
        PetscSFGetGraph

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

    def setGraph(self, nroots : int, local : ndarray, remote : ndarray):
        """Set graph.

        The *nleaves* argument is determined from the size of local and/or remote.
        

        Collective

        Parameters
        -------
        nroots
            Number of root vertices on the current process (these are possible targets for other process to attach leaves)
        local
            Locations of leaves in leafdata buffers, pass `None` for contiguous storage
        remote
            Remote locations of root vertices for each leaf on the current process. Should be ``2*nleaves`` long as (rank, index) pairs.

        See also
        --------
        PetscSFSetGraph

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

    def setRankOrder(self, flag):
        cdef PetscBool bval = asBool(flag)
        CHKERR( PetscSFSetRankOrder(self.sf, bval) )

    #

    def getMulti(self) -> SF:
        cdef SF sf = SF()
        CHKERR( PetscSFGetMultiSF(self.sf, &sf.sf) )
        PetscINCREF(sf.obj)
        return sf

    def createInverse(self) -> SF:
        cdef SF sf = SF()
        CHKERR( PetscSFCreateInverseSF(self.sf, &sf.sf) )
        return sf

    def computeDegree(self) -> int:
        cdef const PetscInt *cdegree = NULL
        cdef PetscInt nroots
        CHKERR( PetscSFComputeDegreeBegin(self.sf, &cdegree) )
        CHKERR( PetscSFComputeDegreeEnd(self.sf, &cdegree) )
        CHKERR( PetscSFGetGraph(self.sf, &nroots, NULL, NULL, NULL) )
        degree = array_i(nroots, cdegree)
        return degree

    def createEmbeddedRootSF(self, selected) -> SF:
        cdef PetscInt nroots = asInt(len(selected))
        cdef PetscInt *cselected = NULL
        selected = iarray_i(selected, &nroots, &cselected)
        cdef SF sf = SF()
        CHKERR( PetscSFCreateEmbeddedRootSF(self.sf, nroots, cselected, &sf.sf) )
        return sf

    def createEmbeddedLeafSF(self, selected) -> SF:
        cdef PetscInt nleaves = asInt(len(selected))
        cdef PetscInt *cselected = NULL
        selected = iarray_i(selected, &nleaves, &cselected)
        cdef SF sf = SF()
        CHKERR( PetscSFCreateEmbeddedLeafSF(self.sf, nleaves, cselected, &sf.sf) )
        return sf

    def createSectionSF(self, Section rootSection, remoteOffsets, Section leafSection) -> SF:
        cdef SF sectionSF = SF()
        cdef PetscInt noffsets = 0
        cdef PetscInt *cremoteOffsets = NULL
        if remoteOffsets is not None:
            remoteOffsets = iarray_i(remoteOffsets, &noffsets, &cremoteOffsets)
        CHKERR( PetscSFCreateSectionSF(self.sf, rootSection.sec, cremoteOffsets,
                                       leafSection.sec, &sectionSF.sf) )
        return sectionSF

    def distributeSection(self, Section rootSection, Section leafSection=None):
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
        cdef SF csf = SF()
        CHKERR( PetscSFCompose(self.sf, sf.sf, &csf.sf))
        return csf

    def bcastBegin(self, unit, ndarray rootdata, ndarray leafdata, op) -> None:
        """Begin pointwise broadcast with root value being reduced to leaf value, to be concluded with call to bcastEnd
        Collective

        Parameters
        ----------
        self
            star forest on which to communicate
        unit
            data type
        rootdata
            buffer to broadcast
        leafdata
            buffer to be reduced with values from each leaf's respective root
        op
            operation to use for reduction

        See also
        --------
        PetscSFBcastBegin

        """

        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFBcastBegin(self.sf, dtype, <const void*>PyArray_DATA(rootdata),
                                  <void*>PyArray_DATA(leafdata), cop) )

    def bcastEnd(self, unit, ndarray rootdata, ndarray leafdata, op) -> None:
        """end a broadcast & reduce operation started with bcastBegin
        Collective

        Parameters
        ----------
        self
            star forest on which to communicate
        unit
            data type
        rootdata
            buffer to broadcast
        leafdata
            buffer to be reduced with values from each leaf's respective root
        op
            operation to use for reduction
            
        See also
        --------
        PetscSFBcastEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFBcastEnd(self.sf, dtype, <const void*>PyArray_DATA(rootdata),
                                <void*>PyArray_DATA(leafdata), cop) )

    def reduceBegin(self, unit, ndarray leafdata, ndarray rootdata, op) -> None:
        """Begin reduction of leafdata into rootdata, to be completed with call to reduceEnd
        Collective

        Parameters
        ----------
        self
            star forest on which to communicate
        unit
            data type
        leafdata
            values to reduce
        rootdata
            result of reduction of values from all leaves of each root
        op
            reduction operation
            
        See also
        --------
        PetscSFReduceBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFReduceBegin(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                   <void*>PyArray_DATA(rootdata), cop) )

    def reduceEnd(self, unit, ndarray leafdata, ndarray rootdata, op) -> None:
        """End a reduction operation started with reduceBegin
        Collective

        Parameters
        ----------
        self
            star forest on which to communicate
        unit
            data type
        leafdata
            values to reduce
        rootdata
            result of reduction of values from all leaves of each root
        op
            reduction operation
            
        See also
        --------
        PetscSFReduceBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFReduceEnd(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                 <void*>PyArray_DATA(rootdata), cop) )

    def scatterBegin(self, unit, ndarray multirootdata, ndarray leafdata) -> None:
        """Begin pointwise scatter operation from multi-roots to leaves, to be completed with scatterEnd

        Parameters
        ----------
        self
            star forest
        unit
            data type
        multirootdata
            root buffer to send to each leaf, one unit of data per leaf
        leafdata
            leaf data to be update with personal data from each respective root
            
        See also
        --------
        PetscSFScatterBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFScatterBegin(self.sf, dtype, <const void*>PyArray_DATA(multirootdata),
                                    <void*>PyArray_DATA(leafdata)) )

    def scatterEnd(self, unit, ndarray multirootdata, ndarray leafdata) -> None:
        """End pointwise scatter operation that was started with scatterBegin
        Collective

        Parameters
        ----------
        self
            star forest
        unit
            data type
        multirootdata
            root buffer to send to each leaf, one unit of data per leaf
        leafdata
            leaf data to be update with personal data from each respective root
            
        See also
        --------
        PetscSFScatterEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFScatterEnd(self.sf, dtype, <const void*>PyArray_DATA(multirootdata),
                                  <void*>PyArray_DATA(leafdata)) )

    def gatherBegin(self, unit, ndarray leafdata, ndarray multirootdata) -> None:
        """Begin pointwise gather of all leaves into multi-roots, to be completed with gatherEnd
        Collective

        Parameters
        ----------
        self
            star forest
        unit
            data type
        leafdata
            leaf data to gather to roots
        multirootdata
            root buffer to gather into, amount of space per root is equal to its degree
            
        See also
        --------
        PetscSFGatherBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFGatherBegin(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                   <void*>PyArray_DATA(multirootdata)) )

    def gatherEnd(self, unit, ndarray leafdata, ndarray multirootdata) -> None:
        """End pointwise gather operation that was started with gatherBegin
        Collective

        Parameters
        ----------
        sf
            star forest
        unit
            data type
        rootdata
            root values to be updated, input state is seen by first process to perform an update
        leafdata
            leaf values to use in reduction
        leafupdate
            state at each leaf's respective root immediately prior to my atomic update
        op
            operation to use for reduction
            
        See also
        --------
        PetscSFGatherEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        CHKERR( PetscSFGatherEnd(self.sf, dtype, <const void*>PyArray_DATA(leafdata),
                                 <void*>PyArray_DATA(multirootdata)) )

    def fetchAndOpBegin(self, unit, rootdata, leafdata, leafupdate, op) -> None:
        """Begin operation that fetches values from root and updates atomically by applying operation using my leaf value, to be completed with fetchAndOpEnd
        Collective

        Parameters
        ----------
        sf
            star forest
        unit
            data type
        rootdata
            root values to be updated, input state is seen by first process to perform an update
        leafdata
            leaf values to use in reduction
        leafupdate
            state at each leaf's respective root immediately prior to my atomic update
        op
            operation to use for reduction
            
        See also
        --------
        PetscSFFetchAndOpBegin

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFFetchAndOpBegin(self.sf, dtype, <void*>PyArray_DATA(rootdata),
                                       <const void*>PyArray_DATA(leafdata),
                                       <void*>PyArray_DATA(leafupdate), cop) )

    def fetchAndOpEnd(self, unit, rootdata, leafdata, leafupdate, op) -> None:
        """End operation started in matching call to fetchAndOpBegin to fetch values from roots and update atomically by applying operation using my leaf value
        Collective

        Parameters
        ----------
        self
            star forest on which to communicate
        rootdata
            buffer to broadcast
        leafdata
            buffer to be reduced with values from each leaf's respective root
        op
            operation to use for reduction
            
        See also
        --------
        PetscSFFetchAndOpEnd

        """
        cdef MPI_Datatype dtype = mpi4py_Datatype_Get(unit)
        cdef MPI_Op cop = mpi4py_Op_Get(op)
        CHKERR( PetscSFFetchAndOpEnd(self.sf, dtype, <void*>PyArray_DATA(rootdata),
                                     <const void*>PyArray_DATA(leafdata),
                                     <void*>PyArray_DATA(leafupdate), cop) )

# --------------------------------------------------------------------

del SFType

# --------------------------------------------------------------------
