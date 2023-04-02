# --------------------------------------------------------------------

cdef class DMComposite(DM):

    def create(self, comm=None):
        """Creates a DMCOMPOSITE, used to generate “composite” vectors made up of several subvectors.

        Collective.

        Parameters
        ----------
        comm - the processors that will share the global vector
        Output Parameter
        packer - the DMCOMPOSITE object

        See also
        --------
        petsc.DMCompositeCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDM newdm = NULL
        CHKERR( DMCompositeCreate(ccomm, &newdm) )
        PetscCLEAR(self.obj); self.dm = newdm
        return self

    def addDM(self, DM dm, *args):
        """aDD a DM vector to a DMCOMPOSITE.

        Collective.
        Add DM to composite.

        Parameters
        ----------
        dmc - the DMCOMPOSITE object
        dm - the DM object

        See also
        --------
        petsc.DMCompositeAddDM

        """
        CHKERR( DMCompositeAddDM(self.dm, dm.dm) )
        cdef object item
        for item in args:
            dm = <DM?> item
            CHKERR( DMCompositeAddDM(self.dm, dm.dm) )

    def getNumber(self):
        """Gets the number of DM objects in the DMCOMPOSITE representation.

        Get number of sub-DMs contained in the `DMComposite`.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.DMCompositeGetNumberDM

        """
        cdef PetscInt n = 0
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        return toInt(n)
    getNumberDM = getNumber

    def getEntries(self):
        """Gets the DM for each entry in a DMCOMPOSITE.

        Get tuple of sub-DMs contained in the `DMComposite`.

        Not collective.

        Parameters
        ----------
        dms - array of sufficient length (see DMCompositeGetNumberDM()) to hold the individual DM

        See also
        --------
        petsc.DMCompositeGetEntriesArray

        """
        cdef PetscInt i, n = 0
        cdef PetscDM *cdms = NULL
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        cdef object tmp = oarray_p(empty_p(n), NULL, <void**>&cdms)
        CHKERR( DMCompositeGetEntriesArray(self.dm, cdms) )
        cdef DM entry = None
        cdef list entries = []
        for i from 0 <= i < n:
            entry = subtype_DM(cdms[i])()
            entry.dm = cdms[i]
            PetscINCREF(entry.obj)
            entries.append(entry)
        return tuple(entries)

    def scatter(self, Vec gvec, lvecs):
        """Scatters from a global packed vector into its individual local vectors.

        Scatter coupled global vector into split local vectors.

        Collective.

        Parameters
        ----------
        gvec - the global vector
        lvecs - array of local vectors, NULL for any that are not needed
        Note
        This is a non-variadic alternative to DMCompositeScatter()


        See also
        --------
        petsc.DMCompositeScatterArray

        """
        cdef PetscInt i, n = 0
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        cdef PetscVec *clvecs = NULL
        cdef object tmp = oarray_p(empty_p(n), NULL, <void**>&clvecs)
        for i from 0 <= i < n:
            clvecs[i] = (<Vec?>lvecs[<Py_ssize_t>i]).vec
        CHKERR( DMCompositeScatterArray(self.dm, gvec.vec, clvecs) )

    def gather(self, Vec gvec, imode, lvecs):
        """Gathers into a global packed vector from its individual local vectors.

        Gather split local vectors into coupled global vector.

        Collective.

        Parameters
        ----------
        gvec - the global vector
        imode - INSERT_VALUES or ADD_VALUES
        lvecs - the individual sequential vectors, NULL for any that are not needed

        See also
        --------
        petsc.DMCompositeGatherArray

        """
        cdef PetscInsertMode cimode = insertmode(imode)
        cdef PetscInt i, n = 0
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        cdef PetscVec *clvecs = NULL
        cdef object tmp = oarray_p(empty_p(n), NULL, <void**>&clvecs)
        for i from 0 <= i < n:
            clvecs[i] = (<Vec?>lvecs[<Py_ssize_t>i]).vec
        CHKERR( DMCompositeGatherArray(self.dm, cimode, gvec.vec, clvecs) )

    def getGlobalISs(self):
        """Gets the index sets for each composed object in a DMCOMPOSITE.

        Collective.

        Parameters
        ----------
        Output Parameter
        is - the array of index sets
        Notes
        The is entries should be destroyed with ISDestroy(), the is array should be freed with PetscFree()

        These could be used to extract a subset of vector entries for a “multi-physics” preconditioner

        Use DMCompositeGetLocalISs() for index sets in the packed local numbering, and DMCompositeGetISLocalToGlobalMappings() for to map local sub-DM (including ghost) indices to packed global indices.

        See also
        --------
        petsc.DMCompositeGetGlobalISs

        """
        cdef PetscInt i, n = 0
        cdef PetscIS *cis = NULL
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        CHKERR( DMCompositeGetGlobalISs(self.dm, &cis) )
        cdef object isets = [ref_IS(cis[i]) for i from 0 <= i < n]
        for i from 0 <= i < n:
            CHKERR( ISDestroy(&cis[i]) )
        CHKERR( PetscFree(cis) )
        return isets

    def getLocalISs(self):
        """Gets index sets for each component of a composite local vector.

        Output Parameter
        is - array of serial index sets for each each component of the DMCOMPOSITE
        Notes
        At present, a composite local vector does not normally exist. This function is used to provide index sets for MatGetLocalSubMatrix(). In the future, the scatters for each entry in the DMCOMPOSITE may be be merged into a single scatter to a composite local vector. The user should not typically need to know which is being done.

        To get the composite global indices at all local points (including ghosts), use DMCompositeGetISLocalToGlobalMappings().

        To get index sets for pieces of the composite global vector, use DMCompositeGetGlobalISs().

        Each returned IS should be destroyed with ISDestroy(), the array should be freed with PetscFree().

        Not collective.

        See also
        --------
        petsc.DMCompositeGetLocalISs

        """
        cdef PetscInt i, n = 0
        cdef PetscIS *cis = NULL
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        CHKERR( DMCompositeGetLocalISs(self.dm, &cis) )
        cdef object isets = [ref_IS(cis[i]) for i from 0 <= i < n]
        for i from 0 <= i < n:
            CHKERR( ISDestroy(&cis[i]) )
        CHKERR( PetscFree(cis) )
        return isets

    def getLGMaps(self):
        """gets an ISLocalToGlobalMapping for each DM in the DMCOMPOSITE, maps to the composite global space.

        Collective.

        Output Parameter
        ltogs - the individual mappings for each packed vector. Note that this includes all the ghost points that individual ghosted DMDA may have.
        Note
        Each entry of ltogs should be destroyed with ISLocalToGlobalMappingDestroy(), the ltogs array should be freed with PetscFree().

        See also
        --------
        petsc.DMCompositeGetISLocalToGlobalMappings

        """
        cdef PetscInt i, n = 0
        cdef PetscLGMap *clgm = NULL
        CHKERR( DMCompositeGetNumberDM(self.dm, &n) )
        CHKERR( DMCompositeGetISLocalToGlobalMappings(self.dm, &clgm) )
        cdef object lgms = [ref_LGMap(clgm[i]) for i from 0 <= i < n]
        for i from 0 <= i < n:
            CHKERR( ISLocalToGlobalMappingDestroy(&clgm[i]) )
        CHKERR( PetscFree(clgm) )
        return lgms

    def getAccess(self, Vec gvec, locs=None):
        """TODO.
        Get access to specified parts of global vector.
        Use via `with` context manager (PEP 343).

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc._DMComposite_access

        """
        return _DMComposite_access(self, gvec, locs)
