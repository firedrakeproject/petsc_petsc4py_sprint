
cdef class DMLabel(Object):

    def __cinit__(self):
        self.obj = <PetscObject*> &self.dmlabel
        self.dmlabel  = NULL

    def destroy(self):
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
        CHKERR( DMLabelDestroy(&self.dmlabel) )
        return self

    def view(self, Viewer viewer=None):
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
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( DMLabelView(self.dmlabel, vwr) )

    def create(self, name, comm=None):
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
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_SELF)
        cdef PetscDMLabel newdmlabel = NULL
        cdef const char *cname = NULL
        name = str2bytes(name, &cname)
        CHKERR( DMLabelCreate(ccomm, cname, &newdmlabel) )
        PetscCLEAR(self.obj); self.dmlabel = newdmlabel
        return self

    def duplicate(self):
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
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelDuplicate(self.dmlabel, &new.dmlabel) )
        return new

    def reset(self):
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
        CHKERR( DMLabelReset(self.dmlabel) )

    def insertIS(self, IS iset, value):
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
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelInsertIS(self.dmlabel, iset.iset, cvalue)  )
        return self

    def setValue(self, point, value):
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
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelSetValue(self.dmlabel, cpoint, cvalue) )

    def getValue(self, point):
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
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = 0
        CHKERR( DMLabelGetValue(self.dmlabel, cpoint, &cvalue) )
        return toInt(cvalue)

    def getDefaultValue(self):
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
        cdef PetscInt cvalue = 0
        CHKERR( DMLabelGetDefaultValue(self.dmlabel, &cvalue) )
        return toInt(cvalue)

    def setDefaultValue(self, value):
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
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelSetDefaultValue(self.dmlabel, cvalue) )

    def clearValue(self, point, value):
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
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelClearValue(self.dmlabel, cpoint, cvalue) )

    def addStratum(self, value):
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
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelAddStratum(self.dmlabel, cvalue) )

    def addStrata(self, strata):
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
        cdef PetscInt *istrata = NULL
        cdef PetscInt numStrata = 0
        fields = iarray_i(strata, &numStrata, &istrata)
        CHKERR( DMLabelAddStrata(self.dmlabel, numStrata, istrata) )

    def addStrataIS(self, IS iset):
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
        CHKERR( DMLabelAddStrataIS(self.dmlabel, iset.iset) )

    def getNumValues(self):
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
        cdef PetscInt numValues = 0
        CHKERR( DMLabelGetNumValues(self.dmlabel, &numValues) )
        return toInt(numValues)

    def getValueIS(self):
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
        cdef IS iset = IS()
        CHKERR( DMLabelGetValueIS(self.dmlabel, &iset.iset) )
        return iset

    def stratumHasPoint(self, value, point):
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
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = asInt(value)
        cdef PetscBool ccontains = PETSC_FALSE
        CHKERR( DMLabelStratumHasPoint(self.dmlabel, cvalue, cpoint, &ccontains) )
        return toBool(ccontains)

    def hasStratum(self, value):
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
        cdef PetscInt cvalue = asInt(value)
        cdef PetscBool cexists = PETSC_FALSE
        CHKERR( DMLabelHasStratum(self.dmlabel, cvalue, &cexists) )
        return toBool(cexists)

    def getStratumSize(self, stratum):
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
        cdef PetscInt cstratum = asInt(stratum)
        cdef PetscInt csize = 0
        CHKERR( DMLabelGetStratumSize(self.dmlabel, stratum, &csize) )
        return toInt(csize)

    def getStratumIS(self, stratum):
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
        cdef PetscInt cstratum = asInt(stratum)
        cdef IS iset = IS()
        CHKERR( DMLabelGetStratumIS(self.dmlabel, cstratum, &iset.iset) )
        return iset

    def setStratumIS(self, stratum, IS iset):
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
        cdef PetscInt cstratum = asInt(stratum)
        CHKERR( DMLabelSetStratumIS(self.dmlabel, cstratum, iset.iset) )

    def clearStratum(self, stratum):
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
        cdef PetscInt cstratum = asInt(stratum)
        CHKERR( DMLabelClearStratum(self.dmlabel, cstratum) )

    def computeIndex(self):
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
        CHKERR( DMLabelComputeIndex(self.dmlabel) )

    def createIndex(self, pStart: int, pEnd: int):
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
        cdef PetscInt cpstart = asInt(pStart), cpend = asInt(pEnd)
        CHKERR( DMLabelCreateIndex(self.dmlabel, cpstart, cpend) )

    def destroyIndex(self):
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
        CHKERR( DMLabelDestroyIndex(self.dmlabel) )

    def hasValue(self, value):
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
        cdef PetscInt cvalue = asInt(value)
        cdef PetscBool cexists = PETSC_FALSE
        CHKERR( DMLabelHasValue(self.dmlabel, cvalue, &cexists) )
        return toBool(cexists)

    def hasPoint(self, point):
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
        cdef PetscInt cpoint = asInt(point)
        cdef PetscBool cexists = PETSC_FALSE
        CHKERR( DMLabelHasPoint(self.dmlabel, cpoint, &cexists) )
        return toBool(cexists)

    def getBounds(self):
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
        cdef PetscInt cpstart = 0, cpend = 0
        CHKERR( DMLabelGetBounds(self.dmlabel, &cpstart, &cpend) )
        return toInt(cpstart), toInt(cpend)

    def filter(self, start, end):
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
        cdef PetscInt cstart = 0, cend = 0
        CHKERR( DMLabelFilter(self.dmlabel, cstart, cend) )

    def permute(self, IS permutation):
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
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelPermute(self.dmlabel, permutation.iset, &new.dmlabel) )
        return new

    def distribute(self, SF sf):
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
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelDistribute(self.dmlabel, sf.sf, &new.dmlabel) )
        return new

    def gather(self, SF sf):
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
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelGather(self.dmlabel, sf.sf, &new.dmlabel) )
        return new

    def convertToSection(self):
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
        cdef Section section = Section()
        cdef IS iset = IS()
        CHKERR( DMLabelConvertToSection(self.dmlabel, &section.sec, &iset.iset) )
        return section, iset

    def getNonEmptyStratumValuesIS(self):
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
        cdef IS iset = IS()
        CHKERR( DMLabelGetNonEmptyStratumValuesIS(self.dmlabel, &iset.iset) )
        return iset
