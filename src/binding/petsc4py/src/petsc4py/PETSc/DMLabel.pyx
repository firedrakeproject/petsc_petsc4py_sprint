
cdef class DMLabel(Object):
    """A subset of the mesh from a `DM`."""
    def __cinit__(self):
        self.obj = <PetscObject*> &self.dmlabel
        self.dmlabel  = NULL

    def destroy(self) -> Self:
        """Destroy the `DMLabel`.

        Collective.

        See also
        --------
        petsc.DMLabelDestroy

        """
        CHKERR( DMLabelDestroy(&self.dmlabel) )
        return self

    def view(self, Viewer viewer=None) -> None:
        """View the label.

        Collective.

        Parameters
        ----------
        viewer
            A `Viewer` to display the graph.

        See also
        --------
        petsc.DMLabelView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( DMLabelView(self.dmlabel, vwr) )

    def create(self, name: str, comm: Comm | None = None) -> Self:
        """Create a DMLabel object, which is a multimap.

        Collective.

        Parameters
        ----------
        name
            The label name.
        comm
            The MPI communicator, usually `COMM_SELF`.

        See also
        --------
        petsc.DMLabelCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_SELF)
        cdef PetscDMLabel newdmlabel = NULL
        cdef const char *cname = NULL
        name = str2bytes(name, &cname)
        CHKERR( DMLabelCreate(ccomm, cname, &newdmlabel) )
        PetscCLEAR(self.obj); self.dmlabel = newdmlabel
        return self

    def duplicate(self) -> DMLabel:
        """Duplicate the `DMLabel`.

        Collective.

        See also
        --------
        petsc.DMLabelDuplicate

        """
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelDuplicate(self.dmlabel, &new.dmlabel) )
        return new

    def reset(self) -> None:
        """Destroy internal data structures in the `DMLabel`.

        Not collective.

        See also
        --------
        petsc.DMLabelReset

        """
        CHKERR( DMLabelReset(self.dmlabel) )

    def insertIS(self, IS iset, value: int) -> Self:
        """Set all points in the `IS` to a value.

        Not collective.

        Parameters
        ----------
        iset
            The point IS.
        value
            The point value.

        See also
        --------
        petsc.DMLabelInsertIS

        """
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelInsertIS(self.dmlabel, iset.iset, cvalue) )
        return self

    def setValue(self, point: int, value: int) -> None:
        """Set the value a label assigns to a point.

        If the value is the same as the label's default value (which is
        initially -1, and can be changed with `setDefaultValue`), this function
        will do nothing.

        Not collective.

        Parameters
        ----------
        point
            The point.
        value
            The point value.

        See also
        --------
        petsc.DMLabelSetValue, getValue, setDefaultValue

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelSetValue(self.dmlabel, cpoint, cvalue) )

    def getValue(self, point: int) -> int:
        """Return the value a label assigns to a point.

        If no value was assigned, a default value will be returned
        The default value, initially -1, can be changed with `setDefaultValue`.

        Not collective.

        Parameters
        ----------
        point
            The point.

        See also
        --------
        petsc.DMLabelGetValue, setValue, setDefaultValue

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = 0
        CHKERR( DMLabelGetValue(self.dmlabel, cpoint, &cvalue) )
        return toInt(cvalue)

    def getDefaultValue(self) -> int:
        """Return the default value returned by `getValue`.

        The default value is returned if a point has not been explicitly given
        a value. When a label is created, it is initialized to -1.

        Not collective.

        See also
        --------
        petsc.DMLabelGetDefaultValue, setDefaultValue

        """
        cdef PetscInt cvalue = 0
        CHKERR( DMLabelGetDefaultValue(self.dmlabel, &cvalue) )
        return toInt(cvalue)

    def setDefaultValue(self, value: int) -> None:
        """Set the default value returned by `getValue`.

        The value is used if a point has not been explicitly given a value.
        When a label is created, the default value is initialized to -1.

        Not collective.

        Parameters
        ----------
        value
            The default value.

        See also
        --------
        petsc.DMLabelSetDefaultValue, getDefaultValue

        """
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelSetDefaultValue(self.dmlabel, cvalue) )

    def clearValue(self, point: int, value: int) -> None:
        """Clear the value a label assigns to a point.

        Not collective.

        Parameters
        ----------
        point
            The point.
        value
            The point value.

        See also
        --------
        petsc.DMLabelClearValue

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelClearValue(self.dmlabel, cpoint, cvalue) )

    def addStratum(self, value: int) -> None:
        """Add a new stratum value in a `DMLabel`.

        Parameters
        ----------
        value
            The stratum value.

        See also
        --------
        petsc.DMLabelAddStratum, addStrata, addStrataIS

        """
        cdef PetscInt cvalue = asInt(value)
        CHKERR( DMLabelAddStratum(self.dmlabel, cvalue) )

    def addStrata(self, strata: Sequence[int]) -> None:
        """Add new stratum values in a `DMLabel`.

        Not collective.

        Parameters
        ----------
        strata
            The stratum values.

        See also
        --------
        petsc.DMLabelAddStrata, addStrataIS, addStratum

        """
        cdef PetscInt *istrata = NULL
        cdef PetscInt numStrata = 0
        fields = iarray_i(strata, &numStrata, &istrata)
        CHKERR( DMLabelAddStrata(self.dmlabel, numStrata, istrata) )

    def addStrataIS(self, IS iset) -> None:
        """Add new stratum values in a `DMLabel`.

        Not collective.

        Parameters
        ----------
        iset
            Index set with stratum values.

        See also
        --------
        petsc.DMLabelAddStrataIS, addStrata, addStratum

        """
        CHKERR( DMLabelAddStrataIS(self.dmlabel, iset.iset) )

    def getNumValues(self) -> int:
        """Return the number of values that the `DMLabel` takes.

        Not collective.

        See also
        --------
        petsc.DMLabelGetNumValues

        """
        cdef PetscInt numValues = 0
        CHKERR( DMLabelGetNumValues(self.dmlabel, &numValues) )
        return toInt(numValues)

    def getValueIS(self) -> IS:
        """Return an `IS` of all values that the `DMLabel` takes.

        Not collective.

        See also
        --------
        petsc.DMLabelGetValueIS

        """
        cdef IS iset = IS()
        CHKERR( DMLabelGetValueIS(self.dmlabel, &iset.iset) )
        return iset

    def stratumHasPoint(self, value: int, point: int) -> bool:
        """Return `True` if the stratum contains a point.

        Not collective.

        Parameters
        ----------
        value
            The stratum value.
        point
            The point.

        See also
        --------
        petsc.DMLabelStratumHasPoint

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cvalue = asInt(value)
        cdef PetscBool ccontains = PETSC_FALSE
        CHKERR( DMLabelStratumHasPoint(self.dmlabel, cvalue, cpoint, &ccontains) )
        return toBool(ccontains)

    def hasStratum(self, value: int) -> bool:
        """Determine whether points exist with the given value.

        Not collective.

        Parameters
        ----------
        value
            The stratum value.

        See also
        --------
        petsc.DMLabelHasStratum

        """
        cdef PetscInt cvalue = asInt(value)
        cdef PetscBool cexists = PETSC_FALSE
        CHKERR( DMLabelHasStratum(self.dmlabel, cvalue, &cexists) )
        return toBool(cexists)

    def getStratumSize(self, stratum: int) -> int:
        """Return the size of a stratum.

        Not collective.

        Parameters
        ----------
        stratum
            The stratum value.

        See also
        --------
        petsc.DMLabelGetStratumSize

        """
        cdef PetscInt cstratum = asInt(stratum)
        cdef PetscInt csize = 0
        CHKERR( DMLabelGetStratumSize(self.dmlabel, cstratum, &csize) )
        return toInt(csize)

    def getStratumIS(self, stratum: int) -> IS:
        """Return an IS with the stratum points.

        Not collective.

        Parameters
        ----------
        stratum
            The stratum value.

        See also
        --------
        petsc.DMLabelGetStratumIS, setStratumIS

        """
        cdef PetscInt cstratum = asInt(stratum)
        cdef IS iset = IS()
        CHKERR( DMLabelGetStratumIS(self.dmlabel, cstratum, &iset.iset) )
        return iset

    def setStratumIS(self, stratum: int, IS iset) -> None:
        """Set the stratum points using an `IS`.

        Not collective.

        Parameters
        ----------
        stratum
            The stratum value.
        iset
            The stratum points.

        See also
        --------
        petsc.DMLabelSetStratumIS, getStratumIS

        """
        cdef PetscInt cstratum = asInt(stratum)
        CHKERR( DMLabelSetStratumIS(self.dmlabel, cstratum, iset.iset) )

    def clearStratum(self, stratum: int) -> None:
        """Remove a stratum.

        Not collective.

        Parameters
        ----------
        stratum
            The stratum value.

        See also
        --------
        petsc.DMLabelClearStratum

        """
        cdef PetscInt cstratum = asInt(stratum)
        CHKERR( DMLabelClearStratum(self.dmlabel, cstratum) )

    def computeIndex(self) -> None:
        """Create an index structure for membership determination.

        Automatically determines the bounds.

        Not collective.

        See also
        --------
        petsc.DMLabelComputeIndex

        """
        CHKERR( DMLabelComputeIndex(self.dmlabel) )

    def createIndex(self, pStart: int, pEnd: int) -> None:
        """Create an index structure for membership determination.

        Not collective.

        Parameters
        ----------
        pStart
            The smallest point.
        pEnd
            The largest point + 1.

        See also
        --------
        petsc.DMLabelCreateIndex, destroyIndex

        """
        cdef PetscInt cpstart = asInt(pStart), cpend = asInt(pEnd)
        CHKERR( DMLabelCreateIndex(self.dmlabel, cpstart, cpend) )

    def destroyIndex(self) -> None:
        """Destroy the index structure.

        Not collective.

        See also
        --------
        petsc.DMLabelDestroyIndex, createIndex

        """
        CHKERR( DMLabelDestroyIndex(self.dmlabel) )

    def hasValue(self, value: int) -> bool:
        """Determine whether a label assigns the value to any point.

        Not collective.

        Parameters
        ----------
        value
            The value.

        See also
        --------
        petsc.DMLabelHasValue, hasPoint

        """
        cdef PetscInt cvalue = asInt(value)
        cdef PetscBool cexists = PETSC_FALSE
        CHKERR( DMLabelHasValue(self.dmlabel, cvalue, &cexists) )
        return toBool(cexists)

    def hasPoint(self, point: int) -> bool:
        """Determine whether a label contains a point.

        The user must call `createIndex` before this function.

        Not collective.

        Parameters
        ----------
        point
            The point.

        See also
        --------
        petsc.DMLabelHasPoint, hasValue

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscBool cexists = PETSC_FALSE
        CHKERR( DMLabelHasPoint(self.dmlabel, cpoint, &cexists) )
        return toBool(cexists)

    def getBounds(self) -> tuple[int, int]:
        """Return the smallest and largest point in the label.

        The returned values are the smallest point and the largest point + 1.

        Not collective.

        See also
        --------
        petsc.DMLabelGetBounds

        """
        cdef PetscInt cpstart = 0, cpend = 0
        CHKERR( DMLabelGetBounds(self.dmlabel, &cpstart, &cpend) )
        return toInt(cpstart), toInt(cpend)

    def filter(self, start: int, end: int) -> None:
        """Remove all points outside of [start, end).

        Not collective.

        Parameters
        ----------
        start
            The first point kept.
        end
            One more than the last point kept.

        See also
        --------
        petsc.DMLabelFilter

        """
        cdef PetscInt cstart = asInt(start), cend = asInt(end)
        CHKERR( DMLabelFilter(self.dmlabel, cstart, cend) )

    def permute(self, IS permutation) -> DMLabel:
        """Create a new label with permuted points.

        Not collective.

        Parameters
        ----------
        permutation
            The point permutation.

        See also
        --------
        petsc.DMLabelPermute

        """
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelPermute(self.dmlabel, permutation.iset, &new.dmlabel) )
        return new

    def distribute(self, SF sf) -> DMLabel:
        """Create a new label pushed forward over the `SF`.

        Collective.

        Parameters
        ----------
        sf
            The map from old to new distribution.

        See also
        --------
        petsc.DMLabelDistribute

        """
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelDistribute(self.dmlabel, sf.sf, &new.dmlabel) )
        return new

    def gather(self, SF sf) -> DMLabel:
        """Gather all label values from leafs into roots.

        This is the inverse operation to `distribute`.

        Collective.

        Parameters
        ----------
        sf
            The `SF` communication map.

        See also
        --------
        petsc.DMLabelGather

        """
        cdef DMLabel new = DMLabel()
        CHKERR( DMLabelGather(self.dmlabel, sf.sf, &new.dmlabel) )
        return new

    def convertToSection(self) -> tuple[Section, IS]:
        """Return a (`Section`, `IS`) tuple that encodes the label.

        Not collective.

        See also
        --------
        petsc.DMLabelConvertToSection

        """
        cdef Section section = Section()
        cdef IS iset = IS()
        CHKERR( DMLabelConvertToSection(self.dmlabel, &section.sec, &iset.iset) )
        return section, iset

    def getNonEmptyStratumValuesIS(self) -> IS:
        """Return an `IS` of all values that the `DMLabel` takes.

        Not collective.

        See also
        --------
        petsc.DMLabelGetNonEmptyStratumValuesIS

        """
        cdef IS iset = IS()
        CHKERR( DMLabelGetNonEmptyStratumValuesIS(self.dmlabel, &iset.iset) )
        return iset
