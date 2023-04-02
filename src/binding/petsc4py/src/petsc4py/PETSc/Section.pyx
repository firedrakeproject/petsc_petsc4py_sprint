# --------------------------------------------------------------------

cdef class Section(Object):
    """TODO."""

    def __cinit__(self):
        self.obj = <PetscObject*> &self.sec
        self.sec  = NULL

    def __dealloc__(self):
        CHKERR( PetscSectionDestroy(&self.sec) )
        self.sec = NULL

    def view(self, Viewer viewer=None, test=1, test=1, test=1, test=1, test=1, test=1):
        """View a PetscSection.

        PetscSectionView(), when viewer is of type PETSCVIEWERHDF5, only saves
        distribution independent data, such as dofs, offsets, constraint dofs,
        and constraint indices. Points that have negative dofs, for instance,
        are not saved as they represent points owned by other processes. Point
        numbering and rank assignment is currently not stored. The saved
        section can be loaded with PetscSectionLoad().

        Collective.

        Parameters
        ----------
        v
                The viewer.

        See also
        --------
        petsc.PetscSectionView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscSectionView(self.sec, vwr) )

    def destroy(self):
        """Free a section object and frees its range if that exists.

        Not collective.

        See also
        --------
        petsc.PetscSectionDestroy

        """
        CHKERR( PetscSectionDestroy(&self.sec) )
        return self

    def create(self, comm=None):
        """Allocate a PetscSection and sets the map contents to the default.

        Typical calling sequence

       PetscSectionCreate(MPI_Comm,PetscSection *);!
       PetscSectionSetNumFields(PetscSection, numFields);
       PetscSectionSetChart(PetscSection,low,high);
       PetscSectionSetDof(PetscSection,point,numdof);
       PetscSectionSetUp(PetscSection);
       PetscSectionGetOffset(PetscSection,point,PetscInt *);
       PetscSectionDestroy(PetscSection);

        The PetscSection object and methods are intended to be used in the PETSc Vec and Mat implementations. The indices returned by the PetscSection are appropriate for the kind of Vec it is associated with. For example, if the vector being indexed is a local vector, we call the section a local section. If the section indexes a global vector, we call it a global section. For parallel vectors, like global vectors, we use negative indices to indicate dofs owned by other processes.

        Collective.

        Parameters
        ----------
        comm - the MPI communicator

        See also
        --------
        petsc.PetscSectionCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscSection newsec = NULL
        CHKERR( PetscSectionCreate(ccomm, &newsec) )
        PetscCLEAR(self.obj); self.sec = newsec
        return self

    def clone(self):
        """Create a shallow (if possible) copy of the PetscSection.

        Collective.

        Parameters
        ----------
        section - the PetscSection
        Output Parameter
        newSection - the copy
        Developer Note
        With standard PETSc terminology this should be called PetscSectionDuplicate()

        See also
        --------
        petsc.PetscSectionClone

        """
        cdef Section sec = <Section>type(self)()
        CHKERR(PetscSectionClone(self.sec, & sec.sec, test, test, test, test, test, test))
        return sec

    def setUp(self):
        """Calculate offsets based upon the number of degrees of freedom for each point.

        Not collective.

        See also
        --------
        petsc.PetscSectionSetUp

        """
        CHKERR( PetscSectionSetUp(self.sec) )

    def reset(self):
        """Free all section data.

        Not collective.

        See also
        --------
        petsc.PetscSectionReset

        """
        CHKERR( PetscSectionReset(self.sec) )

    def getNumFields(self):
        """Return the number of fields in a PetscSection, or 0 if no fields were defined.

        Output Parameter
        numFields - the number of fields defined, or 0 if none were defined

        Not collective.

        See also
        --------
        petsc.PetscSectionGetNumFields

        """
        cdef PetscInt numFields = 0
        CHKERR( PetscSectionGetNumFields(self.sec, &numFields) )
        return toInt(numFields)

    def setNumFields(self,numFields):
        """Set the number of fields in a PetscSection.

        Not collective.

        Parameters
        ----------
        numFields
                The number of fields

        See also
        --------
        petsc.PetscSectionSetNumFields

        """
        cdef PetscInt cnumFields = asInt(numFields)
        CHKERR( PetscSectionSetNumFields(self.sec, cnumFields) )

    def getFieldName(self,field):
        """Return the name of a field in the PetscSection.

        Not collective.

        Parameters
        ----------
        field - the field number
        Output Parameter
        fieldName - the field name
        Note
        Will error if the field number is out of range

        See also
        --------
        petsc.PetscSectionGetFieldName

        """
        cdef PetscInt cfield = asInt(field)
        cdef const char *fieldName = NULL
        CHKERR( PetscSectionGetFieldName(self.sec,cfield,&fieldName) )
        return bytes2str(fieldName)

    def setFieldName(self,field,fieldName):
        """Set the name of a field in the PetscSection.

        Will error if the field number is out of range.

        Not collective.

        Parameters
        ----------
        field - the field number
        fieldName - the field name

        See also
        --------
        petsc.PetscSectionSetFieldName

        """
        cdef PetscInt cfield = asInt(field)
        cdef const char *cname = NULL
        fieldName = str2bytes(fieldName, &cname)
        CHKERR( PetscSectionSetFieldName(self.sec,cfield,cname) )

    def getFieldComponents(self,field):
        """Return the number of field components for the given field.

        Not collective.

        Parameters
        ----------
        field - the field number
        Output Parameter
        numComp - the number of field components
        Developer Note
        This function is misnamed. There is a Num in PetscSectionGetNumFields() but not in this name

        See also
        --------
        petsc.PetscSectionGetFieldComponents

        """
        cdef PetscInt cfield = asInt(field), cnumComp = 0
        CHKERR( PetscSectionGetFieldComponents(self.sec,cfield,&cnumComp) )
        return toInt(cnumComp)

    def setFieldComponents(self,field,numComp):
        """Set the number of field components for the given field.

        Not collective.

        Parameters
        ----------
        field - the field number
        numComp - the number of field components

        See also
        --------
        petsc.PetscSectionSetFieldComponents

        """
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt cnumComp = asInt(numComp)
        CHKERR( PetscSectionSetFieldComponents(self.sec,cfield,cnumComp) )

    def getChart(self):
        """Return the range [pStart, pEnd) in which points (indices) lie for this PetscSection

        Not collective.

        Parameters
        ----------
        s - the PetscSection
        Output Parameters
        pStart - the first point
        pEnd - one past the last point

        See also
        --------
        petsc.PetscSectionGetChart

        """
        cdef PetscInt pStart = 0, pEnd = 0
        CHKERR( PetscSectionGetChart(self.sec, &pStart, &pEnd) )
        return toInt(pStart), toInt(pEnd)

    def setChart(self, pStart, pEnd):
        """Set the range [pStart, pEnd) in which points (indices) lie for this PetscSection.

        Not collective.

        Parameters
        ----------
        pStart - the first point
        pEnd - one past the last point

        See also
        --------
        petsc.PetscSectionSetChart

        """
        cdef PetscInt cStart = asInt(pStart)
        cdef PetscInt cEnd   = asInt(pEnd)
        CHKERR( PetscSectionSetChart(self.sec, cStart, cEnd) )

    def getPermutation(self):
        """Return the permutation of [0, pEnd - pStart) or NULL that was set with PetscSectionSetPermutation().

        Not collective.

        Output Parameter
        perm - The permutation as an IS

        See also
        --------
        petsc.PetscSectionGetPermutation

        """
        cdef IS perm = IS()
        CHKERR( PetscSectionGetPermutation(self.sec, &perm.iset))
        PetscINCREF(perm.obj)
        return perm

    def setPermutation(self, IS perm):
        """Set the permutation for [0, pEnd - pStart).

        Not collective.

        Parameters
        ----------
        perm - the permutation of points

        See also
        --------
        petsc.PetscSectionSetPermutation

        """
        CHKERR( PetscSectionSetPermutation(self.sec, perm.iset))

    def getDof(self,point):
        """Return the number of degrees of freedom associated with a given point.

        Not collective.

        Parameters
        ----------
        point - the point
        Output Parameter
        numDof - the number of dof
        Note
        In a global section, this size will be negative for points not owned by this process.

        See also
        --------
        petsc.PetscSectionGetDof

        """
        cdef PetscInt cpoint = asInt(point), cnumDof = 0
        CHKERR( PetscSectionGetDof(self.sec,cpoint,&cnumDof) )
        return toInt(cnumDof)

    def setDof(self,point,numDof):
        """Set the number of degrees of freedom associated with a given point.

        Not collective.

        Parameters
        ----------
        point - the point
        numDof - the number of dof

        See also
        --------
        petsc.PetscSectionSetDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionSetDof(self.sec,cpoint,cnumDof) )

    def addDof(self,point,numDof):
        """Add to the number of degrees of freedom associated with a given point.

        Not collective.

        Parameters
        ----------
        s - the PetscSection
        point - the point
        numDof - the number of additional dof

        See also
        --------
        petsc.PetscSectionAddDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionAddDof(self.sec,cpoint,cnumDof) )

    def getFieldDof(self,point,field):
        """Return the number of degrees of freedom associated with a field on a given point.

        Not collective.

        Parameters
        ----------
        point - the point
        field - the field
        Output Parameter
        numDof - the number of dof

        See also
        --------
        petsc.PetscSectionGetFieldDof

        """
        cdef PetscInt cpoint = asInt(point), cnumDof = 0
        cdef PetscInt cfield = asInt(field)
        CHKERR( PetscSectionGetFieldDof(self.sec,cpoint,cfield,&cnumDof) )
        return toInt(cnumDof)

    def setFieldDof(self,point,field,numDof):
        """Set the number of degrees of freedom associated with a field on a given point.

        Not collective.

        Parameters
        ----------
        point - the point
        field - the field
        numDof - the number of dof

        See also
        --------
        petsc.PetscSectionSetFieldDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionSetFieldDof(self.sec,cpoint,cfield,cnumDof) )

    def addFieldDof(self,point,field,numDof):
        """Add a number of degrees of freedom associated with a field on a given point.

        Not collective.

        Parameters
        ----------
        s - the PetscSection
        point - the point
        field - the field
        numDof - the number of dof

        See also
        --------
        petsc.PetscSectionAddFieldDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionAddFieldDof(self.sec,cpoint,cfield,cnumDof) )

    def getConstraintDof(self,point):
        """Return the number of constrained degrees of freedom associated with a given point.

        Not collective.

        Parameters
        ----------
        point - the point
        Output Parameter
        numDof - the number of dof which are fixed by constraints

        See also
        --------
        petsc.PetscSectionGetConstraintDof

        """
        cdef PetscInt cpoint = asInt(point), cnumDof = 0
        CHKERR( PetscSectionGetConstraintDof(self.sec,cpoint,&cnumDof) )
        return toInt(cnumDof)

    def setConstraintDof(self,point,numDof):
        """Set the number of constrained degrees of freedom associated with a given point.

        Not collective.

        Parameters
        ----------
        point - the point
        numDof - the number of dof which are fixed by constraints

        See also
        --------
        petsc.PetscSectionSetConstraintDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionSetConstraintDof(self.sec,cpoint,cnumDof) )

    def addConstraintDof(self,point,numDof):
        """Increment the number of constrained degrees of freedom associated with a given point.

        Not collective.

        Parameters
        ----------
        s - the PetscSection
        point - the point
        numDof - the number of additional dof which are fixed by constraints

        See also
        --------
        petsc.PetscSectionAddConstraintDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionAddConstraintDof(self.sec,cpoint,cnumDof) )

    def getFieldConstraintDof(self,point,field):
        """Return the number of constrained degrees of freedom associated with a given field on a point.

        Not collective.

        Parameters
        ----------
        point - the point
        field - the field
        Output Parameter
        numDof - the number of dof which are fixed by constraints

        See also
        --------
        petsc.PetscSectionGetFieldConstraintDof

        """
        cdef PetscInt cpoint = asInt(point), cnumDof = 0
        cdef PetscInt cfield = asInt(field)
        CHKERR( PetscSectionGetFieldConstraintDof(self.sec,cpoint,cfield,&cnumDof) )
        return toInt(cnumDof)

    def setFieldConstraintDof(self,point,field,numDof):
        """Set the number of constrained degrees of freedom associated with a given field on a point.

        Not collective.

        Parameters
        ----------
        point - the point
        field - the field
        numDof - the number of dof which are fixed by constraints

        See also
        --------
        petsc.PetscSectionSetFieldConstraintDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionSetFieldConstraintDof(self.sec,cpoint,cfield,cnumDof) )

    def addFieldConstraintDof(self,point,field,numDof):
        """Increment the number of constrained degrees of freedom associated with a given field on a point.

        Not collective.

        Parameters
        ----------
        s - the PetscSection
        point - the point
        field - the field
        numDof - the number of additional dof which are fixed by constraints

        See also
        --------
        petsc.PetscSectionAddFieldConstraintDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt cnumDof = asInt(numDof)
        CHKERR( PetscSectionAddFieldConstraintDof(self.sec,cpoint,cfield,cnumDof) )

    def getConstraintIndices(self,point):
        """Return the point dof numbers, in [0, dof), which are constrained for a given point.

        Not collective.

        Parameters
        ----------
        s - The PetscSection
        point - The point
        Output Parameter
        indices - The constrained dofs

        See also
        --------
        petsc.PetscSectionGetConstraintIndices
        getConstraintDof,

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt nindex = 0
        cdef const PetscInt *indices = NULL
        CHKERR( PetscSectionGetConstraintDof(self.sec, cpoint, &nindex) )
        CHKERR( PetscSectionGetConstraintIndices(self.sec, cpoint, &indices) )
        return array_i(nindex, indices)

    def setConstraintIndices(self,point,indices):
        """Set the point dof numbers, in [0, dof), which are constrained.

        Not collective.

        Parameters
        ----------
        point - The point
        indices - The constrained dofs

        See also
        --------
        petsc.PetscSectionSetConstraintIndices
        setConstraintDof,

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt nindex = 0
        cdef PetscInt *cindices = NULL
        indices = iarray_i(indices, &nindex, &cindices)
        CHKERR( PetscSectionSetConstraintDof(self.sec,cpoint,nindex) )
        CHKERR( PetscSectionSetConstraintIndices(self.sec,cpoint,cindices) )

    def getFieldConstraintIndices(self,point,field):
        """Return the field dof numbers, in [0, fdof), which are constrained.

        Not collective.

        Parameters
        ----------
        field - The field number
        point - The point
        Output Parameter
        indices - The constrained dofs sorted in ascending order
        Note
        The indices array, which is provided by the caller, must have capacity to hold the number of constrained dofs, e.g., as returned by PetscSectionGetConstraintDof().

        See also
        --------
        petsc.PetscSectionGetFieldConstraintIndices,
        getFieldConstraintDof

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt nindex = 0
        cdef const PetscInt *indices = NULL
        CHKERR( PetscSectionGetFieldConstraintDof(self.sec,cpoint,cfield,&nindex) )
        CHKERR( PetscSectionGetFieldConstraintIndices(self.sec,cpoint,cfield,&indices) )
        return array_i(nindex, indices)

    def setFieldConstraintIndices(self,point,field,indices):
        """Set the field dof numbers, in [0, fdof), which are constrained.

        Not collective.

        Parameters
        ----------
        point - The point
        field - The field number
        indices - The constrained dofs

        See also
        --------
        petsc.PetscSectionSetFieldConstraintIndices,
        setFieldConstraintDof,

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt nindex = 0
        cdef PetscInt *cindices = NULL
        indices = iarray_i(indices, &nindex, &cindices)
        CHKERR( PetscSectionSetFieldConstraintDof(self.sec,cpoint,cfield,nindex) )
        CHKERR( PetscSectionSetFieldConstraintIndices(self.sec,cpoint,cfield,cindices) )

    def getMaxDof(self):
        """Return the maximum number of degrees of freedom on any point in the PetscSection

        Not collective.

        Output Parameter
        maxDof - the maximum dof
        Note
        The returned number is up-to-date without need for PetscSectionSetUp().

        Developer Note
        The returned number is calculated lazily and stashed.

        A call to PetscSectionInvalidateMaxDof_Internal() invalidates the stashed value.

        PetscSectionInvalidateMaxDof_Internal() is called in PetscSectionSetDof(), PetscSectionAddDof() and PetscSectionReset()

        It should also be called every time atlasDof is modified directly.

        See also
        --------
        petsc.PetscSectionGetMaxDof

        """
        cdef PetscInt maxDof = 0
        CHKERR( PetscSectionGetMaxDof(self.sec,&maxDof) )
        return toInt(maxDof)

    def getStorageSize(self):
        """Return the size of an array or local Vec capable of holding all the degrees of freedom defined in a PetscSection.

        Not collective.

        Output Parameter
        size - the size of an array which can hold all the dofs

        See also
        --------
        petsc.PetscSectionGetStorageSize

        """
        cdef PetscInt size = 0
        CHKERR( PetscSectionGetStorageSize(self.sec,&size) )
        return toInt(size)

    def getConstrainedStorageSize(self):
        """Return the size of an array or local Vec capable of holding all unconstrained degrees of freedom in a PetscSection

        Output Parameter
        size - the size of an array which can hold all unconstrained dofs

        Not collective.

        See also
        --------
        petsc.PetscSectionGetConstrainedStorageSize

        """
        cdef PetscInt size = 0
        CHKERR( PetscSectionGetConstrainedStorageSize(self.sec,&size) )
        return toInt(size)

    def getOffset(self,point):
        """Return the offset into an array or Vec for the dof associated with the given point.

        Not collective.

        Parameters
        ----------
        point - the point
        Output Parameter
        offset - the offset
        Note
        In a global section, this offset will be negative for points not owned by this process.

        See also
        --------
        petsc.PetscSectionGetOffset

        """
        cdef PetscInt cpoint = asInt(point), offset = 0
        CHKERR( PetscSectionGetOffset(self.sec,cpoint,&offset) )
        return toInt(offset)

    def setOffset(self,point,offset):
        """Set the offset into an array or Vec for the dof associated with the given point.

        Not collective.

        Parameters
        ----------
        point - the point
        offset - the offset
        Note
        The user usually does not call this function, but uses PetscSectionSetUp()

        See also
        --------
        petsc.PetscSectionSetOffset

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt coffset = asInt(offset)
        CHKERR( PetscSectionSetOffset(self.sec,cpoint,coffset) )

    def getFieldOffset(self,point,field):
        """Return the offset into an array or Vec for the field dof associated with the given point.

        Not collective.

        Parameters
        ----------
        point - the point
        field - the field
        Output Parameter
        offset - the offset
        Note
        In a global section, this offset will be negative for points not owned by this process.

        See also
        --------
        petsc.PetscSectionGetFieldOffset

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt offset = 0
        CHKERR( PetscSectionGetFieldOffset(self.sec,cpoint,cfield,&offset) )
        return toInt(offset)

    def setFieldOffset(self,point,field,offset):
        """Set the offset into an array or Vec for the dof associated with the given field at a point.

        Note
        The user usually does not call this function, but uses PetscSectionSetUp()

        Not collective.

        Parameters
        ----------
        point - the point
        field - the field
        offset - the offset

        See also
        --------
        petsc.PetscSectionSetFieldOffset

        """
        cdef PetscInt cpoint = asInt(point)
        cdef PetscInt cfield = asInt(field)
        cdef PetscInt coffset = asInt(offset)
        CHKERR( PetscSectionSetFieldOffset(self.sec,cpoint,cfield,coffset) )

    def getOffsetRange(self):
        """Return the full range of offsets [start, end) for a PetscSection.

        Not collective.

        Output Parameters
        start - the minimum offset
        end - one more than the maximum offset

        See also
        --------
        petsc.PetscSectionGetOffsetRange

        """
        cdef PetscInt oStart = 0, oEnd = 0
        CHKERR( PetscSectionGetOffsetRange(self.sec,&oStart,&oEnd) )
        return toInt(oStart),toInt(oEnd)

    def createGlobalSection(self, SF sf):
        """Create a section describing the global field layout using the local section and a PetscSF describing the section point overlap.

        Parameters
        ----------
        s - The PetscSection for the local field layout
        sf - The PetscSF describing parallel layout of the section points (leaves are unowned local points)
        includeConstraints - By default this is PETSC_FALSE, meaning that the global field vector will not possess constrained dofs
        localOffsets - If PETSC_TRUE, use local rather than global offsets for the points
        Output Parameter
        gsection - The PetscSection for the global field layout
        Notes
        If we have a set of local sections defining the layout of a set of local vectors, and also a PetscSF to determine which section points are shared and the ownership, we can calculate a global section defining the parallel data layout, and the associated global vector.

        This gives negative sizes and offsets to points not owned by this process

        See also
        --------
        petsc.PetscSectionCreateGlobalSection

        """
        cdef Section gsec = Section()
        CHKERR( PetscSectionCreateGlobalSection(self.sec,sf.sf,PETSC_FALSE,PETSC_FALSE,&gsec.sec) )
        return gsec
