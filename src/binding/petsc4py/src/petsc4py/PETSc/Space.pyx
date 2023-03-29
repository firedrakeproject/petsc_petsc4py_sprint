# --------------------------------------------------------------------

class SpaceType(object):
    POLYNOMIAL = S_(PETSCSPACEPOLYNOMIAL)
    PTRIMMED   = S_(PETSCSPACEPTRIMMED)
    TENSOR     = S_(PETSCSPACETENSOR)
    SUM        = S_(PETSCSPACESUM)
    POINT      = S_(PETSCSPACEPOINT)
    SUBSPACE   = S_(PETSCSPACESUBSPACE)
    WXY        = S_(PETSCSPACEWXY)

# --------------------------------------------------------------------

cdef class Space(Object):

    Type = SpaceType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.space
        self.space  = NULL

    def setUp(self) -> None:
        """Construct data structures for the `Space`.

        Collective on `Space`.

        See also
        --------
        petsc.PetscSpaceSetUp

        """
        CHKERR( PetscSpaceSetUp(self.space) )

    def create(self, comm: Comm | None = None):
        """Creates an empty `Space` object.

        The type can then be set with `setType`.

        Collective.

        Parameters
        ----------
        comm
            The communicator for the `Space` object.

        See also
        --------
        petsc.PetscSpaceCreate

        """

        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscSpace newsp = NULL
        CHKERR( PetscSpaceCreate(ccomm, &newsp) )
        PetscCLEAR(self.obj); self.space = newsp
        return self

    def destroy(self) -> Self:
        """Destroys the `Space` object

        Collective.

        See also
        --------
        petsc.PetscSpaceDestroy

        """
        CHKERR( PetscSpaceDestroy(&self.space) )
        return self

    def view(self, Viewer viewer=None) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscSpaceView(self.space, vwr) )

    def setFromOptions(self) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSetFromOptions

        """
        CHKERR( PetscSpaceSetFromOptions(self.space) )

    def getDimension(self) -> int:
        """Return the dimension of this space, i.e. the number of basis vectors.

        See also
        --------
        petsc.PetscSpaceGetDimension

        """
        cdef PetscInt cdim
        CHKERR( PetscSpaceGetDimension(self.space, &cdim))
        return toInt(cdim)

    def getDegree(self) -> tuple(int, int):
        """Return the polynomial degrees that characterize this space.

        Returns
        -------
        minDegree : int
            The degree of the largest polynomial space contained in the space.
        maxDegree : int
            The degree of the smallest polynomial space containing the space.

        See also
        --------
        petsc.PetscSpaceGetDegree

        """
        cdef PetscInt cdegmax, cdegmin
        CHKERR( PetscSpaceGetDegree(self.space, &cdegmin, &cdegmax))
        return toInt(cdegmin), toInt(cdegmax)

    def setDegree(self, degree: int | None, maxDegree: int | None) -> None:
        """Set the degree of approximation for this space.

        Parameters
        ----------
        degree
            The degree of the largest polynomial space contained in the space.
        maxDegree
            The degree of the largest polynomial space containing the space.
            One of degree and maxDegree can be `PETSC_DETERMINE`. TODO: None?

        See also
        --------
        petsc.PetscSpaceSetDegree

        """
        assert( (degree != None) & (maxDegree != None))
        cdef PetscInt cdegree = PETSC_DETERMINE
        if degree is not None: cdegree = asInt(degree)
        cdef PetscInt cmaxdegree = PETSC_DETERMINE
        if maxDegree is not None: cmaxdegree = asInt(maxDegree)
        CHKERR( PetscSpaceSetDegree(self.space, cdegree, cmaxdegree) )

    def getNumVariables(self) -> int:
        """Return the number of variables for this space.

        See also
        --------
        petsc.PetscSpaceGetNumVariables

        """
        cdef PetscInt cnvars
        CHKERR( PetscSpaceGetNumVariables(self.space, &cnvars))
        return toInt(cnvars)

    def setNumVariables(self, n: int) -> None:
        """Set the number of variables for this space.

        Parameters
        ----------
        n
            n - The number of variables, e.g. ``x``, ``y``, ``z``...

        See also
        --------
        petsc.PetscSpaceSetNumVariables

        """
        cdef PetscInt cn = asInt(n)
        CHKERR( PetscSpaceSetNumVariables(self.space, cn) )

    def getNumComponents(self) -> int:
        """Return the number of components for this space.

        See also
        --------
        petsc.PetscSpaceGetNumComponents

        """
        cdef PetscInt cncomps
        CHKERR( PetscSpaceGetNumComponents(self.space, &cncomps))
        return toInt(cncomps)

    def setNumComponents(self, nc: int) -> None:
        """Set the number of components for this space.

        Parameters
        ----------
        nc
            The number of components.

        See also
        --------
        petsc.PetscSpaceSetNumComponents

        """
        cdef PetscInt cnc = asInt(nc)
        CHKERR( PetscSpaceSetNumComponents(self.space, cnc) )

    #def evaluate(self, points):
    #    cdef PetscInt  cnpoints = 0, cdim=0, cnfuncs=0
    #    cdef PetscReal *cpoints = NULL
    #    cdef PetscReal *B = NULL, *D = NULL, *H = NULL
    #    points = iarray_r(points, &cnpoints,  &cpoints)
    #    # Get the dimension of the space
    #    CHKERR( PetscSpaceGetDimension( self.space, &cnfuncs) )
    #    CHKERR( PetscSpace)
    #    CHKERR( PetscSpaceEvaluate(self.space, cnpoints, &cpoints, &B, &D, &H) )
    #    return array_r(cnpoints*cdim, B), array_r(cnpoints*cnc, D), array_r(, H)

    def getType(self) -> str:
        """Gets the `Space.Type` (as a string) from the object.

        Not collective.

        See also
        --------
        petsc.PetscSpaceGetType

        """
        cdef PetscSpaceType cval = NULL
        CHKERR( PetscSpaceGetType(self.space, &cval) )
        return bytes2str(cval)

    def setType(self, space_type: Space.Type | str) -> Self:
        """Build a particular `Space`.

        Collective.

        Parameters
        ----------
        space_type
            The kind of space.

        See also
        --------
        petsc.PetscSpaceSetType

        """
        cdef PetscSpaceType cval = NULL
        space_type = str2bytes(space_type, &cval)
        CHKERR( PetscSpaceSetType(self.space, cval) )
        return self

    def getSumConcatenate(self) -> bool:
        """Return the concatenate flag for this space.

        A concatenated sum space will have number of components equal to the
        sum of the number of components of all subspaces. A non-concatenated,
        or direct sum space will have the same number of components as its
        subspaces.

        See also
        --------
        petsc.PetscSpaceSumGetConcatenate

        """
        cdef PetscBool concatenate
        CHKERR( PetscSpaceSumGetConcatenate(self.space, &concatenate))
        return toBool(concatenate)

    def setSumConcatenate(self, concatenate: bool) -> None:
        """Set the concatenate flag for this space

        A concatenated sum space will have number of components equal to the
        sum of the number of components of all subspaces. A non-concatenated,
        or direct sum space will have the same number of components as its
        subspaces.

        Parameters
        ----------
        concatenate
            Are subspaces concatenated components (`True`) or direct summands (`False`)?

        See also
        --------
        petsc.PetscSpaceSumSetConcatenate

        """
        cdef PetscBool cconcatenate = asBool(concatenate)
        CHKERR( PetscSpaceSumSetConcatenate(self.space, concatenate))

    def getSumNumSubspaces(self) -> int:
        """Get the number of spaces in the sum.

        See also
        --------
        petsc.PetscSpaceSumGetNumSubspaces

        """
        cdef PetscInt numSumSpaces
        CHKERR( PetscSpaceSumGetNumSubspaces(self.space, &numSumSpaces))
        return toInt(numSumSpaces)

    def getSumSubspace(self, s: int) -> Space:
        """Return a space in the sum.

        Parameters
        ----------
        s
            The space number.

        See also
        --------
        petsc.PetscSpaceSumGetSubspace

        """
        cdef Space subsp = Space()
        cdef PetscInt cs = asInt(s)
        CHKERR( PetscSpaceSumGetSubspace(self.space, s, &subsp.space) )
        return subsp

    def setSumSubspace(self, s: int, Space subsp) -> None:
        """Set a space in the sum.

        Parameters
        ----------
        s
            The space number.
        subsp
            The number of spaces.

        See also
        --------
        petsc.PetscSpaceSumSetSubspace

        """
        cdef PetscInt cs = asInt(s)
        CHKERR( PetscSpaceSumSetSubspace(self.space, cs, subsp.space) )

    def setSumNumSubspaces(self, numSumSpaces: int) -> None:
        """Set the number of spaces in the sum.

        Parameters
        ----------
        numSumSpaces
            The number of spaces.

        See also
        --------
        petsc.PetscSpaceSumSetNumSubspaces

        """
        cdef PetscInt cnumSumSpaces = asInt(numSumSpaces)
        CHKERR( PetscSpaceSumSetNumSubspaces(self.space, cnumSumSpaces) )

    def getTensorNumSubspaces(self) -> int:
        """Return the number of spaces in the tensor product.

        See also
        --------
        petsc.PetscSpaceTensorGetNumSubspaces

        """
        cdef PetscInt cnumTensSpaces = 0
        CHKERR( PetscSpaceTensorGetNumSubspaces(self.space, &cnumTensSpaces) )
        return toInt(cnumTensSpaces)

    def setTensorSubspace(self, s: int, Space subsp) -> None:
        """Set a space in the tensor product.

        Parameters
        ----------
        s
            The space number.
        subsp
            The number of spaces.

        See also
        --------
        petsc.PetscSpaceTensorSetSubspace

        """
        cdef PetscInt cs = asInt(s)
        CHKERR( PetscSpaceTensorSetSubspace(self.space, cs, subsp.space) )

    def getTensorSubspace(self, s: int) -> Space:
        """Get a space in the tensor product.

        Parameters
        ----------
        s
            The space number.

        See also
        --------
        petsc.PetscSpaceTensorGetSubspace

        """
        cdef PetscInt cs = asInt(s)
        cdef Space subsp = Space()
        CHKERR( PetscSpaceTensorGetSubspace(self.space, cs, &subsp.space) )
        return subsp

    def setTensorNumSubspaces(self, numTensSpaces: int) -> None:
        """Set the number of spaces in the tensor product.

        Parameters
        ----------
        numTensSpaces
            The number of spaces.

        See also
        --------
        petsc.PetscSpaceTensorSetNumSubspaces

        """
        cdef PetscInt cnumTensSpaces = asInt(numTensSpaces)
        CHKERR( PetscSpaceTensorSetNumSubspaces(self.space, cnumTensSpaces) )

    def getPolynomialTensor(self) -> bool:
        """Return whether a function space is a space of tensor polynomials.

        Return `True` if a function space is a space of tensor polynomials
        (the space is spanned by polynomials whose degree in each variable is
        bounded by the given order), as opposed to polynomials (the space is
        spanned by polynomials whose total degree—summing over all variables—is
        bounded by the given order).

        See also
        --------
        petsc.PetscSpacePolynomialGetTensor

        """
        cdef PetscBool ctensor
        CHKERR( PetscSpacePolynomialGetTensor(self.space, &ctensor) )
        return toBool(ctensor)

    def setPolynomialTensor(self, tensor: bool) -> None:
        """Set whether a function space is a space of tensor polynomials

        Set to `True` for a function space which is a space of tensor polynomials
        (the space is spanned by polynomials whose degree in each variable is
        bounded by the given order), as opposed to polynomials (the space is
        spanned by polynomials whose total degree—summing over all variables—is
        bounded by the given order).

        Parameters
        ----------
        tensor
            `True` for a tensor polynomial space, `False` for a polynomial space.

        See also
        --------
        petsc.PetscSpacePolynomialSetTensor

        """
        cdef PetscBool ctensor = asBool(tensor)
        CHKERR( PetscSpacePolynomialSetTensor(self.space, ctensor) )

    def setPointPoints(self, Quad quad) -> None:
        """Sets the evaluation points for the space to based on a quad.

        Sets the evaluation points for the space to coincide with the points of a quadrature rule.

        Logically collective.

        Parameters
        ----------
        quad
            The `Quad` defining the points.

        See also
        --------
        petsc.PetscSpacePointSetPoints

        """
        CHKERR( PetscSpacePointSetPoints(self.space, quad.quad))

    def getPointPoints(self) -> Quad:
        """Get the evaluation points for the space as the points of a quad rule.

        Logically collective.

        See also
        --------
        petsc.PetscSpacePointGetPoints

        """
        cdef Quad quad = Quad()
        CHKERR( PetscSpacePointGetPoints(self.space, &quad.quad))
        return quad

    def setPTrimmedFormDegree(self, formDegree: int) -> None:
        """Set the form degree of the trimmed polynomials.

        Parameters
        ----------
        formDegree
            The form degree.

        See also
        --------
        petsc.PetscSpacePTrimmedSetFormDegree

        """
        cdef PetscInt cformDegree = asInt(formDegree)
        CHKERR( PetscSpacePTrimmedSetFormDegree(self.space, cformDegree) )

    def getPTrimmedFormDegree(self) -> int:
        """Return the form degree of the trimmed polynomials.

        See also
        --------
        petsc.PetscSpacePTrimmedGetFormDegree

        """
        cdef PetscInt cformDegree = 0
        CHKERR( PetscSpacePTrimmedGetFormDegree(self.space, &cformDegree) )
        return toInt(cformDegree)

    def viewFromOptions(self, name: str, Object obj=None) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceViewFromOptions

        """
        cdef const char *cname = NULL
        _ = str2bytes(name, &cname)
        cdef PetscObject  cobj = NULL
        if obj is not None: cobj = obj.obj[0]
        CHKERR( PetscSpaceViewFromOptions(self.space, cobj, cname) )

# --------------------------------------------------------------------

class DualSpaceType(object):
    LAGRANGE = S_(PETSCDUALSPACELAGRANGE)
    SIMPLE   = S_(PETSCDUALSPACESIMPLE)
    REFINED  = S_(PETSCDUALSPACEREFINED)
    BDM      = S_(PETSCDUALSPACEBDM)

# --------------------------------------------------------------------

cdef class DualSpace(Object):

    Type = DualSpaceType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.dualspace
        self.dualspace  = NULL

    def setUp(self) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSetUp

        """
        CHKERR( PetscDualSpaceSetUp(self.dualspace) )

    def create(self, comm: Comm | None = None) -> Self:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscDualSpace newdsp = NULL
        CHKERR( PetscDualSpaceCreate(ccomm, &newdsp) )
        PetscCLEAR(self.obj); self.dualspace = newdsp
        return self

    def view(self, Viewer viewer=None) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscDualSpaceView(self.dualspace, vwr) )

    def destroy(self) -> Self:
        """Destroy the `DualSpace` object.

        Collective.

        See also
        --------
        petsc.PetscDualSpaceDestroy

        """
        CHKERR( PetscDualSpaceDestroy(&self.dualspace) )
        return self

    def duplicate(self) -> DualSpace:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceDuplicate

        """
        cdef DualSpace spNew = DualSpace()
        CHKERR( PetscDualSpaceDuplicate(self.dualspace, &spNew.dualspace) )

    def getDM(self) -> DM:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceGetDM

        """
        cdef DM dm = DM()
        CHKERR( PetscDualSpaceGetDM(self.dualspace, &dm.dm) )
        return dm

    def setDM(self, DM dm) -> None:
        """Get the `DM` representing the reference cell.

        Not collective.

        Parameters
        ----------
        dm
            The reference cell.

        See also
        --------
        petsc.PetscDualSpaceSetDM

        """
        CHKERR( PetscDualSpaceSetDM(self.dualspace, dm.dm) )

    def getDimension(self) -> int:
        """Return the dimension of the dual space.

        The dimension of the dual space, i.e. the number of basis functionals.

        Not collective.

        See also
        --------
        petsc.PetscDualSpaceGetDimension

        """
        cdef PetscInt cdim
        CHKERR( PetscDualSpaceGetDimension(self.dualspace, &cdim))
        return toInt(cdim)

    def getNumComponents(self) -> int:
        """Return the number of components for this space.

        See also
        --------
        petsc.PetscDualSpaceGetNumComponents

        """
        cdef PetscInt cncomps
        CHKERR( PetscDualSpaceGetNumComponents(self.dualspace, &cncomps))
        return toInt(cncomps)

    def setNumComponents(self, nc: int) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSetNumComponents

        """
        cdef PetscInt cnc = asInt(nc)
        CHKERR( PetscDualSpaceSetNumComponents(self.dualspace, cnc) )

    def getType(self) -> str:
        """Gets the `DualSpace.Type` name (as a string) from the object.

        Not collective.

        See also
        --------
        petsc.PetscDualSpaceGetType

        """
        cdef PetscDualSpaceType cval = NULL
        CHKERR( PetscDualSpaceGetType(self.dualspace, &cval) )
        return bytes2str(cval)

    def setType(self, dualspace_type: DualSpace.Type | str) -> Self:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSetType

        """
        cdef PetscDualSpaceType cval = NULL
        space_type = str2bytes(dualspace_type, &cval)
        CHKERR( PetscDualSpaceSetType(self.dualspace, cval) )
        return self

    def getOrder(self) -> int:
        """Get the order of the dual space.

        Not collective.

        See also
        --------
        petsc.PetscDualSpaceGetOrder

        """
        cdef PetscInt corder
        CHKERR( PetscDualSpaceGetOrder(self.dualspace, &corder))
        return toInt(corder)

    def setOrder(self, order: int) -> None:
        """Set the order of the dual space.

        Not collective.

        Parameters
        ----------
        order
            The order.

        See also
        --------
        petsc.PetscDualSpaceSetOrder

        """
        cdef PetscInt corder = asInt(order)
        CHKERR( PetscDualSpaceSetOrder(self.dualspace, corder) )

    def getNumDof(self) -> ndarray:
        """Get the number of degrees of freedom for each spatial dimension.

        Not collective.

        TODO: numDof - An array of length dim+1 which holds the number of dofs for each dimension
        why +1?

        See also
        --------
        petsc.PetscDualSpaceGetNumDof

        """
        cdef const PetscInt *cndof = NULL
        cdef PetscInt cdim = 0
        CHKERR( PetscDualSpaceGetDimension(self.dualspace, &cdim) )
        CHKERR( PetscDualSpaceGetNumDof(self.dualspace, &cndof) )
        return array_i(cdim + 1, cndof)

    def getFunctional(self, i: int) -> Quad:
        """Get the i-th basis functional in the dual space.

        Not collective.

        Parameters
        ----------
        i
            The basis number.

        See also
        --------
        petsc.PetscDualSpaceGetFunctional

        """
        cdef PetscInt ci = asInt(i)
        cdef Quad functional = Quad()
        CHKERR( PetscDualSpaceGetFunctional( self.dualspace, ci, &functional.quad) )
        return functional

    def getInteriorDimension(self) -> int:
        """Get the interior dimension of the dual space.

        The interior dimension of the dual space, i.e. the number of basis
        functionals assigned to the interior of the reference domain.

        Not collective.

        See also
        --------
        petsc.PetscDualSpaceGetInteriorDimension

        """
        cdef PetscInt cintdim = 0
        CHKERR( PetscDualSpaceGetInteriorDimension(self.dualspace, &cintdim) )
        return toInt(cintdim)

    def getLagrangeContinuity(self) -> bool:
        """Return the flag for element continuity.

        Not collective.

        See also
        --------
        petsc.PetscDualSpaceLagrangeGetContinuity

        """
        cdef PetscBool ccontinuous = PETSC_FALSE
        CHKERR( PetscDualSpaceLagrangeGetContinuity(self.dualspace, &ccontinuous))
        return toBool(ccontinuous)

    def setLagrangeContinuity(self, continuous: bool) -> None:
        """Indicate whether the element is continuous.

        Not collective.

        Parameters
        ----------
        continuous
            The flag for element continuity.

        See also
        --------
        petsc.PetscDualSpaceLagrangeSetContinuity

        """
        cdef PetscBool ccontinuous = asBool(continuous)
        CHKERR( PetscDualSpaceLagrangeSetContinuity(self.dualspace, ccontinuous))

    def getLagrangeTensor(self) -> bool:
        """Get the tensor nature of the dual space.

        Not collective.

        See also
        --------
        petsc.PetscDualSpaceLagrangeGetTensor

        """
        cdef PetscBool ctensor = PETSC_FALSE
        CHKERR( PetscDualSpaceLagrangeGetTensor(self.dualspace, &ctensor))
        return toBool(ctensor)

    def setLagrangeTensor(self, tensor: bool) -> None:
        """Set the tensor nature of the dual space.

        Not collective.

        Parameters
        ----------
        tensor
            Whether the dual space has tensor layout (vs. simplicial).

        See also
        --------
        petsc.PetscDualSpaceLagrangeSetTensor

        """
        cdef PetscBool ctensor = asBool(tensor)
        CHKERR( PetscDualSpaceLagrangeSetTensor(self.dualspace, ctensor))

    def getLagrangeTrimmed(self) -> bool:
        """Return the trimmed nature of the dual space.

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceLagrangeGetTrimmed

        """
        cdef PetscBool ctrimmed = PETSC_FALSE
        CHKERR( PetscDualSpaceLagrangeGetTrimmed(self.dualspace, &ctrimmed))
        return toBool(ctrimmed)

    def setLagrangeTrimmed(self, trimmed: bool) -> None:
        """Set the trimmed nature of the dual space.

        Not collective.

        Parameters
        ----------
        trimmed
            Whether the dual space represents to dual basis of a trimmed
            polynomial space (e.g. Raviart-Thomas and higher order /
            other form degree variants).

        See also
        --------
        petsc.PetscDualSpaceLagrangeSetTrimmed

        """
        cdef PetscBool ctrimmed = asBool(trimmed)
        CHKERR( PetscDualSpaceLagrangeSetTrimmed(self.dualspace, ctrimmed))

    def viewFromOptions(self, name: str, Object obj=None) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceViewFromOptions

        """
        cdef const char *cname = NULL
        _ = str2bytes(name, &cname)
        cdef PetscObject  cobj = NULL
        if obj is not None: cobj = obj.obj[0]
        CHKERR( PetscDualSpaceViewFromOptions(self.dualspace, cobj, cname) )

    def setSimpleDimension(self, dim: int) -> None:
        """Set the number of functionals in the dual space basis.

        Logically collective.

        Parameters
        ----------
        dim
            The basis dimension.

        See also
        --------
        petsc.PetscDualSpaceSimpleSetDimension

        """
        cdef PetscInt cdim = asInt(dim)
        CHKERR( PetscDualSpaceSimpleSetDimension(self.dualspace, cdim) )

    def setSimpleFunctional(self, func: int, Quad functional) -> None:
        """Set the given basis element for this dual space.

        Not collective.

        Parameters
        ----------
        func
            The basis index.
        functional
            The basis functional.

        See also
        --------
        petsc.PetscDualSpaceSimpleSetFunctional

        """
        cdef PetscInt cfunc = asInt(func)
        CHKERR( PetscDualSpaceSimpleSetFunctional(self.dualspace, cfunc, functional.quad) )

del SpaceType
del DualSpaceType