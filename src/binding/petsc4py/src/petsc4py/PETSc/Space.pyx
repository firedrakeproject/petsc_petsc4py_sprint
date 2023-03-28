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
        CHKERR( PetscSpaceSetUp(self.space) )

    def create(self, comm: Comm | None = None):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceGetDimension

        """
        cdef PetscInt cdim
        CHKERR( PetscSpaceGetDimension(self.space, &cdim))
        return toInt(cdim)

    def getDegree(self) -> tuple(int, int):
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceGetDegree

        """
        cdef PetscInt cdegmax, cdegmin
        CHKERR( PetscSpaceGetDegree(self.space, &cdegmin, &cdegmax))
        return toInt(cdegmin), toInt(cdegmax)

    def setDegree(self, degree: int | None, maxDegree: int | None) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceGetNumVariables

        """
        cdef PetscInt cnvars
        CHKERR( PetscSpaceGetNumVariables(self.space, &cnvars))
        return toInt(cnvars)

    def setNumVariables(self, n: int) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSetNumVariables

        """
        cdef PetscInt cn = asInt(n)
        CHKERR( PetscSpaceSetNumVariables(self.space, cn) )

    def getNumComponents(self) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceGetNumComponents

        """
        cdef PetscInt cncomps
        CHKERR( PetscSpaceGetNumComponents(self.space, &cncomps))
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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceGetType

        """
        cdef PetscSpaceType cval = NULL
        CHKERR( PetscSpaceGetType(self.space, &cval) )
        return bytes2str(cval)

    def setType(self, space_type: Space.Type | str) -> Self:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSetType

        """
        cdef PetscSpaceType cval = NULL
        space_type = str2bytes(space_type, &cval)
        CHKERR( PetscSpaceSetType(self.space, cval) )
        return self

    def getSumConcatenate(self) -> bool:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSumGetConcatenate

        """
        cdef PetscBool concatenate
        CHKERR( PetscSpaceSumGetConcatenate(self.space, &concatenate))
        return toBool(concatenate)

    def setSumConcatenate(self, concatenate: bool) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSumSetConcatenate

        """
        cdef PetscBool cconcatenate = asBool(concatenate)
        CHKERR( PetscSpaceSumSetConcatenate(self.space, concatenate))

    def getSumNumSubspaces(self) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSumGetNumSubspaces

        """
        cdef PetscInt numSumSpaces
        CHKERR( PetscSpaceSumGetNumSubspaces(self.space, &numSumSpaces))
        return toInt(numSumSpaces)

    def getSumSubspace(self, s: int) -> Space:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSumGetSubspace

        """
        cdef Space subsp = Space()
        cdef PetscInt cs = asInt(s)
        CHKERR( PetscSpaceSumGetSubspace(self.space, s, &subsp.space) )
        return subsp

    def setSumSubspace(self, s:int, Space subsp) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSumSetSubspace

        """
        cdef PetscInt cs = asInt(s)
        CHKERR( PetscSpaceSumSetSubspace(self.space, cs, subsp.space) )

    def setSumNumSubspaces(self, numSumSpaces: int) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceSumSetNumSubspaces

        """
        cdef PetscInt cnumSumSpaces = asInt(numSumSpaces)
        CHKERR( PetscSpaceSumSetNumSubspaces(self.space, cnumSumSpaces) )

    def getTensorNumSubspaces(self) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceTensorGetNumSubspaces

        """
        cdef PetscInt cnumTensSpaces = 0
        CHKERR( PetscSpaceTensorGetNumSubspaces(self.space, &cnumTensSpaces) )
        return toInt(cnumTensSpaces)

    def setTensorSubspace(self, s: int, Space subsp) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceTensorSetSubspace

        """
        cdef PetscInt cs = asInt(s)
        CHKERR( PetscSpaceTensorSetSubspace(self.space, cs, subsp.space) )

    def getTensorSubspace(self, s: int) -> Space:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceTensorGetSubspace

        """
        cdef PetscInt cs = asInt(s)
        cdef Space subsp = Space()
        CHKERR( PetscSpaceTensorGetSubspace(self.space, cs, &subsp.space) )
        return subsp

    def setTensorNumSubspaces(self, numTensSpaces: int) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpaceTensorSetNumSubspaces

        """
        cdef PetscInt cnumTensSpaces = asInt(numTensSpaces)
        CHKERR( PetscSpaceTensorSetNumSubspaces(self.space, cnumTensSpaces) )

    def getPolynomialTensor(self) -> bool:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpacePolynomialGetTensor

        """
        cdef PetscBool ctensor
        CHKERR( PetscSpacePolynomialGetTensor(self.space, &ctensor) )
        return toBool(ctensor)

    def setPolynomialTensor(self, tensor: bool) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpacePolynomialSetTensor

        """
        cdef PetscBool ctensor = asBool(tensor)
        CHKERR( PetscSpacePolynomialSetTensor(self.space, ctensor) )

    def setPointPoints(self, Quad quad) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpacePointSetPoints

        """
        CHKERR( PetscSpacePointSetPoints(self.space, quad.quad))

    def getPointPoints(self) -> Quad:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpacePointGetPoints

        """
        cdef Quad quad = Quad()
        CHKERR( PetscSpacePointGetPoints(self.space, &quad.quad))
        return quad

    def setPTrimmedFormDegree(self, formDegree: int) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscSpacePTrimmedSetFormDegree

        """
        cdef PetscInt cformDegree = asInt(formDegree)
        CHKERR( PetscSpacePTrimmedSetFormDegree(self.space, cformDegree) )

    def getPTrimmedFormDegree(self) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSetDM

        """
        CHKERR( PetscDualSpaceSetDM(self.dualspace, dm.dm) )

    def getDimension(self) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceGetDimension

        """
        cdef PetscInt cdim
        CHKERR( PetscDualSpaceGetDimension(self.dualspace, &cdim))
        return toInt(cdim)

    def getNumComponents(self) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceGetOrder

        """
        cdef PetscInt corder
        CHKERR( PetscDualSpaceGetOrder(self.dualspace, &corder))
        return toInt(corder)

    def setOrder(self, order: int) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSetOrder

        """
        cdef PetscInt corder = asInt(order)
        CHKERR( PetscDualSpaceSetOrder(self.dualspace, corder) )

    def getNumDof(self) -> ndarray:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceGetFunctional

        """
        cdef PetscInt ci = asInt(i)
        cdef Quad functional = Quad()
        CHKERR( PetscDualSpaceGetFunctional( self.dualspace, ci, &functional.quad) )
        return functional

    def getInteriorDimension(self, intdim: int) -> int:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceGetInteriorDimension

        """
        cdef PetscInt cintdim = asInt(intdim)
        CHKERR( PetscDualSpaceGetInteriorDimension(self.dualspace, &cintdim) )
        return toInt(cintdim)

    def getLagrangeContinuity(self) -> bool:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceLagrangeGetContinuity

        """
        cdef PetscBool ccontinuous = PETSC_FALSE
        CHKERR( PetscDualSpaceLagrangeGetContinuity(self.dualspace, &ccontinuous))
        return toBool(ccontinuous)

    def setLagrangeContinuity(self, continuous: bool) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceLagrangeSetContinuity

        """
        cdef PetscBool ccontinuous = asBool(continuous)
        CHKERR( PetscDualSpaceLagrangeSetContinuity(self.dualspace, ccontinuous))

    def getLagrangeTensor(self) -> bool:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceLagrangeGetTensor

        """
        cdef PetscBool ctensor = PETSC_FALSE
        CHKERR( PetscDualSpaceLagrangeGetTensor(self.dualspace, &ctensor))
        return toBool(ctensor)

    def setLagrangeTensor(self, tensor: bool) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceLagrangeSetTensor

        """
        cdef PetscBool ctensor = asBool(tensor)
        CHKERR( PetscDualSpaceLagrangeSetTensor(self.dualspace, ctensor))

    def getLagrangeTrimmed(self) -> bool:
        """TODO

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

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
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSimpleSetDimension

        """
        cdef PetscInt cdim = asInt(dim)
        CHKERR( PetscDualSpaceSimpleSetDimension(self.dualspace, cdim) )

    def setSimpleFunctional(self, func: int, Quad functional) -> None:
        """TODO

        Not collective.

        Parameters
        ----------
        TODO
            TODO.

        See also
        --------
        petsc.PetscDualSpaceSimpleSetFunctional

        """
        cdef PetscInt cfunc = asInt(func)
        CHKERR( PetscDualSpaceSimpleSetFunctional(self.dualspace, cfunc, functional.quad) )

del SpaceType
del DualSpaceType