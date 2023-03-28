# --------------------------------------------------------------------

class FEType(object):
    BASIC     = S_(PETSCFEBASIC)
    OPENCL    = S_(PETSCFEOPENCL)
    COMPOSITE = S_(PETSCFECOMPOSITE)

# --------------------------------------------------------------------

cdef class FE(Object):

    Type = FEType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.fe
        self.fe = NULL

    def view(self, Viewer viewer=None) -> None:
        """View a `FE` object

        Collective.

        Parameters
        ----------
        viewer
            A `Viewer` to display the graph.

        See also
        --------
        petsc.PetscFEView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscFEView(self.fe, vwr) )

    def destroy(self) -> Self:
        """Destroys the `FE` object

        Collective.

        See also
        --------
        petsc.PetscFEDestroy

        """
        CHKERR( PetscFEDestroy(&self.fe) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Creates an empty `FE` object.

        The type can then be set with `setType`.

        Collective.

        Parameters
        ----------
        comm
            The communicator for the `FE` object


        See also
        --------
        petsc.PetscFECreate, setType

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscFE newfe = NULL
        CHKERR( PetscFECreate(ccomm, &newfe) )
        PetscCLEAR(self.obj); self.fe = newfe
        return self

    def createDefault(self, dim: int, nc: int, isSimplex: bool, qorder: int, prefix: str = None, comm: Comm | None = None) -> Self:
        """Create a `FE` for basic FEM computation.

        Collective.

        Parameters
        ----------
        dim
            The spatial dimension.
        nc
            The number of components.
        isSimplex
            Flag for simplex reference cell, otherwise its a tensor product.
        qorder
            The quadrature order or PETSC_DETERMINE (TODO: in params)to use `Space`
            polynomial degree.
        prefix
            The options prefix, or None.
        comm
            The MPI communicator.

        See also
        --------
        petsc.PetscFECreateDefault

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscFE newfe = NULL
        cdef PetscInt cdim = asInt(dim)
        cdef PetscInt cnc = asInt(nc)
        cdef PetscInt cqorder = asInt(qorder)
        cdef PetscBool cisSimplex = asBool(isSimplex)
        cdef const char *cprefix = NULL
        if prefix:
             prefix = str2bytes(prefix, &cprefix)
        CHKERR( PetscFECreateDefault(ccomm, cdim, cnc, cisSimplex, cprefix, cqorder, &newfe))
        PetscCLEAR(self.obj); self.fe = newfe
        return self

    def createLagrange(self, dim: int, nc: int, isSimplex: bool, k: int, qorder: int, comm=None) -> Self:
        """Create a `FE` for the basic Lagrange space of degree k.

        Collective.

        Parameters
        ----------
        dim
            The spatial dimension.
        nc
            The number of components.
        isSimplex
            Flag for simplex reference cell, otherwise its a tensor product.
        k
            The degree k of the space.
        qorder
            The quadrature order or `PETSC_DETERMINE` TODO: to use `Space`
            polynomial degree.
        comm
            The MPI comm.

        See also
        --------
        petsc.PetscFECreateLagrange

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscFE newfe = NULL
        cdef PetscInt cdim = asInt(dim)
        cdef PetscInt cnc = asInt(nc)
        cdef PetscInt ck = asInt(k)
        cdef PetscInt cqorder = asInt(qorder)
        cdef PetscBool cisSimplex = asBool(isSimplex)
        CHKERR( PetscFECreateLagrange(ccomm, cdim, cnc, cisSimplex, ck, cqorder, &newfe))
        PetscCLEAR(self.obj); self.fe = newfe
        return self

    def getQuadrature(self) -> Quad:
        """Return the `Quad` used to calculate inner products.

        Not collective.

        See also
        --------
        petsc.PetscFEGetQuadrature

        """
        cdef Quad quad = Quad()
        CHKERR( PetscFEGetQuadrature(self.fe, &quad.quad) )
        return quad

    def getDimension(self) -> int:
        """Return the dimension of the finite element space on a cell.

        Not collective.

        See also
        --------
        petsc.PetscFEGetDimension

        """
        cdef PetscInt cdim = 0
        CHKERR( PetscFEGetDimension(self.fe, &cdim) )
        return toInt(cdim)

    def getSpatialDimension(self) -> int:
        """Return the spatial dimension of the element.

        Not collective.

        See also
        --------
        petsc.PetscFEGetSpatialDimension

        """
        cdef PetscInt csdim = 0
        CHKERR( PetscFEGetSpatialDimension(self.fe, &csdim) )
        return toInt(csdim)

    def getNumComponents(self) -> int:
        """Return the number of components in the element.

        Not collective.

        See also
        --------
        petsc.PetscFEGetNumComponents

        """
        cdef PetscInt comp = 0
        CHKERR( PetscFEGetNumComponents(self.fe, &comp) )
        return toInt(comp)

    def setNumComponents(self, comp: int) -> None:
        """Sets the number of field components in the element.

        Not collective.

        Parameters
        ----------
        comp
            The number of field components.

        See also
        --------
        petsc.PetscFESetNumComponents

        """
        cdef PetscInt ccomp = asInt(comp)
        CHKERR( PetscFESetNumComponents(self.fe, comp) )

    def getNumDof(self) -> ndarray:
        """Return the number of dofs.

        Return the number of dofs (dual basis vectors) associated to mesh
        points on the reference cell of a given dimension.

        Not collective.

        See also
        --------
        petsc.PetscFEGetDimension, petsc.PetscFEGetNumDof

        """
        cdef const PetscInt *numDof = NULL
        cdef PetscInt cdim = 0
        CHKERR( PetscFEGetDimension(self.fe, &cdim) )
        CHKERR( PetscFEGetNumDof(self.fe, &numDof) )
        return array_i(cdim, numDof)

    def getTileSizes(self) -> tuple(int, int, int, int):
        """Return the tile sizes for evaluation.

        Not collective.

        Returns
        -------
        blockSize : int
            The number of elements in a block.
        numBlocks : int
            The number of blocks in a batch.
        batchSize : int
            The number of elements in a batch.
        numBatches : int
            The number of batches in a chunk.

        See also
        --------
        petsc.PetscFEGetTileSizes

        """
        cdef PetscInt blockSize = 0, numBlocks = 0
        cdef PetscInt batchSize = 0, numBatches = 0
        CHKERR( PetscFEGetTileSizes(self.fe, &blockSize, &numBlocks, &batchSize, &numBatches) )
        return toInt(blockSize), toInt(numBlocks), toInt(batchSize), toInt(numBatches)

    def setTileSizes(self, blockSize: int, numBlocks: int, batchSize: int, numBatches: int) -> None:
        """Sets the tile sizes for evaluation.

        Not collective.

        Parameters
        ----------
        blockSize
            The number of elements in a block.
        numBlocks
            The number of blocks in a batch.
        batchSize
            The number of elements in a batch.
        numBatches
            The number of batches in a chunk.

        See also
        --------
        petsc.PetscFESetTileSizes

        """
        cdef PetscInt cblockSize = asInt(blockSize), cnumBlocks = asInt(numBlocks)
        cdef PetscInt cbatchSize = asInt(batchSize), cnumBatches = asInt(numBatches)
        CHKERR( PetscFESetTileSizes(self.fe, blockSize, numBlocks, batchSize, numBatches) )

    def getFaceQuadrature(self) -> Quad:
        """Return the `Quad` used to calculate inner products on faces.

        Not collective.

        See also
        --------
        petsc.PetscFEGetFaceQuadrature

        """
        cdef Quad quad = Quad()
        CHKERR( PetscFEGetFaceQuadrature(self.fe, &quad.quad) )
        return quad

    def setQuadrature(self, Quad quad) -> Self:
        """Sets the `Quad` used to calculate inner products.

        Not collective.

        Parameters
        ----------
        quad
            The `Quad` object.

        See also
        --------
        petsc.PetscFESetQuadrature

        """
        CHKERR( PetscFESetQuadrature(self.fe, quad.quad) )
        return self

    def setFaceQuadrature(self, Quad quad) -> Quad:
        """Set the `Quad` used to calculate inner products on faces.

        Not collective.

        Parameters
        ----------
        quad
            The `Quad` object

        See also
        --------
        petsc.PetscFESetFaceQuadrature

        """
        CHKERR( PetscFESetFaceQuadrature(self.fe, quad.quad) )
        return self

    def setType(self, fe_type) -> Self:
        """Build a particular `FE`.

        Collective.

        Parameters
        ----------
        fe_type
            The kind of FEM space. TODO: params

        See also
        --------
        petsc.PetscFESetType

        """
        cdef PetscFEType cval = NULL
        fe_type = str2bytes(fe_type, &cval)
        CHKERR( PetscFESetType(self.fe, cval) )
        return self

    def getBasisSpace(self) -> Space:
        """Return the `Space` used for the approximation of the solution for the `FE`.

        Not collective.

        See also
        --------
        petsc.PetscFEGetBasisSpace

        """
        cdef Space sp = Space()
        CHKERR( PetscFEGetBasisSpace(self.fe, &sp.space ) )
        return sp

    def setBasisSpace(self, Space sp) -> None:
        """Sets the `Space` used for the approximation of the solution.

        Not collective.

        Parameters
        ----------
        sp
            The `Space` object.

        See also
        --------
        petsc.PetscFESetBasisSpace

        """
        CHKERR( PetscFESetBasisSpace(self.fe, sp.space ) )

    def setFromOptions(self) -> None:
        """TODO

        Collective.

        Parameters
        ----------
        unit
            Data type.

        See also
        --------
        petsc.PetscFESetFromOptions

        """
        CHKERR( PetscFESetFromOptions(self.fe) )

    def setUp(self) -> None:
        """TODO

        Collective.

        Parameters
        ----------
        unit
            Data type.

        See also
        --------
        petsc.PetscFESetUp

        """
        CHKERR( PetscFESetUp(self.fe) )

    def getDualSpace(self) -> DualSpace:
        """Return the `DualSpace` used to define the inner product for the `FE`.

        Not collective.

        See also
        --------
        petsc.PetscFEGetDualSpace, DualSpace

        """
        cdef DualSpace dspace = DualSpace()
        CHKERR( PetscFEGetDualSpace(self.fe, &dspace.dualspace) )
        return dspace

    def setDualSpace(self, DualSpace dspace) -> None:
        """Sets the `DualSpace` used to define the inner product.

        Not collective.

        Parameters
        ----------
        dspace
            The `DualSpace` object

        See also
        --------
        petsc.PetscFESetDualSpace, DualSpace

        """
        CHKERR( PetscFESetDualSpace(self.fe, dspace.dualspace) )

    def viewFromOptions(self, name: str, Object obj=None) -> None:
        """TODO

        Collective.

        Parameters
        ----------
        unit
            Data type.

        See also
        --------
        petsc.PetscFEViewFromOptions

        """
        cdef const char *cname = NULL
        _ = str2bytes(name, &cname)
        cdef PetscObject  cobj = NULL
        if obj is not None: cobj = obj.obj[0]
        CHKERR( PetscFEViewFromOptions(self.fe, cobj, cname) )

# --------------------------------------------------------------------

del FEType

# --------------------------------------------------------------------
