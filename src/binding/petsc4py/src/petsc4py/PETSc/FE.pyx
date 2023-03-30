# --------------------------------------------------------------------

class FEType(object):
    BASIC     = S_(PETSCFEBASIC)
    OPENCL    = S_(PETSCFEOPENCL)
    COMPOSITE = S_(PETSCFECOMPOSITE)

# --------------------------------------------------------------------

cdef class FE(Object):
    """A PETSc object that manages a finite element space."""

    Type = FEType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.fe
        self.fe = NULL

    def view(self, Viewer viewer=None) -> None:
        """View a `FE` object.

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
        """Destroy the `FE` object.

        Collective.

        See also
        --------
        petsc.PetscFEDestroy

        """
        CHKERR( PetscFEDestroy(&self.fe) )
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """Create an empty `FE` object.

        The type can then be set with `setType`.

        Collective.

        Parameters
        ----------
        comm
            The communicator for the `FE` object.

        See also
        --------
        petsc.PetscFECreate, setType

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscFE newfe = NULL
        CHKERR( PetscFECreate(ccomm, &newfe) )
        PetscCLEAR(self.obj); self.fe = newfe
        return self

    # TODO:
    def createDefault(self, dim: int, nc: int, isSimplex: bool, qorder: int = DETERMINE, prefix: str = None, comm: Comm | None = None) -> Self:
        """Create a `FE` for basic FEM computation.

        Collective.

        Parameters
        ----------
        dim
            The spatial dimension.
        nc
            The number of components.
        isSimplex
            Flag for simplex reference cell, otherwise it's a tensor product.
        qorder
            The quadrature order or `DETERMINE` to use `Space` polynomial
            degree.
        prefix
            The options prefix, or `None`.
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

    def createLagrange(self, dim: int, nc: int, isSimplex: bool, k: int, qorder: int = DETERMINE, comm: Comm | None = None) -> Self:
        """Create a `FE` for the basic Lagrange space of degree k.

        Collective.

        Parameters
        ----------
        dim
            The spatial dimension.
        nc
            The number of components.
        isSimplex
            Flag for simplex reference cell, otherwise it's a tensor product.
        k
            The degree of the space.
        qorder
            The quadrature order or `DETERMINE` to use `Space` polynomial
            degree.
        comm
            The MPI communicator.

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
        petsc.PetscFEGetQuadrature, setQuadrature

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
        petsc.PetscFEGetNumComponents, setNumComponents

        """
        cdef PetscInt comp = 0
        CHKERR( PetscFEGetNumComponents(self.fe, &comp) )
        return toInt(comp)

    def setNumComponents(self, comp: int) -> None:
        """Set the number of field components in the element.

        Not collective.

        Parameters
        ----------
        comp
            The number of field components.

        See also
        --------
        petsc.PetscFESetNumComponents, getNumComponents

        """
        cdef PetscInt ccomp = asInt(comp)
        CHKERR( PetscFESetNumComponents(self.fe, comp) )

    def getNumDof(self) -> ndarray:
        """Return the number of dofs.

        Return the number of dofs (dual basis vectors) associated with mesh
        points on the reference cell of a given dimension.

        Not collective.

        See also
        --------
        petsc.PetscFEGetNumDof

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
        petsc.PetscFEGetTileSizes, setTileSizes

        """
        cdef PetscInt blockSize = 0, numBlocks = 0
        cdef PetscInt batchSize = 0, numBatches = 0
        CHKERR( PetscFEGetTileSizes(self.fe, &blockSize, &numBlocks, &batchSize, &numBatches) )
        return toInt(blockSize), toInt(numBlocks), toInt(batchSize), toInt(numBatches)

    def setTileSizes(self, blockSize: int, numBlocks: int, batchSize: int, numBatches: int) -> None:
        """Set the tile sizes for evaluation.

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
        petsc.PetscFESetTileSizes, getTileSizes

        """
        cdef PetscInt cblockSize = asInt(blockSize), cnumBlocks = asInt(numBlocks)
        cdef PetscInt cbatchSize = asInt(batchSize), cnumBatches = asInt(numBatches)
        CHKERR( PetscFESetTileSizes(self.fe, blockSize, numBlocks, batchSize, numBatches) )

    def getFaceQuadrature(self) -> Quad:
        """Return the `Quad` used to calculate inner products on faces.

        Not collective.

        See also
        --------
        petsc.PetscFEGetFaceQuadrature, setFaceQuadrature

        """
        cdef Quad quad = Quad()
        CHKERR( PetscFEGetFaceQuadrature(self.fe, &quad.quad) )
        return quad

    def setQuadrature(self, Quad quad) -> Self:
        """Set the `Quad` used to calculate inner products.

        Not collective.

        Parameters
        ----------
        quad
            The `Quad` object.

        See also
        --------
        petsc.PetscFESetQuadrature, getQuadrature

        """
        CHKERR( PetscFESetQuadrature(self.fe, quad.quad) )
        return self

    def setFaceQuadrature(self, Quad quad) -> Quad:
        """Set the `Quad` used to calculate inner products on faces.

        Not collective.

        Parameters
        ----------
        quad
            The `Quad` object.

        See also
        --------
        petsc.PetscFESetFaceQuadrature, getFaceQuadrature

        """
        CHKERR( PetscFESetFaceQuadrature(self.fe, quad.quad) )
        return self

    def setType(self, fe_type: Type) -> Self:
        """Build a particular `FE`.

        Collective.

        Parameters
        ----------
        fe_type
            The kind of FEM space.

        See also
        --------
        petsc.PetscFESetType

        """
        cdef PetscFEType cval = NULL
        fe_type = str2bytes(fe_type, &cval)
        CHKERR( PetscFESetType(self.fe, cval) )
        return self

    def getBasisSpace(self) -> Space:
        """Return the `Space` used for the approximation of the `FE` solution.

        Not collective.

        See also
        --------
        petsc.PetscFEGetBasisSpace, setBasisSpace

        """
        cdef Space sp = Space()
        CHKERR( PetscFEGetBasisSpace(self.fe, &sp.space ) )
        return sp

    def setBasisSpace(self, Space sp) -> None:
        """Set the `Space` used for the approximation of the solution.

        Not collective.

        Parameters
        ----------
        sp
            The `Space` object.

        See also
        --------
        petsc.PetscFESetBasisSpace, getBasisSpace

        """
        CHKERR( PetscFESetBasisSpace(self.fe, sp.space ) )

    def setFromOptions(self) -> None:
        """Set parameters in a `FE` from the options database.

        Collective.

        See also
        --------
        petsc.PetscFESetFromOptions, petsc_options

        """
        CHKERR( PetscFESetFromOptions(self.fe) )

    def setUp(self) -> None:
        """Construct data structures for the `FE` after the `Type` has been set.

        Collective.

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
        petsc.PetscFEGetDualSpace, setDualSpace, DualSpace

        """
        cdef DualSpace dspace = DualSpace()
        CHKERR( PetscFEGetDualSpace(self.fe, &dspace.dualspace) )
        return dspace

    def setDualSpace(self, DualSpace dspace) -> None:
        """Set the `DualSpace` used to define the inner product.

        Not collective.

        Parameters
        ----------
        dspace
            The `DualSpace` object.

        See also
        --------
        petsc.PetscFESetDualSpace, getDualSpace, DualSpace

        """
        CHKERR( PetscFESetDualSpace(self.fe, dspace.dualspace) )

    def viewFromOptions(self, name: str, Object obj=None) -> None:
        """View from a `FE` based on values in the options database.

        Collective.

        Parameters
        ----------
        name
            Command line option name.
        obj
            Optional object that provides the options prefix.

        See also
        --------
        petsc.PetscFEViewFromOptions, petsc_options

        """
        cdef const char *cname = NULL
        _ = str2bytes(name, &cname)
        cdef PetscObject  cobj = NULL
        if obj is not None: cobj = obj.obj[0]
        CHKERR( PetscFEViewFromOptions(self.fe, cobj, cname) )

# --------------------------------------------------------------------

del FEType

# --------------------------------------------------------------------
