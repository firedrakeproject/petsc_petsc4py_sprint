# --------------------------------------------------------------------

cdef class Quad(Object):
    """
    Quadrature rule for integration.
    """
    def __cinit__(self):
        self.obj = <PetscObject*> &self.quad
        self.quad = NULL

    def view(self, Viewer viewer=None) -> None:
        """View a `Quad` object

        Collective.

        Parameters
        ----------
        viewer
            A `Viewer` to display the graph.

        See also
        --------
        petsc.PetscQuadratureView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PetscQuadratureView(self.quad, vwr) )

    def create(self, comm: Comm | None = None) -> Self:
        """Create a `Quad` object

        Collective.

        Parameters
        ----------
        comm
            The communicator for the `Quad` object.

        See also
        --------
        petsc.PetscQuadratureCreate

        """
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
        cdef PetscQuadrature newquad = NULL
        CHKERR( PetscQuadratureCreate(ccomm, &newquad) )
        PetscCLEAR(self.obj); self.quad = newquad
        return self

    def duplicate(self) -> Quad:
        """Create a deep copy of the `Quad` object.

        Collective.

        See also
        --------
        petsc.PetscQuadratureDuplicate

        """
        cdef Quad newquad = Quad()
        CHKERR( PetscQuadratureDuplicate(self.quad, &newquad.quad) )
        return newquad

    def destroy(self) -> Self:
        """Destroys the `Quad` object

        Collective.

        See also
        --------
        petsc.PetscQuadratureDestroy

        """
        CHKERR( PetscQuadratureDestroy(&self.quad) )
        return self

    def getData(self) -> tuple(ndarray, ndarray):
        """Return the data defining the `Quad`.

        Not collective.

        Returns
        -------
        points : ndarray
            The quadrature points.
        weights : ndarray
            The quadrature weights.

        See also
        --------
        petsc.PetscQuadratureGetData

        """
        cdef PetscInt cdim = 0
        cdef PetscInt cnc = 0
        cdef PetscInt cnpoints = 0
        cdef const PetscReal *cpoints = NULL
        cdef const PetscReal *cweights = NULL
        CHKERR( PetscQuadratureGetData(self.quad, &cdim, &cnc, &cnpoints, &cpoints, &cweights))
        return array_r(cnpoints*cdim, cpoints), array_r(cnpoints*cnc, cweights)

    def getNumComponents(self) -> int:
        """Return the number of components for functions to be integrated

        Not collective.

        See also
        --------
        petsc.PetscQuadratureGetNumComponents

        """
        cdef PetscInt cnc = 0
        CHKERR( PetscQuadratureGetNumComponents(self.quad, &cnc) )
        return toInt(cnc)

    def setNumComponents(self, nc: int) -> None:
        """Return the number of components for functions to be integrated.

        Not collective.

        Parameters
        ----------
        nc
            The number of components.

        See also
        --------
        petsc.PetscQuadratureSetNumComponents

        """
        cdef PetscInt cnc = asInt(nc)
        CHKERR( PetscQuadratureSetNumComponents(self.quad, cnc) )

    def getOrder(self) -> int:
        """Return the order of the method in the `Quad`.

        Not collective.

        See also
        --------
        petsc.PetscQuadratureGetOrder

        """
        cdef PetscInt corder = 0
        CHKERR( PetscQuadratureGetOrder(self.quad, &corder))
        return toInt(corder)

    def setOrder(self, order: int) -> None:
        """Set the order of the method in the `Quad`.

        Not collective.

        Parameters
        ----------
        order
            The order of the quadrature, i.e. the highest degree polynomial
            that is exactly integrated.

        See also
        --------
        petsc.PetscQuadratureSetOrder

        """
        cdef PetscInt corder = asInt(order)
        CHKERR( PetscQuadratureSetOrder(self.quad, corder))


# --------------------------------------------------------------------
