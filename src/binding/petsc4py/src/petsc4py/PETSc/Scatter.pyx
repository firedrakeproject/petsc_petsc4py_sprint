# --------------------------------------------------------------------

class ScatterType(object):
   """Scatter type.

   See Also
   --------
   petsc.VecScatterType

   """
   BASIC      = S_(PETSCSFBASIC)
   NEIGHBOR   = S_(PETSCSFNEIGHBOR)
   ALLGATHERV = S_(PETSCSFALLGATHERV)
   ALLGATHER  = S_(PETSCSFALLGATHER)
   GATHERV    = S_(PETSCSFGATHERV)
   GATHER     = S_(PETSCSFGATHER)
   ALLTOALL   = S_(PETSCSFALLTOALL)
   WINDOW     = S_(PETSCSFWINDOW)

# --------------------------------------------------------------------

cdef class Scatter(Object):
    """Scatter object.

    The object used to perform data movement between vectors.

    See Also
    --------
    Vec, SF, petsc.VecScatter

    """


    Type = ScatterType
    Mode = ScatterMode

    #

    def __cinit__(self):
        self.obj = <PetscObject*> &self.sct
        self.sct = NULL

    def __call__(self, x, y, addv=None, mode=None):
        self.scatter(x, y, addv, mode)

    #

    def view(self, Viewer viewer=None) -> None:
        """View the scatter.

        Collective.

        Parameters
        ----------
        viewer
          A `Viewer` instance or `None` for the default viewer.

        See Also
        --------
        petsc.VecScatterView

        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( VecScatterView(self.sct, vwr) )

    def destroy(self) -> None:
        """Destroy the scatter.

        Collective.

        See Also
        --------
        petsc.VecScatterDestroy

        """
        CHKERR( VecScatterDestroy(&self.sct) )
        return self

    def create(
        self,
        Vec vec_from,
        IS is_from,
        Vec vec_to,
        IS is_to,
    ) -> Self:
        """Create a scatter object.

        Collective.

        Parameters
        ----------
        vec_from
          A representative vector from which to scatter the data.
        is_from
          The indices of ``vec_from`` to scatter. If `None`, uses all indices.
        vec_to
          A representative vector to which scatter the data.
        is_from
          The indices of ``vec_to`` where to receive. If `None`, uses all indices.

        Examples
        --------
        The scatter object can be used to repeatedly perform data movement and it can be seen as the PETSc equivalent of NumPy-like indexing and slicing, with support for parallel communications:

        >>> v1 = PETSc.Vec().createWithArray([1, 2, 3])
        >>> v2 = PETSc.Vec().createWithArray([0, 0, 0])
        >>> sct = PETSc.Scatter().create(v1,None,v2,None)
        >>> sct.scatter(v1,v2) # v2[:] = v1[:]
        >>> sct.scatter(v2,v1,mode=PETSC.Scatter.Mode.REVERSE) # v1[:] = v2[:]

        >>> v1 = PETSc.Vec().createWithArray([1, 2, 3, 4])
        >>> v2 = PETSc.Vec().createWithArray([0, 0])
        >>> is1 = PETSc.IS().createStride(2, 3, -2)
        >>> sct = PETSc.Scatter().create(v1,is1,v2,None)
        >>> sct.scatter(v1,v2) # v2[:] = v1[3:0:-2]
        >>> sct.scatter(v2,v1,mode=PETSC.Scatter.Mode.REVERSE) # v1[3:0:-2] = v2[:]

        See Also
        --------
        IS, petsc.VecScatterCreate

        """
        cdef PetscIS cisfrom = NULL, cisto = NULL
        if is_from is not None: cisfrom = is_from.iset
        if is_to   is not None: cisto   = is_to.iset
        cdef PetscScatter newsct = NULL
        CHKERR( VecScatterCreate(
                vec_from.vec, cisfrom, vec_to.vec, cisto, &newsct) )
        PetscCLEAR(self.obj); self.sct = newsct
        return self

    def setType(self, scatter_type: Type | str) -> None:
        """Set the type of the scatter.

        Logically collective.

        See Also
        --------
        getType, petsc.VecScatterSetType

        """
        cdef PetscScatterType cval = NULL
        vec_type = str2bytes(scatter_type, &cval)
        CHKERR( VecScatterSetType(self.sct, cval) )

    def getType(self) -> str:
        """Return the type of the scatter.

        Not collective.

        See Also
        --------
        setType, petsc.VecScatterGetType

        """
        cdef PetscScatterType cval = NULL
        CHKERR( VecScatterGetType(self.sct, &cval) )
        return bytes2str(cval)

    def setFromOptions(self) -> None:
        """Configure the scatter from the options database.

        Collective.

        See Also
        --------
        petsc_options, petsc.VecScatterSetFromOptions

        """
        CHKERR( VecScatterSetFromOptions(self.sct) )

    def setUp(self) -> None:
        """Set up the internal data structures for using the scatter.

        Collective.

        See Also
        --------
        petsc.VecScatterSetUp

        """
        CHKERR( VecScatterSetUp(self.sct) )
        return self

    def copy(self) -> Scatter:
        """Return a copy of the scatter"""
        cdef Scatter scatter = Scatter()
        CHKERR( VecScatterCopy(self.sct, &scatter.sct) )
        return scatter

    @classmethod
    def toAll(cls, Vec vec) -> Scatter:
        """Create a scatter that communicates a vector to all sharing processes.

        Collective.

        Parameters
        ----------
        vec
          The vector to scatter from. The resulting scatter will have the same communicator.

        See Also
        --------
        toZero, petsc.VecScatterCreateToAll

        """
        cdef Scatter scatter = Scatter()
        cdef Vec ovec = Vec()
        CHKERR( VecScatterCreateToAll(
            vec.vec, &scatter.sct, &ovec.vec) )
        return (scatter, ovec)

    @classmethod
    def toZero(cls, Vec vec) -> Scatter:
        """Create a scatter that communicates a vector to rank zero of the sharing communicator.

        Collective.

        Parameters
        ----------
        vec
          The vector to scatter from. The resulting scatter will have the same communicator.

        See Also
        --------
        toAll, petsc.VecScatterCreateToZero

        """
        cdef Scatter scatter = Scatter()
        cdef Vec ovec = Vec()
        CHKERR( VecScatterCreateToZero(
            vec.vec, &scatter.sct, &ovec.vec) )
        return (scatter, ovec)
    #

    scatterBegin = begin
    scatterEnd = end

    #

    def begin(
        self,
        Vec vec_from,
        Vec vec_to,
        addv : InsertMode | int | None = None,
        mode: ScatterMode | int | None = None,
    ) -> None:
        """Begin a generalized scatter from one vector into another.

        Collective.

        This call has to be concluded with a call to `end`.
        For additional details on the Parameters, see `scatter`.

        See also
        --------
        create, end, petsc.VecScatterBegin

        """
        cdef PetscInsertMode  caddv = insertmode(addv)
        cdef PetscScatterMode csctm = scattermode(mode)
        CHKERR( VecScatterBegin(self.sct, vec_from.vec, vec_to.vec,
                                caddv, csctm) )

    def end(
        self,
        Vec vec_from,
        Vec vec_to,
        addv : InsertMode | int | None = None,
        mode: ScatterMode | int | None = None,
    ) -> None:
        """Finish a generalized scatter from one vector into another.

        Collective.

        This call has to be preceded by a call to `begin`.
        For additional details on the Parameters, see `scatter`.

        See also
        --------
        create, begin, petsc.VecScatterEnd

        """
        cdef PetscInsertMode  caddv = insertmode(addv)
        cdef PetscScatterMode csctm = scattermode(mode)
        CHKERR( VecScatterEnd(self.sct, vec_from.vec, vec_to.vec,
                              caddv, csctm) )

    def scatter(
        self,
        Vec vec_from,
        Vec vec_to,
        addv : InsertMode | int | None = None,
        mode: ScatterMode | int | None = None,
    ) -> None:
        """Perform a generalized scatter from one vector into another.

        Collective.

        Parameters
        ----------
        vec_from
            The source vector.
        vec_to
            The destination vector.
        addv
            Insertion mode. Possible values are:

            - `InsertMode.INSERT_VALUES` Replace existing entries with new
              values (default).

            - `InsertMode.ADD_VALUES` Add new values to existing one.
        mode
            Scatter mode. Possible values are:

            - `ScatterMode.FORWARD` If ``vec_from`` and ``vec_to`` are compatible with the vectors used to create the scatter.

            - `ScatterMode.REVERSE` If ``vec_from`` and ``vec_to`` are swapped with respect to the vectors used to create the scatter.

        See also
        --------
        create, begin, end, petsc.VecScatterBegin, petsc.VecScatterEnd

        """
        cdef PetscInsertMode  caddv = insertmode(addv)
        cdef PetscScatterMode csctm = scattermode(mode)
        CHKERR( VecScatterBegin(self.sct, vec_from.vec, vec_to.vec,
                                caddv, csctm) )
        CHKERR( VecScatterEnd(self.sct, vec_from.vec, vec_to.vec,
                              caddv, csctm) )

# --------------------------------------------------------------------

del ScatterType

# --------------------------------------------------------------------
