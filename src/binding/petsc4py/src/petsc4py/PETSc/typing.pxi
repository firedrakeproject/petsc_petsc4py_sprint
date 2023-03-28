cdef Any
cdef Self
cdef Union
cdef Literal
cdef Optional
cdef NoReturn

cdef Callable
cdef Hashable
cdef Iterable
cdef Iterator
cdef Sequence
cdef Mapping

cdef PathLike

cdef Scalar
cdef ArrayInt
cdef ArrayReal
cdef ArrayComplex
cdef ArrayScalar

# --- TAO ---

cdef TAOObjectiveFunction
cdef TAOGradientFunction
cdef TAOObjectiveGradientFunction
cdef TAOHessianFunction
cdef TAOUpdateFunction
cdef TAOMonitorFunction
cdef TAOConvergedFunction
cdef TAOJacobianFunction
cdef TAOResidualFunction
cdef TAOJacobianResidualFunction
cdef TAOVariableBoundsFunction
cdef TAOConstraintsFunction

# --- TS ---

cdef TSRHSFunction
cdef TSRHSJacobian
cdef TSIFunction
cdef TSIJacobian
cdef TSI2Function
cdef TSI2Jacobian
cdef TSMonitorFunction
cdef TSEventHandlerFunction
cdef TSPostEventFunction
cdef TSPreStepFunction
cdef TSPostStepFunction
