# Author:  Lisandro Dalcin
# Contact: dalcinl@gmail.com
"""Typing support."""

from __future__ import annotations
from typing import (
    Callable,
)
from numpy.typing import (
    NDArray,
)
from .PETSc import (
    Vec,
    Mat,
    KSP,
    SNES,
    TS,
    TAO,
    DM,
)

__all__ = [
    "Scalar",
    "ArrayInt",
    "ArrayReal",
    "ArrayComplex",
    "ArrayScalar",
    "DMCoarsenHookFunction",
    "DMRestrictHookFunction",
    "KSPRHSFunction",
    "KSPOperatorsFunction",
    "KSPConvergenceTestFunction",
    "KSPMonitorFunction",
    "TSRHSFunction",
    "TSRHSJacobian",
    "TSRHSJacobianP",
    "TSIFunction",
    "TSIJacobian",
    "TSIJacobianP",
    "TSI2Function",
    "TSI2Jacobian",
    "TSI2JacobianP",
    "TSMonitorFunction",
    "TSPreStepFunction",
    "TSPostStepFunction",
    "TSEventHandlerFunction",
    "TSPostEventFunction",
    "TSPreStepFunction",
    "TSPostStepFunction",
    "TAOObjectiveFunction",
    "TAOGradientFunction",
    "TAOObjectiveGradientFunction",
    "TAOHessianFunction",
    "TAOUpdateFunction",
    "TAOMonitorFunction",
    "TAOConvergedFunction",
    "TAOJacobianFunction",
    "TAOResidualFunction",
    "TAOJacobianResidualFunction",
    "TAOVariableBoundsFunction",
    "TAOConstraintsFunction",
]

# --- Sys ---

Scalar = float | complex
"""Scalar type.

Scalars can be either `float` or `complex` (but not both) depending on how
PETSc was configured (``./configure --with-scalar-type=real|complex``).

"""

ArrayInt = NDArray[int]
"""Array of `int`."""

ArrayReal = NDArray[float]
"""Array of `float`."""

ArrayComplex = NDArray[complex]
"""Array of `complex`."""

ArrayScalar = NDArray[Scalar]
"""Array of `Scalar` numbers."""

# --- DM ---

DMCoarsenHookFunction = Callable[[DM, DM], None]
DMRestrictHookFunction = Callable[[DM, Mat, Vec, Mat, DM], None]

# --- KSP ---

KSPRHSFunction = Callable[[KSP, Vec], None]
"""`PETSc.KSP` right hand side function callback."""

KSPOperatorsFunction = Callable[[KSP, Mat, Mat], None]
"""`PETSc.KSP` operators function callback."""

KSPConvergenceTestFunction = Callable[[KSP, int, float], KSP.ConvergedReason]
"""`PETSc.KSP` convergence test callback."""

KSPMonitorFunction = Callable[[KSP, int, float], None]
"""`PETSc.KSP` monitor callback."""

# --- SNES ---

SNESMonitorFunction = Callable[[SNES, int, float], None]
"""`SNES` monitor callback."""

SNESObjFunction = Callable[[SNES, Vec], None]
"""`SNES` objective function callback."""

SNESFunction = Callable[[SNES, Vec, Vec], None]
"""`SNES` residual function callback."""

SNESJacobianFunction = Callable[[SNES, Vec, Mat, Mat], None]
"""`SNES` Jacobian callback."""

SNESGuessFunction = Callable[[SNES, Vec], None]
"""`SNES` initial guess callback."""

SNESUpdateFunction = Callable[[SNES, int], None]
"""`SNES` step update callback."""

SNESLSPreFunction = Callable[[Vec, Vec], None]
"""`SNES` linesearch pre-check update callback."""

SNESNGSFunction = Callable[[SNES, Vec, Vec], None]
"""`SNES` nonlinear Gauss-Seidel callback."""

SNESConvergedFunction = Callable[[SNES, int, tuple[float, float, float]], SNES.ConvergedReason]
"""`SNES` convergence test callback."""

# --- TS ---

TSRHSFunction = Callable[[TS, float, Vec, Vec], None]
"""`TS` right hand side function callback."""

TSRHSJacobian = Callable[[TS, float, Vec, Mat, Mat], None]
"""`TS` right hand side Jacobian callback."""

TSRHSJacobianP = Callable[[TS, float, Vec, Mat], None]
"""`TS` right hand side parameter Jacobian callback."""

TSIFunction = Callable[[TS, float, Vec, Vec, Vec], None]
"""`TS` implicit function callback."""

TSIJacobian = Callable[[TS, float, Vec, Vec, float, Mat, Mat], None]
"""`TS` implicit Jacobian callback."""

TSIJacobianP = Callable[[TS, float, Vec, Vec, float, Mat], None]
"""`TS` implicit parameter Jacobian callback."""

TSI2Function = Callable[[TS, float, Vec, Vec, Vec, Vec], None]
"""`TS` implicit 2nd order function callback."""

TSI2Jacobian = Callable[[TS, float, Vec, Vec, Vec, float, float, Mat, Mat], None]
"""`TS` implicit 2nd order Jacobian callback."""

TSI2JacobianP = Callable[[TS, float, Vec, Vec, Vec, float, float, Mat], None]
"""`TS` implicit 2nd order parameter Jacobian callback."""

TSMonitorFunction = Callable[[TS, int, float, Vec], None]
"""`TS` monitor callback."""

TSPreStepFunction = Callable[[TS], None]
"""`TS` pre-step callback."""

TSPostStepFunction = Callable[[TS], None]
"""`TS` post-step callback."""

TSEventHandlerFunction = Callable[[TS, float, Vec, NDArray[Scalar]], None]
"""`TS` event handler callback."""

TSPostEventFunction = Callable[[TS, NDArray[int], float, Vec, bool], None]
"""`TS` post-event handler callback."""

# --- TAO ---

TAOObjectiveFunction = Callable[[TAO, Vec], float]
"""`TAO` objective function callback."""

TAOGradientFunction = Callable[[TAO, Vec, Vec], None]
"""`TAO` objective gradient callback."""

TAOObjectiveGradientFunction =  Callable[[TAO, Vec, Vec], float]
"""`TAO` objective function and gradient callback."""

TAOHessianFunction = Callable[[TAO, Vec, Mat, Mat], None]
"""`TAO` objective Hessian callback."""

TAOUpdateFunction = Callable[[TAO, int], None]
"""`TAO` update callback."""

TAOMonitorFunction = Callable[[TAO], None]
"""`TAO` monitor callback."""

TAOConvergedFunction = Callable[[TAO], None]
"""`TAO` convergence test callback."""

TAOJacobianFunction = Callable[[TAO, Vec, Mat, Mat], None]
"""`TAO` Jacobian callback."""

TAOResidualFunction = Callable[[TAO, Vec, Vec], None]
"""`TAO` residual callback."""

TAOJacobianResidualFunction = Callable[[TAO, Vec, Mat, Mat], None]
"""`TAO` Jacobian residual callback."""

TAOVariableBoundsFunction = Callable[[TAO, Vec, Vec], None]
"""`TAO` variable bounds callback."""

TAOConstraintsFunction = Callable[[TAO, Vec, Vec], None]
"""`TAO` constraints callback."""
