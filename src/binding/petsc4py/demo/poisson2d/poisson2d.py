# Poisson in 2D
# =============
#
# Solve a constant coefficient Poisson problem on a regular grid
#
# .. math::
#
#     - u_{xx} - u_{yy} = 1 \quad\textsf{in}\quad [0,1]^2\\
#     u = 0 \quad\textsf{on the boundary.}

# This is a naïve, parallel implementation, using :math:`n` interior grid
# points per dimension and a lexicographic ordering of the nodes.

# This code is kept as simple as possible to show the minimal differences with
# the sequential version. However, simplicity comes at a price. Here we use a
# naive decomposition that does not lead to an optimal communication complexity
# for the matrix-vector product. See dmda.c for an example of optimal
# communication complexity and nodes ordering.
#
# This demo is structured as a script to be executed using:
#
# .. code-block:: console
#
#   $ python poisson2d.py
#
# potentially with additional options passed at the end of the command.
#
# At the start of your script, call `petsc4py.init` passing `sys.argv` so that
# command-line arguments to the script are passed through to PETSc.

import sys
import petsc4py

petsc4py.init(sys.argv)

# import PETSc module

from petsc4py import PETSc

# Access options database
OptDB = PETSc.Options()

# Grid size and spacing
# default values are specified within getters
n = OptDB.getInt('n', 5)
h = 1.0 / (n + 1)

# Create sparse matrix
# You can omit the comm argument if your objects leave on PETSc.COMM_WORLD
# but it is a dangerous choice to rely on default values for such
# important arguments
A = PETSc.Mat()
A.create(comm=PETSc.COMM_WORLD)

# Specify global matrix size with tuple
A.setSizes((n * n, n * n))
# the call above implicitly assumes using PETSC_DECIDE for local sizes
# and it is equivalent to
#     A.setSizes(((PETSc.DECIDE, n * n), (PETSc.DECIDE, n * n)))

# Set type and customize the object with command line
A.setType(PETSc.Mat.Type.AIJ)
A.setFromOptions()

# Set preallocation
# this call already handles the AIJ derived types for us
A.setPreallocationNNZ(5)

index_to_grid = lambda r: (r // n, r % n)

rstart, rend = A.getOwnershipRange()
for row in range(rstart, rend):
    i, j = index_to_grid(row)
    A[row, row] = 4.0 / h**2
    if i > 0:
        column = row - n
        A[row, column] = -1.0 / h**2
    if i < n - 1:
        column = row + n
        A[row, column] = -1.0 / h**2
    if j > 0:
        column = row - 1
        A[row, column] = -1.0 / h**2
    if j < n - 1:
        column = row + 1
        A[row, column] = -1.0 / h**2

A.assemblyBegin()
A.assemblyEnd()

A.viewFromOptions('-view_mat')

# Create Krylov solver context
ksp = PETSc.KSP()
ksp.create(comm=A.getComm())
ksp.setType(PETSc.KSP.Type.CG)
ksp.getPC().setType(PETSc.PC.Type.GAMG)
ksp.setOperators(A)
ksp.setFromOptions()

# Solve linear system
x, b = A.createVecs()
b.set(1.0)
ksp.solve(b, x)
x.viewFromOptions('-view_sol')
