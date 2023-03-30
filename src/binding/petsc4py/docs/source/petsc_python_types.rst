.. _petsc_python_types:

PETSc Python types
==================

Here we discuss details about Python-aware PETSc types that can be used within the library.

.. _petsc_python_mat:

PETSc Python matrix type
========================

PETSc provides a convenient way to compute the action of linear operators coded
in Python through the MATPYTHON class.

The implementation is not limited to the action only, but it also exposes
additional methods that can be used within the library. A template class for
the supported methods is given below.

.. literalinclude:: ../../demo/python_types/matpython_protocol.py

The low-level, Cython implementation exposing the methods is contained in

**${PETSC_DIR}/src/binding/petsc4py/src/petsc4py/PETSc/libpetsc4py.pyx**

In the example below, we create an operator the applies the Laplacian operator
on a two-dimensional grid, and use it solve the associated linear system.

.. literalinclude:: ../../demo/python_types/mat.py

.. _petsc_python_ksp:

PETSc Python linear solver type
===============================

.. _petsc_python_snes:

PETSc Python nonlinear solver type
==================================

.. _petsc_python_tao:

PETSc Python optimization solver type
=====================================
