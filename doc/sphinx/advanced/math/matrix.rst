.. _matrix-system:

Matrix classes
==============

The purpose of this small class system is to provide an easy way for
users to design slip interaction matrices.  These are square matrices
constructed with certain typical block patterns.

FlatVector
----------

This is just an interface to a standard, dense vector use to provide
matrix-vector products using the Matrix classes

.. doxygenclass:: neml::FlatVector
   :members:
   :undoc-members:

Matrix
------

Generic superclass of all matrices, could be expanded to handle sparse
matrices if required in the future.

.. doxygenclass:: neml::Matrix
   :members:
   :undoc-members:

SquareMatrix
------------

Dense, square matrix with various helpers for setting up block structures.
Initialization options are:

================  =============================================================
Option            Description
================  =============================================================
zero              The zero matrix
identity          The identity matrix
diagonal          Diagonal matrix where the user gives the diagonal entries
diagonal_blocks   Diagonal matrix where the entries are blocked
block             General block matrix
dense             Dense matrix, user gives all entries in row-major form
================  =============================================================

.. doxygenclass:: neml::SquareMatrix
   :members:
   :undoc-members:
