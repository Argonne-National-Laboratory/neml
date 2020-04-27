.. _math-helpers:

Mathematical helpers
====================

NEML contains many helper functions to do tensor operations on 
Mandel-notation representations of tensors stored as pointers as well as common linear algebra
operations like solving systems of equations and inverting matrices.
These helpers are fairly self documenting.
The interfaces are shown below.

Additionally, the :ref:`crystal-plasticity` module uses a new system for
tensor operations using objects, rather than raw pointers.
The new
system considerably simplifies many common mathematical operations.
The long-term plan for NEML is to switch the macroscale plasticity modules
over to this new approach.

Tensor and rotation objects
---------------------------

.. toctree::
   :maxdepth: 1

   math/rotations
   math/tensors

Generic matrix system
---------------------

.. toctree::
   :maxdepth: 1

   math/matrix

Pointer array math functions
----------------------------

.. doxygenfile:: math/nemlmath.cxx
