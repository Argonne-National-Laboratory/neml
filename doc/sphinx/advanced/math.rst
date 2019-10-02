Mathematical helpers
====================

NEML contains many helper functions to do tensor operations on 
Mandel-notation representations of tensors as well as common linear algebra
operations like solve systems of equations and invert matrices.
These helpers are fairly self documenting.
The interfaces are shown below.

Additionally, the :ref:`crystal-plasticity` module uses a new system for
tensor operations using objects, rather than raw pointers.
These objects are documented in the :ref:`crystal-plasticity` section.
The new
system considerably simplifies many common mathematical operations.
The long-term plan for NEML is to switch the macroscale plasticity modules
over to this new approach.

.. doxygenfile:: math/nemlmath.cxx
