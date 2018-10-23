Kinematic hardening
===================

Overview
--------

This object provides the interface for all simple kinematic hardening models.
These models provide a :math:`\mathbf{q}` consisting of a single 
backstress :math:`\mathbf{X}`, implemented as a length 6
Mandel vector representing a full 2nd order tensor.
The function maps between this backstress and similar backstrain 
:math:`\bm{\chi}`, likewise
implemented with a length 6 Mandel vector representing a full 2nd
order tensor.

The interface is

.. math::
   \mathbf{X}, \frac{\partial\mathbf{X}}{\partial \bm{\chi}}
   \leftarrow
   \mathcal{K}\left(\bm{\chi}, T\right).

Implementations
---------------

.. toctree::
   kin_linear

Class description
-----------------

.. doxygenclass:: neml::KinematicHardeningRule
   :members:
   :undoc-members:
