Simple hardening rules
======================

Overview
--------

These rules are intended for use with associative plasticity, where the
evolution rate of each hardening variable is determined by the map between
the "strain-like" history variable and the "stress-like" internal variable
that feeds into the yield surface.
These models implement the interface

.. math::
   \mathbf{q}, \frac{\partial \mathbf{q}}{\partial \bm{\alpha}} 
   \leftarrow \mathcal{H}\left( \bm{\alpha}, T \right)

where :math:`\bm{\alpha}` are the actual "strain-like" history variables
tracked and maintained by the model and :math:`\mathbf{q}` are the stress-like
internal variables provided to the yield surface.

Implementations
---------------

.. toctree::
   simple/isotropic
   simple/kinematic
   simple/combined

Class description
-----------------

.. doxygenclass:: neml::HardeningRule
   :members:
   :undoc-members:
