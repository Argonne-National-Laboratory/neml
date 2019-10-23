Isotropic hardening
===================

Overview
--------

This object provides the interface for all simple isotropic hardening 
models.
These models all provide a :math:`\mathbf{q}` vector consisting of a
single entry, :math:`Q`, the isotropic hardening variable.
The function maps between a single "strain-like" variable :math:`\alpha`
and this scalar isotropic hardening variable.
The interface is then

.. math::
   Q, \frac{\partial Q}{\partial \alpha}
   \leftarrow
   \mathcal{I}\left(\alpha, T \right).

Implementations
---------------

.. toctree::
   iso_linear
   iso_interpolated
   iso_voce
   iso_power
   iso_combined

Class description
-----------------

.. doxygenclass:: neml::IsotropicHardeningRule
   :members:
   :undoc-members:
