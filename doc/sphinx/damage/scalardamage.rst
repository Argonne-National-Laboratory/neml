Scalar damage models
====================

Overview
--------

This object defines a scalar damage model with the interface

.. math::
   \omega_{n+1}, \frac{\partial \omega_{n+1}}{\partial \omega}, \frac{\partial \omega_{n+1}}{\partial \bm{\sigma}}, \frac{\partial \omega_{n+1}}{\partial \bm{\varepsilon}} \leftarrow \mathcal{W}\left( \omega_{n+1}, \omega{n}, \bm{\sigma}_{n+1}, \bm{\sigma}_{n} \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n}, T_{n+1}, T_{n}, t_{n+1}, t_{n}, \right)

where :math:`\omega_{n+1}` is the current value of the scalar damage parameter.

Implementations
---------------

.. toctree::
   scalarrate
   combined
   standard
   work

Class description
-----------------

.. doxygenclass:: neml::ScalarDamage
   :members:
   :undoc-members:
