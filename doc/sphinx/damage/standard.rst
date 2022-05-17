Standard damage 
===============

Overview
--------

This object implements a "standard" damage model where the single
damage variable varies only with the scalar effective inelastic 
strain rate.
This simplifies the damage update to

.. math::
   \omega_{n+1} = \omega_n + w\left(\bm{\sigma}_{n+1}, \omega_{n+1}\right) 
      \Delta \varepsilon_{eff}^{in}.

A separate interface defines the damage update function :math:`w`.

Implementations
---------------

.. toctree::
   powerlaw
   exponential

Class description
-----------------

.. doxygenclass:: neml::StandardScalarDamage
   :members:
   :undoc-members:
