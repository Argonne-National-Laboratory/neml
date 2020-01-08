Scalar creep models
===================

Overview
--------

Scalar creep models provide the interface

.. math::
   \dot{\varepsilon}^{cr}
   \leftarrow
   \mathcal{C}\left(\sigma_{eq}, \varepsilon_{eq}, t, T \right).

They return a scalar creep rate as a function of effective stress, effective
strain, time, and temperature.

Implementations
---------------

.. toctree::

   s_power
   s_normpower
   s_km
   s_norton
   s_mukherjee
   s_generic
   s_blackburn
   s_swindeman

Class description
-----------------

.. doxygenclass:: neml::ScalarCreepRule
