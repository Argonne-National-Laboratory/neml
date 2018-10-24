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
   s_km
   s_norton
   s_mukherjee

Class description
-----------------

.. doxygenclass:: neml::ScalarCreepRule
