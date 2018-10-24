Power law creep
===============

Overview
--------

This object implements power law creep

.. math::
   \dot{\varepsilon}^{cr} = A \sigma_{eq}^n

for temperature dependent parameters :math:`A` and :math:`n`.

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
A          Interpolate               Prefactor                               No
n          Interpolate               Exponent                                No
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::PowerLawCreep
   :members:
   :undoc-members:
