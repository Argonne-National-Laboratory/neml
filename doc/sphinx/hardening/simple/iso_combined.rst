Combined isotropic hardening
============================

Overview
--------

This object combines several isotropic hardening models into a single effective
value of the isotropic hardening variable.
The isotropic hardening variable is defined as

.. math::
   Q=\sum_{i=1}^{n}Q_{i}\left(\alpha\right).

That is, it calls each individual isotropic hardening function in a list and sums
the results.

The model requires a single history variable (:math:`\alpha`)
and maps to a single hardening variable (:math:`Q`).

.. WARNING::
   All of the NEML yield surfaces assume the opposite of the standard
   sign convention for isotropic and kinematic hardening.
   The hardening model is expected to return a negative value of the
   isotropic hardening stress and a negative value of the backstress.

Parameters
----------

========== =================================== ======================================= =======
Parameter  Object type                         Description                             Default
========== =================================== ======================================= =======
rules      std::vector<IsotropicHardeningRule> List of hardening models                No
========== =================================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::CombinedIsotropicHardeningRule
   :members:
   :undoc-members:
