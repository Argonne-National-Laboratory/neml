Linear kinematic hardening
==========================

Overview
--------

This object implements simple linear kinematic hardening.
The backstress is defined as

.. math::
   \mathbf{X} = -H \bm{\chi}.

The model requires a length 6 history variable (:math:`\bm{\chi}`)
and maps to a length 6 backstress (:math:`\mathbf{X}`).

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
H          Interpolate                         Linear hardening constant               No
========== =================================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::LinearKinematicHardeningRule
   :members:
   :undoc-members:
