Interpolated isotropic hardening
================================

Overview
--------

This object provides interpolated isotropic hardening as a function of
the equivalent plastic strain.
The isotropic hardening variable is defined as

.. math::
   Q = -w(\alpha)

for some interpolation function :math:`w`.

The model requires a single history variable (:math:`\alpha`)
and maps to a single hardening variable (:math:`Q`).

.. WARNING::
   All of the NEML yield surfaces assume the opposite of the standard
   sign convention for isotropic and kinematic hardening.
   The hardening model is expected to return a negative value of the
   isotropic hardening stress and a negative value of the backstress.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``flow``, :cpp:class:`neml::Interpolate`, Interpolated flow stress, No

Class description
-----------------

.. doxygenclass:: neml::InterpolatedIsotropicHardeningRule
   :members:
   :undoc-members:
