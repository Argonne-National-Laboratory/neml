Power law isotropic hardening
=============================

Overview
--------

This object provides power law isotropic hardening.
The isotropic hardening variable is defined as

.. math::
   Q = -\sigma_0 - A \alpha^n.

The model requires a single history variable (:math:`\alpha`), which is
the equivalent plastic strain,
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

   ``s0``, :cpp:class:`neml::Interpolate`, Initial yield stress, No
   ``A``, :cpp:class:`neml::Interpolate`, Prefactor, No
   ``n``, :cpp:class:`neml::Interpolate`, Work hardening exponent, No

Class description
-----------------

.. doxygenclass:: neml::PowerLawIsotropicHardeningRule
   :members:
   :undoc-members:
