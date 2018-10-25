Norton-Bailey creep
===================

Overview
--------

This implements the Norton-Bailey creep model

.. math::
   \varepsilon^{cr} = A \sigma_{eff}^n t^m

implemented in the strain-hardening formulation

.. math::
   \dot{\varepsilon}^{cr} = m A^\frac{1}{m} \sigma^\frac{n}{m} \varepsilon_{eff}^\frac{m-1}{m}.


Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``A``, :cpp:class:`neml::Interpolate`, Prefactor, No
   ``m``, :cpp:class:`neml::Interpolate`, Stress exponent, No
   ``n``, :cpp:class:`neml::Interpolate`, Time exponent, No

Class description
-----------------

.. doxygenclass:: neml::NortonBaileyCreep
   :members:
   :undoc-members:
