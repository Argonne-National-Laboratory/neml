Blackburn minimum creep rate model
==================================

Overview
--------

This object implements the minimum creep rate model given in [B1972]_:

.. math::
   \dot{\varepsilon}^{cr} = A \sinh{\frac{\beta \sigma_{eq}}{n}}^n \exp{\frac{-Q}{RT}}

for temperature dependent parameters :math:`A`, :math:`\beta`, and :math:`\beta` and constants :math:`R` and :math:`Q`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``A``, :cpp:class:`neml::Interpolate`, Prefactor, No
   ``n``, :cpp:class:`neml::Interpolate`, Exponent, No
   ``beta``, :cpp:class:`neml::Interpolate`, Exponential factor, No
   ``R``, :c:type:`double`, Gas constant, No
   ``Q``, :c:type:`double`, Activation energy, No

Class description
-----------------

.. doxygenclass:: neml::BlackburnMinimumCreep
   :members:
   :undoc-members:
