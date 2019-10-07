Swindeman minimum creep rate model
==================================

Overview
--------

This object implements the minimum creep rate model given in [S1999]_:

.. math::
   \dot{\varepsilon}^{cr} = C \sigma_{eq}^n \exp{V \sigma_{eq}} \exp{\frac{-Q}{T}}

for constants :math:`C`, :math:`n`, :math:`V`, and :math:`Q`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``C``, :c:type:`double`, Prefactor, No
   ``n``, :c:type:`double`, Exponent, No
   ``V``, :c:type:`double`, Exponential factor, No
   ``Q``, :c:type:`double`, Activation energy, No

Class description
-----------------

.. doxygenclass:: neml::SwindemanMinimumCreep
   :members:
   :undoc-members:
