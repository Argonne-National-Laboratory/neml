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

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``A``, :cpp:class:`neml::Interpolate`, Prefactor, No
   ``n``, :cpp:class:`neml::Interpolate`, Exponent, No

Class description
-----------------

.. doxygenclass:: neml::PowerLawCreep
   :members:
   :undoc-members:
