Power law creep
===============

Overview
--------

This object implements a better-normalized version of power law creep

.. math::
   \dot{\varepsilon}^{cr} = \left(\frac{\sigma_{eq}}{\sigma_0}\right)^n

for temperature dependent parameters :math:`\sigma_0` and :math:`n`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``s0``, :cpp:class:`neml::Interpolate`, Normalization stress, No
   ``n``, :cpp:class:`neml::Interpolate`, Exponent, No

Class description
-----------------

.. doxygenclass:: neml::NormalizedPowerLawCreep
   :members:
   :undoc-members:
