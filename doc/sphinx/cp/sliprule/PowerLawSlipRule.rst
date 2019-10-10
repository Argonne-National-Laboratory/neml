PowerLawSlipRule
================

Overview
--------

This implements the standard power law model for the slip rate:

.. math::
   \dot{\gamma}_{g,i} = \dot{\gamma}_0 \left| \frac{\tau_{g,i}}{\bar{\tau}_{g,i}} \right|^{\left(n-1\right)} \frac{\tau_{g,i}}{\bar{\tau}_{g,i}}

where :math:`\dot{\gamma}_0`, the reference slip rate, and :math:`n`, the rate sensitivity, are temperature-dependent parameters.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``strength``, :cpp:class:`neml::SlipHardening`, Slip hardening definition, No
   ``gamma0``, :cpp:class:`neml::Interpolate`, Reference slip rate, No
   ``n``, :cpp:class:`neml::Interpolate`, Rate sensitivity, No

Class description
-----------------

.. doxygenclass:: neml::PowerLawSlipRule
