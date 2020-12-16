FASlipHardening
===============

Overview
--------

This model provides an empirical Frederick-Armstrong type hardening response on each individual slip system.
The hardening evolution on each system is independent of all other systems 
and each slip system can have its own parameters

.. math::
   \dot{\bar{\tau}}_{k} = k_k \left(\dot{\gamma}_k - \frac{\bar{\tau}_k}{\tau_{sat,k}} \left| \dot{\gamma}_k \right| \right)

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``k``, :cpp:class:`std::vector<neml::Interpolate>`, Hardening prefactors, N
   ``sat``, :cpp:class:`std::vector<neml::Interpolate>`, Saturated hardening values, N 

Class description
-----------------

.. doxygenclass:: neml::FASlipHardening
   :members:
   :undoc-members:
