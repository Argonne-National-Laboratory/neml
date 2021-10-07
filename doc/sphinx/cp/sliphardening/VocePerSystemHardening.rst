VocePerSystemHardening
======================

Overview
--------

This model provides a Voce hardening response on each individual slip system.
The hardening evolution on each system is independent of all other systems 
and each slip system can have its own Voce parameters:

.. math::
   \dot{\bar{\tau}}_{k} = k_k \left(1 - \frac{\bar{\tau}_k - \tau_{0,k}}{\tau_{sat,k} - \tau_{0,k}} \right)^{m}

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``k``, :cpp:class:`std::vector<neml::Interpolate>`, Hardening prefactors, N
   ``saturation``, :cpp:class:`std::vector<neml::Interpolate>`, Saturated hardening values, N 
   ``m``, :cpp:class:`std::vector<neml::Interpolate>`, Voce exponents, N
   ``initial``, :code:`std::vector<double>`, Initial strengths, N

Class description
-----------------

.. doxygenclass:: neml::VocePerSystemHardening
   :members:
   :undoc-members:
