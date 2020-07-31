FixedStrengthHardening
======================

Overview
--------

This class implements the simplest hardening rule -- all slip systems have a constant strength

.. math::
   \bar{\tau}_{g,i}=\tau_{g,i}

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``tau_sat``, :c:type:`std::vector<`:cpp:class:`neml::Interpolate`:c:type:`>`, Unrolled slip system hardening vector, N

Class description
-----------------

.. doxygenclass:: neml::FixedStrengthHardening
   :members:
   :undoc-members:
