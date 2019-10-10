VoceSlipHardening
=================

Overview
--------

This class implements a Voce hardening model with

.. math::
   \tau_0 = \tau_s

   f = b \left(\tau_{sat}-\tilde{\tau}\right)

where :math:`\tau_s`, :math:`b`, and :math:`\tau_{sat}` all temperature-dependent material properties

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``tau_sat``, :cpp:class:`neml::Interpolate`, Saturation strength, N
   ``b``, :cpp:class:`neml::Interpolate`, Rate parameter, N
   ``tau_0``, :cpp:class:`neml::Interpolate`, Static strength, N

Class description
-----------------

.. doxygenclass:: neml::VoceSlipHardening
   :members:
