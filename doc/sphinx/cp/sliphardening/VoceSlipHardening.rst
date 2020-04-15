VoceSlipHardening
=================

Overview
--------

This class implements a Voce hardening model with

.. math::
   \tau_0 = \tau_0

   f = b \left(\tau_{sat}-\tilde{\tau}\right)

   \tau_{nye} = k \sqrt{\left\Vert \bm{\alpha}\right\Vert _{F}}

where :math:`\tau_s`, :math:`b`, :math:`\tau_{sat}`, and :math:`k` all temperature-dependent material properties.

The value of :math:`k` defaults to zero, which means by default the model is
independent of the Nye tensor.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``tau_sat``, :cpp:class:`neml::Interpolate`, Saturation strength, N
   ``b``, :cpp:class:`neml::Interpolate`, Rate parameter, N
   ``tau_0``, :cpp:class:`neml::Interpolate`, Static strength, N
   ``k``, :cpp:class:`neml::Interpolate`, Nye hardening constant, 0

Class description
-----------------

.. doxygenclass:: neml::VoceSlipHardening
   :members:
   :undoc-members:
