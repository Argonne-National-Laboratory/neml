LinearSlipHardening
===================

Overview
--------

This class implements linear hardening with

.. math::
   \tau_0 = \tau_0

   f = k_1

   \tau_{nye} = k_2 \left\Vert \bm{\alpha}\right\Vert _{F}

where :math:`k_1` and :math:`k_2` are temperature-dependent material properties.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``k_1``, :cpp:class:`neml::Interpolate`, Slip hardening coefficient, N
   ``k_2``, :cpp:class:`neml::Interpolate`, Nye hardening coefficient, N

Class description
-----------------

.. doxygenclass:: neml::LinearSlipHardening
   :members:
   :undoc-members:
