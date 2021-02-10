Walker flow rule
================

Overview
---------

This class implements the full Walker viscoplastic flow rule
described in :doc:`../walker`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``eps0``, :cpp:class:`neml::Interpolate`, Reference strain rate, No
   ``softening``, :cpp:class:`neml::SofteningModeling`, Softening model, No
   ``scaling``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No
   ``n``, :cpp:class:`neml::Interpolate`, Rate sensitivity, No
   ``k``, :cpp:class:`neml::Interpolate`, Constant part of drag stress, No
   ``m``, :cpp:class:`neml::Interpolate`, Drag stress exponent, No
   ``R``, :cpp:class:`neml::IsotropicHardening`, Isotropic hardening model, No
   ``D``, :cpp:class:`neml::DragStress`, Drag stress model, No
   ``X``, :cpp:class:`neml::KinematicHardening`, Kinematic hardening model, No

Class description
-----------------

.. doxygenclass:: neml::WalkerFlowRule
   :members:
   :undoc-members:
