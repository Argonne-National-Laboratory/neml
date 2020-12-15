Walker isotropic hardening models
=================================

Overview
--------

Models providing the isotropic hardening contribution in the :doc:`../walker`.

Base class
----------

.. doxygenclass:: neml::IsotropicHardening
   :members:
   :undoc-members:

Constant isotropic hardening
----------------------------

Isotropic hardening fixed to zero.

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``scale``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No scaling

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::ConstantIsotropicHardening
   :members:
   :undoc-members:

Walker's specific model
-----------------------

The specific model described in :doc:`../walker`.

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``r0``, :cpp:class:`neml::Interpolate`, Hardening rate prefactor, No
   ``Rinf``, :cpp:class:`neml::Interpolate`, Final value of isotropic hardening, No
   ``r1``, :cpp:class:`neml::Interpolate`, Static recovery prefactor, No
   ``r2``, :cpp:class:`neml::Interpolate`, Static recovery exponent, No
   ``scale``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No scaling

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::WalkerIsotropicHardening
   :members:
   :undoc-members:
