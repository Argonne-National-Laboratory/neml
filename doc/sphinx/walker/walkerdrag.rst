Walker drag stress models
=========================

Overview
--------

Models providing the drag stress contribution in the :doc:`../walker`.

Base class
----------

.. doxygenclass:: neml::DragStress
   :members:
   :undoc-members:

Constant drag stress
--------------------

Fixed, non-evolving drag stress.

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8
   
   ``value``, :code:`double`, Value of the drag stress, No
   ``scale``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No scaling

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::ConstantDragStress
   :members:
   :undoc-members:

Walker's specific model
-----------------------

The specific drag stress evolution model described in :doc:`../walker`.

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``do``, :cpp:class:`neml::Interpolate`, Hardening prefactor, No
   ``d1``, :cpp:class:`neml::Interpolate`, Recovery prefactor, No
   ``d2``, :cpp:class:`neml::Interpolate`, Recovery exponent, No
   ``D_xi``, :cpp:class:`neml::Interpolate`, Saturated drag stress, No
   ``D_0``, :code:`double`, Initial drag stress, No
   ``softening``, :cpp:class:`neml::SofteningModel`, Softening model, No
   ``scale``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No scaling

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::WalkerDragStress
   :members:
   :undoc-members:
