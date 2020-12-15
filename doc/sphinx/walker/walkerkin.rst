Walker kinematic hardening models
=================================

Overview
--------

Models providing the kinematic hardening contribution in the :doc:`../walker`.

Base class
----------

.. doxygenclass:: neml::KinematicHardening
   :members:
   :undoc-members:

Frederick-Armstrong hardening
-----------------------------

The Frederick-Armstrong [FA2007]_ model, implemented in the Walker subsystem:

.. math::
   \dot{X} = c \dot{\bm{\varepsilon}}_{vp} - g \bm{X}

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8
   
   ``c``, :cpp:class:`neml::Interpolate`, Hardening constant, No
   ``g``, :cpp:class:`neml::Interpolate`, Dynamic recover constant, No
   ``scale``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No scaling

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::FAKinematicHardening
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

   ``c0``, :cpp:class:`neml::Interpolate`, Constant hardening parameter, No
   ``c1``, :cpp:class:`neml::Interpolate`, Hardening evolution prefactor, No
   ``c2``, :cpp:class:`neml::Interpolate`, Hardening evolution exponent, No
   ``l0``, :cpp:class:`neml::Interpolate`, Dynamic recovery exponential rate, No
   ``l1``, :cpp:class:`neml::Interpolate`, Dynamic recovery evolution prefactor, No
   ``l``, :cpp:class:`neml::Interpolate`, Constant dynamic recovery coeficient, No
   ``b0``, :cpp:class:`neml::Interpolate`, Recovery direction constant, No
   ``x0``, :cpp:class:`neml::Interpolate`, Static recovery prefactor, No
   ``x1``, :cpp:class:`neml::Interpolate`, Static recovery exponent, No
   ``softening``, :cpp:class:`neml::SofteningModel`, Softening model, No
   ``scale``, :cpp:class:`neml::ThermalScaling`, Thermal scaling model, No scaling

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::WalkerKinematicHardening
   :members:
   :undoc-members:
