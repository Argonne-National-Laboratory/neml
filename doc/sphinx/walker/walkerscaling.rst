Walker thermal scaling models
=============================

Overview
--------

Scaling models providing the :math:`\chi` term in the :doc:`../walker`.

Base class
----------

.. doxygenclass:: neml::ThermalScaling
   :members:
   :undoc-members:

Walker's specific model
-----------------------

The specific scaling model described in :doc:`../walker`.

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``Q``, :cpp:class:`neml::Interpolate`, Activation energy, No
   ``R``, :code:`double`, Universal gas constant, No
   ``Tref``, :code:`double`, Reference temperature, No

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::ArrheniusThermalScaling
   :members:
   :undoc-members:
