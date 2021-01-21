Walker softening models
=======================

Overview
--------

Softening models providing the :math:`\Phi` term in the :doc:`../walker`.

Base class
----------

.. doxygenclass:: neml::SofteningModel
   :members:
   :undoc-members:

Walker's specific model
-----------------------

The specific softening model described in :doc:`../walker`.

Parameters
^^^^^^^^^^

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``phi0``, :cpp:class:`neml::Interpolate`, Softening prefactor, No
   ``phi1``, :cpp:class:`neml::Interpolate`, Softening exponent, No

Class description
^^^^^^^^^^^^^^^^^

.. doxygenclass:: neml::WalkerSofteningModel
   :members:
   :undoc-members:
