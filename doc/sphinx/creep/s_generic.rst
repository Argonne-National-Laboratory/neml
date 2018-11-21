Generic creep law
=================

Overview
--------

This object implements the generic creep law

.. math::
   \dot{\varepsilon}^{cr} = \exp\left(f\left(\log \sigma_{eq}\right) \right)

where :math:`f` is a generic function implemented as an 
:doc:`interpolate <../interpolate>` object.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``cfn``, :cpp:class:`neml::Interpolate`, Creep rate function, No

Class Description
-----------------

.. doxygenclass:: neml::GenericCreep
   :members:
   :undoc-members:
