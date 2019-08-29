Maximum of several effective stresses
=====================================

Overview
--------

The effective stress the maximum of several effective stresses:

.. math::
   \sigma_e = \max\left(\sigma_e^{(1)}, \sigma_e^{(2)}, ... \right)

where each individual effective stress is given by an :ref:`effective-stress` object.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``measures``, :c:type:`std::vector<`:cpp:class:`neml::EffectiveStress`:c:type:`>`, List of effective stress objects, No

Class description
-----------------

.. doxygenclass:: neml::MaxSeveralEffectiveStress
   :members:
   :undoc-members:
