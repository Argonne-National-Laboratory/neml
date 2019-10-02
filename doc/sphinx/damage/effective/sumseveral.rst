Sum of several effective stresses
=================================

Overview
--------

The effective stress is the weighted sum of several effective stresses:

.. math::
   \sigma_e = \sum_{i}w^{(i)}\sigma_e^{(i)}

where each individual effective stress is given by an :ref:`effective-stress` object.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``measures``, :c:type:`std::vector<`:cpp:class:`neml::EffectiveStress`:c:type:`>`, List of effective stress objects, No
   ``weights``, :c:type:`std::vector<`:c:type:`double`:c:type:`>`, List of weights, No

Class description
-----------------

.. doxygenclass:: neml::SumSeveralEffectiveStress
   :members:
   :undoc-members:
