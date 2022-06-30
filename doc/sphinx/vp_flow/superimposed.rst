Superimposed viscoplastic flow rule
===================================

Overview
--------

This model sums the contribution of multiple individual viscoplastic hardening
rules according to the relations:

.. math::

   \dot{\gamma} = \sum_{i=1}^{n_{models}} \dot{\gamma}_i

   \mathbf{g}_{...} = \frac{1}{\dot{\gamma}} \sum_{i=1}^{n_{models}} \dot{\gamma}_i
      \mathbf{g}_{...,i}

   \mathbf{h}_{...} = \bigoplus_{i=1}^{n_{models}} \mathbf{h}_{...,i}

where :math:`...` here is all of :math:`\gamma`, :math:`t`, and :math:`T` and
:math:`\oplus` represents concatenation.  The internal variables is then
the set of all variables maintained by each individual flow rule (and these
variables cannot overlap).

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``flow_rules``, :code:`std::vector<`:cpp:class:`neml::ViscoplasticFlowRule`:code:`>`, List of individual flow rules, No

Class description
-----------------

.. doxygenclass:: neml::SuperimposedViscoPlasticFlowRule
   :members:
   :undoc-members:
