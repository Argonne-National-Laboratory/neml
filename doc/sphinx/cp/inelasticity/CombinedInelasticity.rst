CombinedInelasticity
====================

Overview
--------

This class sums the plastic deformation rates from several InelasticModel objects
into a single response.  Mathematically then the plastic deformation becomes

.. math::
   \mathbf{D}^p = \sum_{i=1}^{n_{models}} \mathbf{D}^p_i

   \mathbf{W}^p = \sum_{i=1}^{n_{models}} \mathbf{W}^p_i

where the :math:`i` subscripts indicate the response of each individual model.

The history variables and associated history rates and derivatives from all
the individual models are `concatenated` together.  That is, the final
set of history variables includes the history variables of all the 
individual models.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``models``, :c:type:`std::vector<`:cpp:class:`neml::InelasticModel`:c:type:`>`, Individual models, No

Class description
-----------------

.. doxygenclass:: neml::CombinedInelasticity
   :members:
