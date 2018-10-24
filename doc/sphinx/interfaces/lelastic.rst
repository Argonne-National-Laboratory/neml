Linear elasticity
=================

Overview
--------

This object implements a simple linear elastic stress update of the type

.. math::
   \bm{\sigma}_{n+1} = \mathbf{\mathfrak{C}}_{n+1} : \bm{\varepsilon}_{n+1}

where :math:`\mathbf{\mathfrak{C}}` is a temperature-dependent elasticity
tensor.
The model does not maintain any  history variables.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``emodel``, :cpp:class:`neml::LinearElasticModel`, "Temperature dependent elastic constants", "No"
   ``alpha``, :cpp:class:`neml::Interpolate`, "Temperature dependent instantaneous CTE", ``0.0``

Class description
-----------------

.. doxygenclass:: neml::SmallStrainElasticity
   :members:
   :undoc-members:
