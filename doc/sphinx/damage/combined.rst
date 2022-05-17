Combined scalar damage models
=============================

Overview
--------

This object combines several different scalar damage models by applying them
all to the same base material model.
It combines the damage models additively, so that the total damage at
any time step is

.. math::
   \omega_{n+1} = \omega_n + \sum_{i=1}^{n}\Delta\omega_{i}

where :math:`\Delta\omega_i` is the increment in damage from the 
ith damage model.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``models``, :code:`std::vector<`:cpp:class:`neml::NEMLScalarDamagedModel_sd`:code:`>`, List of damage models to apply, No

Class description
-----------------

.. doxygenclass:: neml::CombinedDamage
   :members:
   :undoc-members:
