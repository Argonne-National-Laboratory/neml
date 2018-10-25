Blank elastic model
===================

Overview
--------

This model is only used a signal that the object it is assigned to should
take its elastic constants from the next object up the hierarchy.
For example, in a meta-model like the :doc:`../interfaces/creep_plasticity` 
model this would mean the base material model and the creep model should
both take elastic constants from the parent CreepPlasticity object.

Parameters
----------

None

Class description
-----------------

.. doxygenclass:: neml::BlankElasticModel
   :members:
   :undoc-members:
