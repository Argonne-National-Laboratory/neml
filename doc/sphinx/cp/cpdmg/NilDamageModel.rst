NilDamageModel
==============

Overview
--------

This model is mainly for testing the :any:`neml::DamagedStandardKinematicModel` formulation.
It returns the identify for the damage projection operator:

.. math::

   \mathbf{P} = \mathbf{I}

and does not maintain any internal variables.

Parameters
----------

None

Class description
-----------------

.. doxygenclass:: neml::NilDamageModel
   :members:
   :undoc-members:

