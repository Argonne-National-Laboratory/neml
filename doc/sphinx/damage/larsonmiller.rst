Larson Miller Damage
====================

Overview
--------

This model implements the damage model

.. math::
   \dot{\omega} = \frac{1}{t_R\left(\sigma_{eff} \left(1-\omega\right) \right)}

where :math:`\sigma_{eff}` is a modular effective stress, defined by a :ref:`effective-stress` object, and :math:`t_R` is a time-to-rupture Larson-Miller relation provided
by a :ref:`larson-miller` object.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``lmr``, :cpp:class:`neml::LarsonMillerCorrelation`, Parameter, No
   ``estress``,:cpp:class:`neml::EffectiveStress`, Effective stress, No

Class description
-----------------

.. doxygenclass:: neml::LarsonMillerCreepDamage
   :members:
   :undoc-members:
