Modular creep damage
====================

Overview
--------

This object implements the classical Hayhurst-Leckie-Rabotnov-Kachanov creep damage model [HL1977]_.  Specifically, it exactly replicates Eq. 2.6 in that paper, so that the implementation exactly replicates the analytic expressions in Eq. 2.5.
The damage update is given by 

.. math::
   \dot{\omega} = \left(\frac{\sigma_{eff}}{A}\right)^\xi 
      \left(1 - \omega\right)^{\xi-\phi}

where :math:`\sigma_{eff}` is a modular effective stress, defined by a :ref:`effective-stress` object.

.. toctree::
   effective/effective

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``A``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``xi``, :cpp:class:`neml::Interpolate`, Stress exponent, No
   ``phi``, :cpp:class:`neml::Interpolate`, Damage exponent, No
   ``estress``,:cpp:class:`neml::EffectiveStress`, Effective stress, No

Class description
-----------------

.. doxygenclass:: neml::ModularCreepDamage
   :members:
   :undoc-members:
