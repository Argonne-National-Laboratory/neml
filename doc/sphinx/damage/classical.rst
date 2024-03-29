Classical creep damage
======================

Overview
--------

This object implements the classical Hayhurst-Leckie-Rabotnov-Kachanov creep damage model [HL1977]_.
The damage update is given by 

.. math::
   \dot{\omega} = \left(\frac{\sigma_{eff}}{A}\right)^\xi 
      \left(1 - \omega\right)^{-\phi}

   \sigma_{eff} = \sqrt{\frac{3}{2} \operatorname{dev}\left(\bm{\sigma}\right):\operatorname{dev}\left(\bm{\sigma}\right)}.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``A``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``xi``, :cpp:class:`neml::Interpolate`, Stress exponent, No
   ``phi``, :cpp:class:`neml::Interpolate`, Damage exponent, No

Class description
-----------------

.. doxygenclass:: neml::ClassicalCreepDamage
   :members:
   :undoc-members:
