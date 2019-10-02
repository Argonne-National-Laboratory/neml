Classical creep damage
======================

Overview
--------

This object implements the classical Hayhurst-Leckie-Rabotnov-Kachanov creep damage model [HL1977]_.
The damage update is given by 

.. math::
   \omega_{n+1} = \omega_{n} + \left(\frac{\sigma_{eff}}{A}\right)^\xi 
      \left(1 - \omega_{n+1}\right)^{-\phi} \Delta t_{n+1}

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
   ``base``, :cpp:class:`neml::NEMLModel_sd`, Base material model, No
   ``alpha``, :cpp:class:`neml::Interpolate`, Thermal expansion coefficient, ``0.0``
   ``tol``, :c:type:`double`, Solver tolerance, ``1.0e-8``
   ``miter``, :c:type:`int`, Maximum solver iterations, ``50``
   ``verbose``, :c:type:`bool`, Verbosity flag, ``false``

Class description
-----------------

.. doxygenclass:: neml::ClassicalCreepDamageModel_sd
   :members:
   :undoc-members:
