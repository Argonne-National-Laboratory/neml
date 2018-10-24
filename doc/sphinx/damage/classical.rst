Classical creep damage
======================

Overview
--------

This object implements the classical Hayhurst-Leckie-Rabotnov-Kachanov creep damage model [HL1977]_.
The damage update is given by 

.. math::
   \omega_{n+1} = \omega_{n} + \left(\frac{\sigma_{eff}}{A}\right)^\xi 
      \left(1 - \omega_{n+1}\right)^-\phi \Delta t_{n+1}

   \sigma_{eff} = \sqrt{\frac{3}{2} \operatorname{dev}\left(\bm{\sigma}\right):\operatorname{dev}\left(\bm{\sigma}\right)}.

Parameters
----------

========== ====================================== ======================================= =======
Parameter  Object type                            Description                             Default
========== ====================================== ======================================= =======
elastic    LinearElasticModel                     Elasticity model                        No
A          Interpolate                            Parameter                               No
xi         Interpolate                            Stress exponent                         No
phi        Interpolate                            Damage exponent                         No
base       NEMLModel_sd                           Base material model                     No
alpha      Interpolate                            Thermal expansion coefficient           0.0
tol        double                                 Solver tolerance                        1.0e-8
miter      int                                    Maximum solver iterations               50
verbose    bool                                   Verbosity flag                          false
========== ====================================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::ClassicalCreepDamageModel_sd
   :members:
   :undoc-members:
