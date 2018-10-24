Power law damage
================

Overview
--------

This object implements a "standard" damage model proportional to a power law in stress and directly
to the effective inelastic strain.
The damage function is

.. math::

   w = A \sigma_{eff}^n

   \sigma_{eff} = \sqrt{\frac{3}{2} \operatorname{dev}\left(\bm{\sigma}\right):\operatorname{dev}\left(\bm{\sigma}\right)}.

The standard damage model multiplies this function by the inelastic
strain rate in computing the damage update.

Parameters
----------

========== ====================================== ======================================= =======
Parameter  Object type                            Description                             Default
========== ====================================== ======================================= =======
elastic    LinearElasticModel                     Elasticity model                        No
A          Interpolate                            Prefactor                               No
a          Interpolate                            Stress exponent                         No
base       NEMLModel_sd                           Base material model                     No
alpha      Interpolate                            Thermal expansion coefficient           0.0
tol        double                                 Solver tolerance                        1.0e-8
miter      int                                    Maximum solver iterations               50
verbose    bool                                   Verbosity flag                          false
========== ====================================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::NEMLPowerLawDamagedModel_sd
   :members:
   :undoc-members:
