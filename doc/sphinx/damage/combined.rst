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

========== ====================================== ======================================= =======
Parameter  Object type                            Description                             Default
========== ====================================== ======================================= =======
elastic    LinearElasticModel                     Elasticity model                        No
models     std::vector<NEMLScalarDamagedModel_sd> List of damage models to apply          No
base       NEMLModel_sd                           Base material model                     No
alpha      Interpolate                            Thermal expansion coefficient           0.0
tol        double                                 Solver tolerance                        1.0e-8
miter      int                                    Maximum solver iterations               50
verbose    bool                                   Verbosity flag                          false
========== ====================================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::CombinedDamageModel_sd
   :members:
   :undoc-members:
