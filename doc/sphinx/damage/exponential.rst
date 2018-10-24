Exponential damage
==================

Overview
--------

This object implements a "standard" damage model proportional to the dissipated
inelastic energy.
The damage function is

.. math::

   w = \frac{\left(\omega + k_0\right)^{a_f}}{W_0} \sigma_{eff}

   \sigma_{eff} = \sqrt{\frac{3}{2} \operatorname{dev}\left(\bm{\sigma}\right):\operatorname{dev}\left(\bm{\sigma}\right)}.

The standard damage model multiplies this function by the inelastic
strain rate in computing the damage update.
Because the dissipation rate is equal to :math:`\sigma_{eff} \dot{\varepsilon}_{eff}^{in}` this model actually increases damage in proportion to the dissipation.

Parameters
----------

========== ====================================== ======================================= =======
Parameter  Object type                            Description                             Default
========== ====================================== ======================================= =======
elastic    LinearElasticModel                     Elasticity model                        No
W0         Interpolate                            Parameter                               No
k0         Interpolate                            Parameter                               No
af         Interpolate                            Parameter                               No
base       NEMLModel_sd                           Base material model                     No
alpha      Interpolate                            Thermal expansion coefficient           0.0
tol        double                                 Solver tolerance                        1.0e-8
miter      int                                    Maximum solver iterations               50
verbose    bool                                   Verbosity flag                          false
========== ====================================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::NEMLExponentialWorkDamagedModel_sd
   :members:
   :undoc-members:
