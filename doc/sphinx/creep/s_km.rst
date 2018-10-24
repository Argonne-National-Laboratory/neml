Regime switching Kocks-Mecking creep
====================================

Overview
--------

This model implements a creep law based on the Kocks-Mecking normalized
activation energy [KM2003]_.
The basic creep rate law is

.. math::
   F = -\frac{\mu b^3}{kT} 

   \dot{\varepsilon}^{cr} = \dot{\varepsilon}_0 \exp\left( B F \right) \left(\frac{\sigma_{eq}}{\mu} \right)^{A F}.

Here :math:`\dot{\varepsilon}_0`, :math:`A`, :math:`B` are parameters, 
:math:`\mu` is the temperature dependent shear modulus, :math:`k` is the Boltzmann constant, 
and :math:`b` is a Burgers vector length.
The model constants can be fit using a Kocks-Mecking diagram.

The formulation can switch between different Kocks-Mecking models as a
function of normalized stress :math:`s = \frac{\sigma_{eq}}{\mu}`.
The model first computes the normalized stress.
It then consults a table of :math:`n` normalized stress cutoffs, :math:`c_i`.
If :math:`s \le c_1` then the model applies the strain rate equation
above using constants :math:`A_1` and :math:`B_1`.  
If :math:`c_i < s \le c_{i+1}` then the model uses constants :math:`A_{i+1}` and :math:`B_{i+1}`. 
If :math:`s > c_n` then the model uses :math:`A_n` and :math:`B_n`.

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
cuts       std::vector<double>       Normalized stress cutoffs               No
A          std::vector<double>       Corresponding A constants               No
B          std::vector<double>       Corresponding B constants               No
kboltz     double                    Boltzmann constant                      No
b          double                    Burgers vector                          No
eps0       double                    Reference strain rate                   No
emodel     LinearElasticModel        Elastic model (for shear modulus)       No
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::RegionKMCreep
   :members:
   :undoc-members:
