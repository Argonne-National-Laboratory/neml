Mukherjee creep
===============

Overview
--------

This implements the Mukherjee creep model [BMD1969]_

.. math::
   \dot{\varepsilon}^{cr} = A D_0 e^\frac{Q}{RT} \frac{\mu b}{k T} \left( \frac{\sigma_{eq}}{\mu}\right)^n

with :math:`A`, :math:`D_0`, :math:`Q`, and :math:`n` parameters and 
:math:`R` the gas constant, :math:`T` absolute temperature, :math:`\mu`
the temperature dependent shear modulus, :math:`b` a Burgers vector length, and
:math:`k` the Boltzmann constant.

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
emodel     LinearElasticModel        Elasticity model (for shear modulus)    No
A          double                    Prefactor                               No
n          double                    Stress exponent                         No
D0         double                    Zero temperature lattice diffusivity    No
Q          double                    Activation energy for diffusivity       No
b          double                    Burgers vector                          No
k          double                    Boltzmann constant                      No
R          double                    Gas constant                            No
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::MukherjeeCreep
   :members:
   :undoc-members:

