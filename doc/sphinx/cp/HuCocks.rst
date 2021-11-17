Hu and Cocks model for 316H SS
==============================

This subsystem implements a single crystal model for 316H stainless steel
developed by Hu, Cocks, and coworkers [HC2020]_.  The model includes a specific
flow rule as well as a complex hardening model designed to capture the key
features of dislocation forest hardening, precipitation hardening
caused by the carbide and Laves phases, and solid solution strengthening caused
by Mo, C, and Cr in the alloy.  For clarity, the implementation
breaks the complex hardening model into several smaller objects.

ArrheniusSlipRule
-----------------

Overview
""""""""

This class implements a thermally-activated slip rule, where the slip rate
is defined by

.. math::
   \dot{\gamma}_i = \dot{\gamma}_0 \exp \left[ -\frac{\Delta F_0}{kT} \left(1 - \left| \frac{\tau_i}{\tau_{CRSS,i}} \right|^A \right)^B \right] \operatorname{sign}\left(\tau_i \right) 

where :math:`\dot{\gamma}_i` is the slip rate on each system, :math:`\tau_i` is the resolved shear, :math:`\tau_{CRSS,i}` is the slip system strength,

.. math::
   \Delta F_0 = \alpha_0 G_0 b^3

and the remaining terms are parameters, described below.

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``resistance``, :cpp:class:`neml::SlipHardening`, Slip resistance, N
   ``A``, :code:`double`, Energy barrier shape :math:`A`, N
   ``B``, :code:`double`, Energy barrier shape :math:`B`, N
   ``b``, :code:`double`, Burgers vector :math:`b`, N
   ``g0``, :code:`double`, Reference shear rate :math:`\dot{\gamma}_0`, N
   ``G0``, :code:`double`, Activation energy :math:`G_0`, N
   ``k``, :code:`double`, Boltzmann constant :math:`k`, ``1.3806485e-23``

Class description
"""""""""""""""""

.. doxygenclass:: neml::ArrheniusSlipRule
   :members:
   :undoc-members:

Hardening model
---------------

The Hu and Cocks model combines dislocation forest hardening, precipitation
hardening, and solid solution strengthening.  The specific model
described in the original work considers two precipitation reactions:

1. C and Cr form Cr\ :sub:`23`\ C\ :sub:`6`\ carbides
2. Mo forms Fe\ :sub:`2`\ Mo

Precipitation strengthening and solid solution strengthening are coupled:
as the chemical species contributing the to strengthen phase come out of
solution that decreases the solid solution strengthening while increasing the
precipitation strengthen.  The precipitation reactions undergo two phases of
evolution: diffusion-controlled growth drawing the required chemical species
out of solution followed by Ostwald ripening where the larger precipitates 
cannibalize the smaller precipitates. 

The implementation in NEML treats these reactions are general -- the user
can specify as many chemical species feeding into as many precipitation
reactions as needed.  The only current restriction is that two precipitation 
reactions cannot compete for the same chemical species in solution.  This
behavior could be added to the implementation if needed.

Each precipitation-solid solution strengthening group can then be treated
separately.  Similarly, dislocation hardening can be treated separately
from the precipitation/solid solution strengthening.  The implementation
then uses three objects to form the slip system strengths:

1. A `DislocationSpacingHardening`_ object to manage dislocation hardening
2. An `HuCocksPrecipitationModel`_ object for each precipitation reactor
3. A `HuCocksHardening`_ object to sum the contributions of each mechanism on the slip system strength.

For full details of the implementation in NEML see [VM2021]_.

HuCocksHardening
----------------

Overview
""""""""
This object sums the contributions of the individual precipitation reactions
and the dislocation hardening into a single slip system strength.  
The equation it implements is 

.. math::
   \tau_{CRSS,i} = \sqrt{\tau_{d,i}^2 + \tau_p^2} + \tau_s

where :math:`\tau_{d,i}` is the dislocation hardening strength, provided by
the `DislocationSpacingHardening`_ model.

:math:`\tau_p` is the 
total precipitation hardening given by

.. math::
   \tau_p = \frac{\alpha_p G b}{L_p}

with :math:`\alpha_p` an interaction coefficient, :math:`G` the shear modulus,
:math:`b` the Burgers vector, and :math:`L_p` is 

.. math::
   L_p = \sqrt{\frac{1}{\sum_i 2 r_i N_i}}.

with :math:`r_i` and :math:`N_i` internal variables defined by the individual
`HuCocksPrecipitationModel`_ precipitation reaction models.

:math:`\tau_s` is then the total solid solution strengthening hardening given by

.. math::
   \tau_s = \frac{\alpha_s G b}{L_s}

with :math:`\alpha_s` an interaction coefficient and using

.. math::
    L_s = \sqrt{\frac{1}{b\sum_j \frac{c_j}{v_{m,i}}}}

with :math:`c_j` the chemical concentrations contributing to each precipitation reaction and :math:`v_{m,i}` the corresponding molecular volumes.  The chemical concentrations are again defined by the individual `HuCocksPrecipitationModel`_ models. 

Parameters
""""""""""
.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``dmodel``, :cpp:class:`neml::SlipHardening`, Dislocation hardening model, N
   ``pmodels``, :cpp:class:`std::vector<neml::HuCocksPrecipitationModel>`, Precipitation hardening models, N
   ``ap``, :code:`double`, Precipitation hardening interaction coefficient :math:`\alpha_p`, N
   ``ac``, :code:`double`, Solid solution interaction coefficient :math:`\alpha_c`, N
   ``b``, :code:`double`, Burgers vector :math:`b`, N
   ``G``, :cpp:class:`neml::Interpolate`, Shear modulus, N

Class description
"""""""""""""""""
.. doxygenclass:: neml::HuCocksHardening
   :members:
   :undoc-members:

HuCocksPrecipitationModel
-------------------------

Overview
""""""""
This class manages a `single` precipitation reaction.  The model considers the evolution of a single precipitate through
two phases:

1. Nucleation and growth by diffusion of the underlying chemical species out of the solid solution, followed by
2. Ostwald ripening after the solid solution concentrations reach their equilibrium values

During the first phase, the precipitates draw the chemical species out of the solution, affecting the solid solution
strengthening provided by those elements.

The model tracks a precipitation reaction with three internal variables: `math:`f` the volume fraction, :math:`N` the
number volume density, and  :math:`r` the average precipitate radius.  One of these three variables is redundant, but
the model evolves all three (consistency) to improve the numerical stability of the as a whole.  The chemical
concentrations in solution underlying the precipitation reaction can be determined given these three variables describing
the precipitates.

The model applies a different ODE to evolve the internal variables in each of the regimes.  The two regimes are split
by the chemical concentration of the critical species in solution.  The model is in the growth regime when

.. math::
   c_j < c_{eq,j}

for `all` species contributing to the reaction.  Conversely, the model is in the ripening regime when

.. math::
   c_j = c_{eq,j}

for that critical species.

Growth regime
"""""""""""""
In the growth regime the chemical concentrations evolve as

.. math::
   c_j = \frac{c_{0,j} - f c_{p,j}}{1-f}

where :math:`c_j` is the concentration in solution for species :math:`j`, :math:`c_{0,j}` is the initial solution
concentration for that species, and :math:`c_{p.j}` is the chemical concentration of the species in the precipitate.
The model stays in the growth regime until the first species contributing the precipitation reaction reaches the 
solution equilibrium concentration :math:`c_{eq,j}`.

In this regime

.. math::
   \dot{f}_{growth} = \frac{4}{3}\pi \left(\dot{N} r^3 + 3 N r^2 \dot{r} \right)

.. math::
   \dot{r}_{growth} = \frac{D}{r} \frac{c_j - c_{eq,j}}{c_{p,j} - c_{eq,j}} + \frac{\dot{N}_{growth}}{N} \left( r_c - r \right)

with :math:`G_v` the Gibb's free energy driving the reaction

.. math::
   G_v = -\frac{kT}{v_m} \ln \frac{c_{eff}}{c_{eff,eq}}

using

.. math::
    c_{eff} = \prod_{j} c_j.nit

Diffusion of the `slowest species` controls the reaction rate, with 

.. math::
   D = D_0 \exp\left( \frac{-Q_0}{RT} \right)

where :math:`D_0` is the diffusivity at absolute zero, :math:`Q_0` the activation energy, and :math:`R` the universal gas constant.

Finally, 

.. math::
    r_c = -2 \frac{\chi}{G_v}

with :math:`\chi` the interface energy.

For the nucleation rate:

.. math::
    \dot{N}_{growth} = N_0 Z \beta \exp\left(-\frac{G^*}{kT} \right)

with

.. math::
   G^* = \frac{16 \pi \chi^3}{3 G_v^2}

and 

.. math::
   Z \beta = \frac{2 v_m D c_j}{a_m^4} \sqrt{\frac{\chi}{kT}}

where :math:`a_m` is the relevant lattice parameter.

In the growth regime the nucleation rate is positive.

Ripening regime
"""""""""""""""

In the ripening regime the solution chemical concentrations are frozen and do not change.  The critical species, the
element which first reaches the equilibrium concentration, is frozen at that solution equilibrium concentration
:math:`c_{j,eq}` and the other species remained fixed at the final concentrations from the growth phase.

In this regime:

.. math::
   \dot{f}_{ripening} = 0

.. math::
   \dot{r}_{ripening} = \frac{M}{3r^2}

.. math::
   \dot{N}_{ripening} = -\frac{3N}{r} \dot{r}_{ripening}

with 

.. math::
   M = C_f \frac{8 \chi V_m D c_j}{9 R T}

where :math:`C_f` is a coarsening factor and :math:`V_m = N_a v_m` is the molar volume (:math:`N_a` Avagadro's number).

In the ripening regime the nucleation rate is negative, the radius growth rate is positive, and there is no
net growth in volume fraction.

Switching mechanisms
""""""""""""""""""""

A hard switch between the growth and ripening regimes produces an unstable numerical model.  Instead the NEML implementation
mixes the two rates using a sigmoid function:

.. math::
   \dot{r} = f(c_j) \dot{r}_{growth} + (1-f(c_j)) \dot{r}_{ripening}

and

.. math::
   \dot{N} = f(c_j) \dot{N}_{growth} + (1-f(c_j)) \dot{N}_{ripening}

where

.. math::
   f\left(c_{j}\right)=\begin{cases}
         \frac{c_{j}-c_{0,j}}{c_{eq,j}-c_{0,j}} & c_{j}\le c_{eq,j}\\
         1 & c_{j}>c_{eq,j}
   \end{cases}

With this setup the volume fraction evolution equation naturally trends towards zero in the ripening regime.

Additionally, the equations are scaled to equalize the magnitude of the internal variables, again to help with
numerical performance

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``c0``, :cpp:class:`std::vector<neml::Interpolate>`, Initial solution concentration of each species, N
   ``cp``, :cpp:class:`std::vector<neml::Interpolate>`, Precipitate concentration of each species, N
   ``ceq``, :cpp:class:`std::vector<neml::Interpolate>`, Equilibrium solution concentration of each species, N
   ``am``, :code:`double`, Lattice parameter, N
   ``N0``, :code:`double`, Nucleation site density :math:`N_0`, N
   ``Vm``, :code:`double`, Molar volume, N
   ``chi``, :code:`double`, Surface energy, N
   ``D0``, :code:`double`, Reference diffusivity, N
   ``Q0``, :code:`double`, Diffusion activation energy, N
   ``Cf``, :cpp:class:`neml::Interpolate`, Coarsening factor, N
   ``kboltz``, :code:`double`, Boltzmann constant, ``1.3806485e-23``
   ``R``, :code:`double`, Gas constant, ``8.31462``
   ``Na``, :code:`double`, Avagadro's number, ``6.02e23``
   ``rate``, :code:`size_t`, Index of rate-limiting chemical species, ``0``
   ``f_init``, :code:`double`, Initial volume fraction, ``4.18879e-16``
   ``r_init``, :code:`double`, Initial radius, ``1e-9``
   ``N_init``, :code:`double`, Initial number density, ``1e11``
   ``fs``, :code:`double`, Scaling factor on volume fraction, ``0.1``
   ``rs``, :code:`double`, Scaling factor on radius, ``1e-9``
   ``Ns``, :code:`double`, Scaling factor on number density, ``1e12``

Class description
"""""""""""""""""
.. doxygenclass:: neml::HuCocksPrecipitationModel
   :members:
   :undoc-members:


DislocationSpacingHardening
---------------------------

Overview
""""""""
This model can be used as part of the `HuCocksHardening`_ model or as a stand-alone model for forest dislocation
hardening.  The model maintains a single, scalar dislocation density for each slip system, parameterized
as a mean obstacle spacing.  The slip system strength is given in terms of these obstacle spacings as

.. math::
   \tau_{d,i} = \frac{\alpha_d G b}{L_{d,i}}

with :math:`\alpha_d` an interaction coefficient, :math:`G` the shear modulus, :math:`b` the Burgers vector, and
:math:`L_{d,i}` the dislocation obstacle spacing on system :math:`i`.  This obstacle spacing evolves as

.. math::
   \dot{L}_{d,i} = -L_{d,i}^3 \left( J_1 \left| \dot{\gamma}_i \right| + J_2 \sum_{j \ne i} \left| \dot{\gamma}_j \right|   \right) + \frac{K}{L_{d,i}^3}

with :math:`J_1` the self hardening coefficient, :math:`J_2` the latent hardening coefficient, and :math:`K` a temperature dependent parameter describing dislocation recovery.

Parameters
""""""""""
.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``J1``, :cpp:class:`neml::Interpolate`, Self hardening coefficient, N
   ``J2``, :cpp:class:`neml::Interpolate`, Latent hardening coefficient, N
   ``K``, :cpp:class:`neml::Interpolate`, Recovery coefficient, N
   ``L0``, :code:`double`, Initial obstacle spacing, N
   ``a``, :code:`double`, Interaction coefficient :math:`\alpha_d`, N
   ``b``, :code:`double`, Burgers vector, N
   ``G``, :cpp:class:`neml::Interpolate`, Shear modulus, N
   ``L``, :cpp:class:`neml::lattice`, Lattice to extract number of systems, N
   ``varprefix``, :code:`std::string`, Prefix of internal variables, :code:`"spacing"`

Class description
"""""""""""""""""
.. doxygenclass:: neml::DislocationSpacingHardening
   :members:
   :undoc-members:
