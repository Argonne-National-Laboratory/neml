Regime switching model
======================

Overview
--------

This model selects a stress and history update function from a input 
list of :doc:`NEMLModel_sd` objects based on a normalized activation energy.
This activation energy is a function of strain rate and temperature and
the form of the normalized energy comes from the work of Kocks and Mecking
[KM2003]_.

The input to the metamodel is a list of material models and a corresponding
list of activation energy cutoffs :math:`\left[g_1, g_2, \dots, g_n \right]`.
When the caller requests a stress and history update this metamodel 
first calculates the normalized activation energy

.. math::
   g = \frac{k T}{\mu b^3} \ln\frac{\dot{\varepsilon}_0}{\dot{\varepsilon}}.

Here :math:`k` is the Boltzmann constant, :math:`T` the current value of
temperature, :math:`\mu` the current, temperature-dependent shear modulus,
:math:`b` a Burgers vector, :math:`\dot{\varepsilon}_0` a reference strain 
rate, and :math:`\dot{\varepsilon}` the current equivalent strain rate,
computed as

.. math::
   \dot{\varepsilon} = \frac{\varepsilon_{n+1} - \varepsilon_{n}}{t_{n+1} - t_{n}}

with

.. math::
   \varepsilon_{n+1} = \sqrt{\frac{2}{3}\bm{\varepsilon}_{n+1}:\bm{\varepsilon}_{n+1}}

and similarly for state *n*.

If :math:`g<g_1` the metamodel selects the first model in the input list,
if :math:`g \ge g_n` the metamodel selects the last model in the input, and
for values in between the model selects model *i* such that 
:math:`g_{i-1} < g \ge g_{i}`.

The metamodel dispatches calls for the history evolution, algorithmic 
tangent, and energy likewise.

The Kocks-Mecking metamodel requires each model in the input list to use
compatible history variables.
That is, all the possible base model options must use the same history 
variables.
This means the metamodel maintains only one copy of the history variables,
which are common for all the base models.
The history evolution is therefore consistent no matter which base model
is selected for a particular stress update.

.. WARNING::
   The KMRegimeModel metamodel does not check for consitency of history 
   for the base models, beyond checking to make sure that each model 
   uses the same number of history variables.

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
elastic    LinearElasticModel        Temperature dependent elastic constants No
models     std::vector<NEMLModel_sd> Vector of base models                   No
gs         std::vector<double>       Corresponding vector of energies        No
kboltz     double                    Boltzmann's constant                    No
b          double                    Burger's vector length                  No
eps0       double                    Reference strain rate                   No
alpha      Interpolate               Temperature dependent instantaneous CTE 0.0
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::KMRegimeModel
   :members:
   :undoc-members:
