Rate independent plasticity
===========================

Overview
--------

This class implements rate independent plasticity described by:

   The elastic trial state:

   .. math::

      \bm{\varepsilon}^{p}_{tr} = \bm{\varepsilon}^{p}_n

      \bm{\sigma}_{tr} = \mathbf{\mathfrak{C}}_{n+1} : 
         \left( \bm{\varepsilon}_{n+1} - \bm{\varepsilon}_{tr}^p  \right)

      \bm{\alpha}_{tr} = \bm{\alpha}_{n}

   The plastic correction:

   .. math::
      \bm{\sigma}_{n+1} = \mathbf{\mathfrak{C}}_{n+1} : 
         \left( \bm{\varepsilon}_{n+1} - \bm{\varepsilon}_{n+1}^p \right)

      \bm{\varepsilon}_{n+1}^p = 
         \begin{cases}
            \bm{\varepsilon}^{p}_{tr} & f\left(\bm{\sigma}_{tr},\bm{\alpha}_tr\right)\le0\\
            \bm{\varepsilon}^{p}_{tr}+\mathbf{g}_{n+1}\left( \bm{\sigma}_{n+1}, \bm{\alpha}_{n+1}, T_{n+1} \right)\Delta\gamma_{n+1} & f\left(\bm{\sigma}_{tr},\bm{\alpha}_tr\right)>0
         \end{cases}

      \bm{\alpha}_{n+1} = 
         \begin{cases}
            \bm{\alpha}_{tr} & f\left(\bm{\sigma}_{tr},\bm{\alpha}_tr\right)\le0\\
            \bm{\alpha}_{tr}+\mathbf{h}_{n+1}\left( \bm{\sigma}_{n+1}, \bm{\alpha}_{n+1}, T_{n+1} \right)\Delta\gamma_{n+1} & f\left(\bm{\sigma}_{tr},\bm{\alpha}_tr\right)>0
         \end{cases}

   Solving for :math:`\Delta \gamma_{n+1}` such that

   .. math::
      f\left(\bm{\sigma}_{n+1}, \bm{\alpha}_{n+1} \right) = 0

In these equations :math:`f` is a yield function, :math:`\mathbf{g}_{n+1}` is
a flow function, evaluated at the next state, and :math:`\mathbf{h}_{n+1}` is 
the rate of evolution for the history variables, evaluated at the next
state.
NEML integrates all three of these functions into a RateIndependentFlowRule
interface.

If the step is plastic the stress update is solved through fully-implicit 
backward Euler integration.
The algorithmic tangent is then computed using an implicit function scheme.
The work and energy are integrated with a trapezoid rule from the final values
of stress and plastic strain.

This model maintains a vector of history variables defined by the
model's RateIndependentFlowRule interface.

At the end of the step the model (optionally) checks to ensure the step
met the Kuhn-Tucker conditions

.. math::
   \Delta \gamma_{n+1} \ge 0
   f\left(\bm{\sigma}_{n+1}, \bm{\alpha}_{n+1} \right) \le 0
   \Delta \gamma_{n+1} f\left(\bm{\sigma}_{n+1}, \bm{\alpha}_{n+1} \right) = 0. 

Parameters
----------

========== ======================= ======================================= =======
Parameter  Object type             Description                             Default
========== ======================= ======================================= =======
elastic    LinearElasticModel      Temperature dependent elastic constants No
surface    RateIndependentFlowRule Flow rule interface                     No
alpha      Interpolate             Temperature dependent instantaneous CTE 0.0
tol        double                  Integration tolerance                   1.0e-8
miter      int                     Maximum number of integration iters     50
verbose    bool                    Print lots of convergence info          F
kttol      double                  Tolerance on the Kuhn-Tucker conditions 1.0e-2
check_kt   bool                    Flag to actually check KT               F
========== ======================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::SmallStrainRateIndependentPlasticity
   :members:
   :undoc-members:
