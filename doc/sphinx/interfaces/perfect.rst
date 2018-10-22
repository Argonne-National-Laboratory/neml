Perfect plasticity
==================

Overview
--------

This class implements rate independent perfect plasticity described by:

   The elastic trial state:

   .. math::

      \bm{\varepsilon}^{p}_{tr} = \bm{\varepsilon}^{p}_n

      \bm{\sigma}_{tr} = \mathbf{\mathfrak{C}}_{n+1} : 
         \left( \bm{\varepsilon}_{n+1} - \bm{\varepsilon}_{tr}^p  \right)

   The plastic correction:

   .. math::
      \bm{\sigma}_{n+1} = \mathbf{\mathfrak{C}}_{n+1} : 
         \left( \bm{\varepsilon}_{n+1} - \bm{\varepsilon}_{n+1}^p \right)

      \bm{\varepsilon}_{n+1}^p = 
         \begin{cases}
            \bm{\varepsilon}^{p}_{tr} & f\left(\bm{\sigma}_{tr}\right)\le0\\
            \bm{\varepsilon}^{p}_{tr}+\frac{\partial f_{n+1}}{\partial\bm{\sigma}_{n+1}}\Delta\gamma_{n+1} & f\left(\bm{\sigma}_{tr}\right)>0
         \end{cases}

   Solving for :math:`\Delta \gamma_{n+1}` such that

   .. math::
      f\left(\bm{\sigma}_{n+1} \right) = 0

In these equations :math:`f` is a yield function, parameterized by the yield
stress :math:`\sigma_0`.

If the step is plastic the stress update is solved through fully-implicit 
backward Euler integration.
The algorithmic tangent is then computed using an implicit function scheme.
The work and energy are integrated with a trapezoid rule from the final values
of stress and plastic strain.

This model does not maintain any history variables.

Parameters
----------

========== ===================== ======================================= =======
Parameter  Object type           Description                             Default
========== ===================== ======================================= =======
elastic    LinearElasticModel    Temperature dependent elastic constants No
surface    YieldSurface          The yield surface                       No
ys         Interpolate           The yield stress as a function of T     No
alpha      Interpolate           Temperature dependent instantaneous CTE 0.0
tol        double                Integration tolerance                   1.0e-8
miter      int                   Maximum number of integration iters     50
verbose    bool                  Print lots of convergence info          F
max_divide int                   Maximum number of adaptive subdivisions 8
========== ===================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::SmallStrainPerfectPlasticity
   :members:
   :undoc-members:
