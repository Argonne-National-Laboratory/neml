.. _integration:

Integrating models
==================

A generic material model has two jobs:

1. Integrate the stress and some set of internal variables from time step :math:`n` to time step :math:`n+1`, given the strain increment (and, for some models, the time increment and temperature increment) as input.
2. Provide the calling routine (often the stress update at a Gauss point in a finite element code) the algorithmic tangent -- the derivative of the stress at time :math:`n+1` with respect to the strain increment.

In alternate formulations the model might be stress-based -- the calling routine provides the stress and expects the strain -- or in multiphysics simulations the calling routine might require a similar algorithmic tangent with respect to the temperature or some other variable.  The following derivation can be adjusted to accommodate either condition.

NEML uses a single algorithm to accomplish these tasks for the:

- :doc:`Perfect plasticity <interfaces/perfect>`
- :doc:`Rate indepement plasticity <interfaces/rate_independent>`
- :doc:`General rate dependent <interfaces/general_integrator>`

models.  The same algorithm could be used, pointlessly, for the :doc:`linear elasticity <interfaces/lelastic>` update.  These classes encompass all flavors of small strain constitutive models commonly used in structural simulations.

The remainder of the models, including the :ref:`damage <damage>` models, :doc:`creep + plasticity <interfaces/creep_plasticity>`, and :doc:`regime switching <interfaces/km_regime>` models use specialized integration algorithms.  The common feature for these model is that they take as input one of the base material models listed above and modify that basic stress update formulation.

The standard material model update formulation in NEML can be described as the implicit time integration of the stress rate and a set of generic internal varaibles:

.. math::

   \dot{\bm{x}}=\left[\begin{array}{c}
   \dot{\bm{\sigma}}\left(\bm{\sigma},\bm{h},t,\bm{\varepsilon},T\right)\\
   \dot{\bm{h}}\left(\bm{\sigma},\bm{h},t,\bm{\varepsilon},T\right)
   \end{array}\right]

As suggested by this equation, the implementation divides the variables in the stress and update functions into two sets: the implicit variables :math:`\bm{x}` (i.e. the stress and generic set of history variables) and the explicit variables :math:`\bm{y}` (i.e. the strain if we are only interested in the classical algorithmic tangent).  For a backward Euler implementation of the rate equations we can formulate the generic residual equation:

.. math::
   \bm{R} = \bm{x}_{n+1} - \bm{x}_n - \dot{\bm{x}} \left(\bm{x} ; \bm{y} \right) \Delta t

i.e. the material point integration only updates the implicit variables and leaves the explicit variables fixed.  This generic equation can be solved with Newton's method.  That requires the Jacobian:

.. math::
   \bm{J} = \bm{I} - \frac{\partial \dot{\bm{x}}}{\partial \bm{x}}.

The stress update algorithms for the all of the standard models can be cast into this form.  Each model must then implement the correct residual and jacobian.  The remainder of the integration can use a common algorithm.

The complete, generalized algorithmic tangent is the total derivative of the implicit variables :math:`\bm{x}` with respect to the explicit variables :math:`\bm{y}`.  The classical algorithmic tangent, required by FE solvers, is simply the :math:`\frac{d \bm{\sigma}}{d \bm{\varepsilon}}` block of this generic derivative.  Note that multiphysics codes might be interested in other blocks of this "generalized" tangent matrix, for example the :math:`\frac{d \bm{\sigma}}{d T}` block.

Calculating the generalized tangent matrix is a straightforward application of the implicit function theorem.  After successfully solving the residual equation we have the condition:

.. math::
   d \bm{R} = \bm{0} = \frac{\partial \bm{R}}{\partial \bm{x}} : d \bm{x} + \frac{\partial \bm{R}}{\partial \bm{y}} : d \bm{y}

which implies that

.. math::
   \frac{d \bm{x}}{d \bm{y}} = - \left( \frac{\partial \bm{R}}{\partial \bm{x}} \right)^{-1} : \frac{\partial \bm{R}}{\partial \bm{y}}

The partial derivative :math:`\bm{J} = \frac{\partial \bm{R}}{\partial \bm{x}}` is the jacobian of the local constitutive update.  The other partial of the residual with respect to the explicit variables is not required during the local stress update and must be determined only for the tangent calculation.  The final expression for the generalized algorithmic tangent is then:

.. math::
   \bm{T} = - \bm{J}^{-1} \cdot \bm{E}

with :math:`\bm{T}` the generalized tangent and

.. math::
   \bm{E} = \frac{\partial \bm{R}}{\partial \bm{y}}

Again we are only typically interested in certain blocks of :math:`\bm{T}`.

However, calculating the whole generalized tangent is useful when doing adaptive substepping.  A typical scheme might subdivide the input values of the explicit variables :math:`\bm{y}` into several smaller steps called substeps.  Here we are typically interested in the tangent matrix over the whole step and not the individual tangents for each substep.  As demonstrated by [PRH2001]_ there is a recursive formula for calculating the complete adaptive tangent.  Consider the substepping scheme defined by:

.. math::
   \bm{y}^{i+1} = \bm{y}^{i} + \alpha_{i+1} \left(\bm{y}_{n+1} - \bm{y}_n\right)

The recursive formula for the tangent :math:`\bm{T}^{i+1}`  covering the complete step from :math:`\bm{y}_n` to the current substep is

.. math::
   \bm{T}^{i+1} = \bm{J}^{-1}_{i+1} \cdot \left(\alpha_{i+1} \bm{E}_{i+1} + \bm{T}^{i} \right)

where the partial derivatives :math:`\bm{J}_{i+1}` and :math:`\bm{E}_{i+1}` are for the current substep.  Applying this recursion relation though each substep produces the consistent tangent for the whole step.

Note this algorithm depends on propagating the whole generalized consistent tangent, not just the derivative of the stress with respect to the strain.  This is because the history variables also evolve throughout the substepping.  However, as described again in [PRH2001]_ some optimizations are possible.  Only minor `columns` of :math:`\bm{T}` pertaining to the strain :math:`\bm{\varepsilon}` are required for standard FE codes and so the recursive relation can be restricted to the approach subblocks of the generalized tangent.  Additionally, some types of internal variables, notably the plastic multiplier for rate independent plasticity models, do not propagate from substep to substep but instead reset with each substep.  The minor `rows` for these sorts of internal variables can be omitted from the recursive propagation.  Currently NEML does not make either optimization.
