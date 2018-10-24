Solver interface
================

A common task in NEML is solving a nonlinear system of equations.
For example, because NEML uses implicit time integration it must solve a
nonlinear system each time it integrates a material model defined with a 
rate form.
NEML provides an opaque mechanism for using different nonlinear equation
solvers throughout the code.
It does this by having objects that need to solve nonlinear equations
solvable by inheriting from :cpp:class:`neml::Solvable`.

.. doxygenclass:: neml::Solvable
   :members:

This class requires an object implement three virtual methods:

   1. ``nparams``: Return the number of variables in the nonlinear system to be solved.
   2. ``init_x``: Given a vector of length ``nparams`` and a :cpp:class:`neml::TrialState` object setup an initial guess to start the nonlinear solution iterations.
   3. ``RJ``: Given the current guess at the solution ``x`` (length ``nparams``) and the :cpp:class:`neml::TrialState` object return the residual equations (``R``, length ``nparams``) and the Jacobian of the residual equations with respect to the variables (``J``, ``nparams`` :math:`\times` ``nparams``).

A :cpp:class:`neml::TrialState` is a completely generic object that contains any information
beyond the current guess at the solution the class will need to construct
an initial guess and to calculate the residual and the Jacobian.
Essentially, it contains any variables that are held fixed during the
nonlinear solution process.
While :cpp:class:`neml::TrialState` objects all inherit from a base class, this is currently 
entirely cosmetic.

.. doxygenclass:: neml::TrialState
   :members:

Essentially they are C++ structs holding field data.
Below is an example TrialState for the :cpp:class:`neml::SmallStrainPerfectPlasticity` object.

.. code-block:: c++
   
   /// Small strain perfect plasticity trial state
   //  Store data the solver needs and can be passed into solution interface
   class SSPPTrialState : public TrialState {
    public:
     double ys, T;    // Yield stress, temperature
     double ee_n[6];  // Previous value of elastic strain
     double s_n[6];   // Previous stress
     double s_tr[6];  // Elastic trial stress
     double e_np1[6]; // Next strain
     double e_n[6];   // Previous strain
     double S[36];    // Compliance
     double C[36];    // Stiffness
   };

All this is a ``struct`` containing the information required to setup the nonlinear residual
equations.
Note this object does not contain the current value of stress or history.
This information is contained (and updated) in the solution vector ``x``.

Currently, NEML has two options for solving nonlinear equations.
NEML contains a built-in implementation of the Newton-Raphson method or NEML 
can use the `NOX <https://trilinos.org/packages/nox-and-loca/>` solver contained in
the `Trilinos <https://trilinos.org/>` package, developed by Sandia National Laboratories.
The solver is configured at build time, using the CMake configuration.
