General rate dependent models
=============================

Overview
--------

This class is a generic interface for integrating rate-dependent models.
It solves the nonlinear equations described by

.. math::
   \bm{\sigma}_{n+1} = \bm{\sigma}_{n}+\dot{\bm{\sigma}}\left(\bm{\sigma}_{n+1},\bm{\alpha}_{n+1},\dot{\bm{\varepsilon}}_{n+1},T_{n+1},\dot{T}_{n+1},t_{n+1}\right)\Delta t_{n+1}

   \bm{\alpha}_{n+1} = \bm{\alpha}_{n}+\dot{\bm{\alpha}}\left(\bm{\sigma}_{n+1},\bm{\alpha}_{n+1},\dot{\bm{\varepsilon}}_{n+1},T_{n+1},\dot{T}_{n+1},t_{n+1}\right)\Delta t_{n+1}

In these equations :math:`\dot{\bm{\sigma}}` is some generic stress rate law
and :math:`\dot{\bm{\alpha}}` is some generic history evolution law.
These equations are defined by a GeneralFlowRule interface.
The only current purpose of this general integration routine is to integrate 
viscoplastic material models, but it could be used for other purposes in the
future.

The integrator uses fully implicit backward Euler integration for both the
stress and the history.
It returns the algorithmic tangent, computed using the implicit function 
theorem.
The work and energy are integrated with a trapezoid rule from the final values
of stress and inelastic strain.

This model maintains a vector of history variables defined by the
model's GeneralFlowRule interface.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Temperature dependent elastic constants, No
   ``surface``, :cpp:class:`neml::GeneralFlowRule`, Flow rule interface, No
   ``alpha``, :cpp:class:`neml::Interpolate`, Temperature dependent instantaneous CTE, ``0.0``
   ``tol``, :c:type:`double`, Integration tolerance, ``1.0e-8``
   ``miter``, :c:type:`int`, Maximum number of integration iters, ``50``
   ``verbose``, :c:type:`bool`, Print lots of convergence info, ``false``
   ``max_divide``, :c:type:`int`, Max adaptive integration divides, ``8``

Class description
-----------------

.. doxygenclass:: neml::GeneralIntegrator
   :members:
   :undoc-members:
