Creep + plasticity
==================

Overview
--------

A SmallStrainCreepPlasticity model combines any of the previous types of
material models with a rate dependent creep model.
The model solves the simultaneous nonlinear equations

.. math::
   \bm{\sigma}_{n+1}^{base}\left( 
   \bm{\varepsilon}_{n+1}^{base}, \bm{\varepsilon}_{n}^{base},
   T_{n+1}, T_{n},
   t_{n+1}, t_{n},
   \bm{\sigma}_{n}^{base},
   \bm{\alpha}_{n}^{base}
   \right) = 
   \bm{\sigma}_{n+1}^{cr}

   \bm{\varepsilon}_{n+1} = \bm{\varepsilon}_{n+1}^{base} + 
      \bm{\varepsilon}_{n+1}^{cr}\left(\bm{\varepsilon}_{n}^{cr},\bm{\sigma}_{n+1}^{cr},\Delta t_{n+1},T_{n+1}\right).

Here the model with superscript *base* is the base material model and 
the model with superscript *cr* is the creep model (a CreepModel object).
The base model can be any :doc:`NEMLModel_sd` object.
The creep model add no additional history variables, it is solely a function
of stress, temperature, and strain.

Unlike base material models creep models in NEML are configured to return
strain as a function of stress, rather than stress as a function of strain.
Because of this
these equations can be combined into a single nonlinear equation

.. math::
   \bm{\varepsilon}_{n+1} = \bm{\varepsilon}_{n+1}^{base}+\bm{\varepsilon}_{n+1}^{cr}\left(\bm{\varepsilon}_{n}-\bm{\varepsilon}_{n}^{base},\bm{\sigma}_{n+1}\left(\bm{\varepsilon}_{n+1}^{base},\bm{\varepsilon}_{n}^{ep},\boldsymbol{h}_{n},\Delta t_{n+1},T_{n+1}\right),\Delta t_{n+1},T_{n+1}\right)

The implementation solves this nonlinear equation and provides the appropriate
Jacobian using a matrix decomposition formula.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Temperature dependent elastic constants, No
   ``plastic``, :cpp:class:`neml::NEMLModel_sd`, Base material model, No
   ``creep``, :cpp:class:`neml::CreepModel`, Rate dependent creep model, No
   ``alpha``, :cpp:class:`neml::Interpolate`, Temperature dependent instantaneous CTE, ``0.0``
   ``tol``, :c:type:`double`, Integration tolerance, ``1.0e-8``
   ``miter``, :c:type:`int`, Maximum number of integration iters, ``50``
   ``verbose``, :c:type:`bool`, Print lots of convergence info, ``false``
   ``sf``, :c:type:`double`, Scale factor on strain equation, ``1.0e6``

.. NOTE::
   The scale factor is multiplied by a strain residual equation that may involve
   very small values of strain.
   The default value works well for values of nominal strain (i.e. in/in or mm/mm).

Class description
-----------------

.. doxygenclass:: neml::SmallStrainCreepPlasticity
   :members:
   :undoc-members:
