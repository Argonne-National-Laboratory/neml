Model damaged by a scalar variable
==================================

Overview
--------

This object implements a damage model that uses a single damage variable 
to degrade the stress of a base material model.
It implements the stress update function

.. math::
   \bm{\sigma}_{n+1}^\prime = (1 - \omega_{n+1}) 
      \bm{\sigma}\left( 
      \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n},
      T_{n+1}, T_{n},
      t_{n+1}, t_{n},
      \bm{\sigma}_{n} / (1 - \omega_{n}),
      \bm{\alpha}_{n+1},
      \bm{\alpha}_{n},
      u_n, p_n
      \right)

where :math:`\omega` is the damage variable and :math:`\bm{\sigma}` is the 
base material stress update.
A :cpp:class:`neml::ScalarDamage` model provides the definition of :math:`\omega` as well
as the associated derivatives to form the Jacobian.

The damage model maintains the set of history variables from the base 
material plus one additional history variable for the damage.

This class has the option for element extinction, useful in FEA simulations of damage.  If the ``ekill`` option is set to true once the material point reaches a damage threshold of ``dkill`` the constitutive response will be replaced by a linear elastic response with an elastic stiffness of :math:`\mathbf{\mathfrak{C}}/f` where the factor :math:`f` is given by the parameter ``sfact``.

.. note::
   The scalar damage model passes in the modified stress :math:`\bm{\sigma} / (1 - \omega)` to the base stress update model in addition to modifying the stress update formula as shown in the above equation.

.. warning::
   The model also passes the modified stress :math:`\bm{\sigma} / (1-\omega)` to the damage update equation.  That is, the stress passed into these functions is the modified effective stress, not the actual external stress.  This means that the damage equations implemented in NEML vary slightly from the correpsonding literature equations working with the unmodified stress directly.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``base``, :cpp:class:`neml::NEMLModel_sd`, Base material model, No
   ``damage``, :cpp:class:`neml::ScalarDamage`, Damage model, No
   ``alpha``, :cpp:class:`neml::Interpolate`, Thermal expansion coefficient, ``0.0``
   ``tol``, :code:`double`, Solver tolerance, ``1.0e-8``
   ``miter``, :code:`int`, Maximum solver iterations, ``50``
   ``verbose``, :code:`bool`, Verbosity flag, ``false``
   ``ekill``, :code:`bool`, Trigger element death, ``false``
   ``dkill``, :code:`double`, Critical damage threshold, ``0.5``
   ``sfact``, :code:`double`, Stiffness factor for dead element, ``100000``

Scalar damage models
--------------------

.. toctree::
   scalardamage

Class description
-----------------

.. doxygenclass:: neml::NEMLScalarDamagedModel_sd
   :members:
   :undoc-members:

