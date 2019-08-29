Modular creep damage
====================

Overview
--------

This object implements the classical Hayhurst-Leckie-Rabotnov-Kachanov creep damage model [HL1977]_.  Specifically, it exactly replicates Eq. 2.6 in that paper, so that the implementation exactly replicates the analytic expressions in Eq. 2.5.
The damage update is given by 

.. math::
   \omega_{n+1} = \omega_{n} + \left(\frac{\sigma_{eff}}{A}\right)^\xi 
      \left(1 - \omega_{n+1}\right)^{\xi-\phi} \Delta t_{n+1}

where :math:`\sigma_{eff}` is modular, defined by a :ref:`effective-stress` object.

This class has the option for element extinction, useful in FEA simulations of damage.  If the ``ekill`` option is set to true once the material point reaches a damage threshold of ``dkill`` the constitutive response will be replaced by a linear elastic response with an elastic stiffness of :math:`\mathbf{\mathfrak{C}}/f` where the factor :math:`f` is given by the parameter ``sfact``.

.. toctree::
   effective/effective

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``A``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``xi``, :cpp:class:`neml::Interpolate`, Stress exponent, No
   ``phi``, :cpp:class:`neml::Interpolate`, Damage exponent, No
   ``estress``,:cpp:class:`neml::EffectiveStress`, Effective stress, No
   ``base``, :cpp:class:`neml::NEMLModel_sd`, Base material model, No
   ``alpha``, :cpp:class:`neml::Interpolate`, Thermal expansion coefficient, ``0.0``
   ``tol``, :c:type:`double`, Solver tolerance, ``1.0e-8``
   ``miter``, :c:type:`int`, Maximum solver iterations, ``50``
   ``verbose``, :c:type:`bool`, Verbosity flag, ``false``
   ``ekill``, :c:type:`bool`, Trigger element death, ``false``
   ``dkill``, :c:type:`double`, Critical damage threshold, ``0.5``
   ``sfact``, :c:type:`double`, Stiffness factor for dead element, ``100000``

Class description
-----------------

.. doxygenclass:: neml::ModularCreepDamageModel_sd
   :members:
   :undoc-members:
