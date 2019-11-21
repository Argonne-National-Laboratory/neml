Larson Miller Damage
====================

Overview
--------

This model implements the damage model

.. math::
   \omega_{n+1} = \omega{n} + \frac{\Delta t_{n+1}}{t_R\left(\sigma_{eff} \left(1-\omega_{n+1}\right) \right)}

where :math:`\sigma_{eff}` is a modular effective stress, defined by a :ref:`effective-stress` object, and :math:`t_R` is a time-to-rupture Larson-Miller relation provided
by a :ref:`larson-miller` object.

This model implements classical time-fraction damage summation in a finite 
element framework.  It currently does not include graceful element deletion
routines, though those may be added in the future.  This means the simulation
can only progress up to the point where the simulation reaches 
:math:`\omega = 1` at the first material point in the calculation.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``lmr``, :cpp:class:`neml::LarsonMillerCorrelation`, Parameter, No
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

.. doxygenclass:: neml::LarsonMillerCreepDamageModel_sd
   :members:
   :undoc-members:
