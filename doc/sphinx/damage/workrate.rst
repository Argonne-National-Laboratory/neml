Workrate dependent damage model
===============================

Overview
--------

This object implements a "plastic work rate" based damage model where the
damage rate is a function of inelastic work rate.
The object also has a second classical damage term that can account for hold
time effects during cyclic loading.
This damage law, thus, captures the creep-fatigue interaction.
The damage rate is

.. math::

   \dot{\omega} = \frac{\left(\dot{W}/Q\right)^m}{f(\dot{W})} + G \frac{\left(\sigma_{eff}/H \right)^\xi}{\left(1 - \omega \right)^\phi}

   \sigma_{eff} = \sqrt{\frac{3}{2} \operatorname{dev}\left(\bm{\sigma}\right):\operatorname{dev}\left(\bm{\sigma}\right)}

   f(\dot{W})   = \left(P \dot{W} + A \right)^n

   \dot{W}      = \frac{\sigma_{eff} \dot{\varepsilon_{eff}^{in}}}{1 - \omega}.

The functional dependence of damage on work rate :math:`f(\dot{W})` is selected
following the ductility exhaustion methodology, and obtaining a critical work
to failure.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``workrate``, :cpp:class:`neml::Interpolate`, work rate function, No
   ``P``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``A``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``n``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``Q``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``m``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``G``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``H``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``xi``, :cpp:class:`neml::Interpolate`, Stress exponent, No
   ``phi``, :cpp:class:`neml::Interpolate`, Damage exponent, No
   ``base``, :cpp:class:`neml::NEMLModel_sd`, Base material model, No
   ``alpha``, :cpp:class:`neml::Interpolate`, Thermal expansion coefficient, ``0.0``
   ``tol``, :c:type:`double`, Solver tolerance, ``1.0e-8``
   ``miter``, :c:type:`int`, Maximum solver iterations, ``50``
   ``verbose``, :c:type:`bool`, Verbosity flag, ``false``

Class description
-----------------

.. doxygenclass:: neml::WorkRateFunctionDamage_sd
   :members:
   :undoc-members:
