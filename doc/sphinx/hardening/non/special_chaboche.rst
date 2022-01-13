ChabocheVoceRecovery: a special variant of the Chaboche model
=============================================================

Overview
--------

This model describes a variant of the :doc:`Chaboche kinematic hardening <chaboche>` with 
:doc:`Voce <../simple/iso_voce>` isotropic hardening.  The only differences compared to using the standard
:doc:`Chaboche <chaboche>` model with Voce isotropic hardening are:

1. The Voce model is reparameterized as described below
2. The Voce model includes static recovery
3. There is a subtle difference in the definition of the static recovery on the Chaboche backstress terms.

Voce reparameterization and recovery
""""""""""""""""""""""""""""""""""""

The Voce model for this hardening option is reparameterized so that the strength is given as 

.. math::
   Q = -\sigma_0 - R

.. WARNING::
   All of the NEML yield surfaces assume the opposite of the standard
   sign convention for isotropic and kinematic hardening.
   The hardening model is expected to return a negative value of the
   isotropic hardening stress and a negative value of the backstress.

.. math::
   \dot{R} = \sqrt{\frac{2}{3}} \theta_0 \left(1 - \frac{R}{R_{max}} \right) \dot{p} + r_1 \left| R_{min} - R \right|^{r_2} \operatorname{sign} \left( R_{min} - R \right)

where :math:`\dot{p}` is the scalar plastic strain rate and the remaining undefined terms are parameters.  

The first term, proportional to the scalar plastic strain rate, is identical to the standard :doc:`Voce <../simple/iso_voce>`
model in NEML, just reparameterized.  The second term provides power law static recovery, which can reduce the value 
of the isotropic strength down to the threshold strength :math:`R_{min}`.

Slight change to Chaboche static recovery
"""""""""""""""""""""""""""""""""""""""""

The  *standard* :doc:`Chaboche kinematic hardening <chaboche>` in NEML uses a static recovery term of

.. math::
   \dot{\mathbf{X}}_i = - \sqrt{\frac{3}{2}} A_i \left\Vert \mathbf{X}_i \right\Vert^{a_i-1} \mathbf{X}_i

*This* model instead uses

.. math::
   \dot{\mathbf{X}}_i = - A_i \left(  \sqrt{\frac{3}{2}} \left\Vert \mathbf{X}_i \right\Vert \right)^{a_i-1} \mathbf{X}_i

The only change is the location of the :math:`\sqrt{\frac{3}{2}}` term.  This change makes the current model 
directly compatible with a 1D formulation.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8
   
   ``s0``, :cpp:class:`neml::Interpolate`, Threshold stress, No
   ``theta0``, :cpp:class:`neml::Interpolate`, Initial isotropic hardening slope, No
   ``Rmax``, :cpp:class:`neml::Interpolate`, Maximum (saturated) isotropic hardening strength, No
   ``Rmin``, :cpp:class:`neml::Interpolate`, Minimum (threshold) isotropic strength for static recovery, No
   ``r1``, :cpp:class:`neml::Interpolate`, Prefactor for isotropic recovery, No
   ``r2``, :cpp:class:`neml::Interpolate`, Exponent for isotropic recovery, No
   ``c``, :code:`std::vector<`:cpp:class:`neml::Interpolate`:code:`>`, Values of constant C for each backstress, No
   ``gmodels``, :code:`std::vector<`:cpp:class:`neml::GammaModel`:code:`>`, The gamma functions for each backstress, No
   ``A``, :code:`std::vector<`:cpp:class:`neml::Interpolate`:code:`>`, The value of A for each backstress, No
   ``a``, :code:`std::vector<`:cpp:class:`neml::Interpolate`:code:`>`, The value of a for each backstress, No
   ``noniso``, :code:`bool`, Include the nonisothermal term?, ``true``

The number of backstresses is set implicitly from the lengths of these vectors.
The model will return an error if they have different lengths.

Class description
-----------------

.. doxygenclass:: neml::ChabocheVoceRecovery
   :members:
   :undoc-members:
