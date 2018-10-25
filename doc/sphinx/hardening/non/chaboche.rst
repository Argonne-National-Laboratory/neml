Chaboche nonassociative hardening
=================================

Overview
--------

The Chaboche model has been developed over a long period by Chaboche and 
coworkers [C2008]_ [C1989a]_ [C1989b]_.
It is one of a few canonical high temperature constitutive models for
metals that accounts for the interaction of creep and kinematic hardening.

The Chaboche model includes an extension of the combined isotropic/kinematic
hardening model originally proposed by Frederick and Armstrong [FA2007]_.
The hardening variables are a standard isotropic hardening variable defined by
a strain-like equivalent inelastic strain :math:`\alpha` mapped to a 
stress-like isotropic hardening variable describing the expansion of the
flow surface (:math:`Q`).
Supplementing this associative isotropic hardening is nonassociative 
kinematic hardening describing the shift in the flow surface in stress space.
This kinematic hardening is often called *the* Chaboche model.
It consists of a total backstress summed from several individual backstress
contributions.

The total history vector :math:`\bm{\alpha}` is 

.. math::
   \bm{\alpha}=\left[\begin{array}{ccccc} \alpha & \boldsymbol{X}_{1} & \boldsymbol{X}_{2} & \ldots & \boldsymbol{X}_{n}\end{array}\right]

where each :math:`\bm{X}_i` is one of :math:`n` individual backstresses.
The model converts this vector of history variables into the stress-like
quantities needed by the yield surface with the map

.. math::
   \mathbf{q}=\left[\begin{array}{cc} Q\left(\alpha\right) & \sum_{i=1}^{n}\mathbf{X}_{i}\end{array}\right]

The map for the isotropic hardening is provided by another object 
implementing the :doc:`isotropic hardening interface <../simple/isotropic>`.

The history evolution equation for the isotropic part of the model is associative and proportional only to the scalar inelastic strain rate 

.. math::
   \mathbf{h}_\gamma^\alpha = \sqrt{\frac{2}{3}}.

The evolution equations for each individual backstress are:

.. math::

   \mathbf{h}_{\gamma}^{X_i} = -\frac{2}{3} C_i \frac{\operatorname{dev}\left(\bm{\sigma} - \mathbf{X}\right)}{\left\Vert \operatorname{dev}\left(\bm{\sigma} - \mathbf{X}\right) \right\Vert} - \sqrt{\frac{2}{3}} \gamma_i\left(\alpha, T \right) \mathbf{X}_i   

   \mathbf{h}_{t}^{X_i} = -\sqrt{\frac{3}{2}} A_i \left\Vert \mathbf{X}_i \right\Vert ^ {a_i - 1} \mathbf{x}_i

   \mathbf{h}_{t}^{x_i} = -\sqrt{\frac{2}{3}} \frac{\mathrm{d}C_i / \mathrm{d}t}{C_i} \mathbf{X}_i

The model parameters for each backstress are then :math:`C_i`, :math:`A_i`, :math:`a_i`, and :math:`\gamma_i` as a function of the equivalent inelastic strain.
Various options for gamma are :ref:`described below <gamma-models>`.
The model maintains :math:`1 + 6n` history variables, where :math:`n` is the number of
backstresses.

The implementation has an option to turn on or off the part of the backstress evolution proportional to the temperature rate.
Note this contribution is in any event zero for isothermal loading.

For one backstress (:math:`n=1`), :math:`A_1 = 0`, :math:`a_1 = 1`, and the 
non-isothermal term turned off the model degenerates to the classical Frederick-Armstrong model [FA2007]_.

.. WARNING::
   All of the NEML yield surfaces assume the opposite of the standard
   sign convention for isotropic and kinematic hardening.
   The hardening model is expected to return a negative value of the
   isotropic hardening stress and a negative value of the backstress.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``iso``, :cpp:class:`neml::IsotropicHardeningRule`, Isotropic hardening, No
   ``c``, :c:type:`std::vector<`:cpp:class:`neml::Interpolate`:c:type:`>`, Values of constant C for each backstress, No
   ``gmodels``, :c:type:`std::vector<`:cpp:class:`neml::GammaModel`:c:type:`>`, The gamma functions for each backstress, No
   ``A``, :c:type:`std::vector<`:cpp:class:`neml::Interpolate`:c:type:`>`, The value of A for each backstress, No
   ``a``, :c:type:`std::vector<`:cpp:class:`neml::Interpolate`:c:type:`>`, The value of a for each backstress, No
   ``noniso``, :c:type:`bool`, Include the nonisothermal term?, ``true``

The number of backstresses is set implicitly from the lengths of these vectors.
The model will return an error if they have different lengths.

Class description
-----------------

.. doxygenclass:: neml::Chaboche
   :members:
   :undoc-members:

.. _gamma-models:

Gamma models
------------

The :math:`\gamma` parameter describes dynamic backstress recovery in the Chaboche model.
This tends to send the backstress to some saturated shift of the yield surface with
increasing inelastic strain.
The Chaboche model allows this dynamic recovery coefficient to vary with the accumulated
inelastic strain.
These objects then define the dynamic recovery parameter with the interface

.. math::
   \gamma, \frac{\partial\gamma}{\partial \alpha} 
   \leftarrow
   \mathcal{G}\left( \alpha, T \right)

Class description
"""""""""""""""""

.. doxygenclass:: neml::GammaModel
   :members:
   :undoc-members:

Constant gamma
^^^^^^^^^^^^^^

This function returns a value of :math:`\gamma` that is independent of inelastic strain.
It still might depend on temperature.  The implementation is

.. math::
   \gamma = C.

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``g``, :cpp:class:`neml::Interpolate`, Value of gamma as a function of T, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::ConstantGamma
   :members:
   :undoc-members:

Saturating gamma
^^^^^^^^^^^^^^^^

This gamma function begins a given value and transitions towards a second, saturated
value as a function of accumulated inelastic strain.
It implements the function

.. math::
   \gamma = \gamma_{s} + \left(\gamma_0 - \gamma_s \right) e^{-\beta \alpha}.

Parameters 
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``g0``, :cpp:class:`neml::Interpolate`, Initial value of gamma, No
   ``gs``, :cpp:class:`neml::Interpolate`, Final value of gamma, No
   ``beta``, :cpp:class:`neml::Interpolate`, Controls the saturation rate, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::SatGamma
   :members:
   :undoc-members:
