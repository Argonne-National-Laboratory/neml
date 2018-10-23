Chaboche nonassociative hardening
=================================

Overview
--------

The Chaboche model has been developed over a long period by Chaboche and 
coworkers .
It is one of a few "canonical" high temperature constitutive models for
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
   \bm{\alpha} = \bm{\alpha}=\left[\begin{array}{ccccc} \alpha & \boldsymbol{X}_{1} & \boldsymbol{X}_{2} & \ldots & \boldsymbol{X}_{n}\end{array}\right]

where each :math:`\bm{X}_i` is one of :math:`n` individual backstresses.
The model converts this vector of history variables into the stress-like
quantities needed by the yield surface with the map

.. math::
   \mathbf{q}=\left[\begin{array}{cc} Q\left(\alpha\right) & \sum_{i=1}^{n}\mathbf{X}_{i}\end{array}\right]

The map for the isotropic hardening is provided by another object 
implementing the :doc:`isotropic hardening interface <../simple/isotropic>`.

Parameters
----------


Class description
-----------------


Gamma models
------------


Constant gamma
^^^^^^^^^^^^^^


Parameters
""""""""""


Class description
"""""""""""""""""


Saturating gamma
^^^^^^^^^^^^^^^^


Parameters 
""""""""""


Class description
"""""""""""""""""


