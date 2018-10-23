Isotropic linear elasticity
===========================

Overview
--------

This interface represents an isotropic linear elastic model.  In Mandel
notation the stiffness tensor is 

.. math::
   \left[\begin{array}{cccccc}
      C_{1111} & C_{1122} & C_{1122} & 0 & 0 & 0\\
      C_{1122} & C_{1111} & C_{1122} & 0 & 0 & 0\\
      C_{1122} & C_{1122} & C_{1111} & 0 & 0 & 0\\
      0 & 0 & 0 & 2C_{1212} & 0 & 0\\
      0 & 0 & 0 & 0 & 2C_{1212} & 0\\
      0 & 0 & 0 & 0 & 0 & 2C_{1212}
   \end{array}\right].

In NEML :math:`C_{1111}`, :math:`C_{1122}`, and :math:`C_{1212}` are determined
by two scalar elasticity constants, some combination of the Young's modulus
:math:`E`, the Poisson's ratio :math:`\nu`, the shear modulus :math:`\mu`, 
and the bulk modulus:math:`K`.
The input to this interfaces is the temperature :math:`T`.
Internally the object first converts the provide constants to the bulk and shear
modulus and then constructs the stiffness tensor as

.. math::

   C_{1111} = \frac{4}{3} G + K

   C_{1122} = K - \frac{2}{3} G

   C_{1212} = G.

The compliance tensor is the matrix inverse of the stiffness tensor in
Mandel notation.
However, the implementation uses an analytic formula based on the scalar
elastic constants to improve performance.
Simple formulas link the bulk and shear moduli to the other scalar
elastic constants.

Parameters
----------

The user provides two modulus values ``m1`` and ``m2`` and two strings
defining which constants are being provided, ``m1_type`` and ``m2_type``.
Valid moduli types are ``"shear"``, ``"bulk"``, ``"youngs"``, and 
``"poissons"``.
The implementation checks to ensure the user provides valid moduli types
and that they provided two unique moduli.
Any combination of two scalar elastic constants fully defines the isotropic
elasticity tensor

Class description
-----------------

.. doxygenclass:: neml::IsotropicLinearElasticModel
   :members:
   :undoc-members:
