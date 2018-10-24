IsoJ2I1
=======

Overview
--------

This object implements the yield function

.. math::
      f\left(\bm{\sigma}, \mathbf{q}, T\right) = 
      J_2\left(\bm{\sigma}\right) + \sqrt{\frac{2}{3}}Q + 
      \operatorname{sign}\left(I_1\left(\bm{\sigma}\right)\right)
      h \left(I_1\left(\bm{\sigma}\right)\right)^l

   J_2\left(\mathbf{Y}\right) = \sqrt{\frac{3}{2}
      \operatorname{dev}\left(\mathbf{Y}\right):
      \operatorname{dev}\left(\mathbf{Y}\right)}

   I_1\left(\mathbf{Y}\right) = \operatorname{tr}\left(\mathbf{Y}\right).

It assumes a "stress-like" history vector of

.. math::
   \mathbf{q}=\left[\begin{array}{c}Q\end{array}\right]

where :math:`Q` is the isotropic hardening stress.

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

   ``h``, :cpp:class:`neml::Interpolate`, Power law prefactor, No
   ``l``, :cpp:class:`neml::Interpolate`, Power law exponent, No

Class description
-----------------

.. doxygenclass:: neml::IsoJ2I1
   :members:
   :undoc-members:
