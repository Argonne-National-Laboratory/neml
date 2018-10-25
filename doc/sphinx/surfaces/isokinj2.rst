IsoKinJ2
========

Overview
--------

This object implements standard :math:`J_2` plasticity.  The yield
function is

.. math::
   f\left(\bm{\sigma}, \mathbf{q}, T\right) = 
      J_2\left(\bm{\sigma} + \mathbf{X}\right) + \sqrt{\frac{2}{3}} Q

   J_2\left(\mathbf{Y}\right) = \sqrt{\frac{3}{2}
      \operatorname{dev}\left(\mathbf{Y}\right):
      \operatorname{dev}\left(\mathbf{Y}\right)}.

It assumes a "stress-like" history vector of

.. math::
   \mathbf{q}=\left[\begin{array}{cc}Q & \mathbf{X}\end{array}\right]

where :math:`Q` is the isotropic hardening stress and 
:math:`\mathbf{X}` is the backstress.

.. WARNING::
   All of the NEML yield surfaces assume the opposite of the standard
   sign convention for isotropic and kinematic hardening.
   The hardening model is expected to return a negative value of the
   isotropic hardening stress and a negative value of the backstress.

Parameters
----------

None.

Class description
-----------------

.. doxygenclass:: neml::IsoKinJ2
   :members:
   :undoc-members:
