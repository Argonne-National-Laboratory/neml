.. _interpolate-functions:

Interpolation functions
=======================

Throughout the code NEML uses interpolation functions to represent
model parameters that depend on temperature or, occasionally, some other
variable.
A NEML interpolation function is a scalar function of a single variable
:math:`f\left( x \right)`.

The interface must also provide the first derivative of the function with
respect that the variable.

These interpolation functions are used throughout NEML in place of scalar
parameters.
Except where noted, currently only for the InterpolatedIsotrpicHardeningRule,
interpolation functions are used to give the parameter as a function of
temperature.
A constant parameter, e.g. one that does not depend on temperature can be
expressed by using a :ref:`constant` object.

Interpolate
-----------

The interface for all interpolate objects.

.. doxygenclass:: neml::Interpolate
   :members:
   :undoc-members:

.. _constant:

ConstantInterpolate
-------------------

Expresses a constant parameter

.. math::
   f\left( x \right) = C.

.. doxygenclass:: neml::ConstantInterpolate
   :members:
   :undoc-members:

PolynomialInterpolate
---------------------

Polynomial interpolation of the type

.. math::
   f\left( x \right) = \sum_{i=1}^{n}c_{i}x^{i-1}.

The polynomial coefficients are given starting with the highest order, 
like in the `numpy <https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.polyval.html>`_ library.

.. doxygenclass:: neml::PolynomialInterpolate
   :members:
   :undoc-members:

.. _linear:

PiecewiseLinearInterpolate
--------------------------

Piecewise linear interpolation where parameters give the *x* and *y* coordinates
of the endpoints of each linear segment.
For :math:`x_{i} \le x \le x_{i+1}`

.. math::
   f\left( x \right) = \frac{y_{i+1} - y_{i}}{x_{i+1} - x_{i}}
      \left(x - x_i \right) + y_i.

For :math:`x < x_1` the function returns :math:`y_1` and for :math:`x > x_n`
the function returns :math:`y_n`.

.. doxygenclass:: neml::PiecewiseLinearInterpolate
   :members:
   :undoc-members:

PiecewiseLogLinearInterpolate
-----------------------------

Exactly like :ref:`linear` except the interpolation in *y* is done in log space.
The input *y* values are given as :math:`y_i = \ln{v_i}` where :math:`v_i`
is the actual value of the function at :math:`x_i`.
For :math:`x_{i} \le x \le x_{i+1}`

.. math::
   f\left( x \right) = \exp{\left(\frac{y_{i+1} - y_{i}}{x_{i+1} - x_{i}}
      \left(x - x_i \right) + y_i\right)}.

For :math:`x < x_1` the function returns :math:`\exp{\left(y_1\right)}` and for :math:`x > x_n` the function returns :math:`\exp{\left(y_n\right)}`.

.. doxygenclass:: neml::PiecewiseLogLinearInterpolate
   :members:
   :undoc-members:


MTSShearInterpolate
-------------------

The shear modulus interpolation used in the Mechanical Threshold Stress [MTS1999]_
flow stress model 

.. math::
   f\left( x \right) = V_0 - \frac{D}{e^{T_0 / x} - 1}.

.. doxygenclass:: neml::MTSShearInterpolate
   :members:
   :undoc-members:

Helper Functions
----------------

.. doxygenfunction:: neml::make_vector

.. doxygenfunction:: neml::eval_vector

.. doxygenfunction:: neml::eval_deriv_vector
