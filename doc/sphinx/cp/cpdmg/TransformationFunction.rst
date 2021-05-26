TransformationFunction
======================

These models are part of the :any:`PlanarDamageModel` system.  They map the stress normal to a slip plane :math:`\sigma_{\bot}^{\left(i\right)}` and an 
internal damage variable, defined with a :any:`SlipPlaneDamage` object, to a 
suitable damage metric ranging from 0 (no damage) to 1 (no stiffness in
the direction).

Superclass description
----------------------

.. doxygenclass:: neml::TransformationFunction
   :members:
   :undoc-members:

Individual models
-----------------

SigmoidTransformation
^^^^^^^^^^^^^^^^^^^^^

Overview
""""""""

This transformation function is independent of the normal stress.  It
is a sigmoid function that maps the damage variable to the
interval :math:`[0,1]`:

.. math::

   T\left(d^{\left(i\right)},\sigma_{\bot}^{\left(i\right)}\right)=\begin{cases}
   \frac{1}{1+\left(\frac{d^{\left(i\right)}}{c-d^{\left(i\right)}}\right)^{-\beta}} & d^{\left(i\right)}<c\\
   1 & d^{\left(i\right)}\ge c
   \end{cases}

The parameter :math:`c^{(i)}` is some critical value of the damage variable
:math:`d^{(i)}` on plane :math:`i` and the exponent :math:`\beta^{(i)}` is a
smoothness parameter where :math:`\beta = 1` represents a linear transition
from 0 to 1 and higher values of :math:`\beta` represent a more abrupt onset of
the effects of damage.

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``c``, :c:type:`double`, Critical damage value, No
   ``beta``, :c:type`double`, Abruptness of damage onset, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::SigmoidTransformation
   :members:
   :undoc-members:

SwitchTransformation
^^^^^^^^^^^^^^^^^^^^

Overview
""""""""

This is a "meta-transformation" that takes another :any:`TransformationFunction` as a parameter and modifies it to return a different result.  
Specifically, this function returns 0 (no damage) if the loading on the 
plane is compressive and returns the base :any:`TransformationFunction` value
if the stress on the plane is tensile.  Mathematically:

.. math::

   T\left(d^{\left(i\right)},\sigma_{\bot}^{\left(i\right)}\right)=\begin{cases}
   \tilde{T}\left(d^{\left(i\right)},\sigma_{\bot}^{\left(i\right)}\right) & \sigma_{\bot}^{\left(i\right)}\ge0\\
   0 & \sigma_{\bot}^{\left(i\right)}<0
   \end{cases}

where :math:`\tilde{T}\left(d^{\left(i\right)},\sigma_{\bot}^{\left(i\right)}\right)` is the base :any:`TransformationFunction`.

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``base``, :cpp:class:`neml::TransformationFunction`, Base function, No

Class Description
"""""""""""""""""

.. doxygenclass:: neml::SwitchTransformation
   :members:
   :undoc-members:
