Walker-Krempl rate sensitivity modification
===========================================

Overview
--------

This class modifies an existing :doc:`viscoplastic` to scale
between a rate sensitive response (given by the base model) and a 
rate insensitive response based on a temperature-dependent
parameter :math:`\lambda`.
The approach was developed by Walker and Krempl [WK1978]_ to modify an 
underlying viscoplastic model to return a rate insensitive response
at lower temperatures.

Given the base model scalar inelastic rate 

.. math::
   \dot{\bm{\varepsilon}}_{vp} = \dot{p} \mathbf{g}

the model modifies the scalar part of the flow rate to

.. math::
   \dot{p}_{mod} = \kappa \dot{p}

keeping the flow direction the same.  The scaling function is

.. math::
   \kappa = 1 - \lambda + \frac{\lambda \sqrt{\frac{2}{3}} \left\Vert \dot{\boldsymbol{e}}\right\Vert }{\dot{\varepsilon}_{ref}}

where :math:`\dot{\varepsilon}_{ref}` is a parameter and 
:math:`\left\Vert \dot{\boldsymbol{e}}\right\Vert` is the total (not inelastic) deviatoric strain rate.

In the limit :math:`\lambda = 0` the model returns the rate sensitive, viscoplastic flow response from the underlying model.
In the limit :math:`\lambda \rightarrow 1` the model returns an equivalent rate
insensitive flow rate, asymptotically approximating an equivalent rate
independent model.  An interpolation function can change the value of
:math:`\lambda` as a function of temperature, for example to transition
from a rate insensitive response at lower temperatures to a rate
sensitive response at higher temperatures.

The implementation applies the same scale factor to all time derivatives in
the base model's history evolution equations, maintaining a consistent
modified time integration scheme for all internal variables.

The user must determine an appropriate value of :math:`\lambda` (which
cannot be set exactly to :math:`\lambda = 1`, as this results in numerical
instability) and an appropriate reference rate :math:`\dot{\varepsilon}_{ref}`.
A value of :math:`\lambda = 0.99` seems to work well and the
reference rate should be set to a value much slower than the actual 
strain rates used by the material model.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``flow``, :cpp:class:`neml::ViscoPlasticFlowRule`, Viscoplastic flow rule interface, No
   ``lambda``, :cpp:class:`neml::Interpolate`, Scaling parameter, No
   ``eps0``, :code:`double`, Reference strain rate, No

Class description
-----------------

.. doxygenclass:: neml::WalkerKremplSwitchRule
   :members:
   :undoc-members:
