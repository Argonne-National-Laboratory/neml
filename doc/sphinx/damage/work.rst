Dissipated work damage
======================

Overview
--------

This object implements a damage model based on a critical value of
dissipated work:

.. math::

   \dot{D}=nD^{\frac{n-1}{n}}\frac{\dot{W}}{W_{crit}\left(\dot{W}\right)}

   \dot{W} = \bm{\sigma}:\dot{\boldsymbol{\varepsilon}}_{inelastic}.

The model has two parameters: :math:`W_{crit}` the critical work to failure, 
as a function of the work rate, and :math:`n`, a parameter controlling 
the onset of the appearance of damage in the material flow stress.
In this implementation the critical work is provided as a NEML interpolate
function, meaning it can have a wide variety of functional forms.  
In principle, the critical work might also depend on temperatures.  
However, at least for one material (Alloy 617) the temperature dependence
is relatively negligible.

If requested the model will first take the log of the work rate
before passing it to the :math:`W_{crit}` function and uses :math:`10^f` of the
returned value.  This means the user is providing the function on a log-log
scale.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``Wcrit``, :cpp:class:`neml::Interpolate`, Critical work, No
   ``n``, :code:`double`, Damage exponent, No
   ``eps``, :code:`double`, Numerical offset, ``1.0e-30``
   ``log``, :code:`bool`, Log transform the work to failure relation, ``false``

Class description
-----------------

.. doxygenclass:: neml::WorkDamage
   :members:
   :undoc-members:
