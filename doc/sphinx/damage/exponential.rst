Exponential damage
==================

Overview
--------

This object implements a "standard" damage model proportional to the dissipated
inelastic energy.
The damage function is

.. math::

   w = \frac{\left(\omega + k_0\right)^{a_f}}{W_0} \sigma_{eff}

   \sigma_{eff} = \sqrt{\frac{3}{2} \operatorname{dev}\left(\bm{\sigma}\right):\operatorname{dev}\left(\bm{\sigma}\right)}.

The standard damage model multiplies this function by the inelastic
strain rate in computing the damage update.
Because the dissipation rate is equal to :math:`\sigma_{eff} \dot{\varepsilon}_{eff}^{in}` this model actually increases damage in proportion to the dissipation.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``elastic``, :cpp:class:`neml::LinearElasticModel`, Elasticity model, No
   ``W0``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``k0``, :cpp:class:`neml::Interpolate`, Parameter, No
   ``af``, :cpp:class:`neml::Interpolate`, Parameter, No

Class description
-----------------

.. doxygenclass:: neml::ExponentialWorkDamage
   :members:
   :undoc-members:
