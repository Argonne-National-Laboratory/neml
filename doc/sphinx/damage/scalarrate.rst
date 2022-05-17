Scalar damage, defined in rate form
===================================

Overview
--------

This object further simplifies a :cpp:class`neml::ScalarDamage` model to provide the updated damage value by integrating the damage rate.  The model then defines the updated damage as

.. math::
   \omega_{n+1} = \omega_{n} + \dot{\omega} \Delta t

where :math:`\dot{\omega}` is the damage rate, defined in a subclass with
the `damage_rate` method.  The associated derivatives of the updated damage
with respect to the current damage, the stress, and the strain and defined
analogously.

Implementations
---------------

.. toctree::
   modular
   classical
   larsonmiller

Class description
-----------------

.. doxygenclass:: neml::ScalarDamageRate
   :members:
   :undoc-members:
