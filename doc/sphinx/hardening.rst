Hardening models
================

NEML uses two types of hardening models: simple models intended for use
with associative plasticity and nonassociative hardening models.
The simple models only need to define the map between the "strain-like"
history variables and the "stress-like" internal variables
that feed into the :doc:`surfaces`.
Nonassociative models need to define both this map and the
hardening evolution rate for each "strain-like" internal variable.

The same set of classes are used to implement isotropic and kinematic hardening.

.. toctree::
   hardening/simple
   hardening/nonassociative
