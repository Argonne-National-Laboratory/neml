.. _effective-stress:

Effective stress
================

Overview
--------

These objects provide a modular effective stress system for defining the stress correlating multiaxial stress states to uniaxial creep rupture data.  Mathematically, this class provides the interface

.. math::
   \sigma_e, \frac{\partial \sigma_e}{\partial \bm{\sigma}} \leftarrow \mathcal{S}\left(\bm{\sigma} \right)

The effective stress, :math:`\sigma_e` should be some scalar function of the stress :math:`\bm{\sigma}` such that for uniaxial load the effective stress maps the uniaxial stress tensor to the value of tensile stress.

Implementations
---------------

.. toctree::
   vonmises
   maxprincipal
   huddleston
   maxseveral
   sumseveral

Class description
-----------------

.. doxygenclass:: neml::EffectiveStress
   :members:
   :undoc-members:
