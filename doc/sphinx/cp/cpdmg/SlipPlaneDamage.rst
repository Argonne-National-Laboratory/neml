SlipPlaneDamage
===============

These objects define an internal variable describing some damage process or
measure of damage along a given slip plane.  Specifically, the object
defines the evolution rate of an internal variable describing the damage
process on that plane along with the partial derivatives needed for an 
implicit integration of the rate equations.

Superclass description
----------------------

.. doxygenclass:: neml::SlipPlaneDamage
   :members:
   :undoc-members:

Individual models
-----------------

WorkPlaneDamage
^^^^^^^^^^^^^^^

Overview
""""""""

This model represents damage as evolving with the total dissipated
inelastic work on all slip systems on the given slip plane.  The
rate form of the damage variable is then

.. math::

   \dot{d}^{\left(i\right)}=\sum_{j\in S_{i}}\tau_{j}\dot{\gamma}_{j}

where the index :math:`i` is the slip plane, the set :math:`S_i` are all 
the slip directions on that slip plane, :math:`\tau_j` is the resolved shear
on a given slip system, and :math:`\dot{\gamma}_j` is the slip rate on the
slip system.

Parameters
""""""""""

None

Class description
"""""""""""""""""

.. doxygenclass:: neml::WorkPlaneDamage
   :members:
   :undoc-members:


