PlanarDamageModel
=================

Overview
--------

This model forms the total damage projection operator by degrading the elasticity tensor along individual slip planes.  For each slip plane the framework
can degrade the elastic stiffness in the normal direction of the plane and in
the shear (parallel) directions independently.  The stiffness degradation
in each of these directions is defined by a :any:`TransformationFunction` which
converts a damage internal variable, defined by a :any:`SlipPlaneDamage`
object, and the stress in the normal direction
to the plane to a suitable damage index ranging from 0 (no damage) to 1 
(complete loss of stiffness in that direction).  The mathematical
definition of the projection operator is then

.. math::

   P_{ijkl}=\prod_{i=1}^{n_{planes}}\delta_{ij}\delta_{kl}-N\left(d^{\left(i\right)},\sigma_{\bot}^{\left(i\right)}\right)N_{ijkl}^{\left(i\right)}-S\left(d^{\left(i\right)},\sigma_{\bot}^{\left(i\right)}\right)S_{ijkl}^{\left(i\right)}

where the product proceeds over each individual slip plane :math:`i` defined
by some normal vector `in the current coordinates` :math:`n_{i}^{(i)}`.
The projection for each plane then consists of the identity, a 
:any:`TransformationFunction` for the normal direction :math:`N`, 
the normal projection operator

.. math::

   N_{ijkl}^{\left(i\right)}=n_{i}^{\left(i\right)}n_{j}^{\left(i\right)}n_{k}^{\left(i\right)}n_{l}^{\left(i\right)}

a :any:`TransformationFunction` for the shear direction :math:`S`, and
the shear projection operator

.. math::

   S_{ijkl}^{\left(i\right)}=\left(\delta_{ik}-n_{i}^{\left(i\right)}n_{k}^{\left(i\right)}\right)n_{j}^{\left(i\right)}n_{l}^{\left(i\right)}.

The stress normal to the plane can be calculated as

.. math::

   \sigma_{\bot}^{\left(i\right)}=\sigma_{ij}n_{i}^{\left(i\right)}n_{j}^{\left(i\right)}.

The available :any:`TransformationFunction` options are described here:

.. toctree::
   :maxdepth: 1

   TransformationFunction

The available :any:`SlipPlaneDamage` functions are described here:

.. toctree::
   :maxdepth: 1

   SlipPlaneDamage

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``damage``, :cpp:class:`neml::SlipPlaneDamage`, The damage model, No
   ``shear_transformation``, :cpp:class:`neml::TransformationFunction`, The shear transformation function, No
   ``normal_transformation``, :cpp:class:`neml::TransformationFunction`, The normal transformation function, No
   ``lattice``, :cpp:class:`neml::Lattice`, The lattice object describing the slip geometry, No

Class description
-----------------

.. doxygenclass:: neml::PlanarDamageModel
   :members:
   :undoc-members:
