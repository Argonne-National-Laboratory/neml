StandardKinematicModel
======================

Overview
--------

This KinematicModel implements the kinematic assumptions described in the
:ref:`overview <cp-formulation>` of the crystal plasticity model.
The stress and rotation rates match those developed in that derivation.
The history rate and variables are left to be generic, as are the exact
form of the plastic deformation rate.  These are defined by the

.. toctree::
   :maxdepth: 1

   InelasticModel

class and associated objects.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``emodel``, :cpp:class:`neml::LinearElasticModel`, Elasticity tensor, No
   ``imodel``, :cpp:class:`neml::InelasticModel`, Definition of the plastic deformation rate, No

Class description
-----------------

.. doxygenclass:: neml::StandardKinematicModel
   :members:
   :undoc-members:
