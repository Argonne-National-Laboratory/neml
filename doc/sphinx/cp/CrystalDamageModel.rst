CrystalDamageModel
==================

Overview
--------

This forms the base class for the implementation of a particular crystal 
plasticity continuum damage model.  The class provides the definition of the
projection operator :math:`\mathbf{P}` as well as the associated internal damage
variables :math:`d_i`.  The class also must define
all the partial derivatives required for an implicit integration of the
projection in the context of the crystal plasticity kinematics and the
set of internal variables.

Implementations
----------------

.. toctree::
   :maxdepth: 1

   cpdmg/NilDamageModel
   cpdmg/PlanarDamageModel

Class description
-----------------

.. doxygenclass:: neml::CrystalDamageModel
   :members:
   :undoc-members:
