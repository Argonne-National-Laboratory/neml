Creep models
============

Overview
--------

Creep models in NEML relate stress and strain to a creep strain rate.
This means they do not maintain a set of internal history variables.
Currently, they are only used in conduction with a standard NEML material
model through the :doc:`creep+plasticity <interfaces/creep_plasticity>` 
metamodel.

NEML creep models fulfill the interface

.. math::
   \dot{\bm{\varepsilon}}^{cr}
   \leftarrow 
   \mathcal{C}\left(\bm{\sigma}, \bm{\varepsilon}, t, T \right).

So the creep strain rate can depend on the stress, the total strain,
time, and temperature.

Implementations
---------------

.. toctree::
   creep/j2_creep

Class description
-----------------

.. doxygenclass:: neml::CreepModel
   :members:
   :undoc-members:
