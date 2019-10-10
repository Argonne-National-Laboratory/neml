InelasticModel
==============

Overview
--------

Inelastic models provide the plastic deformation rate and internal variable
history rate equations to the high level objects, along with associated
partial derivatives.
The mathematical interface is:

.. math::
     \mathbf{D}^p, \frac{\partial \mathbf{D}^p}{\partial \bm{\sigma}}, \frac{\partial \mathbf{D}^p}{\partial \mathbf{h}} \leftarrow \mathcal{D}\left(\bm{\sigma}, \mathbf{h}, T \right)

     \mathbf{W}^p, \frac{\partial \mathbf{W}^p}{\partial \bm{\sigma}}, \frac{\partial \mathbf{W}^p}{\partial \mathbf{h}} \leftarrow \mathcal{D}\left(\bm{\sigma}, \mathbf{h}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, T \right)

where all these quantities are defined in the crystal plasticity :ref:`overview <cp-formulation>`.
Note the implementations also have the crystallographic information
available through the model's :doc:`crystallography/Lattice` object.

Implementations
---------------

.. toctree::
   :maxdepth: 1

   inelasticity/AsaroInelasticity
   inelasticity/NoInelasticity
   inelasticity/PowerLaw
   inelasticity/CombinedInelasticity

Class description
-----------------

.. doxygenclass:: neml::InelasticModel
   :members:
