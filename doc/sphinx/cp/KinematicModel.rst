KinematicModel
==============

Overview
--------

A kinematics model provides the interfaces:

.. math::
   \dot{\bm{\sigma}}, \frac{\partial \dot{\bm{\sigma}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\bm{\sigma}}}{\partial \mathbf{h}}, \frac{\partial \dot{\bm{\sigma}}}{\partial \mathbf{D}}, \frac{\partial \dot{\bm{\sigma}}}{\partial \mathbf{W}} \leftarrow \mathcal{S}\left(\bm{\sigma}, \mathbf{h}, \mathbf{D}, \mathbf{W}, \bm{\alpha}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{D}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{W}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, \mathbf{D}, \mathbf{W}, \bm{\alpha}, T \right)

   \bm{\Omega}^e, \frac{\partial \bm{\Omega}^e}{\partial \bm{\sigma}}, \frac{\partial \bm{\Omega}^e}{\partial \mathbf{h}}, \frac{\partial \bm{\Omega}^e}{\partial \mathbf{D}}, \frac{\partial \bm{\Omega}^e}{\partial \mathbf{W}} \leftarrow \mathcal{W}\left(\bm{\sigma}, \mathbf{h}, \mathbf{D}, \mathbf{W}, \bm{\alpha}, T \right)

This general interface defines the stress, history, and orientation rates used
in integrating the single crystal model.
The interface also allows the user to select which parts of the update are
done in a coupled, implicit integration and which parts are left separate 
for an explicit, uncoupled integration.

Implementations
---------------

Currently there is only one implementation of this class, 
implementing the kinematic assumptions described in the
:ref:`overview <cp-formulation>` of the crystal plasticity model.
Other kinematics could be implemented in the future by 
deriving additional subclasses from this abstract base class.

.. toctree::
   :maxdepth: 1

   StandardKinematicModel

Class description
-----------------

.. doxygenclass:: neml::KinematicModel
   :members:
   :undoc-members:
