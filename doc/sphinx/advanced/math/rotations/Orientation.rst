.. _Orientation:

Orientation
===========

This class represents unit quaternions, i.e. 3D rotations, i.e. members 
of the special orthogonal group :math:`\mathrm{SO}\left(3\right)`.
It provides the ``apply`` method capable of applying a rotation to:

* Vector

* RankTwo

* Symmetric

* Skew

* RankFour

* SymSymR4

classes.

The class has static methods that can be used to create rotations from

* Euler angles

  * Kocks convention

  * Bunge convention

* Axis/angle pairs

* Rodrigues vectors

* Rotation matrices

* Hopf coordinates

* Hyperspherical coordinates

Similar methods can be used to convert the quaternion to these representations
for output.

.. doxygenclass:: neml::Orientation
   :members:
   :undoc-members:
