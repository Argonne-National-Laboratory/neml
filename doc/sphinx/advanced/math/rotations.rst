Rotations
=========

Rotations in NEML are stored as quaternions.  In turn, these are stored as
length-4 flat arrays.

NEML provides classes implementing general quaternions and unit quaternions, which are representations of 3D rotations (i.e. elements of the special orthogonal group).  Currently, the general quaternion class is not used in NEML, it just serves as the base class for the Orientation class, which implements rotations.

.. toctree::
   :maxdepth: 1

   rotations/Quaternion
   rotations/Orientation

As with the :doc:`tensors` classes, NEML can either manage the memory of a
quaternion or use externally managed memory, to allow for efficient
block storage of data.

The quaternion classes provide several helpful unary operators, for example methods for inverting and exponentiating quaternions.  The implementation also provides binary operators for composing quaternions (equivalent to composing rotations for the Orientation unit quaternion class).
It also provides several specialized mathematical operators:

.. doxygenfunction:: neml::random_orientations

.. doxygenfunction:: neml::wexp

.. doxygenfunction:: neml::wlog

.. doxygenfunction:: neml::distance

.. doxygenfunction:: neml::rotate_to

.. doxygenfunction:: neml::rotate_to_family
