.. _history:

History object system
=====================

Overview
--------

The History class stores internal variables of material models. 
Originally, NEML stored internal variables in a flat pointer array.
This meant that the programmer or end user had to manually track 
which variable was located at which index in the array and, for non-scalar
variables, correctly flatten and interpret the data.

The History class helps manage model internal variables.  Fundamentally,
it is a dictionary that returns and stores an internal variable 
associated with a string key naming the actual variable in the model
implementation.  The class also manages types so that it returns the correct
type of object.  For example, the History class could manage a ``"direction"``
internal variable with type ``Vector``.  The user can request the
object return ``"direction"`` and the History object would return 
the correct ``Vector`` type.

Currently, the class is configured to store and return objects of the following types:

* Scalars (i.e. :c:type:`double`)

* :doc:`math/tensors/Vector`

* :doc:`math/tensors/RankTwo`

* :doc:`math/tensors/Symmetric`

* :doc:`math/tensors/Skew`

* :doc:`math/rotations/Orientation`

* :doc:`math/tensors/SymSymR4`

Internally, all of these classes are stored in a single flat array.  
The system can either manage its own memory (freeing it on deletion of the
object) or accept a pointer to externally managed memory.
For example, a finite element analysis program can pass in the material
model history as a large array in memory and NEML can `wrap` that memory
using the History class and give descriptive access to each variable `without`
copying the underlying data.
This allows for seamlessly converting a flat array containing the model
internal variables into an instance of the History class without copying.
It also allows the calling program to block the internal variables of
several material points, potentially allowing the compiler to perform
vector optimizations.

Oftentimes, implementations of material models will need to store the
derivatives of a set of internal variables.  The History class contains
methods for duplicating a given History instance but changing the 
types (and corresponding storage) of the new instance to reflect the
types appropriate for storing the derivative of the original History with
respect to some other object type.

Class description
-----------------

.. doxygenclass:: neml::History
   :members:
