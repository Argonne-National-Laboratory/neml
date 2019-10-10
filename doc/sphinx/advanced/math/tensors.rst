Tensors
=======

The tensor objects provide a method for easily manage scalar, vector, and
tensor mathematical operations.
The classes maintain a representation of the tensor type and provide
common mathematical operations, both single tensors operations and 
binary operations between various types of tensors.

The classes are structured with either maintain their own memory or to
use externally managed memory.  The latter mode is useful in finite element
calculations where the calling program can present a large collection of material points in a blocked array.  If this type of memory management can be used it reduces copying when NEML models are called from external program and allows optimizing compilers to attempt vectorization.

The objects are set up to be transparent as to the type of memory model use to store the data.  Generally, initializing an object with a raw pointer sets the objects to external memory management.  This means that when the object is deleted the memory is not deallocated.  Initializing an object with either a default (no parameter) constructor or a STL container sets the classes to manage their own memory and hence the memory will be freed when the object is deleted.

All the tensor classes derive from a common base class (see below).  NEML implements the following tensor types:

.. toctree::
   :maxdepth: 1

   tensors/Vector
   tensors/RankTwo
   tensors/Symmetric
   tensors/Skew
   tensors/RankFour
   tensors/SymSym
   tensors/SymSkew
   tensors/SkewSym

Binary operators
----------------

The module provides a large variety of binary operators.  All scalar multiplication and division operations are covered as are all rational addition/subtraction operations between tensors.
Implemented tensor products include the outer products between vectors
(:math:`C_{ij} = a_i b_j`) and rank two tensors (:math:`C_{ijkl} = A_{ij} B_{kl}`),
matrix-vector products
(:math:`c_j = A_{ij} b_j`),
matrix-matrix products
(:math:`C_{ij} = A_{ik} B_{kj}`),
rank four/rank two composition
(:math:`C_{ij} = A_{ijkl} B_{kl}`),
rank four/rank four composition
(:math:`C_{ijkl} = A_{ijmn} B_{mnkl}`), 
and a few specialized operators that occur in the crystal plasticity 
integration formula.

These binary operators are implemented for all sensible combinations of specialized tensor types (for example, all RankTwo binary operations can
be completed with any combination of RankTwo, Symmetric, and Skew tensors).

Base class description
----------------------

.. doxygenclass:: neml::Tensor
   :members:
   :undoc-members:
