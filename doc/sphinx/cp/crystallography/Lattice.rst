Lattice
=======

Overview
--------

The Lattice class provides slip system information, given the lattice
vectors and symmetry group describing the crystal system.  The class is
intelligent enough to automatically generate the complete set of 
slip group direction and normal vectors given the direction and planes
in Miller indices.

Subclasses
----------

Lattice subclasses specialize the general form to particular types of 
crystal systems, eliminating the need for the user to explicitly provide
the lattice vectors and symmetry group.

.. toctree::
   :maxdepth: 1

   CubicLattice

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``a1``, :cpp:class:`neml::Vector`, First lattice vector, No
   ``a2``, :cpp:class:`neml::Vector`, Second lattice vector, No
   ``a3``, :cpp:class:`neml::Vector`, Third lattice vector, No
   ``symmetry``, :cpp:class:`neml::SymmetryGroup`, Crystal symmetry group, No
   ``isystems``, :c:type:`list_systems`, Initial list of slip systems, ``{}``

Class description
-----------------

.. doxygenclass:: neml::Lattice
   :members:
