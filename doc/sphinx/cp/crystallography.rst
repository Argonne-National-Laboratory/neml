Crystallography classes
=======================

Overview
--------

These classes provide basic crystallographic information for the 
crystal plasticity system.
There are two key classes:

.. toctree::
   :maxdepth: 1

   crystallography/Lattice
   crystallography/SymmetryGroup

The :doc:`crystallography/SymmetryGroup` class provides the symmetry operators
for a particular point group.
The :doc:`crystallography/Lattice` class provides information
about a particular crystal system, defined by a set of lattice vectors,
a point group, and a set of slip systems.  This information includes
the cartesian slip system normal and direction vectors and the resolved
shear stress calculation needed by other objects in the crystal plasticity 
module.
