GeneralLattice
==============

Overview
--------

A general :doc:`Lattice` class where the user can provide the lattice
directions and the symmetry class.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``a1``, :code:`std::vector<double>`, First lattice vector, No
   ``a2``, :code:`std::vector<double>`, Second lattice vector, No
   ``a3``, :code:`std::vector<double>`, Third lattice vector, No
   ``symmetry_group``, :cpp:class:`neml::SymmetryGroup`, Symmetry group, No
   ``slip_systems``, :code:`list_systems`, Initial list of slip systems, ``{}``
   ``twin_systems``, :code:`list_systems`, Initial list of twin systems, ``{}``

Class description
-----------------

.. doxygenclass:: neml::GeneralLattice
