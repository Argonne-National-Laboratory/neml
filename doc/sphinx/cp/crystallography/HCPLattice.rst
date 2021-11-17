HCPLattice
==========

Overview
--------

This class specializes the generic crystal lattice class for HCP systems.
Additionally, it overrides the base class methods which take Miller indices
as inputs to take Miller-Bravais indices instead.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``a``, :code:`double`, Lattice parameter a, No
   ``c``, :code:`double`, Lattice parameter c, No
   ``slip_systems``, :code:`list_systems`, Initial list of slip systems, ``{}``
   ``twin_systems``, :code:`list_systems`, Initial list of twin systems, ``{}``

Class description
-----------------

.. doxygenclass:: neml::HCPLattice
