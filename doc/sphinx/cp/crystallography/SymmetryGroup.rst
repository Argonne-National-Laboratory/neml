SymmetryGroup
=============

Overview
--------

This class represents the symmetry operations associated with a 
crystallographic point group.  Essentially, it is an interface to a list
of quaternion symmetry operations.  These operations have been hardcoded
and verified using an automated unit test.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``sclass``, :c:type:`string`, Point group in Hermann-Mauguin notation, No

Class description
-----------------

.. doxygenclass:: neml::SymmetryGroup
   :members:
