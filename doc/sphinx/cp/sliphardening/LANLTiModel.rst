LANLTiModel
===========

Overview
--------

Put your description and some math here

Parameters
----------

Fill in the table with your actual required and optional parameters.

If it's a required parameter put "N" in the last column.

If it's an optional parameter put the default value in the last column

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``M``, :cpp:class:`neml::SquareMatrix`, Interaction matrix, N
   ``tau_0``, :code:`std::vector<double>`, Initial strengths, N
   ``absvar``, :code:`bool`, If true use absolute value slip rates, :c:type`true`

Class description
-----------------

.. doxygenclass:: neml::LANLTiModel
   :members:
   :undoc-members:
