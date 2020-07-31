GeneralLinearHardening
======================

Overview
--------

General linear interaction hardening -- each hardening variable evolves as a
linear combination of the slip rates (or absolute value of the slip rates)
of all other slip systems.

.. math::
   \dot{\bar{\tau}}_{k}= \sum_{j=1}^{n_{total}} M_{k,j} \dot{\gamma}_j

or

.. math::
   \dot{\bar{\tau}}_{k}= \sum_{j=1}^{n_{total}} M_{k,j} \left|\dot{\gamma}_{j}\right|

where the indices :math:`k` and :math:`j` are the unrolled indices corresponding
to slip group and system numbers.

The initial model strengths are constants.

The :ref:`matrix-system` classes provide the definition of the interaction matrix.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``M``, :cpp:class:`neml::SquareMatrix`, Interaction matrix, N
   ``tau_0``, :c:type:`std::vector<double>`, Initial strengths, N
   ``absvar``, :c:type:`bool`, If true use absolute value slip rates, :c:type`true`

Class description
-----------------

.. doxygenclass:: neml::GeneralLinearHardening
   :members:
   :undoc-members:
