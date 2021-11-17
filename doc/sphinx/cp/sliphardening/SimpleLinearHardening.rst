SimpleLinearHardening
======================

Overview
--------

A simple linear hardening model, commonly used to test features of more
complex systems.

.. math::
   \tau_{i}=\tau_{0,i}+\sum_{j=1}^{n_{slip}}G_{ij}\gamma_{j}

where the indices :math:`i` and :math:`j` are the unrolled indices corresponding
to slip group and system numbers.

The model maintains the integrated slips as the internal variables.

The :ref:`matrix-system` classes provide the definition of the interaction matrix.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``G``, :cpp:class:`neml::SquareMatrix`, Interaction matrix, N
   ``tau_0``, :code:`std::vector<double>`, Initial strengths, N

Class description
-----------------

.. doxygenclass:: neml::SimpleLinearHardening
   :members:
   :undoc-members:
