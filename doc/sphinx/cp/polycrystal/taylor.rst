TaylorModel
===========

Overview
--------

This polycrystal homogenization model implements the standard Taylor approximation.  Here the individual crystal receives the same deformation information and the resulting stresses are averaged.

Mathematically, each crystal receives the same, macroscopic :math:`\mathbf{D}` and :math:`\mathbf{W}` deformation rate objects

.. math::
   \mathbf{D}_i = \mathbf{D}

   \mathbf{W}_i = \mathbf{W}

and the macroscale stress is:

.. math::
   \bf{\sigma} = \frac{1}{n}\sum_{i=1}^{n_{crystal}}\bm{\sigma}_{i}

The stress updates can be completed in parallel using OpenMP threads.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``model``, :cpp:class:`neml::SingleCrystalModel`, Single crystal update, N
   ``qs``, :c:type:`std::vector<`:cpp:class:`neml::Orientation`:c:type:`>`, Vector of orientations, N
   ``nthreads``, :c:type:`int`, Number of threads to use, 1

Class description
-----------------

.. doxygenclass:: neml::TaylorModel
   :members:
   :undoc-members:
