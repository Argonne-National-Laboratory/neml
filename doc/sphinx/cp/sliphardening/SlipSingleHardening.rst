SlipSingleHardening
===================

Overview
--------

This object degenerates the general SlipHardening model so that all
slip groups and systems share the same scalar slip system strength.
This actual internal variables remain generic.

This produces the interface

.. math::
   \bar{\tau}, \frac{\partial \bar{\tau}}{\partial \mathbf{h}} \leftarrow \mathcal{T}\left(\mathbf{h}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, T \right)

where :math:`\bar{\tau}` is now a single scalar.

Implementations
---------------

.. toctree::
   :maxdepth: 1

   SlipSingleStrengthHardening


Class description
-----------------

.. doxygenclass:: neml::SlipSingleHardening
   :members:
