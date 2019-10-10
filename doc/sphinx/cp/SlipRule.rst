SlipRule
========

Overview
--------

These objects provide a relation between the stress, history, and temperature
and the slip rate on each individual slip system as well as the history evolution.  The interface used is:

.. math::
   \dot{\gamma}_{g,i}, \frac{\partial \dot{\gamma}_{g,i}}{\partial \bm{\sigma}}, \frac{\partial \dot{\gamma}_{g,i}}{\partial \mathbf{h}} \leftarrow \mathcal{G}\left( \bm{\sigma}, \mathbf{h}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, T \right)

where :math:`g` indicates the slip group and :math:`i` indicates the system within the group.


Implementations
---------------

.. toctree::
   :maxdepth: 1

   sliprule/SlipStrengthSlipRule

Class description
-----------------

.. doxygenclass:: neml::SlipRule
