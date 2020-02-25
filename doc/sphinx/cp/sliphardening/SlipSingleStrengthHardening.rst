.. _slipsinglestrengthhardening:

SlipSingleStrengthHardening
===========================

Overview
--------

This object further degenerates the general SlipSingleHardening model to where the (single) slip system strength is equal to some static strength, a contribution from a single internal variable, defined by a scalar rate equation, and a contribution from the Nye tensor.
If this scalar history variable is :math:`\tilde{\tau}` the history map is

.. math::
   \bar{\tau} = \tilde{\tau} + \tau_0 + \tau_{nye}

where :math:`\tau_0` is some static strength that does not evolve with time and :math:`\tau_{nye}` is a
function of the Nye tensor :math:`\bm{\alpha}`.  The Nye tensor contribution defaults to zero, which gives
a classical slip hardening model.

The remaining interface must provide the evolution equation, and associated partial derivatives, of the scalar history variable:

.. math::
   \dot{\tilde{\tau}}, \frac{\partial \dot{\tilde{\tau}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\tilde{\tau}}}{\partial \tilde{\tau}}, \tau_0 \leftarrow \mathcal{H}\left(\bm{\sigma}, \tilde{\tau}, \bm{\alpha}, T \right).


Implementations
---------------

.. toctree::
   :maxdepth: 1

   PlasticSlipHardening


Class description
-----------------

.. doxygenclass:: neml::SlipSingleStrengthHardening
   :members:
   :undoc-members:
