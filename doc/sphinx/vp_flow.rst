Viscoplastic flow rule
======================

Overview
--------

This class provides the interface

.. math::
   \dot{\gamma}, 
   \mathbf{g}_\gamma, \mathbf{g}_T, \mathbf{g}_t,
   \mathbf{h}_\gamma, \mathbf{h}_T, \mathbf{h}_t
   \leftarrow \mathcal{V}
   \left(\bm{\sigma}, \bm{\alpha}, T \right).

*and associated partial derivatives*.
:math:`\dot{\gamma}` is the scalar inelastic strain rate, :math:`\mathbf{g}`
is the flow rule, and :math:`\mathbf{h}` is the hardening law.
The subscripts :math:`\gamma` indicates the part proportional to the 
scalar inelastic strain rate, :math:`T` the part proportional to the 
temperature rate, :math:`t` the part directly contributing towards the
total time derivative.
See :doc:`general_flow/viscoplastic` for the specific definition of 
each quantity.

The base class implementation by default provides zero for the :math:`T` and
the :math:`t` quantities, giving a standard model that evolves only in 
proportion to the inelastic strain rate.
Static recovery or thermo-viscoplasticity requires the definition of the
time parts and temperature parts of the flow rule and/or hardening rule.

Implementations
---------------
.. toctree::
   
   vp_flow/perzyna
   vp_flow/chaboche
   vp_flow/yaguchi

Class description
-----------------

.. doxygenclass:: neml::ViscoPlasticFlowRule
   :members:
   :undoc-members:
