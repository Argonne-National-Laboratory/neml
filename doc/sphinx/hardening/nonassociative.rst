Nonassociative hardening rules
==============================

Overview
--------

This object provides the interface for a generic nonassociative hardening rule.
The interface is:

.. math::
   \mathbf{q}, \mathbf{h}_{\gamma}, \mathbf{h}_{t}, \mathbf{h}_{T}
   \leftarrow
   \mathcal{H}\left(\bm{\alpha}, T \right).

Here :math:`\mathbf{q}` is a map between the history variables 
:math:`\bm{\alpha}` 
and the "stress-like" hardening variables needed for the yield stress.
Usually these will be an isotropic hardening variable :math:`Q` and a 
kinematic backstress :math:`\mathbf{X}`.
The functions :math:`\mathbf{h}_\gamma`, :math:`\mathbf{h}_t`, 
and :math:`\mathbf{h}_T` are the parts of the history evolution rate
that are proportional to the scalar inelastic strain rate, the temperature
rate, and time, respectively.
For more information see :doc:`../vp_flow/chaboche` and 
:doc:`../ri_flow/nonassociative`.

Implementations
---------------

.. toctree::
   non/chaboche
   non/special_chaboche

Class description
-----------------

.. doxygenclass:: neml::NonAssociativeHardening
   :members:
   :undoc-members:
