Rate independent associative flow
=================================

Overview
--------
This interface implements an associative flow rule where both the
flow function and the hardening rule are *associated* to the yield surface
by the functional relations

.. math::
   
   \mathbf{g}\left(\bm{\sigma}, \bm{\alpha}, T\right) = 
      \frac{\partial f}{\partial \bm{\sigma}}\left(\bm{\sigma}, 
         \mathbf{q}\left(\bm{\alpha}\right), T\right)

   \mathbf{h}\left(\bm{\sigma}, \bm{\alpha}, T\right) = 
      \frac{\partial f}{\partial \mathbf{q}}\left(\bm{\sigma}, 
         \mathbf{q}\left(\bm{\alpha}\right), T\right)

These quantities have all been defined previous, except for the
function :math:`\mathbf{q}`.
This function maps the "strain-like" set of history vectors to the
"stress-like" set of internal variables that enter the yield surface [SH1997]_.
These "stress-like" internal variables are most commonly an isotropic 
expansion/contraction of the yield surface and a kinematic backstress
shifting the yield surface in space.

A fully-associative flow rule of this type results in a model with 
favorable theoretical and numerical properties.
For example, these models all obey Drucker's postulate [D1959]_ and will 
have symmetric algorithmic tangents.

The hardening rule and yield surface are both defined with
separate interfaces.

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
surface    YieldSurface              Yield surface interface                 No
hardening  HardeningRule             Hardening rule interface                No
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::RateIndependentAssociativeFlow
   :members:
   :undoc-members:
