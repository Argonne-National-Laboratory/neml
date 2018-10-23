Rate independent nonassociative flow
====================================

Overview
--------
This interface implements a nonassociative flow rule where only the
flow function is associated to the yield surface through the relation

.. math::
   
   \mathbf{g}\left(\bm{\sigma}, \bm{\alpha}, T\right) = 
      \frac{\partial f}{\partial \bm{\sigma}}\left(\bm{\sigma}, 
         \mathbf{q}\left(\bm{\alpha}\right), T\right)

The hardening rule is left as a generic interface.
This type of model is much more common than a fully nonassociative model
where neither the flow rule or the hardening rule is associated to the
yield surface.
For example, classical Frederick-Armstrong hardening falls into this
category [FA2007]_.

Parameters
----------

========== =========================== ======================================= =======
Parameter  Object type                 Description                             Default
========== =========================== ======================================= =======
surface    YieldSurface                Yield surface interface                 No
hardening  NonAssociativeHardeningRule Nonassociative hardening rule interface No
========== =========================== ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::RateIndependentNonAssociativeHardening
   :members:
   :undoc-members:
