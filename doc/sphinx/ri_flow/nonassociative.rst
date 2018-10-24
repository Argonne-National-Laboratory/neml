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

The hardening rule is left as a generic interface providing the 
rate of history variable evolution proportional to the equivalent 
plastic strain rate, :math:`\mathbf{h}_\gamma`.
Note that rate independent models cannot use hardening proportional to
time or temperature rate or else the model will not be rate independent!

This type of model is much more common than a fully nonassociative model
where neither the flow rule or the hardening rule is associated to the
yield surface.
For example, classical Frederick-Armstrong hardening falls into this
category [FA2007]_.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``surface``, :cpp:class:`neml::YieldSurface`, Yield surface interface, No
   ``hardening``, :cpp:class:`neml::NonAssociativeHardening`, Nonassociative hardening rule interface, No

Class description
-----------------

.. doxygenclass:: neml::RateIndependentNonAssociativeHardening
   :members:
   :undoc-members:
