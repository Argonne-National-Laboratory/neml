Rate independent flow rule
==========================

Overview
--------

A rate independent flow rule provides the interface

.. math::
   
   f, \frac{\partial f}{\partial \bm{\sigma}}, 
   \mathbf{g}, \frac{\partial \mathbf{g}}{\partial \bm{\sigma}},
      \frac{\partial \mathbf{g}}{\partial \bm{\alpha}},
   \mathbf{h}, \frac{\partial \mathbf{h}}{\partial \bm{\sigma}},
      \frac{\partial \mathbf{h}}{\partial \bm{\alpha}}
   \leftarrow
   \mathcal{F}\left(\bm{\sigma}, \bm{\alpha}, T \right).
      
:math:`f` is a yield surface, :math:`\mathbf{g}` is a
flow function defining the direction of plastic flow, and :math:`\mathbf{h}` 
is a history evolution function.
This interface is used to define a :doc:`interfaces/rate_independent` model.

The base interface is entirely abstract.
It maintains a set of history variables set by the specific implementation.

Implementations
---------------

.. toctree::

   ri_flow/associative
   ri_flow/nonassociative

Class description
-----------------

.. doxygenclass:: neml::RateIndependentFlowRule
   :members:
   :undoc-members:
