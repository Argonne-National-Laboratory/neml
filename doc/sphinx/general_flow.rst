General flow rule
=================

Overview
--------

A general flow rule provides the interface

.. math::
   
   \dot{\bm{\sigma}}, 
      \frac{\partial\dot{\bm{\sigma}}}{\partial\bm{\sigma}},
      \frac{\partial\dot{\bm{\sigma}}}{\partial\bm{\alpha}},
      \frac{\partial\dot{\bm{\sigma}}}{\partial\dot{\bm{\varepsilon}}},
   \dot{\bm{\alpha}}, 
      \frac{\partial\dot{\bm{\alpha}}}{\partial\bm{\sigma}},
      \frac{\partial\dot{\bm{\alpha}}}{\partial\bm{\alpha}},
      \frac{\partial\dot{\bm{\alpha}}}{\partial\dot{\bm{\varepsilon}}}
   \leftarrow \mathcal{G}
   \left(\bm{\sigma},\bm{\alpha},\dot{\bm{\varepsilon}}, T, \dot{T} \right).

So it provides a generic stress rate and history evolution rate as a function
of the current state.

The base interface is entirely abstract.
It maintains a set of history variables set by the specific implementation.

Implementations
---------------

.. toctree::

   general_flow/viscoplastic
   general_flow/walkerkrempl

Class description
-----------------

.. doxygenclass:: neml::GeneralFlowRule
   :members:
   :undoc-members:
