Scalar damage models
====================

Overview
--------

This object implements a damage model that uses a single damage variable 
to degrade the stress of a base material model.
It implements the stress update function

.. math::
   \bm{\sigma}_{n+1}^\prime = (1 - \omega_{n+1}) 
      \bm{\sigma}\left( 
      \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n},
      T_{n+1}, T_{n},
      t_{n+1}, t_{n},
      \bm{\sigma}_{n},
      \bm{\alpha}_{n},
      \bm{\alpha}_{n},
      u_n, p_n
      \right)

where :math:`\omega` is the damage variable and :math:`\bm{\sigma}` is the 
base material stress update.
This object defers damage evolution to another interface.

The damage model maintains the set of history variables from the base 
material plus one additional history variable for the damage.

Implementations
---------------

.. toctree::
   classical
   modular
   standard

Class description
-----------------

.. doxygenclass:: neml::NEMLScalarDamagedModel_sd
   :members:
   :undoc-members:

