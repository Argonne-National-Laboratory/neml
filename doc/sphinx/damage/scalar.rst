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
      \bm{\sigma}_{n} / (1 - \omega_{n}),
      \bm{\alpha}_{n+1},
      \bm{\alpha}_{n},
      u_n, p_n
      \right)

where :math:`\omega` is the damage variable and :math:`\bm{\sigma}` is the 
base material stress update.
This object defers damage evolution to another interface.

The damage model maintains the set of history variables from the base 
material plus one additional history variable for the damage.

.. note::
   The scalar damage model passes in the modified stress :math:`\bm{\sigma} / (1 - \omega)` to the base stress update model in addition to modifying the stress update formula as shown in the above equation.

.. warning::
   The model also passes the modified stress :math:`\bm{\sigma} / (1-\omega)` to the damage update equation.  That is, the stress passed into these functions is the modified effective stress, not the actual external stress.  This means that the damage equations implemented in NEML vary slightly from the correpsonding literature equations working with the unmodified stress directly.

Implementations
---------------

.. toctree::
   classical
   modular
   standard
   larsonmiller
   work

Class description
-----------------

.. doxygenclass:: neml::NEMLScalarDamagedModel_sd
   :members:
   :undoc-members:

