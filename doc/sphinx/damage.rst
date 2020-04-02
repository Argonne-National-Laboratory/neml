.. _damage:

Damage mechanics
================

Overview
--------

NEML implements damage mechanics with metamodel place on top of an existing 
(small strain) base material model.
The metamodel supplements the base material history variables with a 
vector of damage variables :math:`\bm{\omega}`.
In the most generic case, these damage variables somehow modify the base 
material's stress update.
The damage model must also provide the history evolution equations describing
the evolution of damage.
The interface is then

.. math::
   \bm{\sigma}_{n+1}, \bm{\alpha}_{n+1}, \bm{\omega}_{n+1}, \bm{\mathfrak{A}}_{n+1}, u_{n+1}, p_{n+1} \leftarrow
   \mathcal{D}\left( 
   \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n},
   T_{n+1}, T_{n},
   t_{n+1}, t_{n},
   \bm{\sigma}_{n},
   \bm{\alpha}_{n},
   \bm{\alpha}_{n},
   u_n, p_n
   \right)

which is identical to the interface for the 
:doc:`base small strain material models <interfaces/NEMLModel_sd>`, 
except now with the addition of some dependence on the 
damage variables :math:`\bm{\omega}`.

Implementations
---------------

The generic interface described here could have any number of damage 
variables and could modify the stress with those variables in any way.
The more standard damage model, depending on a single scalar damage 
variable is implemented as a NEMLScalarDamagedModel_sd.

.. toctree::
   damage/scalar
   damage/combined

Class description
-----------------

.. doxygenclass:: neml::NEMLDamagedModel_sd
   :members:
   :undoc-members:

