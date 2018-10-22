NEMLModel_sd
============

Overview
--------

The NEMLModel_sd object provides the interface:

.. math::

   \bm{\sigma}_{n+1}, \bm{\alpha}_{n+1}, \bm{\mathfrak{A}}_{n+1}, u_{n+1}, p_{n+1} \leftarrow
   \mathcal{M}\left( 
   \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n},
   T_{n+1}, T_{n},
   t_{n+1}, t_{n},
   \bm{\sigma}_{n},
   \bm{\alpha}_{n},
   u_n, p_n
   \right).

Here :math:`n` indicates values at the previous time step and :math:`n+1` values
at the next time step.
The quantities are stress (:math:`\bm{\sigma}`), strain (:math:`\bm{\varepsilon}`),
the vector of history variables (:math:`\bm{\alpha}`), strain energy (:math:`u`)
dissipated work (:math:`p`), temperature (:math:`T`), time (:math:`t`), and 
the algorithmic tangent (:math:`\mathbf{\mathfrak{A}}`).
For the small strain interface the appropriate algorithmic tangent is

.. math::
   \mathbf{\mathfrak{A}}_{n+1} = \frac{d \bm{\sigma}_{n+1}}{d \bm{\varepsilon}_{n+1}}.

The following sections describe the basic material model implemented from
this generic interfaces.
Another section of the model details continuum damage models, which also
use this same interface.

Implementations
---------------

.. toctree::

   lelastic
   perfect
   rate_independent
   general_integrator
   creep_plasticity
   km_regime

Parameters
----------

========= ===================== ======================================= =======
Parameter Object type           Description                             Default
========= ===================== ======================================= =======
emodel    LinearElasticModel    Temperature dependent elastic constants No
alpha     Interpolate           Temperature dependent instantaneous CTE 0.0
========= ===================== ======================================= =======

Class documentation
-------------------

.. doxygenclass:: neml::NEMLModel_sd
   :members:
   :undoc-members:
