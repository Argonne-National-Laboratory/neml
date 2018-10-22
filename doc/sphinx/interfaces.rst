Interfaces: types of material models
====================================

At the highest level a NEML model returns the stress-strain response of a
material.
This constitutive response depends on the input strain, temperature, and 
time along with a set of internal variables that the model maintains.
When calling a model from another code the user must maintain the
time series of strains, stresses, temperatures, times, history variables,
and (optionally) total strain energy and inelastic dissipation.
The user passes in the stress, history, and energy quantities at the 
previous time step, the strain, temperature, and time at both the previous
and next time steps, and NEML provides the updated stresses, history, 
and energies as output.

A :doc:`interfaces/classes/NEMLModel` object describes a complete constitutive model of this type.
Currently, NEML provides only a single implementation of this abstract
base class, a :doc:`interfaces/classes/NEMLModel_sd` which provides a suitable interface for
small strain kinematics where the user provides the small strain and 
NEML returns the stress.
The idea of having this abstract base class is that in the future NEML
can be expanded to other types of large strain material model updates.

The :doc:`interfaces/classes/NEMLModel_sd` object provides the interface:

.. math::

   \bm{\sigma}_{n+1}, \mathbf{h}_{n+1}, \bm{\mathfrak{A}}_{n+1}, u_{n+1}, p_{n+1} \leftarrow
   \mathcal{M}\left( 
   \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n},
   T_{n+1}, T_{n},
   t_{n+1}, t_{n},
   \bm{\sigma}_{n},
   \mathbf{h}_{n},
   u_n, p_n
   \right).

Here :math:`n` indicates values at the previous time step and :math:`n+1` values
at the next time step.
The quantities are stress (:math:`\bm{\sigma}`), strain (:math:`\bm{\varepsilon}`),
the vector of history variables (:math:`\mathbf{h}`), strain energy (:math:`u`)
dissipated work (:math:`p`), temperature (:math:`T`), time (:math:`t`), and 
the algorithmic tangent (:math:`\mathbf{\mathfrak{A}}`).
For the small strain interface the appropriate algorithmic tangent is

.. math::
   \mathbf{\mathfrak{A}}_{n+1} = \frac{d \bm{\sigma}_{n+1}}{d \bm{\varepsilon}_{n+1}}.

The following sections describe the basic material model implemented from
this generic interfaces.
Another section of the model details continuum damage models, which also
use this same interface.

.. toctree::
   :maxdepth: 2

   interfaces/lelastic
   interfaces/perfect
   interfaces/rate_independent
   interfaces/general_integrator
   interfaces/creep_plasticity
   interfaces/km_regime

.. toctree::
   :hidden:

   interfaces/classes/NEMLModel
   interfaces/classes/NEMLModel_sd
