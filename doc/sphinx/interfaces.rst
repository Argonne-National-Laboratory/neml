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

A :doc:`interfaces/NEMLModel` object describes a complete constitutive model of this type.
Currently, NEML provides only a single implementation of this abstract
base class, a :doc:`interfaces/NEMLModel_sd` which provides a suitable interface for
small strain kinematics where the user provides the small strain and 
NEML returns the stress.
The idea of having this abstract base class is that in the future NEML
can be expanded to other types of large strain material model updates.

.. toctree::

   interfaces/NEMLModel
