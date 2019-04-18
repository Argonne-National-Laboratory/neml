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

NEML aims to provide a generic interface to external finite element software 
through the  :doc:`interfaces/NEMLModel` class.
This class then provides routines that provide access to the constitutive
response for different types of requested stress updates.

Subclasses of this class implement a certain type of stress update.
For example, :doc:`interfaces/NEMLModel_sd` implements small strain
stress updates natively.
However, these subclasses also provide wrapper interfaces around this
native stress update to accommodate requests for other stress updates.
For example, the :doc:`interfaces/NEMLModel_sd` class provides an incremental
large deformations stress update using an objective stress rate.

.. toctree::
   :maxdepth: 5

   interfaces/NEMLModel
