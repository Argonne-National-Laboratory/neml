Single crystal model postprocessors
===================================

This class provides means to alter the state of a crystal model
when certain conditions are met, either based on the internal state of
the crystal model or in response to a stimulus from an external object
(like a full field solver).
These actions are not integrated into the implicit time integration and so
occur explicitly after successfully stress updates.  This may lead to 
decreased numerical convergence for the time step right after the 
action triggers.

.. toctree::
   :maxdepth: 1

   postprocessors/PTRTwinReorientation
