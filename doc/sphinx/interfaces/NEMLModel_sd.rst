NEMLModel_sd
============

Overview
--------

The NEMLModel_sd object natively implements the :ref:`small strain <strain-interfaces>` stress update interface.

It accommodates the :ref:`large strain incremental <strain-interfaces>` stress update interface using the Treusdell objective stress rate of the form:

.. math::
   \dot{\bm{\sigma}}=\hat{\dot{\bm{\sigma}}}-\bm{\sigma}\cdot\bm{L}^{T}-\bm{L}\cdot\bm{\sigma}+\operatorname{tr}\left(\bm{L}\right)\bm{\sigma}

where :math:`\bm{\sigma}` is the Cauchy stress and :math:`\hat{\dot{\bm{\sigma}}}` is the small strain stress rate implied by the small strain kinematics
update interface.
The update calculates the consistent tangents :math:`\mathbf{\mathfrak{A}}`
:math:`\mathbf{\mathfrak{B}}` exactly and provides a helper routine
to recombine these symmetric and skew parts into the full derivative
with respect to the spatial velocity gradient.

.. caution::
   The current Treusdell objective integration does not advect the
   material history variables.
   This means the integration of material models with vector or tensor
   history variables, such as backstresses, will be inaccurate for
   situations requiring large rotations.
   This limitation will be removed in future version of NEML.

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

Class description
-----------------

.. doxygenclass:: neml::NEMLModel_sd
   :members:
   :undoc-members:
