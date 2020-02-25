PlasticSlipHardening
====================

Overview
--------

This class performs yet another simplification of the SlipSingleStrengthHardening model so that the scalar history variable evolves only as a function of the variable itself, temperature, and the absolute sum of the slip rates on all the systems.  That is

.. math::
   \dot{\tilde{\tau}} = f\left(\tilde{\tau}, T \right) \sum_{g=1}^{n_{groups}}\sum_{i=1}^{n_{slip}}\left|\dot{\gamma}_{g,i}\right|

The interface defines the function :math:`f`, its partial derivative with
respect to the history variable, and the static strength:

.. math::
   f, \frac{\partial f}{\partial \tilde{\tau}}, \tau_0 \leftarrow \mathcal{P}\left(\tilde{\tau}, T \right)


Implementations
---------------

.. toctree::
   :maxdepth: 1

   VoceSlipHardening
   LinearSlipHardening


Class description
-----------------

.. doxygenclass:: neml::PlasticSlipHardening
   :members:
   :undoc-members:
