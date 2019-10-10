SlipHardening
=============

Overview
--------

These objects provide the slip system strength required to complete
the constitutive description when using a SlipStrengthSlipRule.
The implementation assumes that the slip system strengths is a function only
of the set of history variables and temperature.
Additionally, these objects provide the definition of the history evolution rate equations and, ultimately, the definition of the model internal variables.

The interface is then:

.. math::
   \bar{\tau}_{g,i}, \frac{\partial \bar{\tau}_{g,i}}{\partial \mathbf{h}} \leftarrow \mathcal{T}\left(\mathbf{h}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, T \right).

Implementations
---------------

.. toctree::
   :maxdepth: 2

   sliphardening/SlipSingleHardening


Class description
-----------------

.. doxygenclass:: neml::SlipHardening
   :members:
