.. _slip-hardening:

SlipHardening
=============

Overview
--------

These objects provide the slip system strength required to complete
the constitutive description when using a SlipStrengthSlipRule.
The implementation assumes that the slip system strengths are functions only
of the set of history variables and temperature.
Additionally, these objects provide the definition of the history evolution rate equations and, ultimately, the definition of the model internal variables.

The interface is then:

.. math::
   \bar{\tau}_{g,i}, \frac{\partial \bar{\tau}_{g,i}}{\partial \mathbf{h}} \leftarrow \mathcal{T}\left(\mathbf{h}, \bm{\alpha}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, \bm{\alpha}, T \right).

Implementations
---------------

.. toctree::
   :maxdepth: 2

   sliphardening/SlipSingleHardening
   sliphardening/FixedStrengthHardening
   sliphardening/GeneralLinearHardening
   sliphardening/SimpleLinearHardening
   sliphardening/VocePerSystemHardening
   sliphardening/FASlipHardening
   sliphardening/LANLTiModel


Class description
-----------------

.. doxygenclass:: neml::SlipHardening
   :members:
   :undoc-members:
