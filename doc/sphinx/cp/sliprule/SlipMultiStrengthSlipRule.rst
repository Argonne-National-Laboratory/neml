SlipMultiStrengthSlipRule
=========================

Overview
--------

These objects provide a relation between the stress, history, and temperature
and the slip rate on each individual slip system where the slip rate is related to the resolved shear stress on the system

.. math::
   \tau_{g,i} = \bm{\sigma} : \left(\mathbf{d}_{g,i}\otimes\mathbf{n}_{g,i}\right)

where :math:`\mathbf{d}_{g,i}` is the slip direction for group `g`, system `i` `in the current coordinates` and :math:`\mathbf{n}_{g,i}` is similarly the slip system normal.  The interface used is:

.. math::
   \dot{\gamma}_{g,i}, \frac{\partial \dot{\gamma}_{g,i}}{\partial \tau_{g,i}}, \frac{\partial \dot{\gamma}_{g,i}}{\partial \bar{\tau}_{g,i}^{j}} \leftarrow \mathcal{G}\left( \tau_{g,i}, \bar{\tau}_{g,i}, \bm{\alpha}, T \right)

   \dot{\mathbf{h}}, \frac{\partial \dot{\mathbf{h}}}{\partial \bm{\sigma}}, \frac{\partial \dot{\mathbf{h}}}{\partial \mathbf{h}} \leftarrow \mathcal{H}\left(\bm{\sigma}, \mathbf{h}, \bm{\alpha}, T \right)

where :math:`g` indicates the slip group, :math:`i` indicates the system within the group, and :math:`\bar{\tau}_{g,i}^j` is a collection of slip system strengths, defined by SlipHardening models:

.. toctree::
   :maxdepth: 2

   ../SlipHardening

The definition of the history evolution is left to the SlipHardening models.

The difference between this class and :ref:`slip-strength-slip-rule` is that the slip system flow is proportional to multiple slip strength models, for example an isotropic and a kinematic strength, instead of a 
single flow strength model.

Implementations
---------------

.. toctree::
   :maxdepth: 1

   KinematicPowerLawSlipRule

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``strength``, :c:type:`std::vector<`:cpp:class:`neml::SlipHardening`:c:type:`>`, List of slip hardening rules, No

Class description
-----------------

.. doxygenclass:: neml::SlipMultiStrengthSlipRule
   :members:
   :undoc-members:
