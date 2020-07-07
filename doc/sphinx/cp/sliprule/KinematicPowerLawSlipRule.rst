KinematicPowerLawSlipRule
=========================

Overview
--------

This implements multi strength flow rule model of the form:

.. math::
   \dot{\gamma}_{g,i}=\dot{\gamma}_{0}\left\langle \frac{\left|\tau_{g,i}-\bar{\tau}_{g,i}^{back}\right|-\bar{\tau}_{g,i}^{iso}}{\bar{\tau}_{g,i}^{resistance}}\right\rangle ^{n}\operatorname{sign}\left(\tau_{g,i}-\bar{\tau}_{g,i}^{back}\right)

where :math:`\dot{\gamma}_0`, the reference slip rate, and :math:`n`, the rate sensitivity, are temperature-dependent parameters.
:ref:`slip-hardening` models provide the back, isotropic, and flow resistance strengths.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``backstrength``, :cpp:class:`neml::SlipHardening`, Back strength definition, No
   ``isostrength``, :cpp:class:`neml::SlipHardening`, Isotropic strength definition, No
   ``flowresistance``, :cpp:class:`neml::SlipHardening`, Flow resistance definition, No
   ``gamma0``, :cpp:class:`neml::Interpolate`, Reference slip rate, No
   ``n``, :cpp:class:`neml::Interpolate`, Rate sensitivity, No

Class description
-----------------

.. doxygenclass:: neml::KinematicPowerLawSlipRule
   :members:
   :undoc-members:
