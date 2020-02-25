SumSlipSingleStrengthHardening
==============================

Overview
--------

This model sums the contributions of multiple :ref:`slipsinglestrengthhardening` models together, i.e.

.. math::
   \bar{\tau} = \sum_{i=1}^{n_{model}}\tilde{\tau}_{i}+\tau_{0,i}+\tau_{nye,i}

Note then that all parts of each model are summed, even the static strengths.

This class will rename the internal variables of the individual hardening
models to avoid overlap in the History object.

Class description
-----------------

.. doxygenclass:: neml::SumSlipSingleStrengthHardening
   :members:
   :undoc-members:
