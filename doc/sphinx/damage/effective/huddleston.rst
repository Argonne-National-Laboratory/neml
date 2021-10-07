Huddleston effective stress
===========================

Overview
--------

The effective stress is that of Huddleston [H1985]_, defined as

.. math::
   \sigma_e = \sigma_{vm} \exp \left( -b \left( \frac{I_1}{S_s} - 1\right) \right)

with 

.. math::
   \sigma_{vm} = \sqrt{\frac{\left(\sigma_1 - \sigma_2\right)^{2} + \left(\sigma_2 - \sigma_3\right)^{2} + \left(\sigma_3 - \sigma_1\right)^{2}}{2}}

the von Mises stress, 

.. math::
   I_1 = \sigma_1 + \sigma_2 + \sigma_3

the first stress invariant and

.. math::
   S_s = \sqrt{\sigma_{1}^{2} + \sigma_{2}^{2} + \sigma_{3}^{2}} 

all in terms of the maximum principal stresses :math:`\sigma_1`, :math:`\sigma_2`, and :math:`\sigma_3`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``b``, :code:`double`, Huddleston parameter, No

Class description
-----------------

.. doxygenclass:: neml::HuddlestonEffectiveStress
   :members:
   :undoc-members:
