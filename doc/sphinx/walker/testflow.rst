Test wrapped flow rule
======================

Overview
--------

As a test of the :doc:`wrapper` system, this class implements a power law flow rule with hard-coded isotropic hardening:

.. math::
   y = \dot{\varepsilon}_0 \left\langle \frac{\sqrt{3/2} \operatorname{dev} \bm{\sigma} - h}{D} \right\rangle^{n}

with

.. math::
   \dot{h} = K

and :math:`h(0) = \sigma_0`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``eps0``, :code:`double`, Reference strain rate, No
   ``D``, :code:`double`, Drag stress, No
   ``n``, :code:`double`, Rate sensitivity exponent, No
   ``s0``, :code:`double`, Initial hardening strength, No
   ``K``, :code:`double`, Hardening modulus, No

Class description
------------------

.. doxygenclass:: neml::TestFlowRule
   :members:
   :undoc-members:
