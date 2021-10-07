Mukherjee creep
===============

Overview
--------

This implements the Mukherjee creep model [BMD1969]_

.. math::
   \dot{\varepsilon}^{cr} = A D_0 e^\frac{Q}{RT} \frac{\mu b}{k T} \left( \frac{\sigma_{eq}}{\mu}\right)^n

with :math:`A`, :math:`D_0`, :math:`Q`, and :math:`n` parameters and 
:math:`R` the gas constant, :math:`T` absolute temperature, :math:`\mu`
the temperature dependent shear modulus, :math:`b` a Burgers vector length, and
:math:`k` the Boltzmann constant.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``emodel``, :cpp:class:`neml::LinearElasticModel`, Elasticity model (for shear modulus), No
   ``A``, :code:`double`, Prefactor, No
   ``n``, :code:`double`, Stress exponent, No
   ``D0``, :code:`double`, Zero temperature lattice diffusivity, No
   ``Q``, :code:`double`, Activation energy for diffusivity, No
   ``b``, :code:`double`, Burgers vector, No
   ``k``, :code:`double`, Boltzmann constant, No
   ``R``, :code:`double`, Gas constant, No

Class description
-----------------

.. doxygenclass:: neml::MukherjeeCreep
   :members:
   :undoc-members:

