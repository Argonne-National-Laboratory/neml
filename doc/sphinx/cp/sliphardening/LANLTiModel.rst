LANLTiModel
===========

Overview
--------

Forest dislocation density based hardening slip model, used to physically describe of 
dislocation dominant plastic deformation.

.. math::
   \tau_{i}=\tau_{0,i}^{s}+X^{i}b^{i}\mu^{i}\sqrt{\rho_{for}^{i}}

The forest dislocation density evloves with strain and varies with dislication trapping
and annihilation rates of each slip system.

.. math::
   \dot{\rho_{for}^{i}}=k_{1}^{i}\sqrt{\rho_{for}^{i}}-k_{2}^{i}\rho_{for}^{i}

A pseudo-slip hardening model used to describe the resistance of the propagation of twinning.

.. math::
   \tau_{i}=\tau_{0,i}^{t}+\mu{i}\sum_{j=1}^{n_{slip}} C^{ij}b^{i}b^{j}\rho_{for}^{j}

where :math:`i` is the unrolled indices corresponding to slip group.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``tau_0``, :code:`std::vector<double>`, Initial strengths, N
   ``C_st``, :cpp:class:`neml::SquareMatrix`, Twin-slip interaction matrix, N
   ``mu``, :code:`std::vector<double>`, Elastic modulus on the system, N
   ``k1``, :code:`std::vector<double>`, Coefficient for trapping of dislocation, N
   ``k2``, :code:`std::vector<double>`, Coefficient for annihilation of dislocation, N
   ``X``, :code:`double`, Dislocation interaction parameter, 0.9
   ``inivalue``, :code:`double`, Initial values of internal variables, 1.0e-6

Class description
-----------------

.. doxygenclass:: neml::LANLTiModel
   :members:
   :undoc-members:
