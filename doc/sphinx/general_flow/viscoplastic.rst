Viscoplastic general flow rule
==============================

Overview
--------

This class specializes the :doc:`../general_flow` interface to match the
standard definition of a generic viscoplastic flow rule.  It defines 
the stress and hardening rates as

.. math::
   \dot{\bm{\sigma}} = \mathbf{\mathfrak{C}}:\left(\dot{\bm{\varepsilon}}-\mathbf{g}_{\gamma}\dot{\gamma}-\mathbf{g}_{T}\dot{T}-\mathbf{g}_{t}\right)

   \dot{\bm{\alpha}} = \mathbf{h}_{\gamma}\dot{\gamma}+\mathbf{h}_{T}\dot{T}+\mathbf{h}_{t}.

A :doc:`../elasticity` interface defines the stiffness tensor and a viscoplastic
flow rule interface defines the scalar flow rate
:math:`\dot{\gamma}`,
the flow functions 
:math:`\mathbf{g}_\gamma`, :math:`\mathbf{g}_T`, and :math:`\mathbf{g}_t` 
and the hardening functions
:math:`\mathbf{h}_\gamma`, :math:`\mathbf{h}_T`, and :math:`\mathbf{h}_t`.

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
elastic    LinearElasticModel        Elasticity model                        No
flow       ViscoPlasticFlowRule      Viscoplastic flow rule interface        No
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::TVPFlowRule
   :members:
   :undoc-members:
