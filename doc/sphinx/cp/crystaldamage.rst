Crystal plasticity damage models
================================

NEML views crystal plasticity continuum damage as an alternative kinematic
model to the one described in the :ref:`main documentation <cp-formulation>`. 
The model takes a projection operator :math:`\mathbf{P}`, defined below, 
which is a function of some set of internal variables describing damage
evolution, :math:`\mathbf{d}`, the current stress state :math:`\bm{\sigma}`, and
uses the projection to project damage onto the elasticity tensor used
in the crystal plasticity stress update formulation.

In the most general case, consider the base stress updated provided by a model as a function of time :math:`\bm{\sigma}^\prime`.  The damage projection will
alter this stress history as:

.. math::

   \bm{\sigma}=\boldsymbol{P}:\bm{\sigma}^{\prime}

so that the new, modified stress rate accounting for damage becomes

.. math::

   \dot{\bm{\sigma}}=\dot{\boldsymbol{P}}:\bm{\sigma}^{\prime}+\boldsymbol{P}:\dot{\bm{\sigma}}^{\prime}

   \dot{\bm{\sigma}}=\dot{\boldsymbol{P}}:\boldsymbol{P}^{-1}:\bm{\sigma}+\boldsymbol{P}:\dot{\bm{\sigma}}^{\prime}

in the case where the damage evolves slowly compared to the evolution of the
stress we can approximate this as

.. math::

   \dot{\bm{\sigma}}=\boldsymbol{P}:\dot{\bm{\sigma}}^{\prime}

which is the form currently implemented in NEML.

Starting then from the basic stress update kinematics defined in the :ref:`main documentation <cp-formulation>`, the
modified stress rate equation becomes

.. math::

   \dot{\sigma}_{ij}=P_{ijkl}C_{klmn}\left(D_{mn}-D_{mn}^{p}\right)-\sigma_{ik}\Omega_{kj}^{\star}+\Omega_{ik}^{\star}\sigma_{kj}

NEML implements this modified stress update as a new :any:`KinematicModel`.  A
:any:`CrystalDamageModel` provides the definition of the projection operator, 
including the history evolution rate equations for any internal variables used
in defining the projection.

.. toctree::
   :maxdepth: 1

   CrystalDamageModel

DamagedStandardKinematicModel
-----------------------------

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``emodel``, :cpp:class:`neml::LinearElasticModel`, Elasticity tensor, No
   ``imodel``, :cpp:class:`neml::InelasticModel`, Definition of the plastic deformation rate, No
   ``dmodel``, :cpp:class:`neml::CrystalDamageModel`, Definition of the damage projection and associated internal variables, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::DamagedStandardKinematicModel
   :members:
   :undoc-members:
