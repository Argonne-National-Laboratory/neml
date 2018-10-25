Elasticity models
=================

Overview
--------

Currently, NEML only has small strain linear elasticity models.
Such models implement the :ref:`LinearElasticModel <le-class>` interface, defined as:

.. math::
   
   \mathbf{\mathfrak{C}}, \mathbf{\mathfrak{S}}, E, \nu, \mu, K \leftarrow 
   \mathcal{E}\left(T\right)

where :math:`\mathbf{\mathfrak{C}}` is the stiffness tensor,
:math:`\mathbf{\mathfrak{S}}` is the compliance tensor,
:math:`E` is the Young's modulus, :math:`\nu` is the Poisson's ratio,
:math:`\mu` is the shear modulus, and :math:`K` is the bulk modulus.
The material model integration algorithms use the stiffness and compliance
directly.
The interface returns the scalar elastic properties for use in calculating
material properties, for example the shear modulus is used in 
calculating normalized activation energy for the :doc:`../interfaces/km_regime`.


Implementations
---------------

.. toctree::
   elastic/isotropic_elastic
   elastic/dummy_elastic


.. _le-class:

Class description
-----------------

.. doxygenclass:: neml::LinearElasticModel
   :members:
   :undoc-members:
