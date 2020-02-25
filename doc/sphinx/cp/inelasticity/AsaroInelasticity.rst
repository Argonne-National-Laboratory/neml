AsaroInelasticity
=================

Overview
--------

This class implements the standard crystal plasticity plastic deformation
rate presented in [A1983]_.
The plastic deformation rate is the sum of simple shear deformations on 
each individual slip direction and plane:

.. math::
   \mathbf{L}^p = \sum_{g=1}^{n_{group}}\sum_{i=1}^{n_{slip}}\dot{\gamma}_{g,i}\left(\mathbf{d}_{g,i}\otimes\mathbf{n}_{g,i}\right)

where :math:`\mathbf{d}_{g,i}` are the slip directions and :math:`\mathbf{n}_{g,i}` are
the slip normals `in the current frame` :math:`\Omega`.
:math:`\mathbf{D}^p` and :math:`\mathbf{W}^p` are the symmetric and skew parts
of this definition.

The slip rate, :math:`\dot{\gamma}_{g,i}`, is a function of stress, 
temperature, and the internal variables.  A separate class implements 
different models for the slip rate:

.. toctree::
   :maxdepth: 1

   ../SlipRule

The implementation also defers the definition of the internal variables to
the :doc:`../SlipRule` object.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``rule``, :cpp:class:`neml::SlipRule`, Slip rate and history definition, No

Class description
-----------------

.. doxygenclass:: neml::AsaroInelasticity
   :members:
   :undoc-members:
