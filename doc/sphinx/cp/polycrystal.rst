Polycrystal models
==================

Overview
--------

This python file provides polycrystal homogenization models.
Here the goal is to homogenize the response of a large collection of single
crystal orientations (i.e. grains) with individual single crystal plastic
responses into a macroscale, effective plastic response.

Polycrystal homogenization models are simply :ref:`NEMLModel_ldi` objects
with that standard interface.  This means they can be used in any of the
NEML python drivers or in finite element analysis exactly like any other
NEML material model.  They can even be defined (somewhat tediously) in the
NEML XML file and saved for future use.  The only difference is that they
take as parameters a :ref:`single-crystal` object (itself a
:ref:`NEMLModel_ldi` object) and a list of orientations as input, instead
of some set of material parameters.

Implementations
---------------

.. toctree::
   polycrystal/taylor

Class description
-----------------

.. doxygenclass:: neml::PolycrystalModel
   :members:
   :undoc-members:
