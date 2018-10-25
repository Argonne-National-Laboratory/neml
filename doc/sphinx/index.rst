NEML -- Nuclear Engineering Material model Library
==================================================

NEML is a modular library for creating constitutive models for 
structural materials.
While it was originally focused on materials slated for use in
high temperature nuclear reactors, it covers a wide range of
constitutive models used for all types of structural materials.

NEML can be tied into any finite element program that can use a C, C++, 
or Fortran interface.
This distribution provides an Abaqus UMAT driver as an example of how
to accomplish this.
Additionally, NEML comes with a Python interface that can be used to
develop, debug, and test material models.
The neml python module includes helper routines to simulate various 
common experimental conditions using NEML material models.

Go :doc:`here <started>` to started with NEML by compiling it, running examples,
and learn how to link it to your finite element code or 
go :doc:`here <about>` to learn more about the structure of the library.

.. toctree::
   :maxdepth: 1

   about
   started
   tutorial
   interpolate
   interfaces
   elasticity
   ri_flow
   general_flow
   vp_flow
   surfaces
   hardening
   creep
   damage
   python
   advanced
   references

