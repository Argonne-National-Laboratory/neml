PTRTwinReorientation
====================

Overview
--------

This class implements the Predominate Twin Reorientation model [TLK1991]_, which
reorients the crystal according to some twin geometry once the integrated
slip on the twin system(s) reaches a critical value.

Specifically, the model reorients the crystal using the twin transformation
for the twin system that reaches a given twinning fraction, defined as

.. math::
   f_i = \frac{\gamma_i}{s_i}

where :math:`\gamma_i` is the integrated slip on twin system :math:`i` and 
:math:`s_i` is the characteristic shear for that twin system.  The
post process expects the base crystal model to provide the integrated
slip rates on each twin system and the post processor itself maintains
internal variables representing the twin fractions as well as a flag
for if the crystal as a whole underwent the twin transformation.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``threshold``, :cpp:class:`neml::Interpolate`, Twin fraction threshold, N
   ``prefix``, :code:`std::string`, Integrated slip internal variable prefix, :code:`"slip"`

Class description
-----------------

.. doxygenclass:: neml::PTRTwinReorientation
   :members:
   :undoc-members:
