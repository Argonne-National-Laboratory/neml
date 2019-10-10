PowerLaw
========

Overview
--------

This class implements isotropic power-law perfect viscoplasticity with the
form:

.. math::
   \mathbf{D}^p = A \sigma_{eq}^{(n-1)} \mathbf{s}

   \mathbf{W}^p = \mathbf{0}

where

.. math::
   \sigma_{eq} = \sqrt{\frac{3}{2} \mathbf{s} : \mathbf{s}},

   \mathbf{s} = \bm{\sigma} - \frac{1}{3} \operatorname{tr}\left( \bm{\sigma} \right) \mathbf{I},

and :math:`A` and :math:`n` are parameters.  Note this is then the standard
:math:`J_2` flow rule.

The implementation does not require any internal variables.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``A``, :cpp:class:`neml::Interpolate`, Prefactor as a function of temperature, No
   ``n``, :cpp:class:`neml::Interpolate`, Rate sensitivity as a function of temperature, No

Class description
-----------------

.. doxygenclass:: neml::PowerLaw
   :members:
