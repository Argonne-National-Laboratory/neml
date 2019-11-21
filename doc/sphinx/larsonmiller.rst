.. _larson-miller:

Larson Miller correlations
==========================

Overview
--------

This class implements a classical Larson-Miller [LM1952]_ correlation between 
stress, temperature, and rupture time.  These correlations are developed using
a large database of creep test results giving the applied stress, temperature, 
and the resulting time to rupture for a given material for many different
experimental conditions.

The Larson-Miller relation takes the form:

.. math::
   \log_{10}{\sigma_R} = f\left(\mathrm{LMP}\right)

where :math:`\sigma_R` is the time to rupture and :math:`f` is some generic
function (often a polynomial regression or a piecewise polynomial regression)
of the Larson-Miller parameter :math:`\mathrm{LMP}`, defined as

.. math::
   \mathrm{LMP} = T \left(C + \log_{10}{t_R} \right)

where :math:`T` is absolute temperature, :math:`C` is a constant fit to data,
and :math:`t_R` is the time to rupture.

NEML implements Larson-Miller relations by using a generic 
:ref:`interpolate <interpolate-functions>` object.  That is, the 
interpolate object provides the function :math:`f` in the above equation
relating the Larson-Miller parameter to the `log` of the stress to rupture.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``function``, :cpp:class:`neml::Interpolate`, Generic Larson-Miller relation, No
   ``lmr``, :c:type:`double`, Constant C from the Larson-Miller parameter, No
   ``tol``, :c:type:`double`, Solver tolerance, ``1.0e-6``
   ``miter``, :c:type:`int`, Maximum solver iterations, ``20``
   ``verbose``, :c:type:`bool`, Verbosity flag, ``false``

Class description
-----------------

.. doxygenclass:: neml::LarsonMillerRelation
   :members:
   :undoc-members:
