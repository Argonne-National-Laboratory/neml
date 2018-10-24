Perzyna viscoplastic flow rule
==============================

Overview
--------

The Perzyna model is an associative viscoplastic model [P1966]_.
It is defined by:

.. math::

   \dot{\gamma} = g\left(\left\langle f\left(\bm{\sigma}, \mathbf{q}\left(\bm{\alpha}\right), T\right)\right\rangle\right)

   \mathbf{g}_{\gamma} = \frac{\partial f}{\partial \bm{\sigma}} 
      \left( \bm{\sigma}, \mathbf{q}\left(\bm{\alpha}\right), T  \right)

   \mathbf{h}_{\gamma} = \frac{\partial f}{\partial \mathbf{q}} 
      \left( \bm{\sigma}, \mathbf{q}\left(\bm{\alpha}\right), T  \right)

where :math:`f` is a flow surface, a hardening interface provides the
:math:`\mathbf{q}` function, and :math:`g` is a :ref:`rate-function`.
The notation :math:`\left\langle \right\rangle` indicates the 
`Macaulay brackets <https://en.wikipedia.org/wiki/Macaulay_brackets>`_.
The implementation uses the default zero values of the time and temperature rate
contributions and so the model evolves only as a function of the inelastic
strain rate.

The model maintains the set of history variables defined by the hardening
model.


Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``surface``, :cpp:class:`neml::YieldSurface`, Flow surface interface, No
   ``hardening``, :cpp:class:`neml::HardeningRule`, Hardening rule interface, No
   ``g``, :cpp:class:`neml::GFlow`, Rate sensitivity function, No

Class description
-----------------

.. doxygenclass:: neml::PerzynaFlowRule
   :members:
   :undoc-members:

.. _rate-function:

Rate function
-------------

These functions define the rate sensitivity of the Peryna flow model.  They have the form 
:math:`g\left(f, T\right)` where `f` is the current value of the flow surface.

.. doxygenclass:: neml::GFlow
   :members:
   :undoc-members:

Power law
^^^^^^^^^

This rate function implements simple power law

.. math::

   g\left(f, T\right) = \left(\frac{f}{\eta}\right)^n

for some temperature-dependent rate sensitivity exponent :math:`n` and fluidity 
:math:`\eta`.

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``n``, :cpp:class:`neml::Interpolate`, Rate sensitivity exponent, No
   ``eta``, :cpp:class:`neml::Interpolate`, Fluidity, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::GPowerLaw
   :members:
   :undoc-members:
