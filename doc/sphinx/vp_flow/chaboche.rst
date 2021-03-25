Chaboche viscoplastic flow rule
===============================

Overview
--------

The Chaboche viscoplastic flow model uses an associative flow function but
nonassociative hardening.
Traditionally the nonassociative hardening model used is the ChabocheHardening
model, but this implementation can use any nonassociative hardening rule.
The model is defined by:

.. math::

   \dot{\gamma} = \dot{\gamma}_0 \sqrt{\frac{3}{2}} \left\langle \frac{f\left(\bm{\sigma}, \mathbf{q}\left(\bm{\alpha}\right), T\right)}{\sqrt{2/3}\eta\left(\bar{\varepsilon}_{vp}, T\right)}\right\rangle^n

   \mathbf{g}_{\gamma} = \frac{\partial f}{\partial \bm{\sigma}}
      \left( \bm{\sigma}, \mathbf{q}\left(\bm{\alpha}\right), T  \right)

The time and temperature rate contributions to the flow function are zero.
The rate sensitivity exponent is generally temperature dependent; the prefactor
:math:`\dot_{\gamma}_0` can be temperature dependent; the
:ref:`fluidity-model` :math:`\eta` can depend on both temperature and inelastic strain.
The hardening model is defined by a NonassociativeHardening model.
The time and temperature rate contributions can be non-zero.

The model maintains the history variables defined by the hardening model.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 15, 30, 50, 8

   ``surface``, :cpp:class:`neml::YieldSurface`, Flow surface interface, No
   ``hardening``, :cpp:class:`neml::NonAssociativeHardeningRule`, Hardening rule interface, No
   ``fluidity``, :cpp:class:`neml::FluidityModel`, Fluidity definition, No
   ``n``, :cpp:class:`neml::Interpolate`, Rate sensitivity exponent, No
   ``prefactor``, :cpp:class:`neml::Interpolate`, Prefactor, 1.

Class description
-----------------

.. doxygenclass:: neml::ChabocheFlowRule
   :members:
   :undoc-members:

.. _fluidity-model:

Fluidity model
--------------

The general Chaboche model allows the fluidity to vary with the integrated effective inelastic strain

.. math::

   \bar{\varepsilon}_{vp}=\int_{0}^{t}\sqrt{\frac{2}{3}\dot{\bm{\varepsilon}}_{vp}:\dot{\bm{\varepsilon}}_{vp}}dt.

These models then define the fluidity as 

.. math::

   \eta \leftarrow \mathcal{N}\left(\bar{\varepsilon}_{vp}, T \right).

.. doxygenclass:: neml::FluidityModel
   :members:
   :undoc-members:

Constant fluidity
^^^^^^^^^^^^^^^^^

This option keeps the fluidity constant with the effective inelastic strain.  

.. math::

   \eta = c

It can still vary with temperature.

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``eta``, :cpp:class:`neml::Interpolate`, Value of the fluidity, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::ConstantFluidity
   :members:
   :undoc-members:

Saturating fluidity
^^^^^^^^^^^^^^^^^^^

This option evolves the fluidity from some initial value through some increment as an 
exponential function of inelastic strain.
The fluidity eventually saturates to a final, fixed value.

.. math::

   \eta = K_0 + A \left(1 - e^{-b \bar{\varepsilon}_{vp}} \right) 

Parameters
""""""""""

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``K0``, :cpp:class:`neml::Interpolate`, Initial fluidity, No
   ``A``, :cpp:class:`neml::Interpolate`, Saturated fluidity is K0 + A, No
   ``b``, :cpp:class:`neml::Interpolate`, Saturation speed exponent, No

Class description
"""""""""""""""""

.. doxygenclass:: neml::SaturatingFluidity
   :members:
   :undoc-members:
