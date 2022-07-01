Linear vicous flow rule
=======================

Overview
--------

This simple model provides a linear viscous response of the type

.. math::

   \dot{\gamma} = \frac{3}{2} \frac{f\left(\bm{\sigma}, \mathbf{0}, T\right)}{\eta\left(T\right)}

   \mathbf{g}_{\gamma} = \frac{\partial f}{\partial \bm{\sigma}} \left( \bm{\sigma}, \mathbf{0}, T  \right)

where :math:`f` is a flow surface

The model does not maintain internal variables.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``surface``, :cpp:class:`neml::YieldSurface`, Flow surface interface, No
   ``eta``, :cpp:class:`neml::Interpolate`, Drag stress, No

Class description
-----------------

.. doxygenclass:: neml::LinearViscousFlow
   :members:
   :undoc-members:
