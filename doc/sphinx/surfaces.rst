Yield surfaces
==============

Overview
--------

NEML YieldSurfaces provide the interface

.. math::
   f, 
      \frac{\partial f}{\partial \bm{\sigma}}, 
      \frac{\partial f}{\partial \mathbf{q}},
      \frac{\partial^2 f}{\partial \bm{\sigma}^2}, 
      \frac{\partial^2 f}{\partial \mathbf{q}^2},
      \frac{\partial^2 f}{\partial \bm{\sigma} \partial \mathbf{q}}, 
      \frac{\partial^2 f}{\partial \mathbf{q} \partial \bm{\sigma}},
   \leftarrow \mathcal{F}\left( \bm{\sigma}, \mathbf{q}\left( \bm{\alpha} \right), T \right).

A hardening interface provides the map between the "strain-like" history
variables :math:`\bm{\alpha}` and the "stress-like" history variables
:math:`\mathbf{q}`.
NEML uses this interface both as a yield surface in rate independent models
and as a flow surface for rate dependent plasticity.

While not required, the existing yield surfaces including
kinematic hardening expect to receive a 
:math:`\mathbf{q}` vector unrolling the isotropic hardening stress 
followed by a Mandel vector describing the total backstress
(i.e. :math:`\mathbf{q}` has length :math:`1 + 6 = 7`).
For pure isotropic hardening :math:`\mathbf{q}` consists of only 
the isotropic hardening stress (i.e. :math:`\mathbf{q}` has length 
:math:`1`).

.. WARNING::
   YieldSurfaces have no knowledge of the hardening model.
   They expect to recieve an unrolled vector of the appropriate 
   length.
   The only check the implementation can do is to make sure 
   that HardeningRule provides the correct length of the 
   :math:`\mathbf{q}` vector.
   It currently cannot check on the order of the provided "stress-like"
   history variables.

For convenience NEML provides a template class that automatically
converts a yield surface with isotropic and kinematic hardening to a
corresponding yield surface that only has isotropic hardening.
Supposing ``ExampleSurfaceIsoKin`` correctly implements combined isotropic and
kinematic hardening then defining a class ``ExampleSurfaceIso`` that 
inherits from ``IsoFunction`` with template arguments
``<ExampleSurfaceIsoKin, Arg A, Arg B, ...>`` where the arguments are
the arguments required to initialize ``ExampleSurfaceIsoKin`` will automatically
provide the isotropic-only version of the surface.

Implementations
---------------
.. toctree::
   
   surfaces/isokinj2
   surfaces/isoj2
   surfaces/isokinj2i1
   surfaces/isoj2i1

Class description
-----------------

.. doxygenclass:: neml::YieldSurface
   :members:
   :undoc-members:
