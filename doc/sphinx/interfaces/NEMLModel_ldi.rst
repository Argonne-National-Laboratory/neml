.. _NEMLModel_ldi:

NEMLModel_ldi
=============

Overview
--------

This object natively implements the :ref:`large strain incremental <strain-interfaces>`.  It accommodates the :ref:`small strain <strain-interfaces>`
by setting the deformation rate tensor equal to the small strain rate

.. math::
   \mathbf{d}_{n+1} = \bm{\varepsilon}_{n+1}

.. math::
   \mathbf{d}_{n} = \bm{\varepsilon}_{n}

and the vorticity to zero

.. math::
   \mathbf{w}_{n+1} = \mathbf{0}

.. math::
   \mathbf{w}_{n} = \mathbf{0}.

This means the skew part of the algorithmic tangent is zero

.. math::
   \mathbf{\mathfrak{B}}_{n+1} = 0.

Currently, this interface is only natively used by the :ref:`crystal plasticity <crystal-plasticity>` material models.

Class description
-----------------

.. doxygenclass:: neml::NEMLModel_ldi
   :members:
   :undoc-members:
