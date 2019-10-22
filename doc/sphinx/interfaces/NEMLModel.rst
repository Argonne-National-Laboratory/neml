.. _NEMLModel:

NEMLModel
=========

.. _strain-interfaces:

Stress update interfaces
------------------------

NEMLModel is the common interface to the constitutive models contained in NEML.
It currently requires models to provide two types of stress updates.
The first type of update is a small strain kinematics incremental update based on the interface

.. math::

   \bm{\sigma}_{n+1}, \bm{\alpha}_{n+1}, \bm{\mathfrak{A}}_{n+1}, u_{n+1}, p_{n+1} \leftarrow
   \mathcal{M}\left( 
   \bm{\varepsilon}_{n+1}, \bm{\varepsilon}_{n},
   T_{n+1}, T_{n},
   t_{n+1}, t_{n},
   \bm{\sigma}_{n},
   \bm{\alpha}_{n},
   u_n, p_n
   \right).

Here :math:`n` indicates values at the previous time step and :math:`n+1` values
at the next time step.
The quantities are stress (:math:`\bm{\sigma}`), strain (:math:`\bm{\varepsilon}`),
the vector of history variables (:math:`\bm{\alpha}`), strain energy (:math:`u`)
dissipated work (:math:`p`), temperature (:math:`T`), time (:math:`t`), and 
the algorithmic tangent (:math:`\mathbf{\mathfrak{A}}`)

.. math::
   \mathbf{\mathfrak{A}}_{n+1} = \frac{d \bm{\sigma}_{n+1}}{d \bm{\varepsilon}_{n+1}}.

The second update is a large strain kinematics update based on the interface

.. math::

   \bm{\sigma}_{n+1},\bm{\alpha}_{n+1},\bm{\mathfrak{A}}_{n+1},\bm{\mathfrak{B}}_{n+1},u_{n+1},p_{n+1}\leftarrow\mathcal{M}\left(\bm{d}_{n+1},\bm{d}_{n},\bm{w}_{n+1},\bm{w}_{n},T_{n+1},T_{n},t_{n+1},t_{n},\bm{\sigma}_{n},\bm{\alpha}_{n},u_{n},p_{n}\right).

Here :math:`\bm{d}` is the deformation rate tensor (the symmetric part of the spatial velocity gradient), :math:`\bm{w}` is the vorticity (the skew part of the spatial velocity gradient), :math:`\mathbf{\mathfrak{A}}` is

.. math::
   \mathbf{\mathfrak{A}}_{n+1} = \frac{d \bm{\sigma}_{n+1}}{d \bm{d}_{n+1}}

and :math:`\mathbf{\mathfrak{B}}` is

.. math::
   \mathbf{\mathfrak{B}}_{n+1} = \frac{d \bm{\sigma}_{n+1}}{d \bm{w}_{n+1}}

while the other quantities are defined identically to the small strain interface.


Implementations
---------------

.. toctree::

   NEMLModel_sd
   NEMLModel_ldi

Parameters
----------

None

Class description
-----------------

.. doxygenclass:: neml::NEMLModel
   :members:
   :undoc-members:
