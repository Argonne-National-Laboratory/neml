uniaxial
========

NEML material models are strain-controlled and 3D.
Degenerating this response to plane strain is trivial -- simply pass in
zeros in the appropriate strain components.
However, another very useful stress state is strain-controlled uniaxial stress.
This is the stress state in many common experimental tests, for 
example standard tension tests.
The stress state in these conditions is:

.. math::
   
   \bm{\sigma}=
      \left[\begin{array}{ccc}
      \sigma & 0 & 0\\
      0 & 0 & 0\\
      0 & 0 & 0
      \end{array}\right]

where :math:`\sigma` is the unknown uniaxial stress.  The strain state is

.. math::
   
   \bm{\varepsilon} = 
      \left[\begin{array}{ccc}
      \varepsilon & \varepsilon_{12} & \varepsilon_{13}\\
      \varepsilon_{12} & \varepsilon_{22} & \varepsilon_{23}\\
      \varepsilon_{13} & \varepsilon_{23} & \varepsilon_{33}
      \end{array}\right]

where :math:`\varepsilon` is the known, controlled input strain and the
remaining strain components are unknowns.

The neml axisym module solves a system of nonlinear equations to impose
this state of strain and stress on a standard, 3D NEML material model.
The module returns the uniaxial stress, the history, the stored
energy and dissipated work, along with the new, uniaxial algorithmic
tangent

.. math::
   
   \mathbf{A} = \frac{\mathrm{d} \sigma}{\mathrm{d} \varepsilon}
