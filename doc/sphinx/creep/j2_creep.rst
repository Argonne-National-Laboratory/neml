J2 creep
========

Overview
--------

This interface specializes a generic creep rate model to :math:`J_2` flow.
The model creeps in the direction of the deviatoric stress and the scalar
creep rate is defined by an :doc:`additional interface <scalar>`.
The creep rate for these types of models is then

.. math::
   \dot{\bm{\varepsilon}}^{cr} = 
      \dot{\varepsilon}^{cr}\left(\sigma_{eff}, \varepsilon_{eff}, t, T \right)
      \frac{\operatorname{dev}\left(\bm{\sigma}\right)}
      {\left\Vert \operatorname{dev}\left(\bm{\sigma}\right) \right\Vert}

with the scalar creep rate a function of effective stress

.. math::
   \sigma_{eff} = \sqrt{\frac{3}{2} \bm{\sigma} : \bm{\sigma}},

effective strain

.. math::
   \varepsilon_{eff} = \sqrt{\frac{2}{3} \bm{\varepsilon} : \bm{\varepsilon}},

time, and temperature.

Scalar creep models
-------------------

.. toctree::

   scalar

Parameters
----------

========== ========================= ======================================= =======
Parameter  Object type               Description                             Default
========== ========================= ======================================= =======
rule       ScalarCreepRule           Scalar creep model                      No
tol        double                    Nonlinear solver tolerance              1.0e-10
miter      int                       Maximum number of nonlinear iterations  25
verbose    bool                      Verbosity flag                          false
========== ========================= ======================================= =======

Class description
-----------------

.. doxygenclass:: neml::J2CreepModel
   :members:
   :undoc-members:
