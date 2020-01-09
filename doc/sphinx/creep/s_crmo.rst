Minimum creep law for 2.25Cr-1Mo steel
======================================

Overview
--------

This is a quite complicated minimum creep rate law for 2.25Cr-1Mo (Gr 22)
steel, to be documented in a PVP paper at some point in the future.

The creep rate law it implements is:

if :math:`\sigma_{eq}\le 60`:

.. math::
   \dot{\varepsilon}^{cr} = \dot{\varepsilon}_{2}

else if :math:`T\le13.571\sigma_{eq}^{0.68127}-1.8\sigma_{eq}+437.63`

.. math::
   \dot{\varepsilon}^{cr} = \dot{\varepsilon}_{1} 

else

.. math::
   \dot{\varepsilon}^{cr} = \dot{\varepsilon}_{2} 

with

.. math::
   \dot{\varepsilon}_{1}=\frac{10^{6.7475+0.011426\sigma_{eq}+\frac{987.72}{U}\log\sigma_{eq}-\frac{13494}{T}}}{100}

.. math::
   \dot{\varepsilon}_{2}=\frac{10^{11.498-\frac{8.2226U}{T}-\frac{20448}{T}+\frac{5862.4}{T}\log\sigma_{eq}}}{100}

and :math:`U` a parameter interpolated linearly as a function of temperature 
from the table:

=========== =====
Temperature Value
=========== =====
644.15      471
673.15      468
723.15      452
773.15      418
823.15      634
873.15      284
894.15      300
922.15      270
=========== =====

Because the model is fully parameterized in the C++ implement it must
be used with units of MPa, hours, and Kelvin.

Parameters
----------

None, all built into class.

Class description
-----------------

.. doxygenclass:: neml::MinCreep225Cr1MoCreep
   :members:
   :undoc-members:
