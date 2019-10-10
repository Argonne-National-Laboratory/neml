SkewSym
=======

This class represents a rank four tensor with symmetries
:math:`C_{ijkl} = -C_{jikl}` and :math:`C_{ijkl} = C_{ijlk}`.
It is stored as a length 18 flat array in a convention that makes it
the natural way to store the derivative of a symmetric rank 2 tensor with
respect to a skew symmetric rank 2 tensor.

.. doxygenclass:: neml::SkewSym
   :members:
   :undoc-members:
