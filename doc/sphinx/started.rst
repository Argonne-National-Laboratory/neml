Getting started
===============

Installing NEML
---------------

Compiling NEML requires the NEML source, a C++ compiler,
and the `Boost <https://www.boost.org/>`_ library.
Additionally, you will need `BLAS <http://www.netlib.org/blas/>`_ and
`LAPACK <http://www.netlib.org/lapack/>`_.

Testing NEML and running the examples described here also requires
compiling the Python bindings.
This has the additional requirements of a
`Python <https://www.python.org/>`_ 3.x installation, including the
development headers,
the `pybind11 <https://github.com/pybind/pybind11>`_ library (which is included in the NEML
source), and the
`numpy <http://www.numpy.org/>`_,
`scipy <https://www.scipy.org/>`_,
and
`networkx <https://networkx.github.io/>`_ python packages.
The `nose <https://nose.readthedocs.io/en/latest/>`_ python package
can be useful in running the provided tests.

Linking NEML into your finite element package may additionally require a
C or Fortran compiler, depending on what language your finite element software
uses and can link to.

We have successfully installed NEML on Linux, Windows, and Mac OS.
Directions for each operating system follow below.

.. toctree::
   :maxdepth: 1

   install/linux
   install/mac
   install/windows
