Getting started
===============

Installing NEML
---------------

Compiling NEML requires the NEML source, and a C++ compiler.
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

Python package
--------------

An easier way to install NEML if you are only interested in the python bindings
is to use the package uploaded to `pypi <https://pypi.org/>`_.
Currently only a source package is available, but the python install scripts
simplify the process of compiling and linking the package.  You still
need a working compiler, the python development headers, and development
versions of BLAS and LAPACK.  These can be installed on Ubuntu, for example,
with 

.. code-block:: console

   sudo apt-get install python3-dev python3-pip cmake libboost-dev libblas-dev liblapack-dev

After that a python package manger, such as pip, can be used to install NEML

.. code-block:: console

   sudo pip3 install neml
