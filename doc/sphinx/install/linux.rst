Linux
=====

.. _basic-library:

Basic library
"""""""""""""

This instructions assume a clean installation of Ubuntu 18.04.1 LTS.
Installation on other Linux systems will be similar.
The package maintainers have compiled and run the library on CentOS and
Debian without difficulty.  The basic steps remain the same.

First install the prerequisites

.. code-block:: console

   apt-get install build-essential git cmake libboost-dev libblas-dev liblapack-dev

Clone the neml source code

.. code-block:: console

   git clone https://github.com/Argonne-National-Laboratory/neml.git

NEML builds in the source directory.
Enter the NEML directory and configure NEML using CMake:

.. code-block:: console

   cmake .

Useful options might include ``-D CMAKE_BUILD_TYPE="Release"`` if you want
to build an optimized version for production runs.

Then simply make the library:

.. code-block:: console

   make

Python bindings
"""""""""""""""

Building the base library is sufficient to link NEML into external finite
element software
However, all the tests and provided examples use the Python bindings.

To build the bindings you will need a few more prerequisites, in addition
to those mentioned above:

.. code-block:: console

   apt-get install python3-dev python3-networkx python3-numpy python3-scipy python3-matplotlib python3-nose

Configure with CMake to setup the python bindings:

.. code-block:: console

   cmake -D WRAP_PYTHON=ON .

And build the library

.. code-block:: console

   make

You can now run the test suit with

.. code-block:: console

   nosetests

Running examples
""""""""""""""""

Once you have the python bindings you can test your compilation of NEML
using the python tests in the :file:`tests/` directory.
If you installed nose, all the tests can be run from the root :file:`neml`
directory by running :command:`nosetests`.

Assuming the tests passed, you can begin to build material models with NEML.
The manual has a :doc:`../tutorial` giving a brief tutorial on setting up a
material model either with the python bindings or the XML input files
and then running that model using the python drivers for some simple
loadings.
Additional examples can be found in the :file:`examples/` directory.


Linking to external software
""""""""""""""""""""""""""""

The main NEML library (in the :file:`lib/`) directory is all that needs to
be linked to your software to call NEML material models.
You only need to include the :file:`src/neml_interface.h` in order to
load material models from XML datafile and use the resulting C++ object
to call for the material response.

The :file:`util/` directory contains example bindings of NEML into
C++, C, and Fortran codes.  The CMake variable ``-D BUILD_UTILS=ON`` option
compiles these example interfaces.
Turing this option on requires a Fortran and C compiler.
Looking at these examples demonstrates how you can integrate NEML into your
finite element code.

Abaqus UMAT interface
---------------------

The :file:`util/abaqus` directory contains a full UMAT interface
that can be used to tie NEML into `Abaqus <https://www.3ds.com/products-services/simulia/products/abaqus/>`_.
This first requires compiling the :ref:`main NEML library <basic-library>`.
Say the full path to :file:`libneml.so` is :envvar:`${NEMLROOT}/lib/libneml.so`.
You would need to alter your abaqus env file (for example :file:`abaqus_v6.env`) to
*add* the library to the ``link_sl`` command.
For example, if the existing ``link_sl`` is:

.. code-block:: bash

   link_sl = [fortCmd,
              '-cxxlib', '-fPIC', '-threads', '-shared','-Wl,--add-needed',
              '%E', '-Wl,-soname,%U', '-o', '%U', '%F', '%A', '%L', '%B', '-parallel',
              '-Wl,-Bdynamic', '-shared-intel']

then you would alter it to

.. code-block:: bash

   link_sl = [fortCmd,
              '${NEMLROOT}/lib/libneml.so', '-V',
              '-cxxlib', '-fPIC', '-threads', '-shared','-Wl,--add-needed',
              '%E', '-Wl,-soname,%U', '-o', '%U', '%F', '%A', '%L', '%B', '-parallel',
              '-Wl,-Bdynamic', '-shared-intel']

You then need to determine the correct number of ``*DEPVAR`` and the correct
``INITIAL CONDITIONS, TYPE=SOLUTION`` to include in your input file in order to
have Abaqus setup and maintain the correct number of history variables for the
NEML model.
The distribution provides a simple program in the :file:`util/abaqus/` directory to
report this information.
The program, called :file:`report` is compiled if the CMake ``BUILD_UTILS`` option is set.
It requires two command line arguments:

**report**

   .. program:: report

   .. option:: file

      Name of the XML input file

   .. option:: model

      Material model to report on in the XML file

The program will print the correct lines to use in your Abaqus input file for that
NEML material.

You should then copy the XML file containing the model you want to run to the
directory containing the Abaqus input file.
You must rename this XML input file to :file:`neml.xml`.
You should rename the model in that file you want to use in Abaqus to ``abaqus``.
The UMAT is hardcoded to load that material from that filename.

The remaining steps are standard for any UMAT.  You need to request Abaqus call the
UMAT in the input file:

.. code-block:: bash

   *MATERIAL, NAME=CUSTOM

   *USER MATERIAL, CONSTANTS=0, UNSYMM

Remembering to also include the output from :file:`report` to initalize the required
history variables.

Finally, run the UMAT

.. code-block:: bash

   abaqus job=xxxx user=/path/to/neml/util/abaqus/nemlumat.f

ANSYS UMATERIAL interface
-------------------------

Directions are in preparation
