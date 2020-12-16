macOS
======

Basic library
"""""""""""""

These instructions assume a clean installation of macOS Catalina 10.15 but
may work on other version of the operating system.  The instructions
rely on `homebrew <https://brew.sh/>`_ to install the required packages.
The procedure here uses the native LLVM install for compiling the library.

First install the prerequisites:

.. code-block:: console

   brew install cmake openblas superlu

Then clone the neml source code

.. code-block:: console

   git clone https://github.com/Argonne-National-Laboratory/neml.git

move into the neml directory, and configure the library using ``cmake``

.. code-block:: console

   cmake -D CMAKE_BUILD_TYPE=Release -D USE_OPENMP=OFF .

The build type can be switched to ``CMAKE_BUILD_TYPE=Debug`` to compile a debug version of the library.  

Finally build the library with

.. code-block:: console

   make

Python bindings
"""""""""""""""

To compile the python bindings install Python3 using homebrew and then
install the python package dependencies using pip:

.. code-block:: console

   brew install python
   pip3 install --user networkx numpy scipy matplotlib nose

Configure the library in the ``neml`` directory

.. code-block:: console

   cmake -D CMAKE_BUILD_TYPE=Release -D WRAP_PYTHON=ON -D PYTHON_EXECUTABLE=$(python3-config --prefix)/bin/python3.9 -D PYTHON_LIBRARY=$(python3-config --prefix)/lib/libpython3.9.dylib -D PYTHON_INCLUDE_DIR=$(python3-config --prefix)/include/python3.9 -D USE_OPENMP=OFF .

These instructions assume the current homebrew python version is ``python3.9``.  The python directories will need to be changed if you have a different version of python installed.

Compile the library as before

.. code-block:: console

   make

Finally, you can run the automated test suite with

.. code-block:: console

   ~/Library/Python/3.9/bin/nosetests

The full path is required as ``pip`` does not install the ``nosetests`` script
in a ``PATH`` location by default.
