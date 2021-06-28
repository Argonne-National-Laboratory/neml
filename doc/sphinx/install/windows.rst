Windows
=======

We compile and test NEML on Windows using the `MSYS2 <https://www.msys2.org/>`_
system.  This provides a linux-like build environment for Windows.  We 
build and test the library using the `mingw-w64 <http://mingw-w64.org>`_
compiler.  MSYS2 supports other compilers and if needed you may be able to
compile and run NEML using one of the other support compilers.  These directions
assume the use of mingw.

Basic library
"""""""""""""

Go to the `MSYS2 website <https://www.msys2.org/>`_  and follow the instructions
there to install the framework.  Follow the directions all the way 
through updating the package repository and installing mingw with the
command:

.. code-block:: console

   pacman -Su
   pacman -S --needed base-devel mingw-w64-x86_64-toolchain

As the directions state, at this point close the original MSYS2 window,
go to the start menu, and run `"MSYS MinGW 64-bit"`.  You are now
ready to obtain the NEML source code and compile the library

Use the `"MSYS MinGW 64-bit"` terminal you opened to navigate to the 
location you want to download and install NEML to.  Install the required
dependencies with

.. code-block:: console

   pacman -S mingw-w64-x86_64-openblas mingw-w64-x86_64-cmake git

You can now clone the NEML source repository and enter the NEML directory
with

.. code-block:: console

   git clone https://github.com/Argonne-National-Laboratory/neml.git
   cd neml

Finally, configure the base library with CMake and build the library with

.. code-block:: console

   cmake . -G"MSYS Makefiles"
   make

This will build the NEML dynamic library, located in `lib` subdirectory in the 
NEML install.

Python bindings
"""""""""""""""

To go on and build the Python binding first follow the instructions in the
previous section for the base library.  Then install a few additional
dependencies with

.. code-block:: console

   pacman -S mingw-w64-x86_64-python-pip mingw-w64-x86_64-python-numpy mingw-w64-x86_64-python-scipy mingw-w64-x86_64-python-networkx mingw-w64-x86_64-python-matplotlib mingw-w64-x86_64-python-nose

Configure and build the Python bindings

.. code-block:: console

   cmake . -G"MSYS Makefiles" -D"WRAP_PYTHON=ON"
   make

This produces the NEML python package in the `neml` subdirectory.
   
You can run the automatic test suite using `nose` with the command

.. code-block:: console

   nosetests-3.8.exe

To use the `neml` Python package outside the source directory you need to
add the `neml\\neml` subdirectory to both the system `PATH` and `PYTHONPATH`
environment variables.
