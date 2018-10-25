Python bindings and helpers
===========================

If you compile NEML with the python bindings the build will 
produce compile python modules in the neml/ directory.
The names of these modules match the underlying C++ files and 
also match the chapter titles in this manual.
So for example, yield surfaces will be in the module surfaces and so on.
All public methods of all classes are available in python.
The modules are configured to take standard Python lists as input, rather
than numpy arrays.

To help developing, testing, and debugging material models NEML provides
several python drivers in the neml python module.
