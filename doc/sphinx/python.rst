Python bindings and helpers
===========================

If you compile NEML with the python bindings the build will 
produce compiled python modules in the :file:`neml/` directory.
The names of these modules match the underlying C++ files and 
also match the chapter titles in this manual.
So for example, yield surfaces will be in the module surfaces and so on.
All public methods of all classes are available in python.
The modules are configured to take standard Python lists as input, rather
than numpy arrays.

To help developing, testing, and debugging material models NEML provides
several python drivers in the neml python module.
These helpers run material models degenerate NEML's 3D formulation to 1D
uniaxial stress, exercise NEML material models in common experimental
loading situations (uniaxial tension, creep, stress relaxation,
strain and strain controlled cyclic loading, with holds, and strain rate
jump tests), run simulations of arbitrarily connected collections of 
uniaxial bars, where a NEML model gives the constitutive response of
each bar, and runs simplified, 1D models of pressure vessels, again
using a NEML model to define the vessel's constitutive response.

.. toctree::
   python/uniaxial
   python/drivers
   python/arbbar
   python/axisym
