Errors in NEML
==============

Currently, the error system in NEML is sadly mixed.
Most of the library returns error codes which are converted to exceptions 
for the python bindings.
However, the parser and object system uses C++ exceptions which are caught
and turned into error codes before entering the C or Fortran bindings.
We plan on fixing this inconsistency in future release by switching to
exceptions (and catching those exceptions before entering the C or Fortran
bindings).

For now, the types of errors are given by an enum:

.. doxygenenum:: neml::Error

The library provides routines for converting from error codes to descriptive
strings and from error codes to exceptions for the python bindings.

.. doxygenfunction:: neml::string_error

.. doxygenfunction:: neml::py_error
