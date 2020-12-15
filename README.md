# NEML: the Nuclear Engineering material Model Library

[![Documentation Status](https://readthedocs.org/projects/neml/badge/?version=stable)](https://neml.readthedocs.io/en/stable/) ![Test Status](https://github.com/Argonne-National-Laboratory/neml/workflows/tests/badge.svg?branch=master)

[![Documentation Status](https://readthedocs.org/projects/neml/badge/?version=dev)](https://neml.readthedocs.io/en/dev/) ![Test Status](https://github.com/Argonne-National-Laboratory/neml/workflows/tests/badge.svg?branch=dev)

NEML (the Nuclear Engineering Material model Library) is a tool for creating
and running structural material models.
While it was originally developed to model high temperature nuclear reactors,
the tool is general enough to apply to most types of structural materials.

The focus of NEML is on modularity and extensibility.
The library is structured so that adding a new feature to an existing material
model should be as simple as possible and require as little code as possible.

NEML material models are modular -- they are built up from smaller pieces into
a complete model.
For example, a model might piece together a temperature-dependent elasticity
model, a yield surface, a flow rule, and several hardening rules.
Each of these submodels is independent of the other objects
so that, for example, switching from conventional J2 plasticity
to a J2 theory requires only a one line change in an input file,
if the model is already implemented, or a relatively small amount of coding
to add the new yield surface if it has not been implemented.
All of these objects are interchangeable.
For example, the damage, viscoplastic, and rate-independent plasticity
models all use the same yield (flow) surfaces, hardening rules, elasticity
models, and so on.

As part of this philosophy, the library only requires new components
provide a few partial derivatives and NEML uses this information to assemble
the Jacobian needed to do a fully implement, backward Euler integration of the
ordinary differential equations comprising the model form and to provide 
the algorithmic tangent needed to integrate the model into an implicit
finite element framework.

Documentation and tutorials are available [here](https://neml.readthedocs.io).

## License

The library is provided under an MIT license found in the
[LICENSE](LICENSE) file.

The NEML distribution contains a copy of
the [pybind11](https://github.com/pybind/pybind11) header library, which
has its own license contained in the pybind11 subdirectory.

The NEML distribution contains a copy of
the [rapidxml](http://rapidxml.sourceforge.net) header library, which
has its own license contained in the rapidxml subdirectory.

## Python package

The python package version of the software available on PiPy requires:

1. A working C++ compiler.
2. Linkable version of [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/).

This version does not support multithreading or linking with external C/C++/Fortran codes.
