# Building and installing neml

For detailed instructions see the [readthedocs](https://neml.readthedocs.io/en/latest/started.html) documentation.

If you are building only the base library you can build in or out of the source
directory.
For the python bindings we recommend building in-source.
These instructions assuming you are building in the source directory.

## Linux (Ubuntu/APT)

To compile just the base library:

```
apt-get install build-essential cmake libblas-dev liblapack-dev
cmake .
make
```

To compile both the base library and the python bindings:

```
apt-get install build-essential cmake libblas-dev liblapack-dev python-dev python-networkx python-numpy python-scipy python-matplotlib python-nose
cmake -D WRAP_PYTHON=ON .
make
```

and run the tests with

```
nosetests
```

## macOS (homebrew)

For just the base library setup prerequisites

```
brew install boost cmake
```

and then build the library

```
cmake .
make
```

For the library and the python bindings setup prerequisites

```
brew install boost cmake python@2
pip install nose numpy scipy nose matplotlib networkx
```

and then build the library
```
cmake -D WRAP_PYTHON=ON .
make 
```
