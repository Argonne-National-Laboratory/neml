# Building and installing _neml_

## macOS (homebrew)

Install prerequisites and set up environment with

```
brew install libxml++3 libxml2 boost
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/local/opt/libxml2/lib/pkgconfig
```

(the addition to the `PKG_CONFIG_PATH` is necessary to allow cmake to properly detect `libxml++3`)
Next build neml with (in the `neml` directory)

```
mkdir build && cd build
cmake ..
make
```

## Linux (with APT)

Install prerequisites with

```
sudo apt install libxml++2.6-dev libboost1.65-dev
```

(the package versions above are for Ubuntu 18.04.01 LTS, adjust accordingly for newer releases)
Next build neml with (in the `neml` directory)

```
mkdir build && cd build
cmake ..
make
```
