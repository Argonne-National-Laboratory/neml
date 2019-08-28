# Ansys NEML Plugin Installation and Usage

## 1. Basic Installation Guide

Most of the files are needed to compile are contained in this repository. The following steps are a guide of how to run NEML on Ansys on NIMBIX server (same logic applies to any ANSYS user)

#### Install on Linux

clone the into local working directory

```bash
git clone https://github.com/khhsu0724/neml.git
git fetch
git checkout ansys_integration
```

Upload the file to NIMIBX (or execute the following commands if the user is not on NIMBIX)

```bash
cd kp-neml
export ANS_USER_PATH=/home/nimbix/data/kp-neml/util/ansys
# Set working directory
cd util/ansys
```

launch Ansys APDL graphic interface

```bash
ansys192 -g
```

launch Ansys APDL graphic interface on NIMBIX

```bash
/ansys_inc/v192/ansys/bin/ansys192 -g
```

#### Calling NEML library using APDL script

An example of the a uniaxial tensile test is given in **kp-neml/util/ansys/inelastic.inp**. This extension uses usermat.F provided by Ansys to call NEML library. The desired material should be put in **umat.xml** with **ansys** tag. There are some example materials provided in umat.xml. To call the material in APDL, use the following syntax in /PREP7 block:
	
```
TB, USER, 1, 0,
TB, STATE, 1, , 1000
```
It should be noted that the more complex the material model is, the more state variable user should provide to run a simulation (which is 1000 in this case), however this number can be reduced if necessary to free up more memory space. Any material used here should be tested in NEML before running in ANSYS.

## 2. Recompile NEML and usermat.F

The following block is a guide in case any code or is changed in NEML or usermat.f.

To recompile NEML, follow the guide here: [installation guide][1]. Make sure all the correct depedencies are installed.


#### Requirements for compilation

* The user must have intel fortran compiler in order to compile
* ANS_PATH must point to directory where ansys is installed

#### Recompile on Linux

```bash
# In case user made changes on github
git pull
cd kp-neml
cmake -D CMAKE_BUILD_TYPE="Release" .
make
mv lib/libneml.o util/ansys
cd util/ansys
./ANSUSERSHARED
```

The libneml in this repository was compiled using openblas, if a user intended to use different blas library, please change CMakeLists.txt to the desired blas, lapack library and their respective direcotry. Also, it is very important to copy the blas linked library into util/ansys. remove libopenblas.so.0 in util/ansys and copy in your own blas linked library before executing ./ANSUSERSHARED.

[1]: https://neml.readthedocs.io/en/dev/started.html?highlight=release

## 3. Using this in Workbench

To use this in workbench, move the workbench project file (.wbpj) into the /util/ansys directory. In workbench after following the installation guide in part 1 , it is suggested for user to read umat.xml from an absolute directory (Right now it's configured to NIMBIX setting). To specify a element to use NEML library material, edit the model and go to geometry and the desired element and add an APDL command with the following lines:

```bash
/UNIT, SI
MPDEL, ALL, matid
MAT, matid
TB, USER, matid, 0,
TB, STATE, matid, , 50
```

ANSYS workbench can be launched from this command:

```bash
/ansys_inc/v192/Framework/bin/Linux64/runwb2
```
