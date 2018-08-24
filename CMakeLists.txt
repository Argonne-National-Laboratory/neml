CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(NEML LANGUAGES CXX)

### Setup modules ###
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


### Need c++ 11 ###
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)


### Select nonlinear solver ###
set(SOLVER "newton" CACHE STRING "Nonlinear solver to use")
set_property(CACHE SOLVER PROPERTY STRINGS newton nox)
if (${SOLVER} MATCHES "newton")
      # Built in, so no configuration
      add_definitions(-DSOLVER_NEWTON)
      set(SOLVER_LIBRARIES "")
elseif (${SOLVER} MATCHES "nox")
      add_definitions(-DSOLVER_NOX)
      set(NOX_INCLUDE_PATH "$ENV{HOME}/packages/trilinos/include" CACHE STRING "NOX include path")
      set(NOX_LIBRARY_PATH "$ENV{HOME}/packages/trilinos/lib" CACHE STRING "NOX library path")
      INCLUDE_DIRECTORIES(${NOX_INCLUDE_PATH})
      set(SOLVER_LIBRARIES "-L${NOX_LIBRARY_PATH} -lnox -lnoxlapack")
else()
      message("Unknown solver")
endif()


### Configure standard-ish libraries ###
set(BLAS_VENDOR Generic)
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(libxml++ REQUIRED)
INCLUDE_DIRECTORIES(${libxml++_INCLUDE_DIRS})


### PLATFORM AND COMPILER SPECIFIC OPTIONS ###
# Make better debug on Intel
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -traceback -ftrapuv")
      set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
endif()

# Some apple things
if(APPLE)
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -undefined dynamic_lookup")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif(APPLE)
######


### DOCUMENTATION ###
option(BUILD_DOCUMENTATION "Build documentation: manual and doxygen" OFF)
if (BUILD_DOCUMENTATION)
      option(BUILD_DOXYGEN "Build and install Doxygen documentation" ON)
      option(BUILD_MANUAL "Build and install LaTeX manual" ON)
      add_subdirectory(doc)
endif()


### PYTHON BINDINGS ###
option(WRAP_PYTHON "Generate a pybind11 wrapper" OFF)
if (WRAP_PYTHON)
      FIND_PACKAGE(PythonLibs 2.7 REQUIRED)
      INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
endif()

### BASE LIBRARY ###
add_subdirectory(src)

### ABAQUS HELPER ###
option(ABAQUS_UTILS "Generate helpers for Abaqus UMATS" OFF)
if (ABAQUS_UTILS)
      add_subdirectory(util/abaqus)
endif()
