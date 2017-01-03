CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

project(NEML LANGUAGES CXX)

# Setup modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Need c++ 11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# This is risky, but I don't use the preprocessor a lot
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wp,-w")

# Select nonlinear solver
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

# Configure libraries
set(BLA_VENDOR Generic)
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

FIND_PACKAGE(libxml++ REQUIRED)
INCLUDE_DIRECTORIES(${libxml++_INCLUDE_DIRS})

# Option for python
option(WRAP_PYTHON "Generate a pybind11 wrapper" OFF)
if (WRAP_PYTHON)
      FIND_PACKAGE(PythonLibs REQUIRED)
      INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml)
else()
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
      # Headers?
endif()

# Make better debug on Intel
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -traceback -ftrapuv")
endif()

# Recurse
add_subdirectory(src)
