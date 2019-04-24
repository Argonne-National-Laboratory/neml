CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(NEML LANGUAGES CXX)

### Version ###
set(NEML_version 1.1.0)

### Setup modules ###
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


### Need c++ 11 ###
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Select nonlinear solver ###
set(SOLVER "newton" CACHE STRING "Nonlinear solver to use")
set_property(CACHE SOLVER PROPERTY STRINGS newton nox)
if (${SOLVER} MATCHES "newton")
      # Built in, so no configuration
      add_definitions(-DSOLVER_NEWTON)
      set(SOLVER_LIBRARIES "")
elseif (${SOLVER} MATCHES "nox")
      add_definitions(-DSOLVER_NOX)
      set(TRILINOS_PATH "" CACHE STRING "Path to trilinos installation")
      FIND_PACKAGE(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
      INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
      set(SOLVER_LIBRARIES ${Trilinos_LIBRARIES})
else()
      message("Unknown solver")
endif()

### Need Boost (only variant.hpp really) ###
FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})

### Configure standard-ish libraries ###
set(BLAS_VENDOR Generic)
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

INCLUDE_DIRECTORIES(SYSTEM rapidxml)

### PLATFORM AND COMPILER SPECIFIC OPTIONS ###
# Make better debug on Intel
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -traceback -ftrapuv")
      set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wextra -Wall -Wno-unused-parameter -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pedantic -g -Wextra -Wall -Wno-unused-parameter")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wextra -Wall -Wno-unused-parameter -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pedantic -g -Wextra -Wall -Wno-unused-parameter")
else()
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
      set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g")
endif()

# For MacOS
if(APPLE)
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -undefined dynamic_lookup")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif(APPLE)
######


### DOCUMENTATION ###
option(BUILD_DOCUMENTATION "Build documentation: manual and doxygen" OFF)
if (BUILD_DOCUMENTATION)
      add_subdirectory(doc)
endif()


### PYTHON BINDINGS ###
option(WRAP_PYTHON "Generate a pybind11 wrapper" OFF)
if (WRAP_PYTHON)
      FIND_PACKAGE(PythonLibs 3.4 REQUIRED)
      INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})
      INCLUDE_DIRECTORIES(SYSTEM pybind11)
endif()

### BASE LIBRARY ###
add_subdirectory(src)

### ABAQUS HELPER ###
option(BUILD_UTILS "Generate interface examples and helpers for Abaqus UMATS" OFF)
if (BUILD_UTILS)
      add_subdirectory(util)
endif()
