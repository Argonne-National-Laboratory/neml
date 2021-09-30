CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(NEML LANGUAGES CXX)

### Version ###
set(NEML_version 1.3.5)

### Setup modules ###
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

### Need c++ 11 ###
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Default to shared libs ###
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

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

### Limits designed to make bad updates fail gracefully ###
set(STRAIN_RATE_LIMIT "1.0e10" CACHE STRING "Limit on strain rate -- NEML may return an error if a strain rate grows beyond this value")
mark_as_advanced(STRAIN_RATE_LIMIT)
add_definitions(-DNEML_STRAIN_RATE_LIMIT=${STRAIN_RATE_LIMIT})

### Configure standard-ish libraries ###
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

INCLUDE_DIRECTORIES(SYSTEM rapidxml)

### PLATFORM AND COMPILER SPECIFIC OPTIONS ###
# Make better debug on Intel
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -traceback -ftrapuv")
      set(CMAKE_CXX_FLAGS_RELEASE "-O2")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -pg")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wextra -Wall -Wno-unused-parameter -Werror -D_GLIBCXX_DEBUG")
      set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pedantic -Wextra -Wall -Wno-unused-parameter")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -march=native -g -pedantic -Wextra -Wall -Wno-unused-parameter")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wextra -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Werror -Wno-unused-private-field -Wno-self-assign-overloaded")
      set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pedantic -Wextra -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Wno-unused-private-field")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -pedantic -Wextra -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Wno-unused-private-field")
endif()

# For MacOS
if(APPLE)
#      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -undefined dynamic_lookup")
#      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
      set(CMAKE_SKIP_INSTALL_RPATH TRUE)
endif(APPLE)
######

### Offer to use OMP for crystal system ###
option(USE_OPENMP "Compile with OpenMP" ON)
if (USE_OPENMP)
      if (APPLE)
            add_definitions(-DUSE_OMP)
            set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-Xpreprocessor -fopenmp")
            set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Xpreprocessor -fopenmp")
            set(OpenMP_CXX_LIBRARIES "-lomp")
            set(OpenMP_C_LIBRARIES "-lomp")
      else()
            find_package(OpenMP)
            if(NOT OPENMP_FOUND)
                  message(WARNING "OpenMP requested but the compiler does not support it!")
            else()
                  add_definitions(-DUSE_OMP)
                  set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
                  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
            endif()
      endif()
endif()
###

### DOCUMENTATION ###
option(BUILD_DOCUMENTATION "Build documentation: manual and doxygen" OFF)
if (BUILD_DOCUMENTATION)
      add_subdirectory(doc)
endif()


### PYTHON BINDINGS ###
option(WRAP_PYTHON "Generate a pybind11 wrapper" OFF)
if (WRAP_PYTHON)
      #find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
      INCLUDE_DIRECTORIES(SYSTEM pybind11)
      # There is some problem with LTO on mingw
      if (WIN32)
	      set(PYBIND11_LTO_CXX_FLAGS "")
      endif()
endif()

### Check if we can validly build shared libraries ###
if ((NOT BUILD_SHARED_LIBS) AND WRAP_PYTHON)
      message(FATAL_ERROR "Cannot build python wrappers with statically-linked based library") 
endif()

# We are building a python package
if (INSTALL)
      set(MODULE_BASE ${INSTALL})
      set(LIBRARY_BASE ${INSTALL})
      set(CMAKE_INSTALL_RPATH ".")
# We are building the wrappers locally
elseif(WRAP_PYTHON)
      set(MODULE_BASE ${PROJECT_BINARY_DIR}/neml)
      set(LIBRARY_BASE ${PROJECT_BINARY_DIR}/neml)
# We are just building the library itself
else()
      set(LIBRARY_BASE ${PROJECT_BINARY_DIR}/lib)
endif()

### BASE LIBRARY ###
add_subdirectory(src)

### ABAQUS HELPER ###
option(BUILD_UTILS "Generate interface examples and helpers for Abaqus UMATS" OFF)
if (BUILD_UTILS)
      add_subdirectory(util)
endif()
