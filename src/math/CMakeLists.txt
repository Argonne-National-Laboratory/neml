include(../macros.cmake)

cmake_policy(SET CMP0076 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml/math)

target_sources(neml PRIVATE nemlmath.cxx rotations.cxx)

if (WRAP_PYTHON)
      pybind(nemlmath)
      pybind(rotations)
endif()
