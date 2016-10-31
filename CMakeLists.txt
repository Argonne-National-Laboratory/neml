CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

project(NEML LANGUAGES CXX)

# Need c++ 11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Configure libraries


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

# Recurse
add_subdirectory(src)
