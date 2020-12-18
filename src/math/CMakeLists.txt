include(../macros.cmake)

target_sources(neml PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}/nemlmath.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/rotations.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/tensors.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/matrix.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/projections.cxx)

if (WRAP_PYTHON)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${MODULE_BASE}/math)
      if (WIN32)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MODULE_BASE}/math)
      endif()

      pybind(nemlmath)
      pybind(rotations)
      pybind(tensors)
      pybind(matrix)
      pybind(projections)
endif()
