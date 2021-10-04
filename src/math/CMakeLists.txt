include(../macros.cmake)

target_sources(neml PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}/nemlmath.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/rotations.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/tensors.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/matrix.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/projections.cxx)

if (WRAP_PYTHON)
      # Make sure we install all the python files in neml/
      file(GLOB python_files ${CMAKE_SOURCE_DIR}/neml/math/*.py)
      install(FILES ${python_files} DESTINATION neml/math)

      pybind(nemlmath math)
      pybind(rotations math)
      pybind(tensors math)
      pybind(matrix math)
      pybind(projections math)
endif()
