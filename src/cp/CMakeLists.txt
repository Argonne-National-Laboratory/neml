include(../macros.cmake)

target_sources(neml PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}/crystallography.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/kinematics.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/sliprules.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/slipharden.cxx)

if (WRAP_PYTHON)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml/cp)
      pybind(crystallography)
      pybind(kinematics)
      pybind(sliprules)
      pybind(slipharden)
endif()
