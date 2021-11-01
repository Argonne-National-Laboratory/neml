include(../macros.cmake)

target_sources(neml PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}/crystallography.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/inelasticity.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/sliprules.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/slipharden.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/kinematics.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/singlecrystal.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/batch.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/polycrystal.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/crystaldamage.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/hucocks.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/postprocessors.cxx
      )

if (WRAP_PYTHON)
      # Make sure we install all the python files in neml/
      file(GLOB python_files ${CMAKE_SOURCE_DIR}/neml/cp/*.py)
      install(FILES ${python_files} DESTINATION neml/cp)

      pybind(crystallography cp)
      pybind(inelasticity cp)
      pybind(sliprules cp)
      pybind(slipharden cp)
      pybind(kinematics cp)
      pybind(singlecrystal cp)
      pybind(batch cp)
      pybind(polycrystal cp)
      pybind(crystaldamage cp)
      pybind(hucocks cp)
      pybind(postprocessors cp)
endif()
