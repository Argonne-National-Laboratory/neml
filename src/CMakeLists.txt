include(macros.cmake)
include (GenerateExportHeader)

### Build library in lib to enable easy in-source build ###
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
if (WIN32)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
endif()

add_library(neml
      objects.cxx
      models.cxx
      solvers.cxx
      surfaces.cxx
      hardening.cxx
      ri_flow.cxx
      visco_flow.cxx
      general_flow.cxx
      nemlerror.cxx
      elasticity.cxx
      parse.cxx
      cinterface.cxx
      interpolate.cxx
      creep.cxx
      damage.cxx
      history.cxx
      larsonmiller.cxx
      walker.cxx
      block.cxx
      deparse.cxx
      )
add_subdirectory(math)
add_subdirectory(cp)

# Required headers
target_include_directories(neml PRIVATE "../include")

# Link the library and generate the export header
target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
generate_export_header(neml EXPORT_FILE_NAME ${PROJECT_SOURCE_DIR}/include/neml_export.h)
generate_export_header(neml EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/neml_export.h)
set_target_properties(neml PROPERTIES PUBLIC_HEADER 
      "${PROJECT_BINARY_DIR}/include/neml_export.h;../include/neml_interface.h;../include/nemlerror.h;../include/windows.h")

# An annoyance for building the python package
if (PYTHON_PACKAGE)
      # Strictly we don't need this for wheels, but for source builds it's handy
      install(TARGETS neml 
            LIBRARY DESTINATION neml)
else()
      install(TARGETS neml 
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include)
endif()

### python bindings in neml ###
if (WRAP_PYTHON)
      # Make sure we install all the python files in neml/
      file(GLOB python_files ${CMAKE_SOURCE_DIR}/neml/*.py)
      install(FILES ${python_files} DESTINATION neml)

      # Copy the stupid library so windows can easily find it for tests
      add_custom_command(TARGET neml POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:neml> ${PROJECT_BINARY_DIR}/neml) 

      pybind(objects)
      pybind(solvers)
      pybind(elasticity)
      pybind(surfaces)
      pybind(hardening)
      pybind(ri_flow)
      pybind(visco_flow)
      pybind(general_flow)
      pybind(models)
      pybind(parse)
      pybind(interpolate)
      pybind(creep)
      pybind(damage)
      pybind(history)
      pybind(larsonmiller)
      pybind(walker)
      pybind(block)
endif()
