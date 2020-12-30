include(macros.cmake)
include (GenerateExportHeader)

if (WRAP_PYTHON)
      add_subdirectory(../pybind11 ${PROJECT_BINARY_DIR}/neml)
endif()

### BASE LIBRARY IN LIB ###
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_BASE})
if (WIN32)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${LIBRARY_BASE})
endif()
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIBRARY_BASE})

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
      block.cxx
      )
add_subdirectory(math)
add_subdirectory(cp)
target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
generate_export_header(neml)
generate_export_header(neml EXPORT_FILE_NAME ${PROJECT_SOURCE_DIR}/src/neml_export.h)

### python bindings in neml ###
if (WRAP_PYTHON)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${MODULE_BASE})
      if (WIN32)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MODULE_BASE})
      endif()

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
      pybind(block)
endif()
