include(macros.cmake)

### BASE LIBRARY IN LIB ###
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
if (WIN32)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
endif()
add_library(neml SHARED
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
      history.cxx)
add_subdirectory(math)
add_subdirectory(cp)
target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES} ${OpenMP_CXX_LIBRARIES})


### python bindings in neml ###
if (WRAP_PYTHON)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml)
      if (WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml)
      endif()
      # Actual modules
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
endif()
