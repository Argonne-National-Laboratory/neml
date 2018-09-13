### Macro helper for python bindings ###
macro(pybind sname)
      add_library(${sname} SHARED "${sname}_wrap.cxx")
      target_link_libraries(${sname} neml ${PYTHON_LIBRARIES})
      set_target_properties(${sname} PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET ${sname} PROPERTY OUTPUT_NAME "${sname}.so")
            set_property(TARGET ${sname} PROPERTY SUFFIX "")
      endif()
endmacro(pybind)

### BASE LIBRARY IN LIB ###
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
add_library(neml SHARED 
      objects.cxx
      models.cxx 
      nemlmath.cxx 
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
      damage.cxx)
target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})


### python bindings in neml ###
if (WRAP_PYTHON)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml)      

      # Actual modules
      pybind(objects)
      pybind(nemlmath)
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
endif()

