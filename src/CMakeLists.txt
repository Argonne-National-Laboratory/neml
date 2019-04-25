### Macro helper for python bindings ###
macro(pybind file)
      get_filename_component(sname ${file} NAME_WE)
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

set(wrapped_src 
      objects.cxx
      models.cxx 
      nemlmath.cxx 
      solvers.cxx 
      surfaces.cxx 
      hardening.cxx 
      ri_flow.cxx
      visco_flow.cxx 
      general_flow.cxx 
      elasticity.cxx
      parse.cxx
      interpolate.cxx
      creep.cxx
      damage.cxx)
set(not_wrapped_src 
      nemlerror.cxx 
      cinterface.cxx)
set(libsrc ${not_wrapped_src} ${wrapped_src})

add_library(objlib OBJECT ${libsrc})

set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(libneml STATIC $<TARGET_OBJECTS:objlib>)
target_link_libraries(libneml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})

add_library(neml SHARED $<TARGET_OBJECTS:objlib>)
target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})

### python bindings in neml ###
if (WRAP_PYTHON)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/neml)      

      # Actual modules
      foreach(file ${wrapped_src})
            pybind(${file})
      endforeach()
endif()
