### Macro helper for python bindings ###
macro(pybind sname)
      set(wrapper "${sname}_wrap.cxx")
      add_library(${sname} SHARED ${wrapper})
      target_link_libraries(${sname} neml ${PYTHON_LIBRARIES})
      set_target_properties(${sname} PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET ${sname} PROPERTY OUTPUT_NAME "${sname}.so")
            set_property(TARGET ${sname} PROPERTY SUFFIX "")
      endif()
      if (WIN32)

      endif()
endmacro(pybind)
