### Macro helper for python bindings ###
macro(pybind sname)
      # Figure out if we need to put this in a subdirectory
      set (variadic_args ${ARGN})
          
      # Did we get any optional args?
      list(LENGTH variadic_args variadic_count)
      if (${variadic_count} GREATER 0)
            list(GET variadic_args 0 subdir)
            set(output_location neml/${subdir})
      else()
            set(output_location neml)
      endif ()

      # Dump to the right location
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${output_location})

      if (WIN32)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${output_location})
      endif()

      set(wrapper "${sname}_wrap.cxx")
      pybind11_add_module(${sname} SHARED ${wrapper})
      target_include_directories(${sname} PRIVATE "${CMAKE_SOURCE_DIR}/include")      
      target_link_libraries(${sname} PRIVATE neml)

      install(TARGETS ${sname} LIBRARY DESTINATION ${output_location})
endmacro(pybind)
