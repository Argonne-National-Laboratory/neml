if (WRAP_PYTHON)
      # Each file as a separate module
      add_library(volumetric SHARED volumetric_wrap.cxx volumetric.cxx pyhelp.cxx)
      target_link_libraries(volumetric ${PYTHON_LIBRARIES})
      set_target_properties(volumetric PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET volumetric PROPERTY OUTPUT_NAME "volumetric.so")
            set_property(TARGET volumetric PROPERTY SUFFIX "")
      endif()
else()
      # Everything in a big shared library
      add_library(neml SHARED 
            neml.cxx 
            volumetric.cxx 
            deviatoric.cxx)
endif()
