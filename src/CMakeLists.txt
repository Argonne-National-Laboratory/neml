if (WRAP_PYTHON)
      # Each file as a separate module
      # volumetric
      add_library(volumetric SHARED volumetric_wrap.cxx volumetric.cxx pyhelp.cxx)
      target_link_libraries(volumetric ${PYTHON_LIBRARIES})
      set_target_properties(volumetric PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET volumetric PROPERTY OUTPUT_NAME "volumetric.so")
            set_property(TARGET volumetric PROPERTY SUFFIX "")
      endif()

      # deviatoric
      add_library(deviatoric SHARED deviatoric_wrap.cxx deviatoric.cxx pyhelp.cxx)
      target_link_libraries(deviatoric ${PYTHON_LIBRARIES} shear surface)
      set_target_properties(deviatoric PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET deviatoric PROPERTY OUTPUT_NAME "deviatoric.so")
            set_property(TARGET deviatoric PROPERTY SUFFIX "")
      endif()

      # shear
      add_library(shear SHARED shear_wrap.cxx shear.cxx pyhelp.cxx)
      target_link_libraries(shear ${PYTHON_LIBRARIES})
      set_target_properties(shear PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET shear PROPERTY OUTPUT_NAME "shear.so")
            set_property(TARGET shear PROPERTY SUFFIX "")
      endif()

      # surface
      add_library(surface SHARED surface_wrap.cxx surface.cxx pyhelp.cxx)
      target_link_libraries(surface ${PYTHON_LIBRARIES})
      set_target_properties(surface PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET surface PROPERTY OUTPUT_NAME "surface.so")
            set_property(TARGET surface PROPERTY SUFFIX "")
      endif()
else()
      # Everything in a big shared library
      add_library(neml SHARED 
            neml.cxx 
            volumetric.cxx 
            deviatoric.cxx
            shear.cxx
            surface.cxx)
endif()
