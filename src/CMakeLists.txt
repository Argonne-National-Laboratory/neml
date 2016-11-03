if (WRAP_PYTHON)
      # Each file as a separate module
      # math
      add_library(nemlmath SHARED nemlmath_wrap.cxx nemlmath.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(nemlmath PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET nemlmath PROPERTY OUTPUT_NAME "nemlmath.so")
            set_property(TARGET nemlmath PROPERTY SUFFIX "")
      endif()

      # volumetric
      add_library(volumetric SHARED volumetric_wrap.cxx volumetric.cxx pyhelp.cxx)
      target_link_libraries(volumetric ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(volumetric PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET volumetric PROPERTY OUTPUT_NAME "volumetric.so")
            set_property(TARGET volumetric PROPERTY SUFFIX "")
      endif()

      # deviatoric
      add_library(deviatoric SHARED deviatoric_wrap.cxx deviatoric.cxx pyhelp.cxx)
      target_link_libraries(deviatoric nemlmath shear surface hardening ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} shear surface)
      set_target_properties(deviatoric PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET deviatoric PROPERTY OUTPUT_NAME "deviatoric.so")
            set_property(TARGET deviatoric PROPERTY SUFFIX "")
      endif()

      # shear
      add_library(shear SHARED shear_wrap.cxx shear.cxx pyhelp.cxx)
      target_link_libraries(shear ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(shear PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET shear PROPERTY OUTPUT_NAME "shear.so")
            set_property(TARGET shear PROPERTY SUFFIX "")
      endif()

      # surface
      add_library(surface SHARED surface_wrap.cxx surface.cxx pyhelp.cxx)
      target_link_libraries(surface nemlmath  ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(surface PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET surface PROPERTY OUTPUT_NAME "surface.so")
            set_property(TARGET surface PROPERTY SUFFIX "")
      endif()

      # hardening
      add_library(hardening SHARED hardening_wrap.cxx hardening.cxx pyhelp.cxx)
      target_link_libraries(hardening nemlmath  ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(hardening PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET hardening PROPERTY OUTPUT_NAME "hardening.so")
            set_property(TARGET hardening PROPERTY SUFFIX "")
      endif()

      # actual material models
      add_library(neml SHARED neml_wrap.cxx neml.cxx pyhelp.cxx)
      target_link_libraries(neml nemlmath deviatoric volumetric ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(neml PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET neml PROPERTY OUTPUT_NAME "neml.so")
            set_property(TARGET neml PROPERTY SUFFIX "")
      endif()
else()
      # Everything in a big shared library
      add_library(neml SHARED 
            neml.cxx 
            volumetric.cxx 
            deviatoric.cxx
            shear.cxx
            surface.cxx
            hardening.cxx
            nemlmath.cxx)
      target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()
