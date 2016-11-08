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

      # Nonlinear solvers
      # elasticity
      add_library(solvers SHARED solvers_wrap.cxx solvers.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(solvers nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(solvers PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET solvers PROPERTY OUTPUT_NAME "solvers.so")
            set_property(TARGET solvers PROPERTY SUFFIX "")
      endif()

      # elasticity
      add_library(elasticity SHARED elasticity_wrap.cxx elasticity.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(elasticity nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(elasticity PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET elasticity PROPERTY OUTPUT_NAME "elasticity.so")
            set_property(TARGET elasticity PROPERTY SUFFIX "")
      endif()

      # actual material models
      add_library(neml SHARED neml_wrap.cxx neml.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(neml nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(neml PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET neml PROPERTY OUTPUT_NAME "neml.so")
            set_property(TARGET neml PROPERTY SUFFIX "")
      endif()
else()
      # Everything in a big shared library
      add_library(neml SHARED 
            neml.cxx 
            nemlmath.cxx
            solvers.cxx
            nemlerror.cxx
            elasticity.cxx)
      target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()
