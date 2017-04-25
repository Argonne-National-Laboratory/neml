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
      add_library(solvers SHARED solvers_wrap.cxx solvers.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(solvers nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES})
      set_target_properties(solvers PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET solvers PROPERTY OUTPUT_NAME "solvers.so")
            set_property(TARGET solvers PROPERTY SUFFIX "")
      endif()

      # elasticity
      add_library(elasticity SHARED interpolate.cxx interpolate_wrap.cxx elasticity_wrap.cxx elasticity.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(elasticity nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(elasticity PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET elasticity PROPERTY OUTPUT_NAME "elasticity.so")
            set_property(TARGET elasticity PROPERTY SUFFIX "")
      endif()

      # yield surfaces
      add_library(surfaces SHARED surfaces_wrap.cxx surfaces.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(surfaces nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(surfaces PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET surfaces PROPERTY OUTPUT_NAME "surfaces.so")
            set_property(TARGET surfaces PROPERTY SUFFIX "")
      endif()

      # hardening rules
      add_library(hardening SHARED interpolate.cxx interpolate_wrap.cxx hardening_wrap.cxx hardening.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(hardening nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(hardening PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET hardening PROPERTY OUTPUT_NAME "hardening.so")
            set_property(TARGET hardening PROPERTY SUFFIX "")
      endif()

      # rate independent flow rules
      add_library(ri_flow SHARED ri_flow_wrap.cxx ri_flow.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(ri_flow surfaces hardening nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(ri_flow PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET ri_flow PROPERTY OUTPUT_NAME "ri_flow.so")
            set_property(TARGET ri_flow PROPERTY SUFFIX "")
      endif()

      # viscoplastic flow stuff
      add_library(visco_flow SHARED interpolate.cxx interpolate_wrap.cxx visco_flow_wrap.cxx visco_flow.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(visco_flow surfaces hardening nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(visco_flow PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET visco_flow PROPERTY OUTPUT_NAME "visco_flow.so")
            set_property(TARGET visco_flow PROPERTY SUFFIX "")
      endif()

      # general flow stuff
      add_library(general_flow SHARED general_flow_wrap.cxx general_flow.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(general_flow surfaces hardening visco_flow elasticity nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(general_flow PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET general_flow PROPERTY OUTPUT_NAME "general_flow.so")
            set_property(TARGET general_flow PROPERTY SUFFIX "")
      endif()

      # actual material models
      add_library(neml SHARED solvers.cxx solvers_wrap.cxx neml_wrap.cxx neml.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(neml solvers ri_flow visco_flow general_flow nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(neml PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET neml PROPERTY OUTPUT_NAME "neml.so")
            set_property(TARGET neml PROPERTY SUFFIX "")
      endif()

      # parser
      add_library(parse SHARED parse_wrap.cxx parse.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(parse neml ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${libxml++_LIBRARIES})
      set_target_properties(parse PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET parse PROPERTY OUTPUT_NAME "parse.so")
            set_property(TARGET parse PROPERTY SUFFIX "")
      endif()

      # various interpolation schemes
      add_library(interpolate SHARED interpolate_wrap.cxx interpolate.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(interpolate nemlmath ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(interpolate PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET interpolate PROPERTY OUTPUT_NAME "interpolate.so")
            set_property(TARGET interpolate PROPERTY SUFFIX "")
      endif()

      # creep models
      add_library(creep SHARED creep_wrap.cxx creep.cxx nemlerror.cxx pyhelp.cxx)
      target_link_libraries(creep nemlmath interpolate solvers ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set_target_properties(creep PROPERTIES PREFIX "")
      if (APPLE)
            set_property(TARGET interpolate PROPERTY OUTPUT_NAME "creep.so")
            set_property(TARGET interpolate PROPERTY SUFFIX "")
      endif()

else()
      # Everything in a big shared library
      add_library(neml SHARED 
            neml.cxx 
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
            creep.cxx)
      target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})
endif()
