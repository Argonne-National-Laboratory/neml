add_executable(report report.cxx)
target_include_directories(report PRIVATE "../../include")
target_link_libraries(report neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})
install(TARGETS report)
