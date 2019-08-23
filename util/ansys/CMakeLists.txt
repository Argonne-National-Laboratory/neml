include_directories(${PROJECT_BINARY_DIR}/src)
add_executable(report report.cxx)
target_link_libraries(neml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})
