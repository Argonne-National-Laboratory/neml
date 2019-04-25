include_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(report report.cxx)
target_link_libraries(report libneml ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SOLVER_LIBRARIES} ${libxml++_LIBRARIES})
