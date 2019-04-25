include_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(fsimple fsimple.f)
target_link_libraries(fsimple libneml)
