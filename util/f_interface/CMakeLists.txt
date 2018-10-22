include_directories(${PROJECT_BINARY_DIR}/src)
add_executable(fsimple fsimple.f)
target_link_libraries(fsimple neml)
