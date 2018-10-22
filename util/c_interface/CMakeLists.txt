include_directories(${PROJECT_BINARY_DIR}/src)
add_executable(csimple csimple.c)
target_link_libraries(csimple neml)
