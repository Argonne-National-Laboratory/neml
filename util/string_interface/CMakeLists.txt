include_directories(${PROJECT_BINARY_DIR}/src)
add_executable(cxxstring cxxstring.cxx)
target_link_libraries(cxxstring neml)
