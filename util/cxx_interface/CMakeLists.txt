include_directories(${PROJECT_BINARY_DIR}/src)
add_executable(cxxsimple cxxsimple.cxx)
target_link_libraries(cxxsimple neml)
