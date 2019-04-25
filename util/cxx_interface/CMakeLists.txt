include_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(cxxsimple cxxsimple.cxx)
target_link_libraries(cxxsimple libneml)
