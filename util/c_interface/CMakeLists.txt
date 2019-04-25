include_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(csimple csimple.c)
target_link_libraries(csimple libneml)
