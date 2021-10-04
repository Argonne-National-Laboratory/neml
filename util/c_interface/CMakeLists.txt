add_executable(csimple csimple.c)
target_include_directories(csimple PRIVATE "../../include")
target_link_libraries(csimple neml)
