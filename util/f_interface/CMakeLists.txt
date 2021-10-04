add_executable(fsimple fsimple.f)
target_include_directories(cxxsimple PRIVATE "../../include")
target_link_libraries(fsimple neml)
