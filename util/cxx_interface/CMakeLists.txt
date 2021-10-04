add_executable(cxxsimple cxxsimple.cxx)
target_include_directories(cxxsimple PRIVATE "../../include")
target_link_libraries(cxxsimple neml)
