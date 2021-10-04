add_executable(cxxstring cxxstring.cxx)
target_include_directories(cxxstring PRIVATE "../../include")
target_link_libraries(cxxstring neml)
