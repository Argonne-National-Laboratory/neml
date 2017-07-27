find_package(Doxygen REQUIRED)

configure_file(${CMAKE_SOURCE_DIR}/doc/doxygen/Doxyfile.in ${CMAKE_BINARY_DIR}/doc/doxygen/Doxyfile @ONLY)
add_custom_target(doc_doxygen ALL
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/doxygen
      COMMENT "Generate Doxygen docs"
      VERATIM)
