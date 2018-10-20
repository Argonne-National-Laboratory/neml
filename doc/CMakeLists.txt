option(SPHINX_HTML "Build html help with Sphinx" ON)
option(SPHINX_PDF "Build PDF manual with Sphinx" ON)
option(DOXYGEN "Build doxygen documentation" ON)

if (SPHINX_HTML OR SPHINX_PDF)
      find_package(Sphinx REQUIRED)
      set(USE_SPHINX ON)
endif()

if (USE_SPHINX OR DOXYGEN)
      find_package(Doxygen REQUIRED)
      set(DOXYGEN_CONF_FILE ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile.in)
      
      add_custom_target(doc-doxygen ${DOXYGEN_EXECUTABLE} 
            ${DOXYGEN_CONF_FILE})
endif()

if (SPHINX_HTML)
      set(SPHINX_HTML_DIR ${PROJECT_BINARY_DIR}/doc/html)
      set(SPHINX_SRC_DIR ${PROJECT_BINARY_DIR}/doc/sphinx)
      add_custom_target(doc-html ${SPHINX_EXECUTABLE}
            -b html
            ${SPHINX_SRC_DIR}
            ${SPHINX_HTML_DIR}
            )
endif()

if (SPHINX_PDF)
      set(SPHINX_PDF_DIR ${PROJECT_BINARY_DIR}/doc/pdf)
      set(SPHINX_SRC_DIR ${PROJECT_BINARY_DIR}/doc/sphinx)
      add_custom_target(doc-pdf ${SPHINX_EXECUTABLE}
            -b latex
            ${SPHINX_SRC_DIR}
            ${SPHINX_PDF_DIR}
            )
endif()
