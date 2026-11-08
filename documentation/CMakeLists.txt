# API documentation: `cmake --build <build> --target doc`
# writes HTML and LaTeX to <build>/documentation
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(doxygen.cfg.in "${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg" @ONLY)
  # INPUT paths in doxygen.cfg are relative to the source root
  add_custom_target(doc
                    COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg"
                    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                    COMMENT "Generating API documentation with Doxygen"
                    VERBATIM)
endif()
