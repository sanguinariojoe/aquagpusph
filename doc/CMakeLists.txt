# ===================================================== #
# Doxygen configuration file                            #
# ===================================================== #
# Included source directories
SET(DOXYGEN_SOURCE_DIR
	${CMAKE_SOURCE_DIR}/src
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/resources
)
STRING(REGEX REPLACE ";" " " DOXYGEN_INPUT_LIST "${DOXYGEN_SOURCE_DIR}")

# exclude some subdirectories
set(DOXYGEN_EXCLUDE_DIR
)
STRING(REGEX REPLACE ";" " " DOXYGEN_EXCLUDE_LIST "${DOXYGEN_EXCLUDE_DIR}")

# General settings
SET(DOXYGEN_OUTPUT_DIR ${CMAKE_BINARY_DIR}/doc/Doxygen)
SET(DOXYGEN_IMAGE_PATH )
SET(DOXYGEN_LAYOUT_FILE )

# Configuration file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/Doxyfile.in
	${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

# ===================================================== #
# Documentation gen target                              #
# ===================================================== #
ADD_CUSTOM_TARGET(DoxygenDoc ALL
	COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_CURRENT_BINARY_DIR}/Doxygen)
ADD_DEPENDENCIES( DoxygenDoc AQUAgpusph )

# ===================================================== #
# Install                                               #
# ===================================================== #
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/Doxygen
    DESTINATION
        ${CMAKE_INSTALL_DOCDIR}
)
