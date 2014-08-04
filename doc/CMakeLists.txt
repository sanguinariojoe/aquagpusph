# ===================================================== #
# Doxygen configuration file                            #
# ===================================================== #
# Included source directories
SET(DOXYGEN_SOURCE_DIR 
	${CMAKE_SOURCE_DIR}/src
	${CMAKE_SOURCE_DIR}/src/CalcServer
	${CMAKE_SOURCE_DIR}/src/CalcServer/Boundary
	${CMAKE_SOURCE_DIR}/src/CalcServer/Movements
	${CMAKE_SOURCE_DIR}/src/CalcServer/Portal
	${CMAKE_SOURCE_DIR}/src/Input
	${CMAKE_SOURCE_DIR}/src/Tokenizer
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/include/CalcServer
	${CMAKE_SOURCE_DIR}/include/CalcServer/Boundary
	${CMAKE_SOURCE_DIR}/include/CalcServer/Movements
	${CMAKE_SOURCE_DIR}/include/CalcServer/Portal
	${CMAKE_SOURCE_DIR}/include/Input
	${CMAKE_SOURCE_DIR}/include/Tokenizer
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
IF(AQUAGPUSPH_3D)
	ADD_DEPENDENCIES( DoxygenDoc AQUAgpusph )
ELSE(AQUAGPUSPH_3D)
	ADD_DEPENDENCIES( DoxygenDoc AQUAgpusph2D )
ENDIF(AQUAGPUSPH_3D)

# ===================================================== #
# Install                                               #
# ===================================================== #
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/Doxygen
    DESTINATION
        ${CMAKE_INSTALL_DOCDIR}
)
