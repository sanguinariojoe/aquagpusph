# ===================================================== #
# Build                                                 #
# ===================================================== #
ADD_CUSTOM_TARGET(tools ALL
	COMMAND ${PYTHON_EXECUTABLE} setup.py build
)
ADD_DEPENDENCIES( tools AQUAgpusph )

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(PYTHON_INSTALL_PREFIX "")
IF(NOT CMAKE_INSTALL_PREFIX STREQUAL /usr/local)
	SET(PYTHON_INSTALL_PREFIX --prefix=${CMAKE_INSTALL_PREFIX})
ENDIF(NOT CMAKE_INSTALL_PREFIX STREQUAL /usr/local)
INSTALL(
    CODE "EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} setup.py install ${PYTHON_INSTALL_PREFIX}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )"
)

