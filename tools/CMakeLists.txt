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
INSTALL(
    CODE "EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} setup.py install --prefix=${CMAKE_INSTALL_PREFIX}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )"
)

