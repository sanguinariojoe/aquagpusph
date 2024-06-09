# ===================================================== #
# Copy everythin in the CMAKE_BINARY_DIR                #
# ===================================================== #
FILE(GLOB_RECURSE FNAMES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*")

FOREACH(FNAME ${FNAMES})
    IF(${FNAME} MATCHES "CMakeLists.txt")
        CONTINUE()
    ENDIF()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${FNAME}
        ${CMAKE_CURRENT_BINARY_DIR}/${FNAME} COPYONLY)
ENDFOREACH()

# ===================================================== #
# Build                                                 #
# ===================================================== #
ADD_CUSTOM_TARGET(tools ALL
    COMMAND ${Python_EXECUTABLE} setup.py build
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
ADD_DEPENDENCIES(tools AQUAgpusph)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(PYTHON_INSTALL_PREFIX "")
IF(NOT CMAKE_INSTALL_PREFIX STREQUAL /usr/local)
	SET(PYTHON_INSTALL_PREFIX --prefix=${CMAKE_INSTALL_PREFIX})
ENDIF(NOT CMAKE_INSTALL_PREFIX STREQUAL /usr/local)
INSTALL(
    CODE "EXECUTE_PROCESS(
        COMMAND ${Python_EXECUTABLE} setup.py install ${PYTHON_INSTALL_PREFIX}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )"
)

