# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenCLMain.xml.in  
	${CMAKE_CURRENT_BINARY_DIR}/OpenCLMain.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenCLMain.xml.in  
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/OpenCLMain.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
	${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/OpenCLMain.xml
)
SOURCE_GROUP("ResourcesGroup" FILES ${RESOURCES_SRCS})
INSTALL(
    FILES
        ${RESOURCES_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources
)

INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/OpenCL
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources
)
