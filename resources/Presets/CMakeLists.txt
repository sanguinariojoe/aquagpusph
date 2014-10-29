# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfd.xml.in  
	${CMAKE_CURRENT_BINARY_DIR}/cfd.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain.xml.in  
	${CMAKE_CURRENT_BINARY_DIR}/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/motion.xml.in  
	${CMAKE_CURRENT_BINARY_DIR}/motion.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfd.xml.in  
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/cfd.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain.xml.in  
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/motion.xml.in  
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/motion.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
	${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/cfd.xml
	${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml
	${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/motion.xml
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
