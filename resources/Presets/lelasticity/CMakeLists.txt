# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/BI.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIFixed.xml
    ${CMAKE_CURRENT_BINARY_DIR}/BIFixed.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIFree.xml
    ${CMAKE_CURRENT_BINARY_DIR}/BIFree.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BI.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIFixed.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIFixed.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIFree.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIFree.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIFixed.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIFree.xml
)
SOURCE_GROUP("ResourcesGroup" FILES ${RESOURCES_SRCS})
INSTALL(
    FILES
        ${RESOURCES_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets/lelasticity
)
