# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfd.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/cfd.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/energy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/energy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/forces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/forces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/motion.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/motion.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/energy.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/energy.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/forces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/forces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/timing.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/timing.report.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfd.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/cfd.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/energy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/forces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/motion.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/motion.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/energy.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/forces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/timing.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/timing.report.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/cfd.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/motion.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/timing.report.xml
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
