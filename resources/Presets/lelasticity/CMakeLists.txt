# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/BI.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIStress.xml
    ${CMAKE_CURRENT_BINARY_DIR}/BIStress.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIMotion.xml
    ${CMAKE_CURRENT_BINARY_DIR}/BIMotion.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/deltaSPH.xml
    ${CMAKE_CURRENT_BINARY_DIR}/deltaSPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/deltaSPH-simple.xml
    ${CMAKE_CURRENT_BINARY_DIR}/deltaSPH-simple.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/deltaSPH-full.xml
    ${CMAKE_CURRENT_BINARY_DIR}/deltaSPH-full.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/elasticBounce.xml
    ${CMAKE_CURRENT_BINARY_DIR}/elasticBounce.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/MLS.xml
    ${CMAKE_CURRENT_BINARY_DIR}/MLS.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/xSPH.xml
    ${CMAKE_CURRENT_BINARY_DIR}/xSPH.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BI.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIStress.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIStress.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BIMotion.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIMotion.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/deltaSPH.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/deltaSPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/deltaSPH-simple.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/deltaSPH.xml-simple.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/deltaSPH-full.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/deltaSPH-full.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/elasticBounce.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/elasticBounce.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/MLS.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/MLS.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/xSPH.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/xSPH.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIStress.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BIMotion.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/deltaSPH.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/deltaSPH-simple.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/deltaSPH-full.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/elasticBounce.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/MLS.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/xSPH.xml
)
SOURCE_GROUP("ResourcesGroup" FILES ${RESOURCES_SRCS})
INSTALL(
    FILES
        ${RESOURCES_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets/lelasticity
)
