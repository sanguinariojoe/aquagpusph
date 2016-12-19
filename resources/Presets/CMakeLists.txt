# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(AQUAGPUSPH_ROOT_PATH ${CMAKE_BINARY_DIR})
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/basic.xml
    ${CMAKE_CURRENT_BINARY_DIR}/basic.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/cfd.xml
    ${CMAKE_CURRENT_BINARY_DIR}/cfd.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lelasticity.xml
    ${CMAKE_CURRENT_BINARY_DIR}/lelasticity.xml @ONLY)

# Create installable version of the file
SET(AQUAGPUSPH_ROOT_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR})
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/basic.xml  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/basic.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/cfd.xml  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/cfd.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lelasticity.xml  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/lelasticity.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/basic.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/cfd.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/lelasticity.xml
)
SOURCE_GROUP("ResourcesGroup" FILES ${RESOURCES_SRCS})
INSTALL(
    FILES
        ${RESOURCES_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets
)

# ===================================================== #
# Configure and install dependencies                    #
# ===================================================== #
ADD_SUBDIRECTORY(basic)
ADD_SUBDIRECTORY(cfd)
ADD_SUBDIRECTORY(lelasticity)