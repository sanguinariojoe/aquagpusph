# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/densityClamp.xml
    ${CMAKE_CURRENT_BINARY_DIR}/densityClamp.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/domain.xml
    ${CMAKE_CURRENT_BINARY_DIR}/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/id_inverse.xml
    ${CMAKE_CURRENT_BINARY_DIR}/id_inverse.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/cubicSpline.xml
    ${CMAKE_CURRENT_BINARY_DIR}/kernels/cubicSpline.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/gauss.xml
    ${CMAKE_CURRENT_BINARY_DIR}/kernels/gauss.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/spiky.xml
    ${CMAKE_CURRENT_BINARY_DIR}/kernels/spiky.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/wendland.xml
    ${CMAKE_CURRENT_BINARY_DIR}/kernels/wendland.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/MLS.xml
    ${CMAKE_CURRENT_BINARY_DIR}/MLS.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/multiresolution.xml
    ${CMAKE_CURRENT_BINARY_DIR}/multiresolution.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/multiresolution/cube.xml
    ${CMAKE_CURRENT_BINARY_DIR}/multiresolution/cube.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/multiresolution/sphere.xml
    ${CMAKE_CURRENT_BINARY_DIR}/multiresolution/sphere.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/noLocalMem.xml
    ${CMAKE_CURRENT_BINARY_DIR}/noLocalMem.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/setBuffer.xml
    ${CMAKE_CURRENT_BINARY_DIR}/setBuffer.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/performance.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/performance.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/timing.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/timing.report.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/densityClamp.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/densityClamp.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/domain.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/id_inverse.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/id_inverse.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/cubicSpline.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/cubicSpline.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/gauss.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/gauss.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/spiky.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/spiky.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/kernels/wendland.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/wendland.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/MLS.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/MLS.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/multiresolution.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/multiresolution.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/multiresolution/cube.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/multiresolution/cube.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/multiresolution/sphere.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/multiresolution/sphere.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/noLocalMem.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/noLocalMem.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/setBuffer.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/setBuffer.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/performance.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/performance.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/timing.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/timing.report.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/densityClamp.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/id_inverse.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/MLS.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/multiresolution.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/noLocalMem.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/setBuffer.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/performance.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/timing.report.xml
)
SOURCE_GROUP("ResourcesGroup" FILES ${RESOURCES_SRCS})
INSTALL(
    FILES
        ${RESOURCES_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets/basic
)
SET(RESOURCES_KERNELS_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/cubicSpline.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/gauss.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/wendland.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/kernels/spiky.xml
)
SOURCE_GROUP("ResourcesKernelsGroup" FILES ${RESOURCES_KERNELS_SRCS})
INSTALL(
    FILES
        ${RESOURCES_KERNELS_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets/basic/kernels
)
SET(RESOURCES_MULTIRESOLUTION_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/multiresolution/cube.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/multiresolution/sphere.xml
)
SOURCE_GROUP("ResourcesMultiresolutionGroup" FILES ${RESOURCES_MULTIRESOLUTION_SRCS})
INSTALL(
    FILES
        ${RESOURCES_MULTIRESOLUTION_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets/basic/multiresolution
)
