# ===================================================== #
# Configuration files                                   #
# ===================================================== #
# Create local version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/noLocalMem.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/noLocalMem.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BI.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/BI.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BINoSlip.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/BINoSlip.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/GP.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/GP.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/GPFreeSlip.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/GPFreeSlip.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/GPEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/GPEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/elasticBounce.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/elasticBounce.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/inlet.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/inlet.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/inletEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/inletEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/outlet.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/outlet.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/outletEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/outletEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/symmetry.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/symmetry.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/variableTimeStep.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/variableTimeStep.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/domain.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/energy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/energy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/fluidEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/fluidEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/forces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/forces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/pressureForces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/pressureForces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/viscousForces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/viscousForces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/id_inverse.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/id_inverse.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/motion.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/motion.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/energy.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/energy.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/forces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/forces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/pressureForces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/pressureForces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/viscousForces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/viscousForces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/performance.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/performance.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/timing.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/timing.report.xml @ONLY)

# Create installable version of the file
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/noLocalMem.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/noLocalMem.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BI.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BI.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/BINoSlip.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BINoSlip.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/GP.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/GP.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/GPFreeSlip.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/GPFreeSlip.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/GPEnergy.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/GPEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/elasticBounce.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/elasticBounce.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/inlet.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/inlet.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/inletEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/inletEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/outlet.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/outlet.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/outletEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/outletEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/symmetry.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/symmetry.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/variableTimeStep.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/variableTimeStep.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/domain.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/energy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/fluidEnergy.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/fluidEnergy.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/forces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/pressureForces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/pressureForces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/viscousForces.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/viscousForces.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/id_inverse.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/id_inverse.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/motion.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/motion.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/energy.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/forces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/pressureForces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/pressureForces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/viscousForces.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/viscousForces.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/performance.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/performance.report.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/timing.report.xml.in  
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/timing.report.xml @ONLY)

# ===================================================== #
# Install                                               #
# ===================================================== #
SET(RESOURCES_SRCS 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/noLocalMem.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BI.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/BINoSlip.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/GP.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/GPFreeSlip.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/GPEnergy.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/fixedGhostParticles.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/elasticBounce.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/inlet.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/inletEnergy.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/outlet.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/outletEnergy.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/symmetry.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/variableTimeStep.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/domain.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/pressureForces.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/viscousForces.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/id_inverse.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/fluidEnergy.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/motion.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/energy.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/forces.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/pressureForces.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/viscousForces.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/performance.report.xml
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/timing.report.xml
)
SOURCE_GROUP("ResourcesGroup" FILES ${RESOURCES_SRCS})
INSTALL(
    FILES
        ${RESOURCES_SRCS}
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets/cfd
)
