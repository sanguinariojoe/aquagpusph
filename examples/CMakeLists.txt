SET(RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)

IF(AQUAGPUSPH_3D)

# ===================================================== #
# perezrojas_etal_stab_2012 (Local)                     #
# ===================================================== #
# Folders
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake/perezrojas_etal_stab_2012)
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/perezrojas_etal_stab_2012)
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
FILE(MAKE_DIRECTORY ${EXAMPLE_DEST_DIR})
# XML definition files
configure_file(${EXAMPLE_ORIG_DIR}/Move/6DOF.xml
	${EXAMPLE_DEST_DIR}/Move/6DOF.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Fluids.xml
	${EXAMPLE_DEST_DIR}/Fluids.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Main.xml
	${EXAMPLE_DEST_DIR}/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Movements.xml
	${EXAMPLE_DEST_DIR}/Movements.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Settings.xml
	${EXAMPLE_DEST_DIR}/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/SPH.xml
	${EXAMPLE_DEST_DIR}/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Time.xml
	${EXAMPLE_DEST_DIR}/Time.xml @ONLY)
# Other simulation data
configure_file(${EXAMPLE_ORIG_DIR}/doc/plot.py
	${EXAMPLE_DEST_DIR}/doc/plot.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/SOUTOIGLESIAS_BOTIA_SPHERIC_TESTCASE9_TLD.pdf
	${EXAMPLE_DEST_DIR}/doc/SOUTOIGLESIAS_BOTIA_SPHERIC_TESTCASE9_TLD.pdf @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Move/T_1-94_A100mm_water.dat
	${EXAMPLE_DEST_DIR}/Move/T_1-94_A100mm_water.dat @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Move/move.py
	${EXAMPLE_DEST_DIR}/Move/move.py @ONLY)
# Particles generation script
configure_file(${EXAMPLE_ORIG_DIR}/Create.py
	${EXAMPLE_DEST_DIR}/Create.py @ONLY)
# Execution script
configure_file(${EXAMPLE_ORIG_DIR}/run.sh
	${EXAMPLE_DEST_DIR}/run.sh @ONLY)
# Target
ADD_CUSTOM_TARGET(perezrojas_etal_stab_2012 ALL
	COMMAND ${PYTHON_EXECUTABLE} ./Create.py
	COMMAND chmod a+x ./run.sh
	WORKING_DIRECTORY ${EXAMPLE_DEST_DIR}
)
ADD_DEPENDENCIES( perezrojas_etal_stab_2012 AQUAgpusph )

# ===================================================== #
# perezrojas_etal_stab_2012 (Install)                   #
# ===================================================== #
# Folders
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake/perezrojas_etal_stab_2012)
SET(EXAMPLE_AUX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/perezrojas_etal_stab_2012)
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/perezrojas_etal_stab_2012)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
# XML definition files
configure_file(${EXAMPLE_ORIG_DIR}/Move/6DOF.xml
	${EXAMPLE_AUX_DIR}/Move/6DOF.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Fluids.xml
	${EXAMPLE_AUX_DIR}/Fluids.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Main.xml
	${EXAMPLE_AUX_DIR}/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Movements.xml
	${EXAMPLE_AUX_DIR}/Movements.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Settings.xml
	${EXAMPLE_AUX_DIR}/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/SPH.xml
	${EXAMPLE_AUX_DIR}/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Time.xml
	${EXAMPLE_AUX_DIR}/Time.xml @ONLY)
# Other simulation data
configure_file(${EXAMPLE_ORIG_DIR}/doc/plot.py
	${EXAMPLE_AUX_DIR}/doc/plot.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/SOUTOIGLESIAS_BOTIA_SPHERIC_TESTCASE9_TLD.pdf
	${EXAMPLE_AUX_DIR}/doc/SOUTOIGLESIAS_BOTIA_SPHERIC_TESTCASE9_TLD.pdf @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Move/T_1-94_A100mm_water.dat
	${EXAMPLE_AUX_DIR}/Move/T_1-94_A100mm_water.dat @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Move/move.py
	${EXAMPLE_AUX_DIR}/Move/move.py @ONLY)
# Particles generation script
configure_file(${EXAMPLE_ORIG_DIR}/Create.py
	${EXAMPLE_AUX_DIR}/Create.py @ONLY)
# Execution script
configure_file(${EXAMPLE_ORIG_DIR}/run.sh
	${EXAMPLE_AUX_DIR}/run.sh @ONLY)
# Install
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/perezrojas_etal_stab_2012/Fluid.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/Fluids.xml
        ${EXAMPLE_AUX_DIR}/Main.xml
        ${EXAMPLE_AUX_DIR}/Movements.xml
        ${EXAMPLE_AUX_DIR}/Settings.xml
        ${EXAMPLE_AUX_DIR}/SPH.xml
        ${EXAMPLE_AUX_DIR}/Time.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/Create.py
        ${EXAMPLE_AUX_DIR}/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/Move/T_1-94_A100mm_water.dat
        ${EXAMPLE_AUX_DIR}/Move/6DOF.xml
        ${EXAMPLE_AUX_DIR}/Move/move.py
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/doc/SOUTOIGLESIAS_BOTIA_SPHERIC_TESTCASE9_TLD.pdf
        ${EXAMPLE_AUX_DIR}/doc/plot.py
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)

ELSE(AQUAGPUSPH_3D)

# ===================================================== #
# LateralWater_1x_DeLeffe (Local)                       #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe)
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Move/6DOF.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Move/6DOF.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Movements.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Movements.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/run.sh @ONLY)
# Target
ADD_CUSTOM_TARGET(LateralWater_1x_DeLeffe ALL
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_DeLeffe/Create.py
	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/Fluid.dat ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_DeLeffe/Fluid.dat
	COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/Fluid.dat
	COMMAND chmod a+x ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_DeLeffe/run.sh
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
ADD_DEPENDENCIES( LateralWater_1x_DeLeffe AQUAgpusph2D )

# ===================================================== #
# LateralWater_1x_DeLeffe (Install)                     #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/LateralWater_1x_DeLeffe)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Move/6DOF.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Move/6DOF.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Movements.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Movements.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_DeLeffe/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/run.sh @ONLY)
# Install
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Fluid.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/Move/6DOF.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/doc/plot.py
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_DeLeffe/doc/lateral_water_1x.txt
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Fluids.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Main.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Movements.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Sensors.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Settings.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/SPH.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Time.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_DeLeffe/Move/6DOF.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)

# ===================================================== #
# LateralWater_1x_Fix (Local)                       #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix)
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Move/6DOF.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Move/6DOF.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Movements.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Movements.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/run.sh @ONLY)
# Target
ADD_CUSTOM_TARGET(LateralWater_1x_Fix ALL
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_Fix/Create.py
	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/Fluid.dat ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_Fix/Fluid.dat
	COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/Fluid.dat
	COMMAND chmod a+x ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_Fix/run.sh
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
# Since generated particles file is the same, wait for previous
# examples have been built before start working with this.
ADD_DEPENDENCIES( LateralWater_1x_Fix AQUAgpusph2D LateralWater_1x_DeLeffe )

# ===================================================== #
# LateralWater_1x_Fix (Install)                     #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/LateralWater_1x_Fix)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Move/6DOF.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Move/6DOF.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Movements.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Movements.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Fix/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/run.sh @ONLY)
# Install
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Fluid.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/Move/6DOF.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/doc/plot.py
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Fix/doc/lateral_water_1x.txt
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Fluids.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Main.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Movements.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Sensors.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Settings.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/SPH.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Time.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Fix/Move/6DOF.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)

# ===================================================== #
# LateralWater_1x_Ghost (Local)                         #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost)
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Move/6DOF.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Move/6DOF.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/GhostParticles.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/GhostParticles.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Movements.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Movements.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/run.sh @ONLY)
# Target
ADD_CUSTOM_TARGET(LateralWater_1x_Ghost ALL
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_Ghost/Create.py
	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/Fluid.dat ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_Ghost/Fluid.dat
	COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/Fluid.dat
	COMMAND chmod a+x ${CMAKE_CURRENT_BINARY_DIR}/LateralWater_1x_Ghost/run.sh
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
# Since generated particles file is the same, wait for previous
# examples have been built before start working with this.
ADD_DEPENDENCIES( LateralWater_1x_Ghost AQUAgpusph2D LateralWater_1x_DeLeffe LateralWater_1x_Fix )

# ===================================================== #
# LateralWater_1x_Ghost (Install)                       #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/LateralWater_1x_Ghost)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Move/6DOF.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Move/6DOF.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/GhostParticles.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/GhostParticles.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Movements.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Movements.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/LateralWater_1x_Ghost/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/run.sh @ONLY)
# Install
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Fluid.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/Move/6DOF.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/doc/plot.py
        ${CMAKE_CURRENT_SOURCE_DIR}/LateralWater_1x_Ghost/doc/lateral_water_1x.txt
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Fluids.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/GhostParticles.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Main.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Movements.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Sensors.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Settings.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/SPH.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Time.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/LateralWater_1x_Ghost/Move/6DOF.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}/Move
)

# ===================================================== #
# lobovsky_etal_jfs_2014 (Local)                        #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014)
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/run.sh @ONLY)
# Target
ADD_CUSTOM_TARGET(lobovsky_etal_jfs_2014 ALL
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/lobovsky_etal_jfs_2014/Create.py
	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/Fluid.dat ${CMAKE_CURRENT_BINARY_DIR}/lobovsky_etal_jfs_2014/Fluid.dat
	COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/Fluid.dat
	COMMAND chmod a+x ${CMAKE_CURRENT_BINARY_DIR}/lobovsky_etal_jfs_2014/run.sh
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
ADD_DEPENDENCIES( lobovsky_etal_jfs_2014 AQUAgpusph2D LateralWater_1x_DeLeffe LateralWater_1x_Ghost )

# ===================================================== #
# lobovsky_etal_jfs_2014 (Install)                      #
# ===================================================== #
# Folders
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/lobovsky_etal_jfs_2014)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
# XML definition files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Fluids.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Fluids.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Main.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Main.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Sensors.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Sensors.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Settings.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Settings.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/SPH.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/SPH.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/Time.xml
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Time.xml @ONLY)
# Execution script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cMake/lobovsky_etal_jfs_2014/run.sh
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/run.sh @ONLY)
# Install
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/Fluid.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/sensor1.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/sensor2.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/sensor3.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/sensor4.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/sensors.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/performance.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_01_ETSIN_Front_Dressler_Martin_2.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_01_ETSIN_Front_Dressler_Martin_3.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_01_ETSIN_Front_Dressler_Martin_4.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_01_ETSIN_Front_Dressler_Martin_5.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_01_ETSIN_Front_Dressler_Martin_6.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_01_ETSIN_Front_Dressler_Martin_7.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_02_ETSIN_Front_Hu_Koshizuka_3.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_02_ETSIN_Front_Hu_Koshizuka_4.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig12_02_ETSIN_Front_Hu_Koshizuka_5.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig20_filtered_1.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig20_filtered_2.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig20_filtered_3.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig29_filtered_1.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig29_filtered_2.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig29_filtered_3.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig29_filtered_4.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig30_filtered_1.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig30_filtered_2.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig30_filtered_3.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig30_filtered_4.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig30_filtered_5.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig31_filtered_1.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig31_filtered_2.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig31_filtered_3.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/lobovsky_etal_jfs_2014/doc/Fig31_filtered_4.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Fluids.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Main.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Sensors.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Settings.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/SPH.xml
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/Time.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}
)
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp/lobovsky_etal_jfs_2014/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)

ENDIF(AQUAGPUSPH_3D)
