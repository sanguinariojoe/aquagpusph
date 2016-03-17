SET(RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)

# ===================================================== #
# In place configuration                                #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)

configure_file(${EXAMPLE_ORIG_DIR}/templates/Fluids.xml
    ${EXAMPLE_DEST_DIR}/templates/Fluids.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Main.xml
    ${EXAMPLE_DEST_DIR}/templates/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Settings.xml
    ${EXAMPLE_DEST_DIR}/templates/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/SPH.xml
    ${EXAMPLE_DEST_DIR}/templates/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Time.xml
    ${EXAMPLE_DEST_DIR}/templates/Time.xml @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_1
    ${EXAMPLE_DEST_DIR}/doc/d38_1 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_2
    ${EXAMPLE_DEST_DIR}/doc/d38_2 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_3
    ${EXAMPLE_DEST_DIR}/doc/d38_3 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_4
    ${EXAMPLE_DEST_DIR}/doc/d38_4 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_5
    ${EXAMPLE_DEST_DIR}/doc/d38_5 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_6
    ${EXAMPLE_DEST_DIR}/doc/d38_6 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_7
    ${EXAMPLE_DEST_DIR}/doc/d38_7 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_8
    ${EXAMPLE_DEST_DIR}/doc/d38_8 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/snapshots_d_38_mm.JPG
    ${EXAMPLE_DEST_DIR}/doc/snapshots_d_38_mm.JPG @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_t.py
    ${EXAMPLE_DEST_DIR}/plot_t.py @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/Create.py
    ${EXAMPLE_DEST_DIR}/Create.py @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/run.sh
    ${EXAMPLE_DEST_DIR}/run.sh @ONLY)

# ===================================================== #
# Installable version (and tagets)                      #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_AUX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp)
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/2D/spheric_testcase10_waveimpact)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file(${EXAMPLE_ORIG_DIR}/templates/Fluids.xml
    ${EXAMPLE_AUX_DIR}/templates/Fluids.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Main.xml
    ${EXAMPLE_AUX_DIR}/templates/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Settings.xml
    ${EXAMPLE_AUX_DIR}/templates/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/SPH.xml
    ${EXAMPLE_AUX_DIR}/templates/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Time.xml
    ${EXAMPLE_AUX_DIR}/templates/Time.xml @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_1
    ${EXAMPLE_AUX_DIR}/doc/d38_1 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_2
    ${EXAMPLE_AUX_DIR}/doc/d38_2 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_3
    ${EXAMPLE_AUX_DIR}/doc/d38_3 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_4
    ${EXAMPLE_AUX_DIR}/doc/d38_4 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_5
    ${EXAMPLE_AUX_DIR}/doc/d38_5 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_6
    ${EXAMPLE_AUX_DIR}/doc/d38_6 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_7
    ${EXAMPLE_AUX_DIR}/doc/d38_7 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/d38_8
    ${EXAMPLE_AUX_DIR}/doc/d38_8 @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/snapshots_d_38_mm.JPG
    ${EXAMPLE_AUX_DIR}/doc/snapshots_d_38_mm.JPG @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_t.py
    ${EXAMPLE_AUX_DIR}/plot_t.py @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/Create.py
    ${EXAMPLE_AUX_DIR}/Create.py @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/run.sh
    ${EXAMPLE_AUX_DIR}/run.sh @ONLY)

INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/templates/Fluids.xml
        ${EXAMPLE_AUX_DIR}/templates/Main.xml
        ${EXAMPLE_AUX_DIR}/templates/Settings.xml
        ${EXAMPLE_AUX_DIR}/templates/SPH.xml
        ${EXAMPLE_AUX_DIR}/templates/Time.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}/templates/
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/Create.py
        ${EXAMPLE_AUX_DIR}/snapshots_d_38_mm.JPG
        ${EXAMPLE_AUX_DIR}/plot_t.py
        ${EXAMPLE_AUX_DIR}/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/doc/d38_1
        ${EXAMPLE_AUX_DIR}/doc/d38_2
        ${EXAMPLE_AUX_DIR}/doc/d38_3
        ${EXAMPLE_AUX_DIR}/doc/d38_4
        ${EXAMPLE_AUX_DIR}/doc/d38_5
        ${EXAMPLE_AUX_DIR}/doc/d38_6
        ${EXAMPLE_AUX_DIR}/doc/d38_7
        ${EXAMPLE_AUX_DIR}/doc/d38_8
        ${EXAMPLE_AUX_DIR}/doc/snapshots_d_38_mm.JPG
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)
