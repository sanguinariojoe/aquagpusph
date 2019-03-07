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
configure_file(${EXAMPLE_ORIG_DIR}/templates/Sensors.xml
    ${EXAMPLE_DEST_DIR}/templates/Sensors.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Settings.xml
    ${EXAMPLE_DEST_DIR}/templates/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/SPH.xml
    ${EXAMPLE_DEST_DIR}/templates/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Time.xml
    ${EXAMPLE_DEST_DIR}/templates/Time.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/h_sensor.cl
    ${EXAMPLE_DEST_DIR}/templates/h_sensor.cl @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/plot_p.py
    ${EXAMPLE_DEST_DIR}/plot_p.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_h.py
    ${EXAMPLE_DEST_DIR}/plot_h.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_t.py
    ${EXAMPLE_DEST_DIR}/plot_t.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/test_case_2_v1p1.pdf
    ${EXAMPLE_DEST_DIR}/doc/test_case_2_v1p1.pdf COPYONLY)
configure_file(${EXAMPLE_ORIG_DIR}/test_case_2_exp_data.dat
    ${EXAMPLE_DEST_DIR}/test_case_2_exp_data.dat COPYONLY)

configure_file(${EXAMPLE_ORIG_DIR}/Create.py
    ${EXAMPLE_DEST_DIR}/Create.py @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/run.sh
    ${EXAMPLE_DEST_DIR}/run.sh @ONLY)

# ===================================================== #
# Installable version (and tagets)                      #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_AUX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp)
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/3D/spheric_testcase2_dambreak)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file(${EXAMPLE_ORIG_DIR}/templates/Fluids.xml
    ${EXAMPLE_AUX_DIR}/templates/Fluids.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Main.xml
    ${EXAMPLE_AUX_DIR}/templates/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Sensors.xml
    ${EXAMPLE_AUX_DIR}/templates/Sensors.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Settings.xml
    ${EXAMPLE_AUX_DIR}/templates/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/SPH.xml
    ${EXAMPLE_AUX_DIR}/templates/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Time.xml
    ${EXAMPLE_AUX_DIR}/templates/Time.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/h_sensor.cl
    ${EXAMPLE_AUX_DIR}/templates/h_sensor.cl @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/plot_p.py
    ${EXAMPLE_AUX_DIR}/plot_p.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_h.py
    ${EXAMPLE_AUX_DIR}/plot_h.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_t.py
    ${EXAMPLE_AUX_DIR}/plot_t.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/doc/test_case_2_v1p1.pdf
    ${EXAMPLE_AUX_DIR}/doc/test_case_2_v1p1.pdf COPYONLY)
configure_file(${EXAMPLE_ORIG_DIR}/test_case_2_exp_data.dat
    ${EXAMPLE_AUX_DIR}/test_case_2_exp_data.dat COPYONLY)

configure_file(${EXAMPLE_ORIG_DIR}/Create.py
    ${EXAMPLE_AUX_DIR}/Create.py @ONLY)

configure_file(${EXAMPLE_ORIG_DIR}/run.sh
    ${EXAMPLE_AUX_DIR}/run.sh @ONLY)

INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/templates/Fluids.xml
        ${EXAMPLE_AUX_DIR}/templates/Main.xml
        ${EXAMPLE_AUX_DIR}/templates/Sensors.xml
        ${EXAMPLE_AUX_DIR}/templates/Settings.xml
        ${EXAMPLE_AUX_DIR}/templates/SPH.xml
        ${EXAMPLE_AUX_DIR}/templates/Time.xml
        ${EXAMPLE_AUX_DIR}/templates/h_sensor.cl
    DESTINATION
        ${EXAMPLE_DEST_DIR}/templates/
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/Create.py
        ${EXAMPLE_AUX_DIR}/test_case_2_exp_data.dat
        ${EXAMPLE_AUX_DIR}/plot_p.py
        ${EXAMPLE_AUX_DIR}/plot_h.py
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
        ${EXAMPLE_AUX_DIR}/doc/test_case_2_v1p1.pdf
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)
