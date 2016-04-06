SET(RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)

# ===================================================== #
# In place configuration                                #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)

configure_file(${EXAMPLE_ORIG_DIR}/templates/BI.xml
    ${EXAMPLE_DEST_DIR}/templates/BI.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Box.xml
    ${EXAMPLE_DEST_DIR}/templates/Box.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Fluid.xml
    ${EXAMPLE_DEST_DIR}/templates/Fluid.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/GP.xml
    ${EXAMPLE_DEST_DIR}/templates/GP.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Main.xml
    ${EXAMPLE_DEST_DIR}/templates/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Motion.xml
    ${EXAMPLE_DEST_DIR}/templates/Motion.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Settings.xml
    ${EXAMPLE_DEST_DIR}/templates/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/SPH.xml
    ${EXAMPLE_DEST_DIR}/templates/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Square.xml
    ${EXAMPLE_DEST_DIR}/templates/Square.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Time.xml
    ${EXAMPLE_DEST_DIR}/templates/Time.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/plot_f.py
    ${EXAMPLE_DEST_DIR}/templates/plot_f.py @ONLY)

file(COPY ${EXAMPLE_ORIG_DIR}/doc/Force_Re100.dat
    DESTINATION ${EXAMPLE_DEST_DIR}/doc/)
file(COPY ${EXAMPLE_ORIG_DIR}/doc/Motion_Body.dat
    DESTINATION ${EXAMPLE_DEST_DIR}/doc/)

configure_file(${EXAMPLE_ORIG_DIR}/Create.py
    ${EXAMPLE_DEST_DIR}/Create.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Motion.py
    ${EXAMPLE_DEST_DIR}/Motion.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_t.py
    ${EXAMPLE_DEST_DIR}/plot_t.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/run.sh
    ${EXAMPLE_DEST_DIR}/run.sh @ONLY)

# ===================================================== #
# Installable version (and tagets)                      #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_AUX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp)
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/2D/cylinder_inside_channel)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file(${EXAMPLE_ORIG_DIR}/templates/BI.xml
    ${EXAMPLE_AUX_DIR}/templates/BI.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Box.xml
    ${EXAMPLE_AUX_DIR}/templates/Box.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Fluid.xml
    ${EXAMPLE_AUX_DIR}/templates/Fluid.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/GP.xml
    ${EXAMPLE_AUX_DIR}/templates/GP.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Main.xml
    ${EXAMPLE_AUX_DIR}/templates/Main.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Motion.xml
    ${EXAMPLE_AUX_DIR}/templates/Motion.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Settings.xml
    ${EXAMPLE_AUX_DIR}/templates/Settings.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/SPH.xml
    ${EXAMPLE_AUX_DIR}/templates/SPH.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Square.xml
    ${EXAMPLE_AUX_DIR}/templates/Square.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/Time.xml
    ${EXAMPLE_AUX_DIR}/templates/Time.xml @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/templates/plot_f.py
    ${EXAMPLE_AUX_DIR}/templates/plot_f.py @ONLY)

file(COPY ${EXAMPLE_ORIG_DIR}/doc/Force_Re100.dat
    DESTINATION ${EXAMPLE_AUX_DIR}/doc/)
file(COPY ${EXAMPLE_ORIG_DIR}/doc/Motion_Body.dat
    DESTINATION ${EXAMPLE_AUX_DIR}/doc/)

configure_file(${EXAMPLE_ORIG_DIR}/Create.py
    ${EXAMPLE_AUX_DIR}/Create.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/Motion.py
    ${EXAMPLE_AUX_DIR}/Motion.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/plot_t.py
    ${EXAMPLE_AUX_DIR}/plot_t.py @ONLY)
configure_file(${EXAMPLE_ORIG_DIR}/run.sh
    ${EXAMPLE_AUX_DIR}/run.sh @ONLY)

INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/templates/BI.xml
        ${EXAMPLE_AUX_DIR}/templates/Box.xml
        ${EXAMPLE_AUX_DIR}/templates/Fluid.xml
        ${EXAMPLE_AUX_DIR}/templates/GP.xml
        ${EXAMPLE_AUX_DIR}/templates/Main.xml
        ${EXAMPLE_AUX_DIR}/templates/Motion.xml
        ${EXAMPLE_AUX_DIR}/templates/Settings.xml
        ${EXAMPLE_AUX_DIR}/templates/SPH.xml
        ${EXAMPLE_AUX_DIR}/templates/Square.xml
        ${EXAMPLE_AUX_DIR}/templates/Time.xml
        ${EXAMPLE_AUX_DIR}/templates/plot_f.py
    DESTINATION
        ${EXAMPLE_DEST_DIR}/templates/
)
INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/Create.py
        ${EXAMPLE_AUX_DIR}/Motion.py
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
        ${EXAMPLE_AUX_DIR}/doc/Force_Re100.dat
        ${EXAMPLE_AUX_DIR}/doc/Motion_Body.dat
    DESTINATION
        ${EXAMPLE_DEST_DIR}/doc
)