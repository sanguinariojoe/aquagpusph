# ===================================================== #
# Local non-installed version                           #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)

configure_file(${EXAMPLE_ORIG_DIR}/RootPath.xml
    ${EXAMPLE_DEST_DIR}/RootPath.xml @ONLY)

# ===================================================== #
# Installable version (and tagets)                      #
# ===================================================== #
SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(EXAMPLE_AUX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CMakeTmp)
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples/)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file(${EXAMPLE_ORIG_DIR}/RootPath.xml
    ${EXAMPLE_AUX_DIR}/RootPath.xml @ONLY)

INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/RootPath.xml
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)

# ===================================================== #
# Setup & install examples                              #
# ===================================================== #
IF(AQUAGPUSPH_3D)
    ADD_SUBDIRECTORY(3D)
ELSE(AQUAGPUSPH_3D)
    ADD_SUBDIRECTORY(2D)
ENDIF(AQUAGPUSPH_3D)
