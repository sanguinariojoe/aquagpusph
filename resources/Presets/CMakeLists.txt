SET(PRESETS_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
FILE(GLOB_RECURSE FNAMES RELATIVE ${PRESETS_ORIG_DIR} "*")

# ===================================================== #
# Configuration                                         #
# ===================================================== #
SET(AQUAGPUSPH_ROOT_PATH ${CMAKE_BINARY_DIR})
SET(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
SET(RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)
SET(PRESETS_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

FOREACH(FNAME ${FNAMES})
    CONFIGURE_FILE(${PRESETS_ORIG_DIR}/${FNAME}
        ${PRESETS_DEST_DIR}/${FNAME} @ONLY)
ENDFOREACH()

# ===================================================== #
# Installable version (and targets)                     #
# ===================================================== #
SET(AQUAGPUSPH_ROOT_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR})
SET(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
SET(RESOURCES_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
SET(PRESETS_AUX_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp)
SET(PRESETS_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources/Presets)

FOREACH(FNAME ${FNAMES})
    CONFIGURE_FILE(${PRESETS_ORIG_DIR}/${FNAME}
        ${PRESETS_AUX_DIR}/${FNAME} @ONLY)
ENDFOREACH()

INSTALL(
    DIRECTORY
        ${PRESETS_AUX_DIR}/
    DESTINATION
        ${PRESETS_DEST_DIR}
    FILES_MATCHING
    PATTERN "*"
)
