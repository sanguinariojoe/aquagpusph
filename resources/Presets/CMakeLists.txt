set(PRESETS_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
file(GLOB_RECURSE FNAMES RELATIVE ${PRESETS_ORIG_DIR} "*")

# ===================================================== #
# Configuration                                         #
# ===================================================== #
set(AQUAGPUSPH_ROOT_PATH ${CMAKE_BINARY_DIR})
set(RESOURCES_OUTPUT_DIR ${CMAKE_BINARY_DIR}/resources)
set(RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)
set(PRESETS_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

foreach(FNAME ${FNAMES})
    configure_file(${PRESETS_ORIG_DIR}/${FNAME}
        ${PRESETS_DEST_DIR}/${FNAME} @ONLY)
endforeach()

# ===================================================== #
# Installable version (and targets)                     #
# ===================================================== #
set(AQUAGPUSPH_ROOT_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR})
set(RESOURCES_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
set(RESOURCES_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
set(PRESETS_AUX_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp)
set(PRESETS_DEST_DIR ${CMAKE_INSTALL_DATADIR}/resources/Presets)

foreach(FNAME ${FNAMES})
    configure_file(${PRESETS_ORIG_DIR}/${FNAME}
        ${PRESETS_AUX_DIR}/${FNAME} @ONLY)
endforeach()

install(
    DIRECTORY
        ${PRESETS_AUX_DIR}/
    DESTINATION
        ${PRESETS_DEST_DIR}
    FILES_MATCHING
    PATTERN "*"
)
