SET(SCRIPTS_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.)
FILE(GLOB_RECURSE FNAMES RELATIVE ${SCRIPTS_ORIG_DIR} "*")

# ===================================================== #
# Configuration                                         #
# ===================================================== #
SET(SCRIPTS_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

FOREACH(FNAME ${FNAMES})
    CONFIGURE_FILE(${SCRIPTS_ORIG_DIR}/${FNAME}
        ${SCRIPTS_DEST_DIR}/${FNAME} COPYONLY)
ENDFOREACH()

# ===================================================== #
# Installable version (and targets)                     #
# ===================================================== #
SET(SCRIPTS_AUX_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp)
SET(SCRIPTS_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources/Scripts)

FOREACH(FNAME ${FNAMES})
    CONFIGURE_FILE(${SCRIPTS_ORIG_DIR}/${FNAME}
        ${SCRIPTS_AUX_DIR}/${FNAME} COPYONLY)
ENDFOREACH()

INSTALL(
    DIRECTORY
        ${SCRIPTS_AUX_DIR}/
    DESTINATION
        ${SCRIPTS_DEST_DIR}
    FILES_MATCHING
    PATTERN "*"
    PATTERN "CMakeLists.txt" EXCLUDE
)
