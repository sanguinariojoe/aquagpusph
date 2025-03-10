set(PRESETS_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB_RECURSE FNAMES RELATIVE ${PRESETS_ORIG_DIR} "*.xml")
foreach(FNAME ${FNAMES})
    configure_file(${PRESETS_ORIG_DIR}/${FNAME}
        ${CMAKE_CURRENT_BINARY_DIR}/${FNAME} COPYONLY)
endforeach()

install(
    DIRECTORY
        ${PRESETS_ORIG_DIR}/
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources/Presets
    FILES_MATCHING
    PATTERN "*"
)
