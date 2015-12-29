# ===================================================== #
# Install                                               #
# ===================================================== #
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/Scripts
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/resources
)

ADD_SUBDIRECTORY(Presets)
