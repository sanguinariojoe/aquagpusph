ADD_SUBDIRECTORY(3D)
ADD_SUBDIRECTORY(2D)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/aqua_example_utils.py
    ${CMAKE_CURRENT_BINARY_DIR}/aqua_example_utils.py COPYONLY)

install(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/aqua_example_utils.py
    DESTINATION
        ${CMAKE_INSTALL_DATADIR}/examples/
)
