include(${CMAKE_SOURCE_DIR}/cMake/Test.cmake)

# Prepare a copy of AQUAgpusph that we have perfectly controlled
set(BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
if(WIN32)
    set(BINARY_NAME AQUAgpusph.exe)
else()
    set(BINARY_NAME AQUAgpusph)
endif()
add_custom_target(COPY_BIN_TO_TESTS ALL)
add_custom_command(TARGET COPY_BIN_TO_TESTS POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:AQUAgpusph> ${CMAKE_CURRENT_BINARY_DIR}/
)
if(WIN32)
    # We are assuming that the user already placed the dependency libs on
    # $PATH, so we just need the executable itself
    add_custom_command(TARGET COPY_BIN_TO_TESTS POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:aquagpusphlib> ${CMAKE_CURRENT_BINARY_DIR}/
    )
endif()

add_subdirectory(ExternalTool)
add_subdirectory(If)
add_subdirectory(While)
add_subdirectory(QueueLimit)
add_subdirectory(Python)
add_subdirectory(CSV)

add_subdirectory(3D/SetScalar)
add_subdirectory(3D/Assert)
add_subdirectory(3D/Copy)
add_subdirectory(3D/Kernel)
add_subdirectory(3D/Set)
add_subdirectory(3D/Reduction)
add_subdirectory(3D/RadixSort)
add_subdirectory(3D/Sort)
add_subdirectory(3D/LinkList)
IF(AQUAGPUSPH_USE_MPI)
    add_subdirectory(3D/MPI_sync)
ENDIF(AQUAGPUSPH_USE_MPI)
add_subdirectory(3D/MLS)
IF(AQUAGPUSPH_USE_VTK)
    add_subdirectory(3D/SaveVTK)
ENDIF(AQUAGPUSPH_USE_VTK)

add_subdirectory(2D/SetScalar)
add_subdirectory(2D/Assert)
add_subdirectory(2D/Copy)
add_subdirectory(2D/Kernel)
add_subdirectory(2D/Set)
add_subdirectory(2D/Reduction)
add_subdirectory(2D/RadixSort)
add_subdirectory(2D/Sort)
add_subdirectory(2D/LinkList)
IF(AQUAGPUSPH_USE_MPI)
    add_subdirectory(2D/MPI_sync)
    add_subdirectory(2D/MPI_plane)
ENDIF(AQUAGPUSPH_USE_MPI)
add_subdirectory(2D/MLS)
IF(AQUAGPUSPH_USE_VTK)
    add_subdirectory(2D/SaveVTK)
ENDIF(AQUAGPUSPH_USE_VTK)

add_subdirectory(2D/Regressions/midpoint_predictor)
