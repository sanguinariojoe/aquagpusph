include(${AQUAGPUSPH_CMAKE_MODULE_PATH}/example.cmake)

add_library(s_rocks_sim SHARED chronosim.cpp)
target_include_directories(s_rocks_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${CHRONO_INCLUDE_DIRS}
                                             ${CMAKE_BINARY_DIR}
                                             ${CMAKE_SOURCE_DIR}
                                             ${CMAKE_SOURCE_DIR}/aquagpusph/ext)
target_include_directories(s_rocks_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})
target_compile_options(s_rocks_sim PRIVATE ${CHRONO_CXX_FLAGS})
target_link_libraries(s_rocks_sim aquagpusphlib ${DEP_LIBS} ${CHRONO_LIBRARIES})
target_compile_definitions(s_rocks_sim PUBLIC AQUAgpusph_EXPORTS)
set_target_properties(s_rocks_sim PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(WIN32)
    set(EXTTOOL_LIB_NAME s_rocks_sim.dll)
else()
    set(EXTTOOL_LIB_NAME libs_rocks_sim.so)
endif()

if(WIN32)
    add_custom_target(copy_s_rocks_sim_dll_to_example ALL)
    add_custom_command(TARGET copy_s_rocks_sim_dll_to_example POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rocks_sim> ${CMAKE_CURRENT_BINARY_DIR}/${EXTTOOL_LIB_NAME}
    )
    add_dependencies(copy_s_rocks_sim_dll_to_example rocks_sim)
endif()

install(TARGETS s_rocks_sim
  LIBRARY DESTINATION ${EXAMPLE_INSTALL_DIR}
)
