include(${AQUAGPUSPH_CMAKE_MODULE_PATH}/example.cmake)

add_library(ball_inflow_sim SHARED chronosim.cpp)
target_include_directories(ball_inflow_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${CHRONO_INCLUDE_DIRS}
                                             ${CMAKE_BINARY_DIR}
                                             ${CMAKE_SOURCE_DIR}
                                             ${CMAKE_SOURCE_DIR}/aquagpusph/ext)
target_include_directories(ball_inflow_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})
target_compile_options(ball_inflow_sim PRIVATE ${CHRONO_CXX_FLAGS})
target_link_libraries(ball_inflow_sim aquagpusphlib ${DEP_LIBS} ${CHRONO_LIBRARIES})
target_compile_definitions(ball_inflow_sim PUBLIC AQUAgpusph_EXPORTS)
set_target_properties(ball_inflow_sim PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(WIN32)
    set(EXTTOOL_LIB_NAME ball_sim.dll)
else()
    set(EXTTOOL_LIB_NAME libs_ball_inflow_sim.so)
endif()

if(WIN32)
    add_custom_target(copy_ball_sim_dll_to_example ALL)
    add_custom_command(TARGET copy_ball_sim_dll_to_example POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rocks_sim> ${CMAKE_CURRENT_BINARY_DIR}/${EXTTOOL_LIB_NAME}
    )
    add_dependencies(copy_ball_sim_dll_to_example rocks_sim)
endif()

install(TARGETS ball_inflow_sim
  LIBRARY DESTINATION ${EXAMPLE_INSTALL_DIR}
)
