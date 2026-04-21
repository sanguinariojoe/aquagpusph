include(${AQUAGPUSPH_CMAKE_MODULE_PATH}/example.cmake)

add_library(britan_ball_sim_com_2 SHARED chronosim.cpp)
target_include_directories(britan_ball_sim_com_2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${CHRONO_INCLUDE_DIRS}
                                             ${CMAKE_BINARY_DIR}
                                             ${CMAKE_SOURCE_DIR}
                                             ${CMAKE_SOURCE_DIR}/aquagpusph/ext)
target_include_directories(britan_ball_sim_com_2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})
target_compile_options(britan_ball_sim_com_2 PRIVATE ${CHRONO_CXX_FLAGS})
target_link_libraries(britan_ball_sim_com_2 aquagpusphlib ${DEP_LIBS} ${CHRONO_LIBRARIES})
target_compile_definitions(britan_ball_sim_com_2 PUBLIC AQUAgpusph_EXPORTS)
set_target_properties(britan_ball_sim_com_2 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(WIN32)
    set(EXTTOOL_LIB_NAME britan_ball_sim_com_2.dll)
else()
    set(EXTTOOL_LIB_NAME libs_britan_ball_sim_com_2.so)
endif()

if(WIN32)
    add_custom_target(copy_britan_ball_sim_com_2_dll_to_example ALL)
    add_custom_command(TARGET copy_britan_ball_sim_com_2_dll_to_example POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rocks_sim> ${CMAKE_CURRENT_BINARY_DIR}/${EXTTOOL_LIB_NAME}
    )
    add_dependencies(copy_britan_ball_sim__com_2dll_to_example rocks_sim)
endif()

install(TARGETS britan_ball_sim_com_2
  LIBRARY DESTINATION ${EXAMPLE_INSTALL_DIR}
)
