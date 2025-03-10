include(${CMAKE_MODULE_PATH}/example.cmake)

add_library(apollo_capsule_sim SHARED chronosim.cpp)
target_include_directories(apollo_capsule_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})
target_compile_options(apollo_capsule_sim PRIVATE ${CHRONO_CXX_FLAGS})
target_link_libraries(apollo_capsule_sim aquagpusphlib ${DEP_LIBS} ${CHRONO_LIBRARIES})
target_compile_definitions(apollo_capsule_sim PUBLIC AQUAgpusph_EXPORTS)
set_target_properties(apollo_capsule_sim PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(WIN32)
    set(EXTTOOL_LIB_NAME apollo_capsule_sim.dll)
else()
    set(EXTTOOL_LIB_NAME libapollo_capsule_sim.so)
endif()

if(WIN32)
    add_custom_target(copy_apollo_capsule_sim_dll_to_example ALL)
    add_custom_command(TARGET copy_apollo_capsule_sim_dll_to_example POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:apollo_capsule_sim> ${CMAKE_CURRENT_BINARY_DIR}/${EXTTOOL_LIB_NAME}
    )
    add_dependencies(copy_apollo_capsule_sim_dll_to_example apollo_capsule_sim)
endif()

install(TARGETS apollo_capsule_sim
  LIBRARY DESTINATION ${EXAMPLE_INSTALL_DIR}
)
