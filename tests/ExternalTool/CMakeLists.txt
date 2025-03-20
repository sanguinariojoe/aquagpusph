# ===================================================== #
# Compile the external tool                             #
# ===================================================== #
add_library(exttooltest SHARED tool.cpp)
target_include_directories(exttooltest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                               ${CMAKE_BINARY_DIR}
                                               ${CMAKE_SOURCE_DIR})
target_link_libraries(exttooltest aquagpusphlib ${DEP_LIBS})
target_compile_definitions(exttooltest PUBLIC AQUAgpusph_EXPORTS)
set_target_properties(exttooltest PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(WIN32)
    set(EXTTOOL_LIB_NAME exttooltest.dll)
else()
    set(EXTTOOL_LIB_NAME libexttooltest.so)
endif()

# ===================================================== #
# The test                                              #
# ===================================================== #
aquagpusph_conf_test(ExternalTool ${CMAKE_CURRENT_SOURCE_DIR}/cMake ${CMAKE_CURRENT_BINARY_DIR} 2)
set_tests_properties(ExternalTool PROPERTIES DEPENDS exttooltest)
if(WIN32)
    add_custom_target(copy_exttooltest_dll_to_test ALL)
    add_custom_command(TARGET copy_exttooltest_dll_to_test POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:exttooltest> ${CMAKE_CURRENT_BINARY_DIR}/${EXTTOOL_LIB_NAME}
    )
    add_dependencies(copy_exttooltest_dll_to_test exttooltest)
    set_tests_properties(ExternalTool PROPERTIES DEPENDS copy_exttooltest_dll_to_test)
endif()
