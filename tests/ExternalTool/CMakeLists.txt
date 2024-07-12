# ===================================================== #
# Compile the external tool                             #
# ===================================================== #
add_library(exttooltest SHARED tool.cpp)
target_include_directories(exttooltest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(exttooltest aquagpusph ${DEP_LIBS})
set_target_properties(exttooltest PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ===================================================== #
# In place configuration                                #
# ===================================================== #
set(RESOURCES_DIR ${CMAKE_BINARY_DIR}/resources)
set(TEST_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
set(TEST_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
file(GLOB_RECURSE FNAMES RELATIVE ${TEST_ORIG_DIR} "${TEST_ORIG_DIR}/*")

foreach(FNAME ${FNAMES})
    get_filename_component(FEXT ${FNAME} EXT)
    if((${FEXT} MATCHES ".xml") OR (${FEXT} MATCHES ".py") OR (${FEXT} MATCHES ".sh"))
        configure_file(${TEST_ORIG_DIR}/${FNAME}
            ${TEST_DEST_DIR}/${FNAME} @ONLY)
    else()
        configure_file(${TEST_ORIG_DIR}/${FNAME}
            ${TEST_DEST_DIR}/${FNAME} COPYONLY)
    endif()
endforeach()

# ===================================================== #
# The test                                              #
# ===================================================== #
add_test(ExternalTool ${BASH_PROGRAM} ${TEST_DEST_DIR}/run.sh)
set_tests_properties(ExternalTool PROPERTIES DEPENDS exttooltest)
