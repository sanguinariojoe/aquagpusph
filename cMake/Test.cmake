# Configure a test

function (aquagpusph_conf_test NAME TEST_ORIG_DIR TEST_DEST_DIR DIMS)
    set(RESOURCES_DIR ${CMAKE_BINARY_DIR}/resources)
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

    add_test(${NAME} ${BASH_PROGRAM} ${TEST_DEST_DIR}/run.sh)
endfunction ()
