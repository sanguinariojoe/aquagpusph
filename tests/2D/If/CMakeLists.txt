# ===================================================== #
# In place configuration                                #
# ===================================================== #
SET(RESOURCES_DIR ${CMAKE_BINARY_DIR}/resources)
SET(TEST_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(TEST_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
FILE(GLOB_RECURSE FNAMES RELATIVE ${TEST_ORIG_DIR} "*")


FOREACH(FNAME ${FNAMES})
    GET_FILENAME_COMPONENT(FEXT ${FNAME} EXT)
    IF((${FEXT} MATCHES ".xml") OR (${FEXT} MATCHES ".py") OR (${FEXT} MATCHES ".sh"))
        configure_file(${TEST_ORIG_DIR}/${FNAME}
            ${TEST_DEST_DIR}/${FNAME} @ONLY)
    ELSE()
        configure_file(${TEST_ORIG_DIR}/${FNAME}
            ${TEST_DEST_DIR}/${FNAME} COPYONLY)
    ENDIF()
ENDFOREACH()

ADD_TEST (2D/If_test ${BASH_PROGRAM} ${TEST_DEST_DIR}/run.sh)
