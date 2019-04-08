# ===================================================== #
# In place configuration                                #
# ===================================================== #
SET(RESOURCES_DIR ${CMAKE_SOURCE_DIR}/resources)
SET(TEST_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
SET(TEST_DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)

CONFIGURE_FILE(${TEST_ORIG_DIR}/main.xml
    ${TEST_DEST_DIR}/main.xml @ONLY)

ADD_TEST (Sort_test ${BASH_PROGRAM} ${TEST_DEST_DIR}/run.sh)
