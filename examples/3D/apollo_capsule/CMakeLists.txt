SET(EXAMPLE_ORIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cMake)
FILE(GLOB_RECURSE FNAMES RELATIVE ${EXAMPLE_ORIG_DIR} "${EXAMPLE_ORIG_DIR}/*")

# ===================================================== #
# Compile the chrono simulation tool                    #
# ===================================================== #
add_library(apollo_capsule_sim SHARED chronosim.cpp)
target_include_directories(apollo_capsule_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CHRONO_INCLUDE_DIRS})
target_compile_options(apollo_capsule_sim PRIVATE ${CHRONO_CXX_FLAGS})
target_link_libraries(apollo_capsule_sim aquagpusph ${DEP_LIBS} ${CHRONO_LIBRARIES})
set_target_properties(apollo_capsule_sim PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ===================================================== #
# Configuration                                         #
# ===================================================== #
SET(RESOURCES_DIR ${CMAKE_BINARY_DIR}/resources)
SET(EXAMPLE_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(BINARY_DIR ${CMAKE_BINARY_DIR}/bin)

FOREACH(FNAME ${FNAMES})
    GET_FILENAME_COMPONENT(FEXT ${FNAME} EXT)
    IF((${FEXT} MATCHES ".xml") OR (${FEXT} MATCHES ".py") OR (${FEXT} MATCHES ".sh"))
        configure_file(${EXAMPLE_ORIG_DIR}/${FNAME}
            ${EXAMPLE_DEST_DIR}/${FNAME} @ONLY)
    ELSE()
        configure_file(${EXAMPLE_ORIG_DIR}/${FNAME}
            ${EXAMPLE_DEST_DIR}/${FNAME} COPYONLY)
    ENDIF()
ENDFOREACH()

# ===================================================== #
# Installable version (and targets)                     #
# ===================================================== #
FILE(RELATIVE_PATH RUN_REL_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/cMake/run.sh)
GET_FILENAME_COMPONENT(EXAMPLE_REL_PATH ${RUN_REL_PATH} PATH)
SET(RESOURCES_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources)
SET(EXAMPLE_AUX_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp)
SET(EXAMPLE_DEST_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/${EXAMPLE_REL_PATH}/..)
SET(BINARY_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

FOREACH(FNAME ${FNAMES})
    GET_FILENAME_COMPONENT(FEXT ${FNAME} EXT)
    IF((${FEXT} MATCHES ".xml") OR (${FEXT} MATCHES ".py") OR (${FEXT} MATCHES ".sh"))
        configure_file(${EXAMPLE_ORIG_DIR}/${FNAME}
            ${EXAMPLE_AUX_DIR}/${FNAME} @ONLY)
    ELSE()
        configure_file(${EXAMPLE_ORIG_DIR}/${FNAME}
            ${EXAMPLE_AUX_DIR}/${FNAME} COPYONLY)
    ENDIF()
ENDFOREACH()

INSTALL(
    DIRECTORY
        ${EXAMPLE_AUX_DIR}/
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    FILES_MATCHING
    PATTERN "*"
    PATTERN "*.sh" EXCLUDE
)

INSTALL(
    FILES
        ${EXAMPLE_AUX_DIR}/run.sh
    DESTINATION
        ${EXAMPLE_DEST_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE             GROUP_READ
                WORLD_EXECUTE             WORLD_READ             
)
