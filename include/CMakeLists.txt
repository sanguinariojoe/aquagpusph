# ===================================================== #
# Embed OpenCL codes                                    #
# ===================================================== #
ADD_CUSTOM_TARGET(opencl_embed_directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/)
SET(embed_targets opencl_embed_directory)
FOREACH(FNAME LinkList RadixSort Reduction Set UnSort)
    FOREACH(FEXT .cl .hcl)
        ADD_CUSTOM_TARGET(opencl_embed_${FNAME}${FEXT} ALL
            COMMAND echo "/** @file" > ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo " * @brief Hardcoded version of the file CalcServer/${FNAME}${FEXT}.in" >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo " */" >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo "" >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND ${XXD_BIN} -i ${FNAME}${FEXT}.in >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer/
        )
        ADD_DEPENDENCIES(opencl_embed_${FNAME}${FEXT} opencl_embed_directory)
        LIST(APPEND embed_targets opencl_embed_${FNAME}${FEXT})
    ENDFOREACH()
ENDFOREACH()

ADD_CUSTOM_TARGET(opencl_embed ALL)
ADD_DEPENDENCIES(opencl_embed ${embed_targets})

# ===================================================== #
# Install                                               #
# ===================================================== #
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
)

INSTALL(
    DIRECTORY
        ${CMAKE_BINARY_DIR}/include/
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.cl"
    PATTERN "*.hcl"
    PATTERN "CMakeFiles" EXCLUDE
)
