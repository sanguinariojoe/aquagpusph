# ===================================================== #
# Embed OpenCL codes                                    #
# ===================================================== #
ADD_CUSTOM_TARGET(opencl_embed ALL
	COMMAND xxd -i LinkList.hcl.in > LinkList.hcl
	COMMAND xxd -i LinkList.cl.in > LinkList.cl
	COMMAND xxd -i RadixSort.hcl.in > RadixSort.hcl
	COMMAND xxd -i RadixSort.cl.in > RadixSort.cl
	COMMAND xxd -i Reduction.hcl.in > Reduction.hcl
	COMMAND xxd -i Reduction.cl.in > Reduction.cl
	COMMAND xxd -i Set.hcl.in > Set.hcl
	COMMAND xxd -i Set.cl.in > Set.cl
	COMMAND xxd -i UnSort.hcl.in > UnSort.hcl
	COMMAND xxd -i UnSort.cl.in > UnSort.cl
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer/
)

# ===================================================== #
# Install                                               #
# ===================================================== #
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/ArgumentsManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/AuxiliarMethods.h
        ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer.h
        ${CMAKE_CURRENT_SOURCE_DIR}/FileManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ProblemSetup.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ScreenManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Singleton.h
        ${CMAKE_CURRENT_SOURCE_DIR}/sphPrerequisites.h
        ${CMAKE_CURRENT_SOURCE_DIR}/TimeManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Variable.h
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
)

INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
)
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/InputOutput
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
)
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/Tokenizer
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
)
