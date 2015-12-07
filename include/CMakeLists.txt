# ===================================================== #
# Embed OpenCL codes                                    #
# ===================================================== #
ADD_CUSTOM_TARGET(opencl_embed ALL
    COMMAND echo "/** @file" > LinkList.hcl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/LinkList.hcl.in" >> LinkList.hcl
    COMMAND echo " */" >> LinkList.hcl
    COMMAND echo "" >> LinkList.hcl
    COMMAND xxd -i LinkList.hcl.in >> LinkList.hcl
    COMMAND echo "/** @file" > LinkList.cl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/LinkList.cl.in" >> LinkList.cl
    COMMAND echo " */" >> LinkList.cl
    COMMAND echo "" >> LinkList.cl
    COMMAND xxd -i LinkList.cl.in >> LinkList.cl
    COMMAND echo "/** @file" > RadixSort.hcl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/RadixSort.hcl.in" >> RadixSort.hcl
    COMMAND echo " */" >> RadixSort.hcl
    COMMAND echo "" >> RadixSort.hcl
    COMMAND xxd -i RadixSort.hcl.in >> RadixSort.hcl
    COMMAND echo "/** @file" > RadixSort.cl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/RadixSort.cl.in" >> RadixSort.cl
    COMMAND echo " */" >> RadixSort.cl
    COMMAND echo "" >> RadixSort.cl
    COMMAND xxd -i RadixSort.cl.in >> RadixSort.cl
    COMMAND echo "/** @file" > Reduction.hcl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/Reduction.hcl.in" >> Reduction.hcl
    COMMAND echo " */" >> Reduction.hcl
    COMMAND echo "" >> Reduction.hcl
    COMMAND xxd -i Reduction.hcl.in >> Reduction.hcl
    COMMAND echo "/** @file" > Reduction.cl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/Reduction.cl.in" >> Reduction.cl
    COMMAND echo " */" >> Reduction.cl
    COMMAND echo "" >> Reduction.cl
    COMMAND xxd -i Reduction.cl.in >> Reduction.cl
    COMMAND echo "/** @file" > Set.hcl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/Set.hcl.in" >> Set.hcl
    COMMAND echo " */" >> Set.hcl
    COMMAND echo "" >> Set.hcl
    COMMAND xxd -i Set.hcl.in >> Set.hcl
    COMMAND echo "/** @file" > Set.cl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/Set.cl.in" >> Set.cl
    COMMAND echo " */" >> Set.cl
    COMMAND echo "" >> Set.cl
    COMMAND xxd -i Set.cl.in >> Set.cl
    COMMAND echo "/** @file" > UnSort.hcl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/UnSort.hcl.in" >> UnSort.hcl
    COMMAND echo " */" >> UnSort.hcl
    COMMAND echo "" >> UnSort.hcl
    COMMAND xxd -i UnSort.hcl.in >> UnSort.hcl
    COMMAND echo "/** @file" > UnSort.cl
    COMMAND echo " * @brief Hardcoded version of the file CalcServer/UnSort.cl.in" >> UnSort.cl
    COMMAND echo " */" >> UnSort.cl
    COMMAND echo "" >> UnSort.cl
    COMMAND xxd -i UnSort.cl.in >> UnSort.cl
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
