# ===================================================== #
# Sources to compile                                    #
# ===================================================== #
SET(CPP_SRCS
    ArgumentsManager.cpp
    AuxiliarMethods.cpp
    FileManager.cpp
    ProblemSetup.cpp
    TimeManager.cpp
    Variable.cpp
    InputOutput/State.cpp
    InputOutput/Report.cpp
    InputOutput/Logger.cpp
    InputOutput/Particles.cpp
    InputOutput/ASCII.cpp
    InputOutput/FastASCII.cpp
    InputOutput/CSV.cpp
    InputOutput/VTK.cpp
    CalcServer/Assert.cpp
    CalcServer/CalcServer.cpp
    CalcServer/Conditional.cpp
    CalcServer/Copy.cpp
    CalcServer/Kernel.cpp
    CalcServer/LinkList.cpp
    CalcServer/Python.cpp
    CalcServer/RadixSort.cpp
    CalcServer/Reduction.cpp
    CalcServer/Set.cpp
    CalcServer/SetScalar.cpp
    CalcServer/Sort.cpp
    CalcServer/Tool.cpp
    CalcServer/UnSort.cpp
    CalcServer/Reports/Dump.cpp
    CalcServer/Reports/Performance.cpp
    CalcServer/Reports/Report.cpp
    CalcServer/Reports/Screen.cpp
    CalcServer/Reports/SetTabFile.cpp
    CalcServer/Reports/TabFile.cpp
    ${OPTIONAL_SRCS}
)

# ===================================================== #
# Embed OpenCL codes                                    #
# ===================================================== #
ADD_CUSTOM_TARGET(opencl_embed_directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/)
SET(embed_targets opencl_embed_directory)
FOREACH(FNAME LinkList MPISync RadixSort Reduction Set Sort UnSort)
    FOREACH(FEXT .cl .hcl)
        ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo "/** @file" > ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo " * @brief Hardcoded version of the file CalcServer/${FNAME}${FEXT}.in" >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo " */" >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND echo "" >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND ${XXD_BIN} -i ${FNAME}${FEXT}.in >> ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer/
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer/${FNAME}${FEXT}.in)

        ADD_CUSTOM_TARGET(opencl_embed_${FNAME}${FEXT} ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
        )
        ADD_DEPENDENCIES(opencl_embed_${FNAME}${FEXT} opencl_embed_directory)
        LIST(APPEND embed_targets opencl_embed_${FNAME}${FEXT})
    ENDFOREACH()
ENDFOREACH()

ADD_CUSTOM_TARGET(opencl_embed ALL)
ADD_DEPENDENCIES(opencl_embed ${embed_targets})

# ===================================================== #
# library compilation                                   #
# ===================================================== #
SOURCE_GROUP("AQUAgpusph" FILES ${CPP_SRCS})

SET(LibName aquagpusph)
ADD_LIBRARY(${LibName} SHARED ${CPP_SRCS})
TARGET_LINK_LIBRARIES(${LibName} ${DEP_LIBS})
ADD_DEPENDENCIES(${LibName} opencl_embed)

if(MSVC)
    set_target_properties(${LibName} PROPERTIES DEBUG_OUTPUT_NAME "${LibName}D")
    set_target_properties(${LibName} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # dirty hack to avoid Debug/Release subdirectory
    set_target_properties(${LibName} PROPERTIES PREFIX "../")
elseif(MINGW)
    set_target_properties(${LibName} PROPERTIES DEBUG_OUTPUT_NAME "${LibName}D")
    set_target_properties(${LibName} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
else(MSVC)
    set_target_properties(${LibName} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    set_target_properties(${LibName} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH})
endif(MSVC)

if(WIN32)
    INSTALL(TARGETS ${LibName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else(WIN32)
    INSTALL(TARGETS ${LibName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif(WIN32)


# ===================================================== #
# program compilation                                   #
# ===================================================== #
SET(AppName AQUAgpusph)
ADD_EXECUTABLE(${AppName} main.cpp)
TARGET_LINK_LIBRARIES(${AppName} ${LibName} ${DEP_LIBS})

if(MSVC)
    set_target_properties(${AppName} PROPERTIES DEBUG_OUTPUT_NAME "${AppName}D")
    set_target_properties(${AppName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # dirty hack to avoid Debug/Release subdirectory
    set_target_properties(${AppName} PROPERTIES PREFIX "../")
elseif(MINGW)
    set_target_properties(${AppName} PROPERTIES DEBUG_OUTPUT_NAME "${AppName}D")
    set_target_properties(${AppName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
else(MSVC)
    set_target_properties(${AppName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    set_target_properties(${AppName} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH})
endif(MSVC)

if(WIN32)
    INSTALL(TARGETS ${AppName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else(WIN32)
    INSTALL(TARGETS ${AppName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif(WIN32)
