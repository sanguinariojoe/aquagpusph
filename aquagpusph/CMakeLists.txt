# ===================================================== #
# Sources to compile                                    #
# ===================================================== #
set(CPP_SRCS
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
add_custom_target(opencl_embed_directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/)
set(embed_targets "")
foreach(FNAME LinkList MPISync RadixSort Reduction Set Sort UnSort)
    foreach(FEXT .cl .hcl)
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            COMMAND ${XXD_BIN} -i ${FNAME}${FEXT}.in > ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer/
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer/${FNAME}${FEXT}.in)

        add_custom_target(opencl_embed_${FNAME}${FEXT} ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CalcServer/${FNAME}${FEXT}
        )
        add_dependencies(opencl_embed_${FNAME}${FEXT} opencl_embed_directory)
        list(APPEND embed_targets opencl_embed_${FNAME}${FEXT})
    endforeach()
endforeach()

# ===================================================== #
# library compilation                                   #
# ===================================================== #
source_group("AQUAgpusph" FILES ${CPP_SRCS})

# NOTE: The library cannot be called aquagpusph, because VSC is apparently
# case insensitive, so it will overwrite it with the executable target
set(LibName aquagpusphlib)
add_library(${LibName} SHARED ${CPP_SRCS})
target_link_libraries(${LibName} ${DEP_LIBS})
target_compile_definitions(${LibName} PUBLIC AQUAgpusph_EXPORTS)
add_dependencies(${LibName} ${embed_targets})

install(TARGETS ${LibName}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


# ===================================================== #
# program compilation                                   #
# ===================================================== #
set(AppName AQUAgpusph)

set(EXE_SRCS main.cpp)
if(WIN32)
    list(APPEND EXE_SRCS main.rc)
endif()
add_executable(${AppName} ${EXE_SRCS})
target_link_libraries(${AppName} ${LibName} ${DEP_LIBS})

install(TARGETS ${AppName}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
  FILES
    ArgumentsManager.hpp
    AuxiliarMethods.hpp
    FileManager.hpp
    ProblemSetup.hpp
    sphPrerequisites.hpp
    TimeManager.hpp
    Variable.hpp
    ${CMAKE_BINARY_DIR}/config.h
  DESTINATION
    include/aquagpusph
)

install(
  FILES
    InputOutput/Logger.hpp
    InputOutput/State.hpp
  DESTINATION
    include/aquagpusph/InputOutput
)

install(
  FILES
    Tokenizer/Tokenizer.hpp
    Tokenizer/Tokenizer_exprtk.hpp
    Tokenizer/Tokenizer_muparser.hpp
  DESTINATION
    include/aquagpusph/Tokenizer
)

install(
  FILES
    CalcServer/CalcServer.hpp
    CalcServer/Tool.hpp
  DESTINATION
    include/aquagpusph/CalcServer
)
