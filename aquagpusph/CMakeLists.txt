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
)
set(PUBLIC_LIBS
    Python::Python
    Python::NumPy
    XercesC::XercesC
    OpenCL::OpenCL
)
set(PRIVATE_LIBS
)
set(PRIVATE_INCS
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/aquagpusph/ext
)
set(PUBLIC_INCS
)

if(HAVE_MPI)
    list(APPEND PUBLIC_LIBS MPI::MPI_C)
    list(APPEND CPP_SRCS CalcServer/MPISync.cpp)
endif(HAVE_MPI)
if(HAVE_MUPARSER)
    list(APPEND PUBLIC_INCS ${MUPARSER_INCLUDE_DIRS})
    list(APPEND PUBLIC_LIBS ${MUPARSER_LIBRARIES})
    list(APPEND CPP_SRCS Tokenizer/Tokenizer_muparser.cpp)
else(HAVE_MUPARSER)
    list(APPEND CPP_SRCS Tokenizer/Tokenizer_exprtk.cpp)
endif(HAVE_MUPARSER)
if(HAVE_VTK)
    if(AQUAGPUSPH_PACKAGE_IGNORE_VTK_DEPENDENCY)
        list(APPEND PRIVATE_LIBS VTK::CommonCore VTK::IOXML)
    else(AQUAGPUSPH_PACKAGE_IGNORE_VTK_DEPENDENCY)
        list(APPEND PUBLIC_LIBS VTK::CommonCore VTK::IOXML)
    endif(AQUAGPUSPH_PACKAGE_IGNORE_VTK_DEPENDENCY)
endif(HAVE_VTK)
if(nlohmann_json_FOUND)
    list(APPEND PUBLIC_LIBS nlohmann_json::nlohmann_json)
endif(nlohmann_json_FOUND)

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
add_library(aquagpusphlib SHARED ${CPP_SRCS})
target_compile_definitions(aquagpusphlib PUBLIC AQUAgpusph_EXPORTS)
target_link_libraries(aquagpusphlib PUBLIC ${PUBLIC_LIBS}
                                    PRIVATE ${PRIVATE_LIBS})
target_include_directories(aquagpusphlib PUBLIC ${PUBLIC_INCS}
                                         PRIVATE ${PRIVATE_INCS})
add_dependencies(aquagpusphlib ${embed_targets})

install(TARGETS aquagpusphlib
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION include
)

# ===================================================== #
# program compilation                                   #
# ===================================================== #
set(AppName AQUAgpusph)

set(EXE_SRCS main.cpp)
if(WIN32)
    list(APPEND EXE_SRCS main.rc)
endif()
add_executable(AQUAgpusph ${EXE_SRCS})
target_include_directories(AQUAgpusph PUBLIC ${PUBLIC_INCS}
                                      PRIVATE ${PRIVATE_INCS})
target_link_libraries(AQUAgpusph aquagpusphlib)

install(TARGETS AQUAgpusph
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
    ${CMAKE_BINARY_DIR}/aquagpusph/config.h
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
