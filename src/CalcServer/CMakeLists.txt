# ===================================================== #
# Dependencies                                          #
# ===================================================== #
IF(HAVE_MPI)
    SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${MPI_INCLUDE_PATH})
    SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${MPI_CXX_LIBRARIES})
    SET(OPTIONAL_SRCS ${OPTIONAL_SRCS} MPISync.cpp)
ENDIF(HAVE_MPI)
IF(HAVE_NCURSES)
    SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${CURSES_INCLUDE_DIRS})
    SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${CURSES_LIBRARIES})
ENDIF(HAVE_NCURSES)
IF(HAVE_VTK)
    SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${VTK_INCLUDE_DIRS})
    SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${VTK_LIBRARIES})
ENDIF(HAVE_VTK)

# ===================================================== #
# Include & Link                                        #
# ===================================================== #
INCLUDE_DIRECTORIES(
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Python_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIRS}
    ${OpenCL_INCLUDE_DIRS}
    ${XercesC_INCLUDE_DIRS}
    ${CLANG_INCLUDE_DIRS}
    ${OPTIONAL_INCLUDE_PATH}
)

SET(DEP_LIBS 
    ${OpenCL_LIBRARIES}
    ${Python_LIBRARIES}
    ${XercesC_LIBRARIES}
    ${CLANG_LIBRARIES}
    ${OPTIONAL_LIBS}
)

# ===================================================== #
# Sources to compile                                    #
# ===================================================== #
SET(Server_CPP_SRCS
    Assert.cpp
    CalcServer.cpp
    Conditional.cpp
    Copy.cpp
    Kernel.cpp
    LinkList.cpp
    Python.cpp
    RadixSort.cpp
    Reduction.cpp
    Set.cpp
    SetScalar.cpp
    Tool.cpp
    UnSort.cpp
    Reports/Performance.cpp
    Reports/Report.cpp
    Reports/Screen.cpp
    Reports/SetTabFile.cpp
    Reports/TabFile.cpp
    ${OPTIONAL_SRCS}
)

# ===================================================== #
# Target                                                #
# ===================================================== #
SOURCE_GROUP("AQUAgpusphServer" FILES ${Server_CPP_SRCS})

IF(AQUAGPUSPH_3D)
    SET(ServerTagetName AQUAgpusphServer)
ELSE(AQUAGPUSPH_3D)
    SET(ServerTagetName AQUAgpusphServer2D)
ENDIF(AQUAGPUSPH_3D)
SET(OpenCLEmbedTargetName opencl_embed)

ADD_LIBRARY(${ServerTagetName} SHARED ${Server_CPP_SRCS})
TARGET_LINK_LIBRARIES(${ServerTagetName} ${DEP_LIBS})
ADD_DEPENDENCIES(${ServerTagetName} ${OpenCLEmbedTargetName})

# ===================================================== #
# Install                                               #
# ===================================================== #
if(MSVC)
    set_target_properties(${ServerTagetName} PROPERTIES DEBUG_OUTPUT_NAME "${ServerTagetName}D")
    set_target_properties(${ServerTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # dirty hack to avoid Debug/Release subdirectory
    set_target_properties(${ServerTagetName} PROPERTIES PREFIX "../")
elseif(MINGW)
    set_target_properties(${ServerTagetName} PROPERTIES DEBUG_OUTPUT_NAME "${ServerTagetName}D")
    set_target_properties(${ServerTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
else(MSVC)
    set_target_properties(${ServerTagetName} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    set_target_properties(${ServerTagetName} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH})
endif(MSVC)

if(WIN32)
    INSTALL(TARGETS ${ServerTagetName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else(WIN32)
    INSTALL(TARGETS ${ServerTagetName}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif(WIN32)
