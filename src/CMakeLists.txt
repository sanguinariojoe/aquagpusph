# ===================================================== #
# Subcomponents                                         #
# ===================================================== #
add_subdirectory(CalcServer)

# ===================================================== #
# Dependencies                                          #
# ===================================================== #
IF(HAVE_NCURSES)
    SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${CURSES_INCLUDE_PATH})
    SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${CURSES_LIBRARY})
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
    ${PYTHON_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIRS}
    ${OPENCL_INCLUDE_DIRS}
    ${XERCESC_INCLUDE_PATH}
    ${OPTIONAL_INCLUDE_PATH}
)

SET(DEP_LIBS 
    ${OPENCL_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${XERCESC_LIBRARIES}
    ${OPTIONAL_LIBS}
    ${MATHEVAL_LIBRARIES}
)

# ===================================================== #
# Sources to compile client                             #
# ===================================================== #
SET(Client_CPP_SRCS
    ArgumentsManager.cpp
    AuxiliarMethods.cpp
    FileManager.cpp
    InputOutput/State.cpp
    InputOutput/Report.cpp
    InputOutput/Log.cpp
    InputOutput/Particles.cpp
    InputOutput/ASCII.cpp
    InputOutput/VTK.cpp
    ProblemSetup.cpp
    ScreenManager.cpp
    TimeManager.cpp
    Variable.cpp
    Tokenizer/Tokenizer.cpp
)

# ===================================================== #
# Client target                                         #
# ===================================================== #
SOURCE_GROUP("AQUAgpusphClient" FILES ${Client_CPP_SRCS})

IF(AQUAGPUSPH_3D)
    SET(ClientTagetName AQUAgpusphClient)
ELSE(AQUAGPUSPH_3D)
    SET(ClientTagetName AQUAgpusphClient2D)
ENDIF(AQUAGPUSPH_3D)

add_library(${ClientTagetName} SHARED ${Client_CPP_SRCS})

target_link_libraries(${ClientTagetName} ${DEP_LIBS})

# ===================================================== #
# Install client                                        #
# ===================================================== #
if(MSVC)
    set_target_properties(${ClientTagetName} PROPERTIES DEBUG_OUTPUT_NAME "${ClientTagetName}D")
    set_target_properties(${ClientTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # dirty hack to avoid Debug/Release subdirectory
    set_target_properties(${ClientTagetName} PROPERTIES PREFIX "../")
elseif(MINGW)
    set_target_properties(${ClientTagetName} PROPERTIES DEBUG_OUTPUT_NAME "${ClientTagetName}D")
    set_target_properties(${ClientTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
else(MSVC)
    set_target_properties(${ClientTagetName} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    set_target_properties(${ClientTagetName} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH})
endif(MSVC)

if(WIN32)
    INSTALL(TARGETS ${ClientTagetName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else(WIN32)
    INSTALL(TARGETS ${ClientTagetName}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif(WIN32)

# ===================================================== #
# Sources to compile app                                #
# ===================================================== #
SET(APP_CPP_SRCS
    main.cpp
)

# ===================================================== #
# App target                                         #
# ===================================================== #
SOURCE_GROUP("AQUAgpusph" FILES ${APP_CPP_SRCS})

IF(AQUAGPUSPH_3D)
    SET(BinTagetName AQUAgpusph)
    SET(ServerTagetName AQUAgpusphServer)
ELSE(AQUAGPUSPH_3D)
    SET(BinTagetName AQUAgpusph2D)
    SET(ServerTagetName AQUAgpusphServer2D)
ENDIF(AQUAGPUSPH_3D)

add_executable(${BinTagetName} ${APP_CPP_SRCS})

target_link_libraries(${BinTagetName} ${ClientTagetName} ${ServerTagetName} ${DEP_LIBS})

# ===================================================== #
# Install App                                           #
# ===================================================== #
if(MSVC)
    set_target_properties(${BinTagetName} PROPERTIES DEBUG_OUTPUT_NAME "${BinTagetName}D")
    set_target_properties(${BinTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # dirty hack to avoid Debug/Release subdirectory
    set_target_properties(${BinTagetName} PROPERTIES PREFIX "../")
elseif(MINGW)
    set_target_properties(${BinTagetName} PROPERTIES DEBUG_OUTPUT_NAME "${BinTagetName}D")
    set_target_properties(${BinTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
else(MSVC)
    set_target_properties(${BinTagetName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    set_target_properties(${BinTagetName} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH})
endif(MSVC)

if(WIN32)
    INSTALL(TARGETS ${BinTagetName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else(WIN32)
    INSTALL(TARGETS ${BinTagetName}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif(WIN32)
