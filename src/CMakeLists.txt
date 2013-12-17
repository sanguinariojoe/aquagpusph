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
IF(HAVE_H5PART)
	SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${HDF5_INCLUDE_PATH} ${H5PART_INCLUDE_PATH})
	SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${HDF5_LIBRARY} ${H5PART_LIBRARY})
ENDIF(HAVE_H5PART)
IF(HAVE_VTK)
	SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${VTK_INCLUDE_DIRS})
	SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${VTK_LIBRARIES})
ENDIF(HAVE_VTK)
IF(HAVE_MPI)
	SET(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${MPI_INCLUDE_PATH})
	SET(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${MPI_LIBRARIES})
ENDIF(HAVE_MPI)

# ===================================================== #
# Include & Link                                        #
# ===================================================== #
INCLUDE_DIRECTORIES(
	${CMAKE_BINARY_DIR}/include
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_CURRENT_BINARY_DIR}
	${PYTHON_INCLUDE_PATH}
	${OPENCL_INCLUDE_DIRS}
	${XERCESC_INCLUDE_PATH}
	${OPTIONAL_INCLUDE_PATH}
)

SET(DEP_LIBS 
	${OPENCL_LIBRARIES}
	${PYTHON_LIBRARIES}
	${XERCESC_LIBRARIES}
	${OPTIONAL_LIBS}
)

# ===================================================== #
# Sources to compile client                             #
# ===================================================== #
SET(Client_CPP_SRCS
	ArgumentsManager.cpp
	AuxiliarMethods.cpp
	FileManager.cpp
	Fluid.cpp
	Input.cpp
	Output.cpp
	ProblemSetup.cpp
	ScreenManager.cpp
	TimeManager.cpp
	Input/ASCII.cpp
	Input/GiD.cpp
	Input/XML.cpp
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
