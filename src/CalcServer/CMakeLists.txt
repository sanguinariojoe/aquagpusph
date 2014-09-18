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
	${PYTHON_INCLUDE_PATH}
	${OPENCL_INCLUDE_DIRS}
	${XERCESC_INCLUDE_PATH}
	${CLANG_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
	${OPTIONAL_INCLUDE_PATH}
)

SET(DEP_LIBS 
	${OPENCL_LIBRARIES}
	${PYTHON_LIBRARIES}
	${XERCESC_LIBRARIES}
	${CLANG_LIBRARIES}
	${OPTIONAL_LIBS}
)

# ===================================================== #
# Sources to compile                                    #
# ===================================================== #
SET(Server_CPP_SRCS
	CalcServer.cpp
	Copy.cpp
	Kernel.cpp
	LinkList.cpp
	RadixSort.cpp
	Reduction.cpp
	Set.cpp
	Tool.cpp
	UnSort.cpp
)

# ===================================================== #
# Target                                                #
# ===================================================== #
SOURCE_GROUP("AQUAgpusphServer" FILES ${Server_CPP_SRCS})

IF(AQUAGPUSPH_3D)
	SET(ServerTagetName AQUAgpusphServer)
	SET(OpenCLEmbedTargetName opencl_embed)
ELSE(AQUAGPUSPH_3D)
	SET(ServerTagetName AQUAgpusphServer2D)
	SET(OpenCLEmbedTargetName opencl_embed)
ENDIF(AQUAGPUSPH_3D)

add_library(${ServerTagetName} SHARED ${Server_CPP_SRCS})

target_link_libraries(${ServerTagetName} ${DEP_LIBS})
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
