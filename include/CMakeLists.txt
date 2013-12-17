# ===================================================== #
# Install                                               #
# ===================================================== #
INSTALL(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/ArgumentsManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/AuxiliarMethods.h
        ${CMAKE_CURRENT_SOURCE_DIR}/CalcServer.h
        ${CMAKE_CURRENT_SOURCE_DIR}/FileManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Fluid.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Input.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Output.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ProblemSetup.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ScreenManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Singleton.h
        ${CMAKE_CURRENT_SOURCE_DIR}/sphPrerequisites.h
        ${CMAKE_CURRENT_SOURCE_DIR}/TimeManager.h
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
        ${CMAKE_CURRENT_SOURCE_DIR}/Input
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
)
INSTALL(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/Tokenizer
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}
)
