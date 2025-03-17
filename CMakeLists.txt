# ===================================================== #
#                                                       #
#    #    ##   #  #   #                           #     #
#   # #  #  #  #  #  # #                          #     #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #
#                            # #             #          #
#                          ##  #             #          #
#                                                       #
# ===================================================== #
#                                                       #
#   Authors:                                            #
#       - Jose Luis Cercos-Pita                         #
#       - Leo Miguel Gonzalez                           #
#       - Antonio Souto-Iglesias                        #
#                                                       #
# ===================================================== #


# ===================================================== #
# General CMake options                                 #
# ===================================================== #
cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message("ERROR: In-source builds are not allowed")
    message("=======================================")
    message("")
    message("cmake cannot be executed in the same folder were AQUAgpusph has")
    message("been downloaded")
    message("See https://gitlab.com/sanguinariojoe/aquagpusph/-/wikis/install#compiling")
    message("")
    message("To fix this issue, you can remove CMakeCache.txt, create a")
    message("subfolder, and execute CMake from there.")
    message("")
    message("For instance, if you originally executed the following command:")
    message("")
    message("cmake -DCMAKE_BUILD_TYPE=Release .")
    message("")
    message("You can fix the problem typing")
    message("rm CMakeCache.txt")
    message("mkdir -p build && cd build")
    message("cmake -DCMAKE_BUILD_TYPE=Release ..")
    message("")
    message("Note that the single period at the end of the CMake command has")
    message("been replace by two periods (pointing to the parent folder)")
    message("")
    message(FATAL_ERROR "In-source builds are not allowed. Please, run cmake in a different folder.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
project(AQUAgpusph)

set(PACKAGE_NAME "AQUAgpusph")
set(PACKAGE_VERSION_MAJOR 5)
set(PACKAGE_VERSION_MINOR 0)
set(PACKAGE_VERSION_PATCH 4)
set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
set(AQUAGPUSPH_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
set(AQUAGPUSPH_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
set(AQUAGPUSPH_VERSION_PATCH ${PACKAGE_VERSION_PATCH})
set(AQUAGPUSPH_VERSION ${PACKAGE_VERSION})
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

set (CMAKE_CXX_STANDARD 17)

# include local modules
include(AddFileDependencies)

# Setup CTest to build tests
include(CTest)

set(AQUAGPUSPH_CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")
list(APPEND CMAKE_MODULE_PATH "${AQUAGPUSPH_CMAKE_MODULE_PATH}")

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# No Debug/Release output paths
set(DEBUG_MAIN_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(RELEASE_MAIN_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(WIN32)
    set(PLATFORM_CP xcopy /Y /S)
    set(PLATFORM_MK mkdir)
else(WIN32)
    set(PLATFORM_CP cp)
    set(PLATFORM_MK mkdir -p)
endif(WIN32)

# ===================================================== #
# Output directories for install target                 #
# ===================================================== #
if(WIN32)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Installation root directory")
else(WIN32)
    set(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "Installation root directory")
    set(INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif(WIN32)

set(CMAKE_INSTALL_DATADIR share/aquagpusph CACHE PATH "Output directory for data and resource files")
set(CMAKE_INSTALL_DOCDIR share/doc/aquagpusph CACHE PATH "Output directory for documentation and license files")
set(CMAKE_INSTALL_BINDIR bin CACHE PATH "Output directory for binary files")
set(CMAKE_INSTALL_LIBDIR lib CACHE PATH "Output directory for Library files")

# used as compiler defines
set(RESOURCEDIR "${CMAKE_INSTALL_DATADIR}")
set(DOCDIR "${CMAKE_INSTALL_DOCDIR}")

message(STATUS "prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "datadir: ${CMAKE_INSTALL_DATADIR}")
message(STATUS "docdir: ${CMAKE_INSTALL_DOCDIR}")
message(STATUS "bindir: ${CMAKE_INSTALL_BINDIR}")


# ===================================================== #
# All the options for the build process                 #
# ===================================================== #
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)
option(AQUAGPUSPH_USE_MPI "Build AQUAgpusph with MPI support, for multidevice running." ON)
option(AQUAGPUSPH_USE_VTK "Build AQUAgpusph with VTK output format support." ON)
option(AQUAGPUSPH_USE_MUPARSER "Replace ExprTk (Experimental support) by MuParser." YES)
option(AQUAGPUSPH_BUILD_EXAMPLES "Build AQUAgpusph examples" ON)
option(AQUAGPUSPH_BUILD_DOC "Build AQUAgpusph documentation" OFF)
option(AQUAGPUSPH_GPU_PROFILE "Profile the GPU during the runtime (consuming additional resources)" OFF)
option(AQUAGPUSPH_CHRONO_EXAMPLES "If AQUAGPUSPH_BUILD_EXAMPLES is ON, build the examples using PROJECTCHRONO" ON)
mark_as_advanced(
  AQUAGPUSPH_GPU_PROFILE
  AQUAGPUSPH_USE_MUPARSER
  AQUAGPUSPH_CHRONO_EXAMPLES
)
if(WIN32)
    set(AQUAGPUSPH_VENDOR_DLLS_FOLDER true CACHE PATH "${CMAKE_CURRENT_SOURCE_DIR}/")
    mark_as_advanced(AQUAGPUSPH_VENDOR_DLLS_FOLDER)
endif()

# ===================================================== #
# Search the packages                                   #
# ===================================================== #

# xxd
find_package(XXD REQUIRED)
if(NOT XXD_FOUND)
message(FATAL_ERROR "xxd not found, but ${PACKAGE_NAME} requires it. Please install vim!")
endif(NOT XXD_FOUND)

# MPI
if(AQUAGPUSPH_USE_MPI)
    find_package(MPI)
    set(HAVE_MPI TRUE)
endif(AQUAGPUSPH_USE_MPI)

# Python
find_package(Python REQUIRED COMPONENTS
    Interpreter
    Development
    NumPy)

if(NOT Python_FOUND)
message(FATAL_ERROR "Python not found, but ${PACKAGE_NAME} requires it. Please install Python!")
endif(NOT Python_FOUND)

# Xerces-C
find_package(XercesC REQUIRED)

if(NOT XercesC_FOUND)
message(FATAL_ERROR "Xerces-C not found, but ${PACKAGE_NAME} requires it. Please install Xerces-C!")
endif(NOT XercesC_FOUND)

# OpenCL
find_package(OpenCL REQUIRED)

if(NOT OpenCL_FOUND)
message(FATAL_ERROR "OpenCL not found, but ${PACKAGE_NAME} requires it. Please install OpenCL!")
endif(NOT OpenCL_FOUND)

if(OpenCL_VERSION_MAJOR LESS 2)
message(FATAL_ERROR "OpenCL ${OpenCL_VERSION_MAJOR}.${OpenCL_VERSION_MINOR} found, but OpenCL >= 2.0 is required.")
endif(OpenCL_VERSION_MAJOR LESS 2)

# VTK
if(AQUAGPUSPH_USE_VTK)
    find_package(VTK REQUIRED COMPONENTS
        CommonCore
        IOXML)
    if(NOT VTK_FOUND)
        message(FATAL_ERROR "VTK not found, but AQUAGPUSPH_USE_VTK is ON. Install VTK or set AQUAGPUSPH_USE_VTK=OFF")
    endif(NOT VTK_FOUND)
    if (VTK_VERSION VERSION_LESS "8.90.0")
        message(FATAL_ERROR "VTK >= 8.9 is required if AQUAGPUSPH_USE_VTK is ON. Update VTK or set AQUAGPUSPH_USE_VTK=OFF")
    endif (VTK_VERSION VERSION_LESS "8.90.0")
    message(STATUS "Found VTK ${VTK_VERSION}")
    set(HAVE_VTK TRUE)
endif(AQUAGPUSPH_USE_VTK)

# NLOHMANN JSON
# VTK 9.4 is already integrating NLOHMANN JSON. To avoid collisions we better
# use that vendored library
if(AQUAGPUSPH_USE_VTK AND (VTK_VERSION VERSION_GREATER_EQUAL 9.4))
    add_definitions(-DVTK_NLOHMANN_JSON)
else()
    find_package(nlohmann_json REQUIRED)
    if(NOT nlohmann_json_FOUND)
    message(FATAL_ERROR "NLOHMANN JSON not found, but ${PACKAGE_NAME} requires it. Please install NLOHMANN JSON!")
    endif(NOT nlohmann_json_FOUND)
endif()

# muparser
if(AQUAGPUSPH_USE_MUPARSER)
    find_package(MuParser REQUIRED)

    if(NOT MUPARSER_FOUND)
    message(FATAL_ERROR "MuParser not found, but AQUAGPUSPH_USE_MUPARSER is on. Install MuParser or set AQUAGPUSPH_USE_MUPARSER=OFF")
    endif(NOT MUPARSER_FOUND)
    set(HAVE_MUPARSER TRUE)
endif(AQUAGPUSPH_USE_MUPARSER)

if(AQUAGPUSPH_BUILD_EXAMPLES AND AQUAGPUSPH_CHRONO_EXAMPLES)
    find_package(Chrono)

    if(NOT Chrono_FOUND)
        message(FATAL_ERROR "PROJECTCHRONO not found, but AQUAGPUSPH_CHRONO_EXAMPLES is on. Install PROJECTCHRONO or set AQUAGPUSPH_CHRONO_EXAMPLES=OFF")
    endif(NOT Chrono_FOUND)
    set(HAVE_CHRONO TRUE)
endif(AQUAGPUSPH_BUILD_EXAMPLES AND AQUAGPUSPH_CHRONO_EXAMPLES)


# GPU profile
if(AQUAGPUSPH_GPU_PROFILE)
    add_definitions(-DHAVE_GPUPROFILE)
endif(AQUAGPUSPH_GPU_PROFILE)

# Tests building
if(BUILD_TESTING)
    find_program(BASH_PROGRAM bash)
    if(NOT BASH_PROGRAM)
        message(FATAL_ERROR "Bash not found, but BUILD_TESTING has been set. Install Bash or set BUILD_TESTING=OFF")
    endif(NOT BASH_PROGRAM)
    mark_as_advanced(BASH_PROGRAM)
endif()

# ===================================================== #
# Search special packages needed to build doc           #
# ===================================================== #
if(AQUAGPUSPH_BUILD_DOC)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen not found, but AQUAGPUSPH_BUILD_DOC is ON. Install Doxygen or set AQUAGPUSPH_BUILD_DOC=OFF")
    endif(NOT DOXYGEN_FOUND)
    if(NOT DOXYGEN_DOT_FOUND)
        message(FATAL_ERROR "Graphviz not found, but AQUAGPUSPH_BUILD_DOC is ON. Install Graphviz or set AQUAGPUSPH_BUILD_DOC=OFF")
    endif(NOT DOXYGEN_DOT_FOUND)
endif(AQUAGPUSPH_BUILD_DOC)

# ===================================================== #
# config.h stuff                                        #
# ===================================================== #
configure_file(config.h.cmake ${CMAKE_BINARY_DIR}/aquagpusph/config.h)

# ===================================================== #
# Global Compiler and Linker Settings                   #
# ===================================================== #
# check for 64-bit platform
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "Platform is 64-bit")
else(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "Platform is 32-bit")
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DAQUA_DEBUG")
message(STATUS "DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

# ===================================================== #
# Compilation parts                                     #
# ===================================================== #
if(WIN32)
    set(OPTIONAL_SRCS ${OPTIONAL_SRCS} win32/getopt.c
                                       win32/dlfcn.c)
    set(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/aquagpusph/win32)
endif()
if(HAVE_MPI)
    set(OPTIONAL_LIBS ${OPTIONAL_LIBS} MPI::MPI_CXX)
    set(OPTIONAL_SRCS ${OPTIONAL_SRCS} CalcServer/MPISync.cpp)
endif(HAVE_MPI)
if(HAVE_MUPARSER)
    set(OPTIONAL_INCLUDE_PATH ${OPTIONAL_INCLUDE_PATH} ${MUPARSER_INCLUDE_DIRS})
    set(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${MUPARSER_LIBRARIES})
    set(OPTIONAL_SRCS ${OPTIONAL_SRCS} Tokenizer/Tokenizer_muparser.cpp)
else(HAVE_MUPARSER)
    set(OPTIONAL_SRCS ${OPTIONAL_SRCS} Tokenizer/Tokenizer_exprtk.cpp)
endif(HAVE_MUPARSER)
if(HAVE_VTK)
    set(OPTIONAL_LIBS ${OPTIONAL_LIBS} VTK::CommonCore VTK::IOXML)
endif(HAVE_VTK)
if(nlohmann_json_FOUND)
    set(OPTIONAL_LIBS ${OPTIONAL_LIBS} nlohmann_json::nlohmann_json)
endif(nlohmann_json_FOUND)
include_directories(
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/aquagpusph/ext
    ${OPTIONAL_INCLUDE_PATH}
)
set(DEP_LIBS
    Python::Python
    Python::NumPy
    XercesC::XercesC
    OpenCL::OpenCL
    ${OPTIONAL_LIBS}
)
if(NOT WIN32)
    add_compile_options(-Wno-write-strings -Wno-deprecated -Wall -Wextra -Wpedantic)
endif()
add_subdirectory(aquagpusph)
add_subdirectory(resources)
if(AQUAGPUSPH_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif(AQUAGPUSPH_BUILD_EXAMPLES)
if(AQUAGPUSPH_BUILD_DOC)
    add_subdirectory(doc)
endif(AQUAGPUSPH_BUILD_DOC)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# ===================================================== #
# Packaging                                             #
# ===================================================== #
include(InstallRequiredSystemLibraries)

if(WIN32)
    file(GLOB_RECURSE FNAMES RELATIVE ${AQUAGPUSPH_VENDOR_DLLS_FOLDER} "${AQUAGPUSPH_VENDOR_DLLS_FOLDER}/*.dll")
    foreach(FNAME ${FNAMES})
        install(FILES ${AQUAGPUSPH_VENDOR_DLLS_FOLDER}/${FNAME} DESTINATION ${CMAKE_INSTALL_BINDIR})
    endforeach()
endif()

set(CPACK_PACKAGE_VENDOR "Jose Luis Cercos-Pita")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${AQUAGPUSPH_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${AQUAGPUSPH_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${AQUAGPUSPH_VERSION_PATCH}")
if(WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "AQUAgpusph windows installer")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "AQUAgpusph-${AQUAGPUSPH_VERSION_MAJOR}.${AQUAGPUSPH_VERSION_MINOR}")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
    set(CPACK_NSIS_HELP_LINK "https:\\\\\\\\gitlab.com/sanguinariojoe/aquagpusph/-/wikis/home")
    set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\canal.etsin.upm.es/aquagpusph/")
    set(CPACK_NSIS_CONTACT "jlcercos@gmail.com")
    set(CPACK_NSIS_MODIFY_PATH ON)
    file(DOWNLOAD
        https://www.python.org/ftp/python/${Python_VERSION}/python-${Python_VERSION}-amd64.exe
        ${CMAKE_CURRENT_BINARY_DIR}/python-amd64.exe
        SHOW_PROGRESS
    )
    install(
        FILES
            ${CMAKE_CURRENT_BINARY_DIR}/python-amd64.exe
        DESTINATION
            ${CMAKE_INSTALL_BINDIR}
    )
    list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS "  MessageBox MB_YESNO 'Install Python ${Python_VERSION}?' /SD IDYES IDNO endInstallPython")
    list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS "    ExecWait '\\\$INSTDIR\\\\bin\\\\python-amd64.exe'")
    list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS "    ExecWait 'py -${Python_VERSION_MAJOR}.${Python_VERSION_MINOR} -m pip install --upgrade numpy'")
    list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS "  endInstallPython:")
    string (REPLACE ";" "\n" CPACK_NSIS_EXTRA_INSTALL_COMMANDS "${CPACK_NSIS_EXTRA_INSTALL_COMMANDS}")
elseif(UNIX)
    set(CPACK_GENERATOR "STGZ")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Redistributable package")
endif()

include(CPack)

if(WIN32)
    configure_file(${PROJECT_SOURCE_DIR}/cMake/NSIS.definitions.nsh.in
                   ${CMAKE_CURRENT_BINARY_DIR}/cMake/NSIS.definitions.nsh)
endif()


# ===================================================== #
# Show a brief                                          #
# ===================================================== #
message("=====================================================")
message("Target: AQUAgpusph")
message("Binary destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
message("Resources destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources")
if(AQUAGPUSPH_BUILD_EXAMPLES)
    message("Examples destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples")
endif(AQUAGPUSPH_BUILD_EXAMPLES)
if(AQUAGPUSPH_BUILD_DOC)
    message("Documentation destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DOCDIR}")
endif(AQUAGPUSPH_BUILD_DOC)
message("OPTIONS:")
if(AQUAGPUSPH_USE_MPI)
    message("    - With MPI")
else(AQUAGPUSPH_USE_MPI)
    message("    - Without MPI")
endif(AQUAGPUSPH_USE_MPI)
if(AQUAGPUSPH_USE_VTK)
    message("    - With VTK")
else(AQUAGPUSPH_USE_VTK)
    message("    - Without VTK")
endif(AQUAGPUSPH_USE_VTK)
if(AQUAGPUSPH_USE_MUPARSER)
    message("    - With MuParser")
else(AQUAGPUSPH_USE_MUPARSER)
    message("    - With exprtk")
endif(AQUAGPUSPH_USE_MUPARSER)
message("=====================================================")
