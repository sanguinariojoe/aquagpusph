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
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed. Please, run cmake in a different folder.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
project(AQUAgpusph)
set(AQUAGPUSPH_VERSION "4.1.1")

set(PACKAGE_NAME "AQUAgpusph")
set(PACKAGE_VERSION_MAJOR "4")
set(PACKAGE_VERSION_MINOR "1")
set(PACKAGE_VERSION_PATCH "1")
set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

# include local modules
include(AddFileDependencies)

# Setup CTest to build tests
include(CTest)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")

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
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Installation root directory")
else(WIN32)
    SET(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "Installation root directory")
    SET(INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif(WIN32)

SET(CMAKE_INSTALL_DATADIR share/aquagpusph CACHE PATH "Output directory for data and resource files")
SET(CMAKE_INSTALL_INCLUDEDIR include/AQUAgpusph CACHE PATH "Output directory for header files")
SET(CMAKE_INSTALL_DOCDIR share/doc/aquagpusph CACHE PATH "Output directory for documentation and license files")
SET(CMAKE_INSTALL_BINDIR bin CACHE PATH "Output directory for binary files")
SET(CMAKE_INSTALL_LIBDIR lib CACHE PATH "Output directory for library files")

# used as compiler defines
SET(RESOURCEDIR "${CMAKE_INSTALL_DATADIR}")
SET(DOCDIR "${CMAKE_INSTALL_DOCDIR}")

MESSAGE(STATUS "prefix: ${CMAKE_INSTALL_PREFIX}")
MESSAGE(STATUS "datadir: ${CMAKE_INSTALL_DATADIR}")
MESSAGE(STATUS "docdir: ${CMAKE_INSTALL_DOCDIR}")
MESSAGE(STATUS "bindir: ${CMAKE_INSTALL_BINDIR}")
MESSAGE(STATUS "libdir: ${CMAKE_INSTALL_LIBDIR}")


# ===================================================== #
# All the options for the build process                 #
# ===================================================== #
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE)
OPTION(AQUAGPUSPH_3D "Build 3D AQUAgpusph version. 2D and 3D versions can be installed on the same system." ON)
OPTION(AQUAGPUSPH_USE_MPI "Build AQUAgpusph with MPI support, for multidevice running." ON)
OPTION(AQUAGPUSPH_USE_NCURSES "Build AQUAgpusph with Ncurses terminal support." ON)
OPTION(AQUAGPUSPH_USE_VTK "Build AQUAgpusph with VTK output format support (on development)." ON)
OPTION(AQUAGPUSPH_BUILD_TOOLS "Build AQUAgpusph tools" ON)
OPTION(AQUAGPUSPH_BUILD_EXAMPLES "Build AQUAgpusph examples" ON)
OPTION(AQUAGPUSPH_BUILD_DOC "Build AQUAgpusph documentation" OFF)
OPTION(AQUAGPUSPH_GPU_PROFILE "Profile the GPU during the runtime (consuming additional resources)" OFF)
MARK_AS_ADVANCED(
  AQUAGPUSPH_GPU_PROFILE
)

# ===================================================== #
# Search the packages                                   #
# ===================================================== #
# 2D/3D
IF(AQUAGPUSPH_3D)
    ADD_DEFINITIONS(-DHAVE_3D)
ELSE(AQUAGPUSPH_3D)
    ADD_DEFINITIONS(-DHAVE_2D)
ENDIF(AQUAGPUSPH_3D)

# xxd
FIND_PACKAGE(XXD REQUIRED)
IF(NOT XXD_FOUND)
MESSAGE(FATAL_ERROR "xxd not found, but ${PACKAGE_NAME} requires it. Please install vim!")
ENDIF(NOT XXD_FOUND)

# MPI
IF(AQUAGPUSPH_USE_MPI)
    FIND_PACKAGE(MPI)
    IF(NOT MPI_FOUND)
        MESSAGE(FATAL_ERROR "MPI not found, but AQUAGPUSPH_USE_MPI is ON. Install MPI or set AQUAGPUSPH_USE_MPI=OFF")
    ENDIF(NOT MPI_FOUND)

    SET(HAVE_MPI TRUE)
ENDIF(AQUAGPUSPH_USE_MPI)

# Python
SET(PYTHON_ADDITIONAL_VERSIONS "2.5" "2.6" "2.7" "2.8" "2.9" "3.0" "3.1" "3.2" "3.3")
FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(PythonInterp REQUIRED)

IF(NOT PYTHONLIBS_FOUND)
MESSAGE(FATAL_ERROR "Python not found, but ${PACKAGE_NAME} requires it. Please install Python!")
ENDIF(NOT PYTHONLIBS_FOUND)

# Numpy
FIND_PACKAGE(Numpy REQUIRED)
IF(NOT NUMPY_FOUND)
MESSAGE(FATAL_ERROR "numpy not found, but ${PACKAGE_NAME} requires it. Please install numpy!")
ENDIF(NOT NUMPY_FOUND)

# Xerces-C
FIND_PACKAGE(Xerces REQUIRED)

IF(NOT XERCESC_FOUND)
MESSAGE(FATAL_ERROR "Xerces-C not found, but ${PACKAGE_NAME} requires it. Please install Xerces-C!")
ENDIF(NOT XERCESC_FOUND)

# OpenCL
FIND_PACKAGE(OpenCL REQUIRED)

IF(NOT OPENCL_FOUND)
MESSAGE(FATAL_ERROR "OpenCL not found, but ${PACKAGE_NAME} requires it. Please install OpenCL!")
ENDIF(NOT OPENCL_FOUND)

# CLang
FIND_PACKAGE(CLang REQUIRED)

IF(NOT CLANG_FOUND)
MESSAGE(FATAL_ERROR "libCLang not found, but ${PACKAGE_NAME} requires it. Please install libCLang!")
ENDIF(NOT CLANG_FOUND)

# NCURSES
IF(AQUAGPUSPH_USE_NCURSES)
    SET(CURSES_NEED_NCURSES TRUE)
    FIND_PACKAGE(Curses)
    IF(NOT CURSES_FOUND)
        MESSAGE(FATAL_ERROR "Curses not found, but AQUAGPUSPH_USE_NCURSES is ON. Install Curses or set AQUAGPUSPH_USE_NCURSES=OFF")
    ENDIF(NOT CURSES_FOUND)

    SET(HAVE_NCURSES TRUE)
ENDIF(AQUAGPUSPH_USE_NCURSES)

# VTK
IF(AQUAGPUSPH_USE_VTK)
    FIND_PACKAGE(VTK)
    IF(NOT VTK_FOUND)
        MESSAGE(FATAL_ERROR "VTK not found, but AQUAGPUSPH_USE_VTK is ON. Install vtk or set AQUAGPUSPH_USE_VTK=OFF")
    ELSE(NOT VTK_FOUND)
        MESSAGE(STATUS "Found VTK ${VTK_VERSION}")
    ENDIF(NOT VTK_FOUND)

    IF(NOT VTK_LIBRARIES)
        SET(VTK_LIBRARIES vtkCommon vtkHybrid vtkIO vtkFiltering)
    ENDIF(NOT VTK_LIBRARIES)

    SET(HAVE_VTK TRUE)
ENDIF(AQUAGPUSPH_USE_VTK)

# muparser
FIND_PACKAGE(MuParser REQUIRED)

IF(NOT MUPARSER_FOUND)
MESSAGE(FATAL_ERROR "MuParser not found, but ${PACKAGE_NAME} requires it. Please install MuParser!")
ENDIF(NOT MUPARSER_FOUND)

# GPU profile
IF(AQUAGPUSPH_GPU_PROFILE)
    ADD_DEFINITIONS(-DHAVE_GPUPROFILE)
ENDIF(AQUAGPUSPH_GPU_PROFILE)

# Tests building
IF(BUILD_TESTING)
    FIND_PROGRAM(BASH_PROGRAM bash)
    IF(NOT BASH_PROGRAM)
        MESSAGE(FATAL_ERROR "Bash not found, but BUILD_TESTING has been set. Install Bash or set BUILD_TESTING=OFF")
    ENDIF(NOT BASH_PROGRAM)
    MARK_AS_ADVANCED(BASH_PROGRAM)
ENDIF()

# ===================================================== #
# Search special packages needed to build doc           #
# ===================================================== #
IF(AQUAGPUSPH_BUILD_DOC)
    FIND_PACKAGE(Doxygen)
    IF(NOT DOXYGEN_FOUND)
        MESSAGE(FATAL_ERROR "Doxygen not found, but AQUAGPUSPH_BUILD_DOC is ON. Install Doxygen or set AQUAGPUSPH_BUILD_DOC=OFF")
    ENDIF(NOT DOXYGEN_FOUND)
    IF(NOT DOXYGEN_DOT_FOUND)
        MESSAGE(FATAL_ERROR "Graphviz not found, but AQUAGPUSPH_BUILD_DOC is ON. Install Graphviz or set AQUAGPUSPH_BUILD_DOC=OFF")
    ENDIF(NOT DOXYGEN_DOT_FOUND)
ENDIF(AQUAGPUSPH_BUILD_DOC)

# ===================================================== #
# config.h stuff                                        #
# ===================================================== #
include(cMake/ConfigureChecks.cmake)
if(CMAKE_COMPILER_IS_GNUCXX)
    configure_file(config.h.gcc.cmake ${CMAKE_BINARY_DIR}/config.h)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    configure_file(config.h.clang.cmake ${CMAKE_BINARY_DIR}/config.h)
endif()
add_definitions(-DHAVE_CONFIG_H)
add_definitions(-Wno-write-strings)
add_definitions(-Wno-deprecated)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

# CMAKE < 3.1 Change this when minimum CMake version is updated (see
# CMAKE_CXX_STANDARD)
SET(CMAKE_CXX_FLAGS "-std=c++11")
# Apparently -std=c++11 implies defining __STRICT_ANSI__ in GCC, disabling union
# structures in CL/cl_platform.h
add_definitions(-U__STRICT_ANSI__)

# ===================================================== #
# Global Compiler and Linker Settings                   #
# ===================================================== #
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/include)

# check for 64-bit platform
IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
    MESSAGE(STATUS "Platform is 64-bit")
ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
    MESSAGE(STATUS "Platform is 32-bit")
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8)

SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DAQUA_DEBUG")
MESSAGE(STATUS "DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

# ===================================================== #
# Compilation parts                                     #
# ===================================================== #
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(resources)
IF(AQUAGPUSPH_BUILD_TOOLS)
    ADD_SUBDIRECTORY(tools)
ENDIF(AQUAGPUSPH_BUILD_TOOLS)
IF(AQUAGPUSPH_BUILD_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF(AQUAGPUSPH_BUILD_EXAMPLES)
IF(AQUAGPUSPH_BUILD_DOC)
    ADD_SUBDIRECTORY(doc)
ENDIF(AQUAGPUSPH_BUILD_DOC)

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(tests)
ENDIF()

# ===================================================== #
# Show a brief                                          #
# ===================================================== #
MESSAGE("=====================================================")
IF(AQUAGPUSPH_3D)
    MESSAGE("Target: AQUAgpusph")
ELSE(AQUAGPUSPH_3D)
    MESSAGE("Target: AQUAgpusph2D")
ENDIF(AQUAGPUSPH_3D)
MESSAGE("Binary destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
MESSAGE("Libraries destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
MESSAGE("Resources destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/resources")
IF(AQUAGPUSPH_BUILD_EXAMPLES)
    MESSAGE("Examples destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/examples")
ENDIF(AQUAGPUSPH_BUILD_EXAMPLES)
IF(AQUAGPUSPH_BUILD_DOC)
    MESSAGE("Documentation destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DOCDIR}")
ENDIF(AQUAGPUSPH_BUILD_DOC)
MESSAGE("OPTIONS:")
IF(AQUAGPUSPH_USE_MPI)
    MESSAGE("    - With MPI")
ELSE(AQUAGPUSPH_USE_MPI)
    MESSAGE("    - Without MPI")
ENDIF(AQUAGPUSPH_USE_MPI)
IF(AQUAGPUSPH_USE_NCURSES)
    MESSAGE("    - With NCurses")
ELSE(AQUAGPUSPH_USE_NCURSES)
    MESSAGE("    - Without NCurses")
ENDIF(AQUAGPUSPH_USE_NCURSES)
IF(AQUAGPUSPH_USE_VTK)
    MESSAGE("    - With VTK")
ELSE(AQUAGPUSPH_USE_VTK)
    MESSAGE("    - Without VTK")
ENDIF(AQUAGPUSPH_USE_VTK)
MESSAGE("=====================================================")
