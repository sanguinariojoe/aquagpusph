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
#       - Cercos Pita, Jose Luis                        #
#       - Miguel Gonzalez, Leo                          #
#       - Souto Iglesias, Antonio                       #
#                                                       #
# ===================================================== #


# ===================================================== #
# General CMake options                                 #
# ===================================================== #
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)
project(AQUAgpusph)
set(AQUAGPUSPH_VERSION "1.4.05")

set(PACKAGE_NAME "AQUAgpusph")
set(PACKAGE_VERSION_MAJOR "1")
set(PACKAGE_VERSION_MINOR "4")
set(PACKAGE_VERSION_PATCH "06")
set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

# include local modules
include(AddFileDependencies)

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
MESSAGE(STATUS "examplesdir: ${CMAKE_INSTALL_EXAMPLESDIR}")


# ===================================================== #
# All the options for the build process                 #
# ===================================================== #
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE)
OPTION(AQUAGPUSPH_3D "Build 3D AQUAgpusph version. 2D and 3D versions can be installed on the same system." ON)
OPTION(AQUAGPUSPH_ENABLE_MPI "Build AQUAgpusph with MPI support. It is necessary if you have H5Part parallel built version." ON)
OPTION(AQUAGPUSPH_USE_NCURSES "Build AQUAgpusph with Ncurses terminal support." ON)
OPTION(AQUAGPUSPH_USE_H5PART "Build AQUAgpusph with H5Part output format support." ON)
OPTION(AQUAGPUSPH_USE_VTK "Build AQUAgpusph with VTK output format support (on development)." ON)
OPTION(AQUAGPUSPH_USE_TECPLOT "Build AQUAgpusph with Tecplot output format support." ON)
OPTION(AQUAGPUSPH_BUILD_EXAMPLES "Build AQUAgpusph examples" ON)
OPTION(AQUAGPUSPH_BUILD_DOC "Build AQUAgpusph documentation" OFF)


# ===================================================== #
# Search the packages                                   #
# ===================================================== #
# 2D/3D
IF(AQUAGPUSPH_3D)
	ADD_DEFINITIONS(-DHAVE_3D)
ENDIF(AQUAGPUSPH_3D)


# Python
SET(PYTHON_ADDITIONAL_VERSIONS "2.5" "2.6" "2.7" "2.8" "2.9" "3.0" "3.1" "3.2" "3.3")
FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(PythonInterp REQUIRED)

IF(NOT PYTHONLIBS_FOUND)
MESSAGE(FATAL_ERROR "Python not found, but ${PACKAGE_NAME} requires it. Please install Python!")
ENDIF(NOT PYTHONLIBS_FOUND)

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

# NCURSES
IF(AQUAGPUSPH_USE_NCURSES)
	SET(CURSES_NEED_NCURSES TRUE)
	FIND_PACKAGE(Curses)
	IF(NOT CURSES_FOUND)
		MESSAGE(FATAL_ERROR "Curses not found, but AQUAGPUSPH_USE_NCURSES is ON. Install Curses or set AQUAGPUSPH_USE_NCURSES=OFF")
	ENDIF(NOT CURSES_FOUND)

	SET(HAVE_NCURSES TRUE)
ENDIF(AQUAGPUSPH_USE_NCURSES)

# H5PART
IF(AQUAGPUSPH_USE_H5PART)
	FIND_PACKAGE(hdf5)
	IF(NOT HDF5_FOUND)
		MESSAGE(FATAL_ERROR "hdf5 not found, but AQUAGPUSPH_USE_H5PART is ON. Install hdf5 or set AQUAGPUSPH_USE_HDF5=OFF")
	ENDIF(NOT HDF5_FOUND)

	FIND_PACKAGE(H5Part)
	IF(NOT H5PART_FOUND)
		MESSAGE(FATAL_ERROR "H5Part not found, but AQUAGPUSPH_USE_H5PART is ON. Install H5Part or set AQUAGPUSPH_USE_H5PART=OFF")
	ENDIF(NOT H5PART_FOUND)

	SET(HAVE_H5PART TRUE)
ENDIF(AQUAGPUSPH_USE_H5PART)

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

# Tecplot
IF(AQUAGPUSPH_USE_TECPLOT)
	SET(HAVE_TECPLOT TRUE)
ENDIF(AQUAGPUSPH_USE_TECPLOT)

# MPI
IF(AQUAGPUSPH_ENABLE_MPI)
	FIND_PACKAGE(MPI REQUIRED)
	IF(NOT MPI_FOUND)
		MESSAGE(FATAL_ERROR "MPI not found, but AQUAGPUSPH_ENABLE_MPI is ON. Install MPI or set AQUAGPUSPH_ENABLE_MPI=OFF")
	ENDIF(NOT MPI_FOUND)

	SET(CMAKE_C_COMPILER mpicc)
	SET(CMAKE_CXX_COMPILER mpicxx)
	SET(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
	SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
	IF(HAVE_H5PART)
		ADD_DEFINITIONS(-DPARALLEL_IO)
		ADD_DEFINITIONS(-DH5_HAVE_PARALLEL)
		ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
	ENDIF(HAVE_H5PART)

	SET(HAVE_MPI TRUE)
ENDIF(AQUAGPUSPH_ENABLE_MPI)

# Eigen3
FIND_PACKAGE(Eigen3 REQUIRED)

IF(NOT EIGEN3_FOUND)
MESSAGE(FATAL_ERROR "Eigen3 not found, but ${PACKAGE_NAME} requires it. Please install Eigen3!")
ENDIF(NOT EIGEN3_FOUND)

# matheval
FIND_PACKAGE(matheval REQUIRED)

IF(NOT MATHEVAL_FOUND)
MESSAGE(FATAL_ERROR "libmatheval not found, but ${PACKAGE_NAME} requires it. Please install libmatheval!")
ENDIF(NOT MATHEVAL_FOUND)

# GPU profile
IF(CMAKE_BUILD_TYPE EQUAL Debug)
	ADD_DEFINITIONS(-DHAVE_GPUPROFILE)
ENDIF(CMAKE_BUILD_TYPE EQUAL Debug)


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
if(CMAKE_COMPILER_IS_GNUCXX)
    include(cMake/ConfigureChecks.cmake)
    configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
    add_definitions(-DHAVE_CONFIG_H)
    add_definitions(-Wno-write-strings)
    add_definitions(-Wno-deprecated)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
    # get linker errors as soon as possible and not at runtime e.g. for modules
    if(UNIX)
    #    SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    endif(UNIX)
endif(CMAKE_COMPILER_IS_GNUCXX)


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

IF(MSVC)
	# set default compiler settings
	#SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /GF /GY")
	SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DAQUA_DEBUG -DHAVE_GPUPROFILE")
	# set default libs
	SET (CMAKE_C_STANDARD_LIBRARIES "kernel32.lib user32.lib gdi32.lib winspool.lib SHFolder.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib ")
	set (CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES}")
	# set linker flag /nodefaultlib
	set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB")
	SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB")

	# Mark 32 bit executables large address aware so they can use > 2GB address space
	# NOTE: This setting only has an effect on machines with at least 3GB of RAM, although it sets the linker option it doesn't set the the linker switch 'Enable Large Addresses'
	IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
		set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /LARGEADDRESSAWARE")
		SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
	ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
ELSE(MSVC)
	SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DAQUA_DEBUG -DHAVE_GPUPROFILE")
	MESSAGE(STATUS "DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
ENDIF(MSVC)
IF(MINGW)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mthreads")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mthreads -Wl,--export-all-symbols")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -mthreads -Wl,--export-all-symbols")
	LINK_LIBRARIES(-lgdi32)
ENDIF(MINGW)
# 0000661: cmake build on Mac OS: dealing with dylib versus so
IF(APPLE)
	SET(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
ENDIF(APPLE)
IF(UNIX)
	ADD_DEFINITIONS(-DUNIX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
ENDIF(UNIX)

# ===================================================== #
# Compilation parts                                     #
# ===================================================== #
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(resources)
IF(AQUAGPUSPH_BUILD_DOC)
	ADD_SUBDIRECTORY(doc)
ENDIF(AQUAGPUSPH_BUILD_DOC)
IF(AQUAGPUSPH_BUILD_EXAMPLES)
	ADD_SUBDIRECTORY(examples)
ENDIF(AQUAGPUSPH_BUILD_EXAMPLES)

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
IF(AQUAGPUSPH_USE_NCURSES)
	MESSAGE("    - With NCurses")
ELSE(AQUAGPUSPH_USE_NCURSES)
	MESSAGE("    - Without NCurses")
ENDIF(AQUAGPUSPH_USE_NCURSES)
IF(AQUAGPUSPH_USE_H5PART)
	MESSAGE("    - With H5Part")
ELSE(AQUAGPUSPH_USE_H5PART)
	MESSAGE("    - Without H5Part")
ENDIF(AQUAGPUSPH_USE_H5PART)
IF(AQUAGPUSPH_ENABLE_MPI)
	MESSAGE("    - MPI supported")
ELSE(AQUAGPUSPH_ENABLE_MPI)
	MESSAGE("    - MPI unsupported")
ENDIF(AQUAGPUSPH_ENABLE_MPI)
MESSAGE("=====================================================")
