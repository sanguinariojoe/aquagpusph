# Find libclang
#
# It defines the following variables
# CLANG_FOUND - True if libclang has been found.
# CLANG_ROOT - where to find libclang
# CLANG_INCLUDE_DIRS - where to find libclang include files
# CLANG_LIBRARIES - libraries to be linked

INCLUDE (FindPackageHandleStandardArgs)

SET(llvm_config_names llvm-config
  llvm-config-5.0 llvm-config-mp-5.0
  llvm-config-4.0 llvm-config-mp-4.0
  llvm-config-3.7 llvm-config-mp-3.7
  llvm-config-3.6 llvm-config-mp-3.6
  llvm-config-3.5 llvm-config-mp-3.5
  llvm-config-3.4 llvm-config-mp-3.4
  llvm-config-3.3 llvm-config-mp-3.3
  llvm-config-3.2 llvm-config-mp-3.2
  llvm-config-3.1 llvm-config-mp-3.1)
FIND_PROGRAM(LLVM_CONFIG_EXECUTABLE NAMES ${llvm_config_names})

IF(LLVM_CONFIG_EXECUTABLE)
  MESSAGE(STATUS "llvm-config found at: ${LLVM_CONFIG_EXECUTABLE}")
ELSE(LLVM_CONFIG_EXECUTABLE)
  MESSAGE(FATAL_ERROR "Could NOT find llvm-config executable")
ENDIF(LLVM_CONFIG_EXECUTABLE)

EXECUTE_PROCESS(
  COMMAND ${LLVM_CONFIG_EXECUTABLE} --prefix
  OUTPUT_VARIABLE LLVM_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

FIND_PATH(CLANG_ROOT
  NAMES
    include/clang-c/Index.h
  HINTS
    ${LLVM_ROOT}
  DOC
    "libclang root directory")


FIND_PATH(CLANG_INCLUDE_DIRS
  NAMES
    clang-c/Index.h
  HINTS
    ${CLANG_ROOT}/include
  DOC
    "libclang include directory")

FIND_LIBRARY(CLANG_LIBRARIES
  NAMES
    clang
  HINTS
    ${CLANG_ROOT}/lib
  PATH_SUFFIXES
    ${_CLANG_POSSIBLE_LIB_SUFFIXES}
  DOC
    "libclang path")

MARK_AS_ADVANCED(LLVM_ROOT LLVM_CONFIG_EXECUTABLE CLANG_ROOT CLANG_INCLUDE_DIRS
                 CLANG_LIBRARIES)

FIND_PACKAGE_HANDLE_STANDARD_ARGS (CLang REQUIRED_VARS CLANG_ROOT
  CLANG_INCLUDE_DIRS CLANG_LIBRARIES)
