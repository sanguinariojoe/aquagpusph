# Find libclang
#
# It defines the following variables
# CLANG_FOUND - True if libclang has been found.
# CLANG_ROOT - where to find libclang
# CLANG_INCLUDE_DIRS - where to find libclang include files
# CLANG_LIBRARIES - libraries to be linked

INCLUDE (FindPackageHandleStandardArgs)

FIND_PATH (CLANG_ROOT
  NAMES include/clang-c/Index.h
  HINTS
  ${LLVM_ROOT}
  DOC "libclang root directory")


FIND_PATH (CLANG_INCLUDE_DIR
  NAMES  clang-c/Index.h
  HINTS
  ${CLANG_ROOT}/include
  ${LLVM_INCLUDE_DIRS}
  DOC "libclang include directory")

FIND_LIBRARY(CLANG_LIBRARY NAMES clang
  HINTS
  ${CLANG_ROOT}/lib
  ${LLVM_LIBRARY_DIRS}
  PATH_SUFFIXES ${_CLANG_POSSIBLE_LIB_SUFFIXES}
  DOC "libclang")

SET (CLANG_INCLUDE_DIRS ${CLANG_INCLUDE_DIR})
SET (CLANG_LIBRARIES ${CLANG_LIBRARY})

MARK_AS_ADVANCED (CLANG_INCLUDE_DIR CLANG_LIBRARY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS (CLang REQUIRED_VARS CLANG_ROOT
  CLANG_INCLUDE_DIR CLANG_LIBRARY)
