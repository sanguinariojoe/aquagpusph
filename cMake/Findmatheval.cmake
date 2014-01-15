# Locate libmatheval include paths and libraries
# libmatheval can be found at http://www.gnu.org/software/libmatheval
# Written by Jose Luis Cercos Pita

# This module defines
# MATHEVAL_INCLUDE_PATH, where to find matheval.h.
# MATHEVAL_LIBRARIES, the libraries to link against to use matheval.
# MATHEVAL_FOUND, If false, don't try to use matheval.

find_path(MATHEVAL_INCLUDE_DIR matheval.h
  "[HKEY_CURRENT_USER\\software\\libmatheval\\src]"
  "[HKEY_CURRENT_USER\\libmatheval\\src]"
  $ENV{MATHEVALROOT}/src/
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(MATHEVAL_LIBRARIES
  NAMES
    matheval
  PATHS
    "[HKEY_CURRENT_USER\\software\\libmatheval\\lib]"
    "[HKEY_CURRENT_USER\\libmatheval\\lib]"
    $ENV{MATHEVALROOT}/${LIB_DESTINATION}
    /usr/local/${LIB_DESTINATION}
    /usr/${LIB_DESTINATION}
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(matheval MATHEVAL_INCLUDE_DIR MATHEVAL_LIBRARIES)

SET(MATHEVAL_INCLUDE_PATH ${MATHEVAL_INCLUDE_DIR})
SET(MATHEVAL_FOUND ${MATHEVAL_FOUND})
MARK_AS_ADVANCED( MATHEVAL_INCLUDE_DIR MATHEVAL_LIBRARIES )
