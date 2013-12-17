# Locate Xerces-C include paths and libraries
# Xerces-C can be found at http://xml.apache.org/xerces-c/
# Written by Frederic Heem, frederic.heem _at_ telsey.it
# Modified by Jos van den Oever

# This module defines
# XERCESC_INCLUDE_PATH, where to find ptlib.h, etc.
# XERCESC_LIBRARIES, the libraries to link against to use pwlib.
# XERCESC_FOUND, If false, don't try to use pwlib.

FIND_PATH(XERCESC_INCLUDE_DIR xercesc/dom/DOM.hpp
  "[HKEY_CURRENT_USER\\software\\xerces-c\\src]"
  "[HKEY_CURRENT_USER\\xerces-c\\src]"
  $ENV{XERCESCROOT}/src/
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(XERCESC_LIBRARIES
  NAMES
    xerces-c
  PATHS
    "[HKEY_CURRENT_USER\\software\\xerces-c\\lib]"
    "[HKEY_CURRENT_USER\\xerces-c\\lib]"
    $ENV{XERCESCROOT}/${LIB_DESTINATION}
    /usr/local/${LIB_DESTINATION}
    /usr/${LIB_DESTINATION}
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Xerces XERCESC_INCLUDE_DIR XERCESC_LIBRARIES)

SET(XERCESC_INCLUDE_PATH ${XERCESC_INCLUDE_DIR})
SET(XERCESC_FOUND ${XERCES_FOUND})
MARK_AS_ADVANCED( XERCESC_INCLUDE_DIR XERCESC_LIBRARIES )