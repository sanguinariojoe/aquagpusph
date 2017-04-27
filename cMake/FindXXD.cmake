# Find xxd executable
#
# It defines the following variables
# XXD_FOUND - True if xxd has been found.
# XXD_BIN - xxd binary path

INCLUDE (FindPackageHandleStandardArgs)

# Lets numpy say us where we can find the headers
FIND_PROGRAM(XXD_BIN xxd)

IF(NOT XXD_BIN)
    MESSAGE(FATAL_ERROR "Failure finding xxd binary")
ENDIF(NOT XXD_BIN)

MARK_AS_ADVANCED(XXD_BIN)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(XXD REQUIRED_VARS XXD_BIN)
