# - Find the H5Part include file and library
#
#  H5PART_FOUND - system has H5Part
#  H5PART_INCLUDE_DIR - the H5Part include directory
#  H5PART_LIBRARIES - The libraries needed to use H5Part
#  H5PART_LIBRARY - set for backwards compatibility with 2.4 CMake

FIND_LIBRARY(H5PART_LIBRARY NAMES H5Part )

FIND_PATH(H5PART_INCLUDE_PATH NAMES H5Part.h )

# Need to provide the *_LIBRARIES
SET(H5PART_LIBRARIES ${H5PART_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set H5PART_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(H5Part DEFAULT_MSG H5PART_LIBRARY H5PART_INCLUDE_PATH)

MARK_AS_ADVANCED(H5PART_INCLUDE_PATH H5PART_LIBRARY)

 
