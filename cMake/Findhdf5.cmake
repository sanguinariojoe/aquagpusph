# - Find the hdf5 include file and library
#
#  HDF5_FOUND - system has H5Part
#  HDF5_INCLUDE_DIR - the H5Part include directory
#  HDF5_LIBRARIES - The libraries needed to use H5Part
#  HDF5_LIBRARY - set for backwards compatibility with 2.4 CMake

FIND_LIBRARY(HDF5_LIBRARY NAMES hdf5)

FIND_PATH(HDF5_INCLUDE_PATH NAMES hdf5.h)

# Need to provide the *_LIBRARIES
SET(HDF5_LIBRARIES ${HDF5_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set HDF5_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(hdf5 DEFAULT_MSG HDF5_LIBRARY HDF5_INCLUDE_PATH)

MARK_AS_ADVANCED(HDF5_INCLUDE_PATH HDF5_LIBRARY)
