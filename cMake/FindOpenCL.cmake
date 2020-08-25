# Module for locating OpenCL.
#
# Customizable variables:
#   OPENCL_ROOT_DIR
#     Specifies OpenCL's root directory. The find module uses this variable to
#     locate OpenCL. The variable will be filled automatically unless explicitly
#     set using CMake's -D command-line option. Instead of setting a CMake
#     variable, an environment variable called OCLROOT can be used.
#     While locating the root directory, the module will try to detect OpenCL
#     implementations provided by AMD's Accelerated Parallel Processing SDK,
#     NVIDIA's GPU Computing Toolkit and Intel's OpenCL SDK by examining the
#     AMDAPPSDKROOT, CUDA_PATH and INTELOCLSDKROOT environment variables,
#     respectively.
#
# Read-only variables:
#   OPENCL_FOUND
#     Indicates whether OpenCL has been found.
#
#   OPENCL_INCLUDE_DIRS
#     Specifies the OpenCL include directories.
#
#   OPENCL_LIBRARIES
#     Specifies the OpenCL libraries that should be passed to
#     target_link_libararies.
#
#
# Copyright (c) 2012 Sergiu Dotenco
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTOPENCLLAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

INCLUDE (FindPackageHandleStandardArgs)

IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET (_OPENCL_POSSIBLE_LIB_SUFFIXES lib/Win64 lib/x86_64 lib/x64)
ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET (_OPENCL_POSSIBLE_LIB_SUFFIXES lib/Win32 lib/x86)
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8)

LIST(APPEND _OPENCL_POSSIBLE_LIB_SUFFIXES lib/nvidia-current)

FIND_PATH(OPENCL_ROOT_DIR
  NAMES OpenCL/cl.h
        CL/cl.h
        include/CL/cl.h
        include/nvidia-current/CL/cl.h
  HINTS ${CUDA_TOOLKIT_ROOT_DIR}
  PATHS ENV OCLROOT
        ENV AMDAPPSDKROOT
        ENV CUDA_PATH
        ENV INTELOCLSDKROOT
  PATH_SUFFIXES cuda
  DOC "OpenCL root directory")

FIND_PATH(OPENCL_INCLUDE_DIRS
  NAMES OpenCL/cl.h CL/cl.h
  HINTS ${OPENCL_ROOT_DIR}
  PATH_SUFFIXES include include/nvidia-current
  DOC "OpenCL include directory")

FIND_LIBRARY(OPENCL_LIBRARIES
  NAMES OpenCL
  HINTS ${OPENCL_ROOT_DIR}
  PATH_SUFFIXES ${_OPENCL_POSSIBLE_LIB_SUFFIXES}
  DOC "OpenCL library")

IF(OPENCL_INCLUDE_DIRS AND OPENCL_LIBRARIES)
  SET(_OPENCL_VERSION_TEST_SOURCE
"
#if __APPLE__
    #include <OpenCL/cl.h>
#else /* !__APPLE__ */
    #include <CL/cl.h>
#endif /* __APPLE__ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

cl_int parse_version(char* v, cl_uint *major, cl_uint *minor)
{
    char *v_major, *v_minor, *tmp;
    cl_uint i_major, i_minor;
    size_t pos;
    if (strncmp(v, \"OpenCL \", 7) != 0) {
        return CL_INVALID_VALUE;
    }

    v_major = &(v[7]);
    pos = strcspn(v_major, \".\");
    if ((pos == strlen(v_major)) || (pos == 0))
        return CL_INVALID_VALUE;

    tmp = (char*)malloc((pos + 1) * sizeof(char));
    strncpy(tmp, v_major, pos);
    tmp[pos] = '\\0';
    i_major = (unsigned int)atoi(tmp);
    free(tmp);

    v_minor = &(v_major[pos]);
    pos = strcspn(v_minor, \" \");
    if ((pos == strlen(v_minor)) || (pos == 0))
        return CL_INVALID_VALUE;

    tmp = (char*)malloc((pos + 1) * sizeof(char));
    strncpy(tmp, v_minor, pos);
    tmp[pos] = '\\0';
    i_minor = (unsigned int)atoi(tmp);
    free(tmp);

    if ((i_major < *major) || ((i_major == *major) && (i_minor < *minor))) {
        *major = i_major;
        *minor = i_minor;
        return CL_SUCCESS;
    }

    return CL_INVALID_PLATFORM;
}

int main()
{
    char *version = NULL;
    cl_int result;
    cl_uint n_platforms, i;
    cl_platform_id *ids;
    size_t n;
    cl_uint major = UINT_MAX, minor = UINT_MAX;

    result = clGetPlatformIDs(0, NULL, &n_platforms);
    if ((result != CL_SUCCESS) || (n_platforms == 0))
        return EXIT_FAILURE;

    ids = (cl_platform_id*)malloc(n_platforms * sizeof(cl_platform_id));
    result = clGetPlatformIDs(n_platforms, ids, NULL);
    if (result != CL_SUCCESS)
        return EXIT_FAILURE;

    for (i = 0; i < n_platforms; i++) {
        result = clGetPlatformInfo(ids[i], CL_PLATFORM_VERSION, 0, NULL, &n);
        if (result != CL_SUCCESS)
            continue;

        char *v = (char*)malloc(n * sizeof(char));
        if (!v)
            continue;

        result = clGetPlatformInfo(ids[i], CL_PLATFORM_VERSION, n, v, NULL);
        if (result != CL_SUCCESS)
            continue;

        result = parse_version(v, &major, &minor);
        if (result == CL_SUCCESS)
            version = v;
        else
            free(v);
    }
    
    printf(\"%s\", version);
    fflush(stdout);
    free(ids);
    free(version);
    return version != NULL ? EXIT_SUCCESS : EXIT_FAILURE;
}
")

  SET(_OPENCL_VERSION_SOURCE
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/openclversion.c")

  FILE(WRITE ${_OPENCL_VERSION_SOURCE} "${_OPENCL_VERSION_TEST_SOURCE}\n")

  TRY_RUN(_OPENCL_VERSION_RUN_RESULT _OPENCL_VERSION_COMPILE_RESULT
    ${CMAKE_BINARY_DIR} ${_OPENCL_VERSION_SOURCE}
    RUN_OUTPUT_VARIABLE _OPENCL_VERSION_STRING
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${OPENCL_INCLUDE_DIRS}"
                "-DLINK_LIBRARIES:STRING=${OPENCL_LIBRARIES}")

  IF(_OPENCL_VERSION_RUN_RESULT EQUAL 0)
    STRING(REGEX REPLACE "OpenCL[ \t]+([0-9]+)\\.[0-9]+.*" "\\1"
      OPENCL_VERSION_MAJOR "${_OPENCL_VERSION_STRING}")
    STRING(REGEX REPLACE "OpenCL[ \t]+[0-9]+\\.([0-9]+).*" "\\1"
      OPENCL_VERSION_MINOR "${_OPENCL_VERSION_STRING}")

    SET(OPENCL_VERSION_COMPONENTS 2)
    SET(OPENCL_VERSION "${OPENCL_VERSION_MAJOR}.${OPENCL_VERSION_MINOR}")
  ENDIF(_OPENCL_VERSION_RUN_RESULT EQUAL 0)

  IF("${OPENCL_VERSION}" STREQUAL "")
    MESSAGE(WARNING "Cannot determine OpenCL's version")
  ENDIF("${OPENCL_VERSION}" STREQUAL "")
  
ENDIF(OPENCL_INCLUDE_DIRS AND OPENCL_LIBRARIES)

MARK_AS_ADVANCED(OPENCL_ROOT_DIR OPENCL_INCLUDE_DIRS OPENCL_LIBRARIES
 OPENCL_VERSION_MAJOR OPENCL_VERSION_MINOR)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenCL REQUIRED_VARS OPENCL_ROOT_DIR
  OPENCL_INCLUDE_DIRS OPENCL_LIBRARIES
  OPENCL_VERSION_MAJOR OPENCL_VERSION_MINOR)
