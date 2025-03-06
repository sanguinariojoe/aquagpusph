#ifndef CLANG_CONFIG_H
#define CLANG_CONFIG_H

/* Name of package */
#define PACKAGE "AQUAgpusph"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "jlcercos@gmail.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "AQUAgpusph"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "AQUAgpusph ${AQUAGPUSPH_VERSION}"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "AQUAgpusph"

/* Define to the version of this package. */
#define PACKAGE_VERSION "${AQUAGPUSPH_VERSION}"

/* Version number of package */
#define VERSION "${AQUAGPUSPH_VERSION}"

/* ========================================================================== */
/* Dependencies                                                               */
/* ========================================================================== */
/* Lowest OpenCL version to be supported */
#define CL_TARGET_OPENCL_VERSION 200
#define OPENCL_PLATFORM_MAJOR ${OpenCL_VERSION_MAJOR}
#define OPENCL_PLATFORM_MINOR ${OpenCL_VERSION_MINOR}

/* Numpy version */
#define NUMPY_VERSION_MAJOR ${NUMPY_VERSION_MAJOR}
#define NUMPY_VERSION_MINOR ${NUMPY_VERSION_MINOR}

/* MPI */
#cmakedefine HAVE_MPI

/* VTK */
#cmakedefine HAVE_VTK

/* ExprTk or MuParser */
#cmakedefine HAVE_MUPARSER

#endif
