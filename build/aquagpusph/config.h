#ifndef AQUAGPUSPH_CONFIG_H
#define AQUAGPUSPH_CONFIG_H

/* Name of package */
#define PACKAGE "AQUAgpusph"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "jlcercos@gmail.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "AQUAgpusph"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "AQUAgpusph 5.0.4"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "AQUAgpusph"

/* Define to the version of this package. */
#define PACKAGE_VERSION "5.0.4"

/* Version number of package */
#define VERSION "5.0.4"

/* ========================================================================== */
/* Dependencies                                                               */
/* ========================================================================== */
/* Lowest OpenCL version to be supported */
#define CL_TARGET_OPENCL_VERSION 200
#define OPENCL_PLATFORM_MAJOR 3
#define OPENCL_PLATFORM_MINOR 0

/* Numpy version */
#define NUMPY_VERSION_MAJOR 
#define NUMPY_VERSION_MINOR 

/* MPI */
#define HAVE_MPI

/* VTK */
#define HAVE_VTK

/* ExprTk or MuParser */
#define HAVE_MUPARSER

#endif
