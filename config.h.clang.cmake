/* This generated file is for internal use. Do not include it from headers. */

#ifdef CLANG_CONFIG_H
#error config.h can only be included once
#else
#define CLANG_CONFIG_H

/* Bug report URL. */
#define BUG_REPORT_URL "${BUG_REPORT_URL}"

/* Default linker to use. */
#define CLANG_DEFAULT_LINKER "${CLANG_DEFAULT_LINKER}"

/* Default C/ObjC standard to use. */
#cmakedefine CLANG_DEFAULT_STD_C LangStandard::lang_${CLANG_DEFAULT_STD_C}

/* Default C++/ObjC++ standard to use. */
#cmakedefine CLANG_DEFAULT_STD_CXX LangStandard::lang_${CLANG_DEFAULT_STD_CXX}

/* Default C++ stdlib to use. */
#define CLANG_DEFAULT_CXX_STDLIB "${CLANG_DEFAULT_CXX_STDLIB}"

/* Default runtime library to use. */
#define CLANG_DEFAULT_RTLIB "${CLANG_DEFAULT_RTLIB}"

/* Default objcopy to use */
#define CLANG_DEFAULT_OBJCOPY "${CLANG_DEFAULT_OBJCOPY}"

/* Default OpenMP runtime used by -fopenmp. */
#define CLANG_DEFAULT_OPENMP_RUNTIME "${CLANG_DEFAULT_OPENMP_RUNTIME}"

/* Default architecture for OpenMP offloading to Nvidia GPUs. */
#define CLANG_OPENMP_NVPTX_DEFAULT_ARCH "${CLANG_OPENMP_NVPTX_DEFAULT_ARCH}"

/* Multilib suffix for libdir. */
#define CLANG_LIBDIR_SUFFIX "${CLANG_LIBDIR_SUFFIX}"

/* Relative directory for resource files */
#define CLANG_RESOURCE_DIR "${CLANG_RESOURCE_DIR}"

/* Directories clang will search for headers */
#define C_INCLUDE_DIRS "${C_INCLUDE_DIRS}"

/* Directories clang will search for configuration files */
#cmakedefine CLANG_CONFIG_FILE_SYSTEM_DIR "${CLANG_CONFIG_FILE_SYSTEM_DIR}"
#cmakedefine CLANG_CONFIG_FILE_USER_DIR "${CLANG_CONFIG_FILE_USER_DIR}"

/* Default <path> to all compiler invocations for --sysroot=<path>. */
#define DEFAULT_SYSROOT "${DEFAULT_SYSROOT}"

/* Directory where gcc is installed. */
#define GCC_INSTALL_PREFIX "${GCC_INSTALL_PREFIX}"

/* Define if we have libxml2 */
#cmakedefine CLANG_HAVE_LIBXML ${CLANG_HAVE_LIBXML}

/* Define if we have z3 and want to build it */
#cmakedefine CLANG_ANALYZER_WITH_Z3 ${CLANG_ANALYZER_WITH_Z3}

/* Define if we have sys/resource.h (rlimits) */
#cmakedefine CLANG_HAVE_RLIMITS ${CLANG_HAVE_RLIMITS}

/* The LLVM product name and version */
#define BACKEND_PACKAGE_STRING "${BACKEND_PACKAGE_STRING}"

/* Linker version detected at compile time. */
#cmakedefine HOST_LINK_VERSION "${HOST_LINK_VERSION}"

/* pass --build-id to ld */
#cmakedefine ENABLE_LINKER_BUILD_ID

/* enable x86 relax relocations by default */
#cmakedefine01 ENABLE_X86_RELAX_RELOCATIONS

/* Enable the experimental new pass manager by default */
#cmakedefine01 ENABLE_EXPERIMENTAL_NEW_PASS_MANAGER

/* Enable each functionality of modules */
#cmakedefine01 CLANG_ENABLE_ARCMT
#cmakedefine01 CLANG_ENABLE_OBJC_REWRITER
#cmakedefine01 CLANG_ENABLE_STATIC_ANALYZER


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

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Define to 1 to build zipios++ sources with iostream. */
#cmakedefine USE_STD_IOSTREAM

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

/* NCurses */
#cmakedefine HAVE_NCURSES

/* VTK */
#cmakedefine HAVE_VTK


#endif
