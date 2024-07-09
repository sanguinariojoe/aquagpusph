/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @brief Set of definitions and macros related with the implementation.
 */

#ifndef SPHPREREQUISITES_H_INCLUDED
#define SPHPREREQUISITES_H_INCLUDED

/** @def XSTR
 * @brief Just an alias for #STR
 */
#define XSTR(x) STR(x)
/** @def STR
 * @brief An auxiliary macro to can print defined variables in compilation time.
 *
 * For instance:
 * @code{.c}
    #define _VAR 1
    #pragma message "_VAR = " STR(_VAR)
   @endcode
 */
#define STR(x) #x

#include <string>

// CMake configuration file
#include <config.h>

#if __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#ifndef icl
/** @def icl
 * @brief Signed integer number
 */
#define icl cl_int
#endif
#ifndef lcl
/** @def lcl
 * @brief Signed 64bits integer number
 */
#define lcl cl_long
#endif
#ifndef uicl
/** @def uicl
 * @brief Unsigned integer number
 */
#define uicl cl_uint
#endif
#ifndef ulcl
/** @def lcl
 * @brief Unsigned 64bits integer number
 */
#define ulcl cl_ulong
#endif
#ifndef fcl
/** @def fcl
 * @brief Floating point number
 */
#define fcl cl_float
#endif
#ifndef dcl
/** @def fcl
 * @brief Double precission floating point number
 */
#define dcl cl_double
#endif

#ifndef vec2
/** @def vec2
 * @brief Vector of 2 real components.
 */
#define vec2 cl_float2
#endif
#ifndef vec3
/** @def vec3
 * @brief Vector of 3 real components.
 */
#define vec3 cl_float3
#endif
#ifndef vec4
/** @def vec4
 * @brief Vector of 4 real components.
 */
#define vec4 cl_float4
#endif
#ifndef vec8
/** @def vec8
 * @brief Vector of 8 real components.
 */
#define vec8 cl_float8
#endif
#ifndef vec16
/** @def vec16
 * @brief Vector of 16 real components.
 */
#define vec16 cl_float16
#endif

#ifndef matrix2x2
/** @def matrix2x2
 * @brief 2D matrix of float numbers.
 */
#define matrix2x2 cl_float4
#endif
#ifndef matrix4x4
/** @def matrix4x4
 * @brief 3D matrix of float numbers.
 */
#define matrix4x4 cl_float16
#endif

#ifndef dvec2
/** @def dvec2
 * @brief Vector of 2 64bits real components.
 */
#define dvec2 cl_double2
#endif
#ifndef dvec3
/** @def dvec3
 * @brief Vector of 3 64bits real components.
 */
#define dvec3 cl_double3
#endif
#ifndef dvec4
/** @def dvec4
 * @brief Vector of 4 64bits real components.
 */
#define dvec4 cl_double4
#endif
#ifndef dvec8
/** @def dvec8
 * @brief Vector of 8 64bits real components.
 */
#define dvec8 cl_double8
#endif
#ifndef dvec16
/** @def dvec16
 * @brief Vector of 16 64bits real components.
 */
#define dvec16 cl_double16
#endif

#ifndef ivec2
/** @def ivec2
 * @brief Vector of 2 integer components.
 */
#define ivec2 cl_int2
#endif
#ifndef ivec3
/** @def ivec3
 * @brief Vector of 3 integer components.
 */
#define ivec3 cl_int3
#endif
#ifndef ivec4
/** @def ivec4
 * @brief Vector of 4 integer components.
 */
#define ivec4 cl_int4
#endif
#ifndef ivec8
/** @def ivec8
 * @brief Vector of 8 integer components.
 */
#define ivec8 cl_int8
#endif
#ifndef ivec16
/** @def ivec16
 * @brief Vector of 16 integer components.
 */
#define ivec16 cl_int16
#endif

#ifndef lvec2
/** @def lvec2
 * @brief Vector of 2 64bits integer components.
 */
#define lvec2 cl_long2
#endif
#ifndef lvec3
/** @def lvec3
 * @brief Vector of 3 64bits integer components.
 */
#define lvec3 cl_long3
#endif
#ifndef lvec4
/** @def lvec4
 * @brief Vector of 4 64bits integer components.
 */
#define lvec4 cl_long4
#endif
#ifndef lvec8
/** @def lvec8
 * @brief Vector of 8 64bits integer components.
 */
#define lvec8 cl_long8
#endif
#ifndef lvec16
/** @def lvec16
 * @brief Vector of 16 64bits integer components.
 */
#define lvec16 cl_long16
#endif

#ifndef uivec2
/** @def uivec2
 * @brief Vector of 2 unsigned integer components.
 */
#define uivec2 cl_uint2
#endif
#ifndef uivec3
/** @def uivec3
 * @brief Vector of 3 unsigned integer components.
 */
#define uivec3 cl_uint3
#endif
#ifndef uivec4
/** @def uivec4
 * @brief Vector of 4 unsigned integer components.
 */
#define uivec4 cl_uint4
#endif
#ifndef uivec8
/** @def uivec8
 * @brief Vector of 8 unsigned integer components.
 */
#define uivec8 cl_uint8
#endif
#ifndef uivec16
/** @def uivec16
 * @brief Vector of 16 unsigned integer components.
 */
#define uivec16 cl_uint16
#endif

#ifndef ulvec2
/** @def ulvec2
 * @brief Vector of 2 64bits unsigned integer components.
 */
#define ulvec2 cl_ulong2
#endif
#ifndef ulvec3
/** @def ulvec3
 * @brief Vector of 3 64bits unsigned integer components.
 */
#define ulvec3 cl_ulong3
#endif
#ifndef ulvec4
/** @def ulvec4
 * @brief Vector of 4 64bits unsigned integer components.
 */
#define ulvec4 cl_ulong4
#endif
#ifndef ulvec8
/** @def ulvec8
 * @brief Vector of 8 64bits unsigned integer components.
 */
#define ulvec8 cl_ulong8
#endif
#ifndef ulvec16
/** @def ulvec16
 * @brief Vector of 16 64bits unsigned integer components.
 */
#define ulvec16 cl_ulong16
#endif

#ifndef __CL_MIN_LOCALSIZE__
/** @def __CL_MIN_LOCALSIZE__
 * @brief Minimum local work size to execute kernels.
 *
 * The selected value have been tested in the following platforms:
 *   - AMD CPU
 *   - Intel CPU
 *   - AMD GPU
 *   - NVidia GPU
 */
#define __CL_MIN_LOCALSIZE__ 64
#endif
#ifndef __CL_MAX_LOCALSIZE__
/** @def __CL_MAX_LOCALSIZE__
 * @brief Maximum local work size to execute kernels.
 *
 * The selected value have been tested in the following platforms:
 *   - AMD CPU
 *   - Intel CPU
 *   - AMD GPU
 *   - NVidia GPU
 */
#define __CL_MAX_LOCALSIZE__ 1024
#endif

/// Helper string for #methodAndClassName function.
static std::string methodAndClassName_str;

/** @brief Function to extract the class and function from
 * @paramname{__PRETTY_FUNCTION__} macro.
 *
 * The GNU compiler (GCC) macro @paramname{__PRETTY_FUNCTION__} contains a lot
 * of information about the name space, the returning value, and the parameters.
 *
 * In order to print log information this data should be simplified.
 *
 * @param pretty_function GCC @paramname{__PRETTY_FUNCTION__} macro
 * @return Class and function name ("Class::function")
 * @see #__METHOD_CLASS_NAME__
 */
inline const std::string
methodAndClassName(const std::string& pretty_function)
{
	std::string all_name, method_name, class_name;
	size_t begin, end;

	// Filter the name removing the preceding and trailing types data
	end = pretty_function.find("(");
	begin = pretty_function.substr(0, end).rfind(" ") + 1;
	end -= begin;
	all_name = pretty_function.substr(begin, end);

	// Get the method name
	begin = all_name.rfind("::");
	if (begin == std::string::npos) {
		// There are no class name
		methodAndClassName_str = all_name;
		return methodAndClassName_str;
	}
	method_name = all_name.substr(begin + 2, std::string::npos);
	end = begin;
	begin = all_name.substr(0, end).rfind("::");
	if (begin == std::string::npos)
		begin = 0;
	else
		begin += 2;
	end -= begin;
	class_name = all_name.substr(begin, end);
	methodAndClassName_str = class_name + "::" + method_name;
	return methodAndClassName_str;
}

/** @def __METHOD_CLASS_NAME__
 * @brief Returns automatically the current class and function names.
 * @see methodAndClassName()
 * @see #addMessageF
 */
#define __METHOD_CLASS_NAME__ methodAndClassName(__PRETTY_FUNCTION__)

/** @def UNUSED
 * @brief Set a parameter as unused
 */
#ifdef __GNUC__
#  define UNUSED_PARAM __attribute__((__unused__))
#else
#  define UNUSED_PARAM
#endif

#endif // SPHPREREQUISITES_H_INCLUDED
