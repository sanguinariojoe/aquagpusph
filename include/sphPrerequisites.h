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

#include <CL/cl.h>

// Assume the 3D version if 2D has not been set
#ifndef HAVE_2D
    /** @def HAVE_3D
     * @brief AQUAgpusph 3D version.
     *
     * If it is defined the AQUAgpusph 3D version will be built, otherwise the
     * 2D version will be selected.
     */
    #define HAVE_3D
#endif

#ifndef vec2
    /// @brief Vector of 2 real components.
    typedef cl_float2 vec2;
#endif
#ifndef vec3
    /// @brief Vector of 3 real components.
    typedef cl_float3 vec3;
#endif
#ifndef vec4
    /// @brief Vector of 4 real components.
    typedef cl_float4 vec4;
#endif

#ifndef ivec2
    /// @brief Vector of 2 integer components.
    typedef cl_int2 ivec2;
#endif
#ifndef ivec3
    /// @brief Vector of 3 integer components.
    typedef cl_int3 ivec3;
#endif
#ifndef ivec4
    /// @brief Vector of 4 integer components.
    typedef cl_int4 ivec4;
#endif

#ifndef ivec2
    /// @brief Vector of 2 unsigned integer components.
    typedef cl_uint2 uivec2;
#endif
#ifndef ivec3
    /// @brief Vector of 3 unsigned integer components.
    typedef cl_uint3 uivec3;
#endif
#ifndef ivec4
    /// @brief Vector of 4 unsigned integer components.
    typedef cl_uint4 uivec4;
#endif

#ifdef HAVE_3D
    #ifndef vec
        /** @brief Vector of real components.
         *
         * The number of components depends on weather the 2D version or 3D
         * version is compiled:
         *   - 2D = 2 components
         *   - 3D = 4 components
         */
        typedef vec4 vec;
    #endif
    #ifndef ivec
        /** @brief Vector of integer components.
         *
         * The number of components depends on weather the 2D version or 3D
         * version is compiled:
         *   - 2D = 2 components
         *   - 3D = 4 components
         */
        typedef ivec4 ivec;
    #endif
    #ifndef uivec
        /** @brief Vector of unsigned integer components.
         *
         * The number of components depends on weather the 2D version or 3D
         * version is compiled:
         *   - 2D = 2 components
         *   - 3D = 4 components
         */
        typedef uivec4 uivec;
    #endif
    #ifndef matrix
        /** @brief Matrix of float numbers.
         *
         * The number of components depends on weather the 2D version or 3D
         * version is compiled:
         *   - 2D = 2x2 components
         *   - 3D = 4x4 components
         */
        typedef cl_float16 matrix;
    #endif
#else
    #ifndef vec
        typedef vec2 vec;
    #endif
    #ifndef ivec
        typedef ivec2 ivec;
    #endif
    #ifndef uivec
        typedef uivec2 uivec;
    #endif
    #ifndef matrix
        typedef cl_float4 matrix;
    #endif
#endif

#ifdef HAVE_3D
/// Null vector
static const vec VEC_ZERO = {0.f, 0.f, 0.f, 0.f};
/// X axis unit vector
static const vec VEC_X = {1.f, 0.f, 0.f, 0.f};
/// Y axis unit vector
static const vec VEC_Y = {0.f, 1.f, 0.f, 0.f};
/// Z axis unit vector
static const vec VEC_Z = {0.f, 0.f, 1.f, 0.f};
#else
/// Null vector
static const vec VEC_ZERO = {0.f, 0.f};
/// X axis unit vector
static const vec VEC_X = {1.f, 0.f};
/// Y axis unit vector
static const vec VEC_Y = {0.f, 1.f};
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

#ifndef __ERROR_SHOW_TIME__
    #ifndef HAVE_NCURSES
        /** @def __ERROR_SHOW_TIME__
         * @brief Time let to the user to notice and read the error.
         *
         * 0 seconds if ncurses is disabled, and 3 seconds otherwise.
         */
        #define __ERROR_SHOW_TIME__ 0u
    #else
        /** @def __ERROR_SHOW_TIME__
         * @brief Time let to the user to notice and read the error.
         *
         * 0 seconds if ncurses is disabled, and 3 seconds otherwise.
         */
        #define __ERROR_SHOW_TIME__ 3u
    #endif
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
inline const std::string methodAndClassName(const std::string& pretty_function)
{
    std::string all_name, method_name, class_name;
    size_t begin, end;

    // Filter the name removing the preceding and trailing types data
    end    = pretty_function.find("(");
    begin  = pretty_function.substr(0,end).rfind(" ") + 1;
    end   -= begin;
    all_name = pretty_function.substr(begin,end);

    // Get the method name
    begin = all_name.rfind("::");
    if (begin == std::string::npos){
        // There are no class name
        methodAndClassName_str = all_name;
        return methodAndClassName_str;
    }
    method_name = all_name.substr(begin+2, std::string::npos);
    end = begin;
    begin = all_name.substr(0,end).rfind("::");
    if(begin == std::string::npos)
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

/** \class mat sphPrerequisites.h sphPrerequisites.h
 * @brief Matrix of real components.
 *
 * The number of components depends on weather the 2D version or 3D version is
 * compiled:
 *   - 2D = 2 \f$ \times \f$ 2 matrix
 *   - 3D = 4 \f$ \times \f$ 4 matrix
 * @note This matrix class have not errors control, so a bad item request may
 * result in a segmentation fault.
 */
class mat
{
public:
    /** @brief Constructor.
     *
     * The matrix will be initialized as a zeroes one.
     */
    mat(){
        row[0] = VEC_ZERO;
        row[1] = VEC_ZERO;
        #ifdef HAVE_3D
            row[0] = VEC_ZERO;
            row[1] = VEC_ZERO;
        #endif
    }

    /** @brief Subscript operator to return a matrix row.
     * @param index Row index.
     * @return Row vector.
     */
    vec const& operator[] (const unsigned &index) const{
        return row[index];
    }
    /** Subscript operator to return a matrix row.
     * @param index Row index.
     * @return Row vector.
     */
    vec& operator[] (const unsigned &index){
        return row[index];
    }

    /** Matrix-Vector Inner product operator.
     * @param V Vector to operate.
     * @return Resulting vector.
     */
    vec operator* (const vec &V){
        unsigned int i,j;
        vec R;
        #ifdef HAVE_3D
            R.x = row[0].x*V.x + row[0].y*V.y + row[0].z*V.z + row[0].w*V.w;
            R.y = row[1].x*V.x + row[1].y*V.y + row[1].z*V.z + row[1].w*V.w;
            R.z = row[2].x*V.x + row[2].y*V.y + row[2].z*V.z + row[2].w*V.w;
            R.w = row[3].x*V.x + row[3].y*V.y + row[3].z*V.z + row[3].w*V.w;
        #else
            R.x = row[0].x*V.x + row[0].y*V.y;
            R.y = row[1].x*V.x + row[1].y*V.y;
        #endif
        return R;
    }

    /** Assignment operator.
     * @param M Matrix to copy.
     * @return Copied matrix.
     */
    mat& operator= (const mat &M) {
        if (this == &M) {   // Same object, simply return it
            return *this;
        }
        #ifdef HAVE_3D
            row[0].x=M[0].x; row[0].y=M[0].y; row[0].z=M[0].z; row[0].w=M[0].w;
            row[1].x=M[1].x; row[1].y=M[1].y; row[1].z=M[1].z; row[1].w=M[1].w;
            row[2].x=M[2].x; row[2].y=M[2].y; row[2].z=M[2].z; row[2].w=M[2].w;
            row[3].x=M[3].x; row[3].y=M[3].y; row[3].z=M[3].z; row[3].w=M[3].w;
        #else
            row[0].x=M[0].x; row[0].y=M[0].y;
            row[1].x=M[1].x; row[1].y=M[1].y;
        #endif
        return *this;
    }
private:
    /// Matrix rows
    #ifdef HAVE_3D
        vec row[4];
    #else
        vec row[2];
    #endif
};

#endif // SPHPREREQUISITES_H_INCLUDED
