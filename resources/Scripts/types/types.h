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
 * @brief Generic types definition file.
 *
 * This file is redirecting to either 2D.hcl or 3D.hcl, depending on #HAVE_3D
 */

#ifndef INFINITY
    #define INFINITY FLT_MAX
#endif

#define vec2 float2
#define vec3 float3
#define vec4 float4
#define vec8 float8
#define vec16 float16
#define dvec2 double2
#define dvec3 double3
#define dvec4 double4
#define dvec8 double8
#define dvec16 double16
#define ivec2 int2
#define ivec3 int3
#define ivec4 int4
#define ivec8 int8
#define ivec16 int16
#define lvec2 long2
#define lvec3 long3
#define lvec4 long4
#define lvec8 long8
#define lvec16 long16
#define uivec2 uint2
#define uivec3 uint3
#define uivec4 uint4
#define uivec8 uint8
#define uivec16 uint16
#define ulvec2 ulong2
#define ulvec3 ulong3
#define ulvec4 ulong4
#define ulvec8 ulong8
#define ulvec16 ulong16
#define svec2 usize2
#define svec3 usize3
#define svec4 usize4
#define svec8 usize8
#define svec16 usize16
#define ssvec2 ssize2
#define ssvec3 ssize3
#define ssvec4 ssize4
#define ssvec8 ssize8
#define ssvec16 ssize16

/** @brief Helper function for #CONVERT
 *
 * The helper is required because the preprocessor is only recursively expanding
 * macros if the definition is not affected by # nor ## string operators.
 * Then, this inner function is concatenating the unexpanded words, while
 * #CONVERT is effectively expanding the type name.
 */
#define _CONVERT(TYPE) convert_ ## TYPE

/** @brief Conversor between complex types.
 *
 * In OpenCL, to convert between complex types the functions convert_TYPEN
 * should be used. Otherwise casting errors will be received.
 *
 * This definition provides a convenient function to become used with the
 * overloaded types vec, ivec and uivec.
 *
 * For instance, to convert a vec variable, v, to an ivec variable, you can
 * call CONVERT(ivec, v);
 */
#define CONVERT(TYPE, v) _CONVERT(TYPE)(v)

/** @brief Utility to can redefine the cell of the particle to be computed.
 * 
 * It can be used for mirrrored particles, which are temporary associated to a
 * different cell.
 *
 * @see #BEGIN_LOOP_OVER_NEIGHS
 */
#define C_I() const usize c_i = icell[i]

#ifdef HAVE_3D
    #include "resources/Scripts/types/3D.h"
#else
    #include "resources/Scripts/types/2D.h"
#endif

/** @brief Macro to easily add the parameters to run #BEGIN_NEIGHS macro,
 * interacting with the local set of particles, i.e. the particles handled by
 * this process.
 * 
 * @see #BEGIN_NEIGHS
 * @note The number of particles, N, is not included
 */
#define LINKLIST_LOCAL_PARAMS                                                  \
    const __global usize * icell,                                              \
    const __global usize * ihoc,                                               \
    svec4 n_cells

/** @brief Macro to easily add the parameters to run #BEGIN_NEIGHS macro,
 * interacting with the local set of particles, i.e. the particles handled by
 * this process.
 * 
 * @see #BEGIN_NEIGHS
 * @note The number of particles, N, is not included
 */
#define LINKLIST_REMOTE_PARAMS                                                 \
    const __global usize * icell,                                              \
    const __global usize * mpi_icell,                                          \
    const __global usize * mpi_ihoc,                                           \
    svec4 n_cells

/** @brief Loop over the particle-by-particle neighbour chains
 *
 * All the code between this macro and END_FOR_NEIGHS will be executed for
 * all the neighbours.
 *
 * The resulting neighs will be automatically identified by the unsigned
 * integer variable j. You are always entitled to discard a neighbour by
 * executing \code{.c} continue \endcode
 *
 * The following variables will be declared, and therefore cannot be redeclared
 * within the loop scope:
 *   - __jhoc_id: Index of the JHOC entry to be read
 *   - j: Index of the neighbour particle
 *
 * @param NPARTS Number of particles (usually \code{.c} N \endcode)
 * @param JHOC Array of head of chains (usually \code{.c} jhoc \endcode)
 * @see #END_FOR_NEIGHS
 * @note This macro is created to replace the old #BEGIN_NEIGHS, which was
 * performing suboptimally
 */
#define FOR_NEIGHS(NPARTS, JHOC)                                               \
    for(unsigned int row = 0; row < NNC; row++) {                              \
        const unsigned int __jhoc_id = i + row * NPARTS                        \
        for(unsigned int j = JHOC[__jhoc_id].x; j < JHOC[__jhoc_id].y; j++) {

/** @brief End of the loop over the neighs to compute the interactions.
 * 
 * @see #FOR_NEIGHS
 * @note This macro is created to replace the old #END_NEIGHS, which was
 * performing suboptimally
 */
#define END_FOR_NEIGHS()                                                       \
        }                                                                      \
    }

/** @brief Null #vec, i.e. filled with zero components.
 */
#define VEC2_ZERO ((float2)(0.f,0.f))
/** @brief Null #vec, i.e. filled with zero components.
 */
#define VEC3_ZERO ((float3)(0.f,0.f,0.f))
/** @brief Null #vec, i.e. filled with zero components.
 */
#define VEC4_ZERO ((float4)(0.f,0.f,0.f,0.f))
/** @brief Null #vec, i.e. filled with zero components.
 */
#define VEC8_ZERO ((float8)(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f))
/** @brief Null #vec, i.e. filled with zero components.
 */
#define VEC16_ZERO ((float16)(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f))

#include "resources/Scripts/types/reductions.h"
