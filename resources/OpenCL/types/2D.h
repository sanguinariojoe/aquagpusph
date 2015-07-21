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
 * @brief Type definitions for the OpenCL kernels (2D version).
 */

#define unit unsigned int
#define vec2 float2
#define vec3 float3
#define vec4 float4
#define ivec2 int2
#define ivec3 int3
#define ivec4 int4
#define uivec2 uint2
#define uivec3 uint3
#define uivec4 uint4

/** @def vec
 * @brief Vector of real components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 4 components
 * 
 * This type should be used for input arguments, but for the local variables
 * maybe you can consider using vec_xyz (which has just 3 components in 3D)
 */
#define vec float2

/** @def ivec
 * @brief Vector of integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 4 components
 * 
 * This type should be used for input arguments, but for the local variables
 * maybe you can consider using ivec_xyz (which has just 3 components in 3D)
 */
#define ivec int2

/** @def ivec
 * @brief Vector of unsigned integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 4 components
 * 
 * This type should be used for input arguments, but for the local variables
 * maybe you can consider using uivec_xyz (which has just 3 components in 3D)
 */
#define uivec uint2

/** @def VEC_ZERO
 * @brief Null #vec, i.e. filled with zero components.
 */
#define VEC_ZERO (float2)(0.f,0.f)
/** @def VEC_ONE
 * @brief Ones #vec, i.e. filled with one components.
 */
#define VEC_ONE (float2)(1.f, 1.f)
/** @def VEC_ALL_ONE
 * @brief VEC_ONE.
 */
#define VEC_ALL_ONE VEC_ONE
/** @def VEC_INFINITY
 * @brief Infinity #vec, i.e. filled with infinity components.
 */
#define VEC_INFINITY (float2)(INFINITY, INFINITY)
/** @def VEC_ALL_INFINITY
 * @brief VEC_INFINITY.
 */
#define VEC_ALL_INFINITY VEC_INFINITY
/** @def VEC_NEG_INFINITY
 * @brief -Infinity #vec, i.e. filled with -infinity components.
 */
#define VEC_NEG_INFINITY (-VEC_INFINITY)
/** @def VEC_ALL_NEG_INFINITY
 * @brief VEC_NEG_INFINITY.
 */
#define VEC_ALL_NEG_INFINITY (-VEC_ALL_INFINITY)

/** @def vec_xyz
 * @brief Vector of real components with the minimum number of components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 3 components
 *
 * This type can be used for the local variables to reduce the VGPRs.
 */
#define vec_xyz vec2

/** @def ivec_xyz
 * @brief Vector of integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 3 components
 *
 * This type can be used for the local variables to reduce the VGPRs.
 */
#define ivec_xyz ivec2

/** @def uivec_xyz
 * @brief Vector of unsigned integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 3 components
 *
 * This type can be used for the local variables to reduce the VGPRs.
 */
#define uivec_xyz uivec2

/** @def XYZ
 * @brief Convenient access to the vector components.
 * 
 * It is useful to be used with #vec_xyz, #ivec_xyz and #uivec_xyz type:
 *   - 2D = .xy
 *   - 3D = .xyz
 */
#define XYZ xy

/** @def C_I
 * @brief Utility to can redefine the cell of the particle to be  computed.
 * 
 * It can be used for mirrrored particles, which are temporary associated to a
 * different cell.
 *
 * @see BEGIN_LOOP_OVER_NEIGHS
 */
#define C_I() const uint c_i = icell[i]

/** @def BEGIN_LOOP_OVER_NEIGHS
 * @brief Loop over the neighs to compute the interactions.
 * 
 * All the code between this macro and END_LOOP_OVER_NEIGHS will be executed for
 * all the neighbours.
 *
 * To use this macro, the main particle (the one which the interactions are
 * intended to be computed) should be identified by an unsigned integer variable
 * i, while the resulting neighs will be automatically identified by the
 * unsigned integer variable j. To discard a neighbour particle, remember
 * calling \code{.c}j++\endcode before \code{.c}continue\endcode
 *
 * The following variables will be declared, and therefore cannot be used
 * elsewhere:
 *   - c_i: The cell where the particle i is placed
 *   - ci: Index of the cell of the neighbour particle j, in the x direction
 *   - cj: Index of the cell of the neighbour particle j, in the x direction
 *   - ck: Index of the cell of the neighbour particle j, in the x direction
 *   - c_j: Index of the cell of the neighbour particle j
 *   - j: Index of the neighbour particle.
 *
 * @see END_LOOP_OVER_NEIGHS
 */
#define BEGIN_LOOP_OVER_NEIGHS()                                               \
    C_I();                                                                     \
    for(int ci = -1; ci <= 1; ci++) {                                          \
        for(int cj = -1; cj <= 1; cj++) {                                      \
            const uint c_j = c_i +                                             \
                             ci +                                              \
                             cj * n_cells.x;                                   \
            uint j = ihoc[c_j];                                                \
            while((j < N) && (icell[j] == c_j)) {

/** @def END_LOOP_OVER_NEIGHS
 * @brief End of the loop over the neighs to compute the interactions.
 * 
 * @see BEGIN_LOOP_OVER_NEIGHS
 */
#define END_LOOP_OVER_NEIGHS()                                                 \
                j++;                                                           \
            }                                                                  \
        }                                                                      \
    }