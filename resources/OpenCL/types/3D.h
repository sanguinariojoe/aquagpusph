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
 * @brief Type definitions for the OpenCL kernels (3D version).
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
 */
#define vec float4

/** @def ivec
 * @brief Vector of integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 4 components
 */
#define ivec int4

/** @def ivec
 * @brief Vector of unsigned integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 4 components
 */
#define uivec uint4

/** @def VEC_ZERO
 * @brief Null #vec, i.e. filled with zero components.
 */
#define VEC_ZERO (float4)(0.f,0.f,0.f,0.f)

#define VEC_ONE (float4)(1.f, 1.f, 1.f, 0.f)
#define VEC_ALL_ONE (float4)(1.f, 1.f, 1.f, 1.f)
#define VEC_INFINITY (float4)(INFINITY, INFINITY, INFINITY, 0.f)
#define VEC_ALL_INFINITY (float4)(INFINITY, INFINITY, INFINITY, INFINITY)
#define VEC_NEG_INFINITY (-VEC_INFINITY)
#define VEC_ALL_NEG_INFINITY (-VEC_ALL_INFINITY)

