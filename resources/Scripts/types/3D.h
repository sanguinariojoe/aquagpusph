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

#define vec float4
#define ivec int4
#define uivec uint4
#define matrix float16

/** @def VEC_ZERO
 * @brief Null #vec, i.e. filled with zero components.
 */
#define VEC_ZERO ((float4)(0.f,0.f,0.f,0.f))
/** @def VEC_ONE
 * @brief Ones #vec, i.e. filled with one components (except the w component).
 */
#define VEC_ONE ((float4)(1.f, 1.f, 1.f, 0.f))
/** @def VEC_ALL_ONE
 * @brief Ones #vec, i.e. filled with one components.
 */
#define VEC_ALL_ONE ((float4)(1.f, 1.f, 1.f, 1.f))
/** @def VEC_INFINITY
 * @brief Infinity #vec, i.e. filled with infinity components (except the w
 * component).
 */
#define VEC_INFINITY ((float4)(INFINITY, INFINITY, INFINITY, 0.f))
/** @def VEC_ALL_INFINITY
 * @brief Infinity #vec, i.e. filled with infinity components.
 */
#define VEC_ALL_INFINITY ((float4)(INFINITY, INFINITY, INFINITY, INFINITY))
/** @def VEC_NEG_INFINITY
 * @brief -Infinity #vec, i.e. filled with -infinity components (except the w
 * component).
 */
#define VEC_NEG_INFINITY (-VEC_INFINITY)
/** @def VEC_ALL_NEG_INFINITY
 * @brief -Infinity #vec, i.e. filled with -infinity components.
 */
#define VEC_ALL_NEG_INFINITY (-VEC_ALL_INFINITY)

/** @def MAT_ZERO
 * @brief Null #matrix, i.e. filled with zero components.
 */
#define MAT_ZERO ((float16)(0.f, 0.f, 0.f, 0.f,                                \
						    0.f, 0.f, 0.f, 0.f,                                \
						    0.f, 0.f, 0.f, 0.f,                                \
						    0.f, 0.f, 0.f, 0.f))
/** @def MAT_ONE
 * @brief Ones #matrix, i.e. filled with one components, except the last row and
 * column.
 */
#define MAT_ONE ((float16)(1.f, 1.f, 1.f, 0.f,                                 \
				 		   1.f, 1.f, 1.f, 0.f,                                 \
				 		   1.f, 1.f, 1.f, 0.f,                                 \
						   0.f, 0.f, 0.f, 0.f))
/** @def MAT_ALL_ONE
 * @brief Ones #matrix, i.e. filled with one components.
 */
#define MAT_ALL_ONE ((float16)(1.f, 1.f, 1.f, 1.f,                             \
						       1.f, 1.f, 1.f, 1.f,                             \
						       1.f, 1.f, 1.f, 1.f,                             \
						       1.f, 1.f, 1.f, 1.f))
/** @def MAT_EYE
 * @brief Eye #matrix , except the south-east component, which is filled with a
 * zero.
 *
 * \f$ m_{ii} = 1 \leftrightarrow i \neq 4;
 *     m_{ii} = 0 \leftrightarrow i = 4;
 *     m_{ij} = 0 \leftrightarrow i \neq j \f$
 */
#define MAT_EYE ((float16)(1.f, 0.f, 0.f, 0.f,                                 \
						   0.f, 1.f, 0.f, 0.f,                                 \
						   0.f, 0.f, 1.f, 0.f,                                 \
						   0.f, 0.f, 0.f, 0.f))
/** @def MAT_ALL_EYE
 * @brief Eye #matrix
 *
 * \f$ m_{ii} = 1; m_{ij} = 1 \leftrightarrow i \neq j \f$
 */
#define MAT_ALL_EYE ((float16)(1.f, 0.f, 0.f, 0.f,                             \
				     		   0.f, 1.f, 0.f, 0.f,                             \
					    	   0.f, 0.f, 1.f, 0.f,                             \
						       0.f, 0.f, 0.f, 1.f))

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
#define vec_xyz vec3

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
#define ivec_xyz ivec3

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
#define uivec_xyz uivec3

/** @def XYZ
 * @brief Convenient access to the vector components.
 * 
 * It is useful to be used with #vec_xyz, #ivec_xyz and #uivec_xyz type:
 *   - 2D = .xy
 *   - 3D = .xyz
 */
#define XYZ xyz

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
            for(int ck = -1; ck <= 1; ck++) {                                  \
                const uint c_j = c_i +                                         \
                                 ci +                                          \
                                 cj * n_cells.x +                              \
                                 ck * n_cells.x * n_cells.y;                   \
                uint j = ihoc[c_j];                                            \
                while((j < N) && (icell[j] == c_j)) {

/** @def END_LOOP_OVER_NEIGHS
 * @brief End of the loop over the neighs to compute the interactions.
 * 
 * @see BEGIN_LOOP_OVER_NEIGHS
 */
#define END_LOOP_OVER_NEIGHS()                                                 \
                    j++;                                                       \
                }                                                              \
            }                                                                  \
        }                                                                      \
    }

/** @def MATRIX_DOT
 * @brief Multiply a matrix by a vector (inner product)
 *
 * @note The vector should have 3 components, not 4.
 */
#define MATRIX_DOT(_M, _V)                                                     \
	((float4)(dot(_M.s012, _V),                                                \
              dot(_M.s456, _V),                                                \
              dot(_M.s89A, _V),                                                \
              0.f))

/** @def MATRIX_DOT_ALL
 * @brief Multiply a matrix by a vector (inner product)
 */
#define MATRIX_DOT_ALL(_M, _V)                                                 \
    ((float4)(dot(_M.s0123, _V),                                               \
              dot(_M.s4567, _V),                                               \
              dot(_M.s89AB, _V),                                               \
              dot(_M.sCDEF, _V)))

/** @def MATRIX_MUL
 * @brief Multiply a matrix by a matrix (inner product)
 *
 * @note The last row and column of each matrix will be ignored. To perform a
 * complete inner product use #MATRIX_MUL_ALL
 */
#define MATRIX_MUL(_M1, _M2) ((float16)(                                            \
    dot(_M1.s012, _M2.s048), dot(_M1.s012, _M2.s159), dot(_M1.s012, _M2.s26A), 0.f, \
    dot(_M1.s456, _M2.s048), dot(_M1.s456, _M2.s159), dot(_M1.s456, _M2.s26A), 0.f, \
    dot(_M1.s89A, _M2.s048), dot(_M1.s89A, _M2.s159), dot(_M1.s89A, _M2.s26A), 0.f, \
    0.f, 0.f, 0.f, 0.f))

/** @def MATRIX_MUL_ALL
 * @brief Multiply a matrix by a matrix (inner product)
 *
 * @note For performance purposes, using #MATRIX_MUL instead of this operator is
 * strongly recommended.
 */
#define MATRIX_MUL_ALL(_M1, _M2) ((float16)(                                                                    \
    dot(_M1.s0123, _M2.s048C), dot(_M1.s0123, _M2.s159D), dot(_M1.s0123, _M2.s26AE), dot(_M1.s0123, _M2.s37BF), \
    dot(_M1.s4567, _M2.s048C), dot(_M1.s4567, _M2.s159D), dot(_M1.s4567, _M2.s26AE), dot(_M1.s4567, _M2.s37BF), \
    dot(_M1.s89AB, _M2.s048C), dot(_M1.s89AB, _M2.s159D), dot(_M1.s89AB, _M2.s26AE), dot(_M1.s89AB, _M2.s37BF), \
    dot(_M1.sCDEF, _M2.s048C), dot(_M1.sCDEF, _M2.s159D), dot(_M1.sCDEF, _M2.s26AE), dot(_M1.sCDEF, _M2.s37BF)))

/** @def MATRIX_TRANSPOSE
 * @brief Transpose a matrix
 */
#define TRANSPOSE s048C159D26AE37BF

/** @def DIAG
 * @brief The matrix diagonal (as vector)
 */
#define DIAG s05A

/** @def MATRIX_FROM_DIAG
 * @brief Build up a matrix from the diagonal information (as vector)
 * @note The component w of the vector is ignored (and 0.f is used instead)
 */
#define MATRIX_FROM_DIAG(_V)                                                   \
	((float16)(_V.x,  0.f,  0.f,  0.f,                                         \
                0.f, _V.y,  0.f,  0.f,                                         \
                0.f,  0.f, _V.z,  0.f,                                         \
                0.f,  0.f,  0.f,  0.f))

/** @brief Perform the outer product of two vectors
 *
 * @param v1 Left operand vector (of 3 components)
 * @param v2 Right operand vector (of 3 components)
 */
matrix outer(vec3 v1, vec3 v2)
{
	matrix m = MAT_ZERO;
	m.s012 = v1.x * v2;
	m.s456 = v1.y * v2;
	m.s89A = v1.z * v2;
	return m;
}
