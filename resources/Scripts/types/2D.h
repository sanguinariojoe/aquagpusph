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

#ifndef INFINITY
    #define INFINITY FLT_MAX
#endif

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

#define vec float2
#define ivec int2
#define uivec uint2
#define matrix float4

/** @brief Helper function for CONVERT()
 *
 * The helper is required because the preprocessor is only recursively expanding
 * macros if the definition is not affected by # nor ## string operators.
 * Then, this inner function is concatenating the unexpanded words, while
 * CONVERT() is effectively expanding the type name.
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
 * For instance, to convert a vec variable, v, to an ivec variable, you can call
 * CONVERT(ivec, v);
 */
#define CONVERT(TYPE, v) _CONVERT(TYPE)(v)

/** @brief Null #vec, i.e. filled with zero components.
 */
#define VEC_ZERO ((float2)(0.f,0.f))
/** @brief Ones #vec, i.e. filled with one components.
 */
#define VEC_ONE ((float2)(1.f, 1.f))
/** @brief VEC_ONE.
 */
#define VEC_ALL_ONE VEC_ONE
/** @brief Infinity #vec, i.e. filled with infinity components.
 */
#define VEC_INFINITY ((float2)(INFINITY, INFINITY))
/** @brief VEC_INFINITY.
 */
#define VEC_ALL_INFINITY VEC_INFINITY
/** @brief -Infinity #vec, i.e. filled with -infinity components.
 */
#define VEC_NEG_INFINITY (-VEC_INFINITY)
/** @brief VEC_NEG_INFINITY.
 */
#define VEC_ALL_NEG_INFINITY (-VEC_ALL_INFINITY)

/** @brief Null #matrix, i.e. filled with zero components.
 */
#define MAT_ZERO ((float4)(0.f, 0.f,                                           \
                           0.f, 0.f))
/** @brief Ones #matrix, i.e. filled with one components.
 */
#define MAT_ONE ((float4)(1.f, 1.f,                                            \
                          1.f, 1.f))
/** @brief MAT_ONE
 */
#define MAT_ALL_ONE MAT_ONE
/** @brief Eye #matrix
 *
 * \f$ m_{ii} = 1; m_{ij} = 1 \leftrightarrow i \neq j \f$
 */
#define MAT_EYE ((float4)(1.f, 0.f,                                            \
                             0.f, 1.f))
/** @brief MAT_EYE
 */
#define MAT_ALL_EYE MAT_EYE

/** @brief Vector of real components with the minimum number of components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 3 components
 *
 * This type can be used for the local variables to reduce the VGPRs.
 */
#define vec_xyz vec2

/** @brief Vector of integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 3 components
 *
 * This type can be used for the local variables to reduce the VGPRs.
 */
#define ivec_xyz ivec2

/** @brief Vector of unsigned integer components.
 *
 * The number of components depends on weather the 2D version or 3D
 * version is compiled:
 *   - 2D = 2 components
 *   - 3D = 3 components
 *
 * This type can be used for the local variables to reduce the VGPRs.
 */
#define uivec_xyz uivec2

/** @brief Convenient access to the vector components.
 * 
 * It is useful to be used with #vec_xyz, #ivec_xyz and #uivec_xyz type:
 *   - 2D = .xy
 *   - 3D = .xyz
 */
#define XYZ xy

/** @brief Utility to can redefine the cell of the particle to be  computed.
 * 
 * It can be used for mirrrored particles, which are temporary associated to a
 * different cell.
 *
 * @see BEGIN_LOOP_OVER_NEIGHS
 */
#define C_I() const uint c_i = icell[i]

/** @brief Loop over the neighs to compute the interactions.
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

/** @brief End of the loop over the neighs to compute the interactions.
 * 
 * @see BEGIN_LOOP_OVER_NEIGHS
 */
#define END_LOOP_OVER_NEIGHS()                                                 \
                j++;                                                           \
            }                                                                  \
        }                                                                      \
    }

/** @brief Neglect a neighbour particle at the neighbours interaction loop
 * 
 * @see BEGIN_LOOP_OVER_NEIGHS
 */
#define NEGLECT_NEIGH()                                                        \
                    j++; continue;

/** @brief Multiply a matrix by a vector (inner product)
 */
#define MATRIX_DOT(_M, _V)                                                     \
    ((float2)(dot(_M.s01, _V),                                                 \
              dot(_M.s23, _V)))

/** @brief #MATRIX_DOT
 */
#define MATRIX_DOT_ALL MATRIX_DOT

/** @brief Multiply a matrix by a matrix (inner product)
 */
#define MATRIX_MUL(_M1, _M2)                                                   \
    ((float4)(dot(_M1.s01, _M2.s02), dot(_M1.s01, _M2.s13),                    \
              dot(_M1.s23, _M2.s02), dot(_M1.s23, _M2.s13)))

/** @brief #MATRIX_MUL
 */
#define MATRIX_MUL_ALL MATRIX_MUL

/** @brief Transpose a matrix
 */
#define TRANSPOSE s0213

/** @brief The matrix diagonal (as vector)
 */
#define DIAG s03

/** @brief Build up a matrix from the diagonal information (as vector)
 */
#define MATRIX_FROM_DIAG(_V)                                                   \
    ((float4)(_V.x,  0.f,                                                      \
               0.f, _V.y))

/** @brief Trace of the matrix
 *
 * i.e. The sum of the diagonal elements of the matrix.
 */
#define MATRIX_TRACE(_M) (_M.s0 + _M.s3)

/** @brief Perform the outer product of two vectors
 *
 * @param v1 Left operand vector
 * @param v2 Right operand vector
 * @return Outer product matrix
 */
matrix outer(const vec v1, const vec v2)
{
    matrix m;
    m.s01 = v1.x * v2;
    m.s23 = v1.y * v2;
    return m;
}

/** @brief Determinant of a matrix
 *
 * @param m Matrix to invert
 * @return Determinant
 */
float det(const matrix m)
{
    return m.s0 * m.s3 - m.s1 * m.s2;
}

/** @brief Inverse of a matrix
 *
 * @param m Matrix to invert
 * @return Inverse of the matrix
 * @remarks If the input matrix m is singular, nan values matrix will be
 * returned. Consider using #MATRIX_INV pseudo-inverse instead
 */
matrix inv(const matrix m)
{
    const float d = 1.f / det(m);
    return ((matrix)( m.s3, -m.s1,
                     -m.s2,  m.s0)) * d;
}

/** @brief Pseudo-inverse of a matrix
 *
 * The SVD Moore-Penrose method is applied:
 * \f[ A^{\dag} = \left( A^T A \right)^{-1} A^T \f]
 */
#define MATRIX_INV(_M)                                                         \
    MATRIX_MUL(inv(MATRIX_MUL(_M.TRANSPOSE, _M)), _M.TRANSPOSE)
