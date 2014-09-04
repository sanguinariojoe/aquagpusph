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

#ifndef vec
    /** @def vec
     * @brief Vector of real components.
     *
     * The number of components depends on weather the 2D version or 3D
     * version is compiled:
     *   - 2D = 2 components
     *   - 3D = 4 components
     */
    #define vec float2
#endif

#ifndef ivec
    /** @def ivec
     * @brief Vector of integer components.
     *
     * The number of components depends on weather the 2D version or 3D
     * version is compiled:
     *   - 2D = 2 components
     *   - 3D = 4 components
     */
    #define ivec int2
#endif

#ifndef uivec
    /** @def ivec
     * @brief Vector of unsigned integer components.
     *
     * The number of components depends on weather the 2D version or 3D
     * version is compiled:
     *   - 2D = 2 components
     *   - 3D = 4 components
     */
    #define uivec uint2
#endif

#ifndef VEC_ZERO
    /** @def VEC_ZERO
     * @brief Null #vec, i.e. filled with zero components.
     */
    #define VEC_ZERO (float2)(0.f,0.f)
#endif

