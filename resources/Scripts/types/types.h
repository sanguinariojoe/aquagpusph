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
#define C_I() const uint c_i = icell[i]

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
    const __global uint * icell,                                               \
    const __global uint * ihoc,                                                \
    uivec4 n_cells
