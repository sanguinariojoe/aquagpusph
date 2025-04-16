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

/**
 * \addtogroup ideal_gas
 * @{
 */

/** @file
 * @brief Sort the internal energy by the cell indexes
 *
 * This is an extension of resources/Scripts/basic/Sort.cl
 */

#include "resources/Scripts/types/types.h"

/** @brief Sort the internal energy.
 *
 * @param eint_in Unsorted internal energy \f$ e \f$.
 * @param eint Sorted internal energy \f$ e \f$.
 * @param deintdt Unsorted internal energy rate of change
 * \f$ \frac{d e}{d t} \f$.
 * @param deintdt_in Sorted internal energy rate of change
 * \f$ \frac{d e}{d t} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void entry(__global float* z,
                    __global float* y_H2,
                    const __global float* dy_H2dt,
                    __global float* y_O2,
                    const __global float* dy_O2dt,
                    __global float* y_N2,
                    const __global float* dy_N2dt,
                    __global float* y_H2O,
                    const __global float* dy_H2Odt,
                    const __global float* z_in,
                    const __global float* y_H2_in,
                    __global float* dy_H2dt_in,
                    const __global float* y_O2_in,
                    __global float* dy_O2dt_in,
                    const __global float* y_N2_in,
                    __global float* dy_N2dt_in,
                    const __global float* y_H2O_in,
                    __global float* dy_H2Odt_in,
                    const __global usize *id_sorted,
                    usize N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    const usize i_out = id_sorted[i];

    z[i_out] = z_in[i];    

    y_H2[i_out] = y_H2_in[i];
    dy_H2dt_in[i_out] = dy_H2dt[i];

    y_O2[i_out] = y_O2_in[i];
    dy_O2dt_in[i_out] = dy_O2dt[i];

    y_N2[i_out] = y_N2_in[i];
    dy_N2dt_in[i_out] = dy_N2dt[i];

    y_H2O[i_out] = y_H2O_in[i];
    dy_H2Odt_in[i_out] = dy_H2Odt[i];
}

/*
 * @}
 */
