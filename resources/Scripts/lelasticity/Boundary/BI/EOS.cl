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

/** @addtogroup lela
 * @{
 */

/** @file
 * @brief Inverse EOS equation application.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Inverse EOS to get the density from the pressure value.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param p Pressure \f$ p \f$.
 * @param rho Density \f$ \rho \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global float* p,
                    __global float* rho,
                    __constant float* refd,
                    uint N,
                    float cs,
                    float p0)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3){
        return;
    }
    rho[i] = refd[iset[i]] + (p[i] - p0) / (cs * cs);
}

/*
 * @}
 */