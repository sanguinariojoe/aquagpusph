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
 * @brief Boundary integral term computation.
 * (See Aqua::CalcServer::Boundary::DeLeffe for details)
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Shepard factor computation.
 *
 * \f[ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f]
 *
 * The shepard renormalization factor is applied for several purposes:
 *   - To interpolate values
 *   - To recover the consistency with the Boundary Integrals formulation
 *   - Debugging
 *
 * In the shepard factor computation the fluid extension particles are not taken
 * into account.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param tangent Tangent \f$ \mathbf{t} \f$.
 * @param binormal Binormal \f$ \mathbf{b} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void compute(const __global int* imove,
                      const __global vec* r,
                      const __global vec* normal,
                      const __global vec* tangent,
                      const __global vec* binormal,
                      const __global float* m,
                      __global float* shepard,
                      // Link-list data
                      const __global uint *icell,
                      const __global uint *ihoc,
                      // Simulation data
                      uint N,
                      uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] < -3) || (imove[i] > 1))
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _SHEPARD_ shepard[i]
    #else
        #define _SHEPARD_ shepard_l[it]
        __local float shepard_l[LOCAL_MEM_SIZE];
    #endif

    _SHEPARD_ = 1.f;
    bool self_added = false;

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != -3){
            j++;
            continue;
        }

        if(i == j){
            // Boundary element trying to interact with itself
            if(!self_added){
                self_added = true;
                _SHEPARD_ -= 0.5f;
            }
            j++;
            continue;
        }
        
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        const vec_xyz n_j = normal[j].XYZ;
        const vec_xyz t_j = tangent[j].XYZ;
        const vec_xyz b_j = binormal[j].XYZ;
        const float r_n = dot(r_ij, n_j);
        const float area_j = m[j];

        // Special case of particles/sensors lying at the boundary
        if((r_n > -1E-8f * H) &&
           (r_n < 1E-8f * H))
        {
            if(!self_added) {
                const float r_t = fabs(dot(r_ij, t_j));
                const float r_b = fabs(dot(r_ij, b_j));
                const float dr = 0.55f * pow(area_j, 1.f / (DIMS - 1.f));
                if ((r_t <= dr) && (r_b <= dr)) {
                    self_added = true;
                    _SHEPARD_ -= 0.5f;
                }
            }
            j++;
            continue;
        }

        {
            const float r_t = fabs(dot(r_ij, t_j));
            const float r_b = fabs(dot(r_ij, b_j));
            _SHEPARD_ += r_n * CONW * kernelS_P(q) * area_j +
                         kernelS_D(fabs(r_n), r_t, r_b, area_j);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        shepard[i] = _SHEPARD_;
    #endif
}


/** @brief Renormalize the differential operators.
 *
 * The main drawback of the boundary integrals formulation is the requirement
 * of the renormalization of the computed differentiqal operators, which is
 * destroying several conservation properties.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @see Boundary/BI/Interactions.cl
 */
__kernel void apply(const __global int* imove,
                    const __global float* shepard,
                    __global vec* grad_p,
                    __global vec* lap_u,
                    __global float* div_u,
                    uint N,
                    float cs)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != 1){
        return;
    }

    float shepard_i = shepard[i];

    grad_p[i] /= shepard_i;
    lap_u[i] /= shepard_i;
    // Following line was causing problems at the free surface
    div_u[i] /= shepard_i;
}
