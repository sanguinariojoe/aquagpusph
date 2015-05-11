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
 * @brief Fixed ghost particles fields interpolation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "../../types/2D.h"
    #include "../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../types/3D.h"
    #include "../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Fixed ghost particles fields interpolation.
 *
 * The dummy particle will be virtually placed in the mirrored position with
 * respect to the associated one, interpolating the fields in such new position.
 *
 * Density = SSM
 * Pressure = SSM + Hydrostatic contribution
 * Normal velocity component = SSM
 * Tangential velocity component = SSM
 *
 * Mirroring particle id is denominated iref.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param associations Mirroring particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param refd Density of reference \f$ \rho_0 \f$ (one per set of particles)
 * @param icell Cell where each particle is located.
 * @param mirrored_cell Cell where each mirrored particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void main(const __global uint* iset,
                   const __global int* imove,
                   const __global uint *associations,
                   const __global vec* r,
                   const __global vec* normal,
                   const __global float* m,
                   __global float* rho,
                   __global float* p,
                   __global vec* u,
                   __constant float* refd,
                   // Link-list data
                   const __global uint *icell,
                   const __global uint *mirrored_cell,
                   const __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells,
                   float p0,
                   vec g)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -1)
        return;
    const uint iref = associations[i];
    if(iref >= N)
        return;
    
    const vec_xyz r_iref = r[iref].XYZ;
    const vec_xyz n_iref = normal[iref].XYZ;

    const vec_xyz dr_i = dot(r_iref - r[i].XYZ, n_iref) * n_iref;
    const vec_xyz r_i = r[i].XYZ + 2.f * dr_i;

    const uint c_i = mirrored_cell[i];
    const float refd_i = refd[iset[i]];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _RHO_ rho[i]
        #define _P_ p[i]
        #define _U_ u[i].XYZ
    #else
        #define _RHO_ rho_l[it]
        #define _P_ p_l[it]
        #define _U_ u_l[it]
        __local float rho_l[LOCAL_MEM_SIZE];
        __local float p_l[LOCAL_MEM_SIZE];
        __local vec_xyz u_l[LOCAL_MEM_SIZE];
    #endif
    _RHO_ = 0.f;
    _P_ = 0.f;
    _U_ = VEC_ZERO;
    float shepard = 0.f;

    // Loop over neighs
    // ================
    for(int ci = -1; ci <= 1; ci++) {
        for(int cj = -1; cj <= 1; cj++) {
            #ifdef HAVE_3D
            for(int ck = -1; ck <= 1; ck++) {
            #else
            const int ck = 0; {
            #endif
                const uint c_j = c_i +
                                ci +
                                cj * n_cells.x +
                                ck * n_cells.x * n_cells.y;
                uint j = ihoc[c_j];
                while((j < N) && (icell[j] == c_j)) {
                    if(imove[j] <= 0){
                        j++;
                        continue;
                    }
                    const vec_xyz r_ij = r[j].XYZ - r_i;
                    const float q = fast_length(r_ij) / H;
                    if(q >= SUPPORT)
                    {
                        j++;
                        continue;
                    }

                    {
                        #include "Interpolation.hcl"
                    }
                    j++;
                }
            }
        }
    }

    if(shepard < 1.0E-6f){
        _RHO_ = refd_i;
        _P_ = p0;
        _U_ = VEC_ZERO;
    }
    else{
        _RHO_ /= shepard;
        _P_ /= shepard;
        _U_ /= shepard;

        _P_ -= refd_i * 2.f * dot(g, dr_i);
    }

    #ifdef LOCAL_MEM_SIZE
        rho[i] = _RHO_;
        p[i] = _P_;
        u[i].XYZ = _U_;
    #endif
}
