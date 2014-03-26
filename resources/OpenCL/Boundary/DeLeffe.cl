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

/* Kernel to use. Take care with support distance (sep), that may vary.
 */
#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
#endif

#ifndef HAVE_3D
    #ifndef NEIGH_CELLS
        /** @def NEIGH_CELLS Number of neighbour cells. In 2D case 8,
         * and the main cells must be computed, but in 3D 27 cells,
         * must be computed.
         */ 
        #define NEIGH_CELLS 9
    #endif
#else
    #ifndef NEIGH_CELLS
        #define NEIGH_CELLS 27
    #endif
#endif

#ifndef M_PI
    #define M_PI 3.14159265359f
#endif
#ifndef iM_PI
    #define iM_PI 0.318309886f
#endif

#ifndef uint
    #define uint unsigned int
#endif

/** Set new values for the boundary elements.
 * @param ifluid Fluid identifiers.
 * @param imove Moving flags.
 * @param f Accelerations.
 * @param drdt Densities variation rate.
 * @param press Pressures.
 * @param dens Densities.
 * @param shepard Shepard terms (0th correction).
 * @param refd Reference density.
 * @param gamma EOS gamma exponent.
 * @param N Total number of vertices.
 * @param cs Speed of sound.
 */
__kernel void Elements(__global int* ifluid, __global int* imove,
                       __global vec* f, __global float* drdt,
                       __global float* press, __global float* dens,
                       __global float* shepard, __constant float* refd,
                       __constant float* gamma,  uint N, float cs)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    /** To set the boundary elements pressure and density values two approaches
     * may be considred (just one of them must be activated):
     *   - _ALL_INTERPOLATED_: Proposed by Ferrand, the density and pressure are interpolated from the neighbour particles.
     *   - _DENSITY_BASED_: Just the density is interpolated, computing the pressure using the EOS.
     *   - _PRESSURE_BASED_: The pressure value is interpolated, and the density is computed the EOS.
     * Pressure field is corrected (in all the cases) by the hydrostatic field.
     * Since density is not corrected, if you are using weakly compressible
     * (where the fluid pressure and density are related) consideer using
     * the _PRESSURE_BASED_ algorithm.
     */
    #define _PRESSURE_BASED_

    // ---- | ------------------------ | ----
    // ---- V ---- Your code here ---- V ----

    if(imove[i]>=0)
        return;

    float shepard_i = shepard[i];
    if(shepard_i < 0.01f)
        shepard_i = 1.f;
    #ifdef _ALL_INTERPOLATED_
        dens[i] = drdt[i] / shepard_i;
        press[i] = (f[i].x + f[i].y) / shepard_i;
    #elif defined _DENSITY_BASED_
        dens[i] = drdt[i] / shepard_i;
        // Batchelor 1967
        const float rdenf = refd[ifluid[i]];
        const float gammf = gamma[ifluid[i]];
        const float ddenf = dens[i] / rdenf;
        const float prb = cs * cs * rdenf / gammf;
        press[i] = max(prb * (pow(ddenf, gammf) - 1.f),
                       f[i].x / shepard_i)
                   + f[i].y/shepard_i;
    #elif defined _PRESSURE_BASED_
        press[i] = max(0.f, (f[i].x + f[i].y) / shepard_i);
        // Reversed Batchelor 1967
        const float rdenf = refd[ifluid[i]];
        const float gammf = gamma[ifluid[i]];
        const float iprb = gammf / (cs * cs * rdenf);
        dens[i] = rdenf * pow(1.f + iprb * press[i], 1.f / gammf);
    #else
        #error "Unknow boundary elements field computation algorithm"
    #endif
    // Rinitializate output variables
    f[i] = VEC_ZERO;
    drdt[i] = 0.f;

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----

}

/** Performs the boundary effect over particles based on the boundary
 * integrals.
 * @param ifluid Fluid identifiers.
 * @param imove Moving flags.
 * @param pos Positions.
 * @param normal Normals.
 * @param v Velocities.
 * @param dens Densities.
 * @param press Pressures.
 * @param mass Masses.
 * @param viscdyn Dynamic viscosity (one per fluid)
 * @param f Accelerations.
 * @param drdt Density variation rates.
 * @param icell Cell where each particle is located.
 * @param ihoc Head particle of chain for each cell.
 * @param N Number of particles.
 * @param lvec Number of cells
 */
__kernel void Boundary(__global int* ifluid, __global int* imove,
                       __global vec* pos, __global vec* normal,
                       __global vec* v, __global float* dens,
                       __global float* press, __global float* mass,
                       __constant float* viscdyn,
                       __global vec* f, __global float* drdt,
                        // Link-list data
                        __global uint *icell, __global uint *ihoc,
                        // Simulation data
                        uint N, uivec lvec)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    // ---- | ------------------------ | ----
    // ---- V ---- Your code here ---- V ----

    const uint c_i = icell[i];
    const vec pos_i = pos[i];
    const vec v_i = v[i];
    const float press_i = press[i];
    const float dens_i = dens[i];
    const float prfac_i = press_i / (dens_i * dens_i);
    #ifdef __NO_SLIP__
        const float viscdyn_i = viscdyn[ifluid[i]];
    #endif

    #ifndef HAVE_3D
        const float conw = 1.f/(h*h);
    #else
        const float conw = 1.f/(h*h*h);
    #endif

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _F_ f[j]
        #define _DRDT_ drdt[j]
    #else
        #define _F_ f_l[it]
        #define _DRDT_ drdt_l[it]
        __local vec f_l[LOCAL_MEM_SIZE];
        __local float drdt_l[LOCAL_MEM_SIZE];
        _F_ = f[i];
        _DRDT_ = drdt[i];
    #endif

    // Loop over neighbour particles
    // =============================
    {
        uint j;
        // Home cell, starting from the next particle
        // ==========================================
        j = i + 1;
        while((j < N) && (icell[j] == c_i) ) {
            #include "DeLeffe.hcl"
            j++;
        }

        // Neighbour cells
        // ===============
        for(uint cell = 1; cell < NEIGH_CELLS; cell++) {
            uint c_j;
            switch(cell) {
                case 0: c_j = c_i + 0; break;
                case 1: c_j = c_i + 1; break;
                case 2: c_j = c_i - 1; break;
                case 3: c_j = c_i + lvec.x; break;
                case 4: c_j = c_i + lvec.x + 1; break;
                case 5: c_j = c_i + lvec.x - 1; break;
                case 6: c_j = c_i - lvec.x; break;
                case 7: c_j = c_i - lvec.x + 1; break;
                case 8: c_j = c_i - lvec.x - 1; break;
                #ifdef HAVE_3D
                    case 9 : c_j = c_i + 0          - lvec.x*lvec.y; break;
                    case 10: c_j = c_i + 1          - lvec.x*lvec.y; break;
                    case 11: c_j = c_i - 1          - lvec.x*lvec.y; break;
                    case 12: c_j = c_i + lvec.x     - lvec.x*lvec.y; break;
                    case 13: c_j = c_i + lvec.x + 1 - lvec.x*lvec.y; break;
                    case 14: c_j = c_i + lvec.x - 1 - lvec.x*lvec.y; break;
                    case 15: c_j = c_i - lvec.x     - lvec.x*lvec.y; break;
                    case 16: c_j = c_i - lvec.x + 1 - lvec.x*lvec.y; break;
                    case 17: c_j = c_i - lvec.x - 1 - lvec.x*lvec.y; break;

                    case 18: c_j = c_i + 0          + lvec.x*lvec.y; break;
                    case 19: c_j = c_i + 1          + lvec.x*lvec.y; break;
                    case 20: c_j = c_i - 1          + lvec.x*lvec.y; break;
                    case 21: c_j = c_i + lvec.x     + lvec.x*lvec.y; break;
                    case 22: c_j = c_i + lvec.x + 1 + lvec.x*lvec.y; break;
                    case 23: c_j = c_i + lvec.x - 1 + lvec.x*lvec.y; break;
                    case 24: c_j = c_i - lvec.x     + lvec.x*lvec.y; break;
                    case 25: c_j = c_i - lvec.x + 1 + lvec.x*lvec.y; break;
                    case 26: c_j = c_i - lvec.x - 1 + lvec.x*lvec.y; break;
                #endif
            }

            j = ihoc[c_j];
            while((j < N) && (icell[j] == c_j)) {
                #include "DeLeffe.hcl"
                j++;
            }            
        }
        // Home cell, starting from the head of chain
        // ==========================================
        j = ihoc[c_i];
        while(j < i) {
            #include "DeLeffe.hcl"
            j++;
        }
    }

    #ifdef LOCAL_MEM_SIZE
        f[i] = _F_;
        drdt[i] = _DRDT_;
    #endif

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----

}

