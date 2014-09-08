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

#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
#endif

#ifndef HAVE_3D
    #ifndef NEIGH_CELLS
		/** @def NEIGH_CELLS
         * @brief Number of neigh cells.
         *
         * In 2D cases 9 cells must be computed, while in 3D simulations 27
         * cells must be computed.
		 */ 
        #define NEIGH_CELLS 9
    #endif
#else
    #ifndef NEIGH_CELLS
		/** @def NEIGH_CELLS
         * @brief Number of neigh cells.
         *
         * In 2D cases 9 cells must be computed, while in 3D simulations 27
         * cells must be computed.
		 */ 
        #define NEIGH_CELLS 27
    #endif
#endif

#ifndef uint
	/** @def uint
     * @brief Short alias for unsigned integer type.
	 */ 
    #define uint unsigned int
#endif

/** @brief Set the fields for the boundary elements using the data interpolated
 * by Aqua::CalcServer::Rates.
 * @param ifluid Fluid index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drdt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param press Pressure \f$ p \f$.
 * @param dens Density \f$ \rho \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param gamma Gamma EOS exponent \f$ \gamma \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @see Aqua::CalcServer::Boundary::DeLeffe
 * @see Aqua::CalcServer::Rates
 */
__kernel void Elements(__global int* ifluid,
                       __global int* imove,
                       __global vec* dvdt,
                       __global float* drdt,
                       __global float* press,
                       __global float* dens,
                       __global float* shepard,
                       __constant float* refd,
                       __constant float* gamma,
                       uint N,
                       float cs)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    /** To set the boundary elements pressure and density values two approaches
     * may be considred (just one of them must be activated):
     *   - _ALL_INTERPOLATED_: Proposed by Ferrand, the density and pressure
     *     are interpolated from the neighbour particles.
     *   - _DENSITY_BASED_: Just the density is interpolated, computing the
     *     pressure using the EOS.
     *   - _PRESSURE_BASED_: The pressure value is interpolated, and the
     *     density is computed the EOS.
     *
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
        press[i] = (dvdt[i].x + dvdt[i].y) / shepard_i;
    #elif defined _DENSITY_BASED_
        dens[i] = drdt[i] / shepard_i;
        // Batchelor 1967
        const float rdenf = refd[ifluid[i]];
        const float gammf = gamma[ifluid[i]];
        const float ddenf = dens[i] / rdenf;
        const float prb = cs * cs * rdenf / gammf;
        press[i] = max(prb * (pow(ddenf, gammf) - 1.f),
                       dvdt[i].x / shepard_i)
                   + dvdt[i].y/shepard_i;
    #elif defined _PRESSURE_BASED_
        press[i] = max(0.f, (dvdt[i].x + dvdt[i].y) / shepard_i);
        // Reversed Batchelor 1967
        const float rdenf = refd[ifluid[i]];
        const float gammf = gamma[ifluid[i]];
        const float iprb = gammf / (cs * cs * rdenf);
        dens[i] = rdenf * pow(1.f + iprb * press[i], 1.f / gammf);
    #else
        #error "Unknow boundary elements field computation algorithm"
    #endif
    // Rinitializate output variables
    dvdt[i] = VEC_ZERO;
    drdt[i] = 0.f;

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----

}

/** @brief Performs the boundary effect on the fluid particles.
 * @param ifluid Fluid index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param dens Density \f$ \rho \f$.
 * @param press Pressure \f$ p \f$.
 * @param mass Mass \f$ m \f$.
 * @param viscdyn Dynamic viscosity \f$ \nu \f$ (one per fluid)
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drdt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param lvec Number of cells in each direction
 */
__kernel void Boundary(__global int* ifluid, __global int* imove,
                       __global vec* pos, __global vec* normal,
                       __global vec* v, __global float* dens,
                       __global float* press, __global float* mass,
                       __constant float* viscdyn,
                       __global vec* dvdt, __global float* drdt,
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
        #define _F_ dvdt[j]
        #define _DRDT_ drdt[j]
    #else
        #define _F_ dvdt_l[it]
        #define _DRDT_ drdt_l[it]
        __local vec dvdt_l[LOCAL_MEM_SIZE];
        __local float drdt_l[LOCAL_MEM_SIZE];
        _F_ = dvdt[i];
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
        dvdt[i] = _F_;
        drdt[i] = _DRDT_;
    #endif

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----

}

