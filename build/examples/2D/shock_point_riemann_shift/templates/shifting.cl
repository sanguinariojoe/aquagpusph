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

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "../../../resources/Scripts/types/types.h"
#include "../../../resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Gradient of shepard factor computation.
 *
 * \f[ \nabla \gamma(\mathbf{x}) = \int_{\Omega}
 *     \nabla W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f]
 *
 * The gradient of the shepard renormalization factor is applied for identifying
 * the particles located at a free surface (the ones with a value of the gradient
 * different to zero:
 *
 * In the gradient of the shepard factor computation the fluid extension
 * particles are not taken into account.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param dshepard Gradient of the Shepard term.
 * \f$ \nabla \gamma(\mathbf{x}) = \int_{\Omega}
 *    \nabla W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */

__kernel void relU(__global float* relu,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* u,
                    vec g,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells)

{
    // find position in global arrays
	const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

	const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;

// Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _RELU_ relu[i]
    #else
        #define _RELU_ relu_l[it]
        __local float relu_l[LOCAL_MEM_SIZE];
        _RELU_ = 0.f;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j){
            j++;
            continue;
        }
        if(imove[j] != 1){
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
        {
			_RELU_ = max(_RELU_, fabs(length((u[j].XYZ - u[i].XYZ))));
		}
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        relu[i] = _RELU_;
    #endif
}


__kernel void drfluid(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* du_shift,
					const __global vec* dshepard,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
					float dr,
					float h,
					float cs,
					float Uref,
					float courant,
					float dt,
                    uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _DUSHIFT_ du_shift[i].XYZ
    #else
        #define _DUSHIFT_ du_shift_l[it]
        __local vec_xyz du_shift_l[LOCAL_MEM_SIZE];
        _DUSHIFT_ = VEC_ZERO.XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != 1){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
		
		const float ctant = -Uref * (2.f * h);
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
			const float f_ij = kernelF(q) * CONF * m[j];

		//	_DUSHIFT_.XYZ += ctant * r_ij * f_ij / rho[j];
			_DUSHIFT_ += ctant * ( 1.f + 0.2f * pow(( kernelW(q) ) , 4.f ) ) * r_ij * f_ij / rho[j];
		//	_DUSHIFT_.XYZ += ctant * ( 1.f + 0.2f * pow(( kernelW(q) / kernelW(0.f) ) , 4.f ) ) * r_ij * f_ij / rho[j];

        //    _DUSHIFT_.XYZ -= ctant * ( 1.f + 0.2f * pow(( kernelW(q) / kernelW(1.f/ SUPPORT) ) , 4.f ) ) * r_ij * kernelF(q) * H * CONF * m[j] / rho[j];
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
		du_shift[i].XYZ = _DUSHIFT_ ;
    #endif
}

__kernel void drbound(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* du_shift,
					const __global vec* dshepard,
					const __global float* shepard,
					const __global vec* normal,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
					float dr,
					float h,
					float cs,
					float Uref,
					float courant,
					float dt,
                    uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _DUSHIFT_ du_shift[i].XYZ
    #else
        #define _DUSHIFT_ du_shift_l[it]
        __local vec_xyz du_shift_l[LOCAL_MEM_SIZE];
        _DUSHIFT_ = du_shift[i].XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != -3){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
		
		const float ctant = -Uref * (2.f * h); // Umax = Ma * cs = 0.1 * cs
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
			const float w_ij = kernelW(q) * CONW * m[j];
			const vec_xyz n_j = normal[j].XYZ;
		 //   _DUSHIFT_.XYZ += ctant * w_ij * normal[j];
		    _DUSHIFT_ += ctant * ( 1.f + 0.2f * pow(( kernelW(q) ) , 4.f ) ) * w_ij * n_j;
		 //   _DUSHIFT_.XYZ += ctant * ( 1.f + 0.2f * pow(( kernelW(q) / kernelW(0.f) ) , 4.f ) ) * w_ij * normal[j];
		 //   _DUSHIFT_.XYZ += ctant * ( 1.f + 0.2f * pow(( kernelW(q) / kernelW(dr/H) ) , 4.f ) ) * w_ij * normal[j];
         //   _DUSHIFT_.XYZ -= ctant * ( 1.f + 0.2f * pow(( kernelW(q) / kernelW(dr/h) ) , 4.f ) ) * kernelW(q) * H * CONW * m[j] * normal[j];
        }
    }END_LOOP_OVER_NEIGHS()

	const float shepard_i = shepard[i];

    #ifdef LOCAL_MEM_SIZE
		du_shift[i].XYZ = _DUSHIFT_ / shepard_i;
    #endif
}

__kernel void du(__global int* imove,
                    __global unsigned int* iset,
                    __global vec* du_shift,
					float dr,
					float Uref,
					float cs,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != 1)
        return;

    const float norm_du = length(du_shift[i]);
	const float M = min(0.25f * Uref, norm_du);
	const vec_xyz du_shift_i = du_shift[i].XYZ;

	du_shift[i].XYZ = M * fast_normalize(du_shift_i);
}

__kernel void correct(__global int* imove,
                    __global unsigned int* iset,
                    __global vec* r,
                    __global vec* r_in,
                    __global vec* u,
                    __global vec* u_in,
					const __global vec* du_shift,
					const float dt,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != 1)
        return;

	float DT = dt;

   // r_in[i] += DT * du_shift[i];
    // r[i] += DT * du_shift[i];
   // u[i] += du_shift[i];
}

