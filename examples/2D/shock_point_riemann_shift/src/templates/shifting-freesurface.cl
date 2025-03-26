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

__kernel void threshold(const __global int* imove,
                    const __global vec* r,
					__global int* kappa,
					const __global vec* inormal,
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
    if(imove[i] != 1)
        return;

	const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _KAPPA_ kappa[i]
    #else
        #define _KAPPA_ kappa_l[it]
        __local int kappa_l[LOCAL_MEM_SIZE];
        _KAPPA_ = 1;
    #endif

	vec_xyz n_i = inormal[i].XYZ;

    BEGIN_LOOP_OVER_NEIGHS(){
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
		const float norma = dot(n_i, inormal[j].XYZ);

		if(acos(fabs(norma)) > 0.2611f)
        {
            _KAPPA_ = 0;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
		kappa[i] = _KAPPA_;
    #endif
}

__kernel void entry(const __global int* imove,
                    const __global vec* r,
					const __global float* lambda_bi,
					const __global float* lambda,
					const __global vec* tangent,
					const __global int* kappa,
					const __global vec* u,
					__global vec* inormal,
					const __global vec* dshepard,
					const __global unsigned int* frees, 
					const __global unsigned int* region, 
					__global vec* du_shift,
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

	if(region[i] == 0)
	 	return;

	const vec_xyz n_i = inormal[i].XYZ;    

    #ifndef HAVE_3D
	vec tau_i = VEC_ZERO;
	tau_i.x = inormal[i].y;
	tau_i.y = -inormal[i].x;
	#else
	vec tau_i = tangent[i];
	#endif

		if(lambda[i] < 0.4f)
		{
			du_shift[i] = VEC_ZERO;
			return;
		}
		

		if(lambda[i] >= 0.4f) {
		//	if (dot(inormal[i], du_shift[i]) > 0.f){

				const matrix out = outer(inormal[i].XYZ, inormal[i].XYZ);
				const matrix I = MAT_EYE - out;
				const vec R = MATRIX_DOT( I, du_shift[i]);
				du_shift[i] = R;
				// du_shift[i] = kappa[i] * R;	
		//	}

		}   


}
