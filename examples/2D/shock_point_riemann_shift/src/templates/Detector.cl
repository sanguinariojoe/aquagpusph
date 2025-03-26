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

/** @addtogroup basic
 * @{
 */

/** @file
 * @brief delta-SPH methods, including the correction terms
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "../../../resources/Scripts/types/types.h"
#include "../../../resources/Scripts/KernelFunctions/Kernel.h"

/** @brief MLS based correction term.
 *
 * Here the MLS renormalization is applied to the correction term
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param mls Kernel MLS transformation matrix \f$ L \f$.
 * @param lap_p_corr Correction term for the Morris Laplacian formula.
 * @param N Number of particles.
 */
__kernel void lambda(const __global int* imove,
                   const __global matrix* mls_fluid,
                   __global float* lambda,
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
	if (imove[i] != 1)
		return;
    #ifndef HAVE_3D

	float a = mls_fluid[i].s0;
	float b = mls_fluid[i].s1;
	float c = mls_fluid[i].s2;
	float d = mls_fluid[i].s3;

	const float lambda1 = 1.f / ( ( 0.5f * (a + d) + 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) ) );
	const float lambda2 = 1.f / ( ( 0.5f * (a + d) - 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) ) );
	
    lambda[i] = fmin(lambda1,lambda2);

	#else

	const float c1 = MATRIX_TRACE(mls_fluid[i]) * MATRIX_TRACE(mls_fluid[i]);
	const float16 c2 = MATRIX_MUL(mls_fluid[i], mls_fluid[i]);

	const float a = - MATRIX_TRACE(mls_fluid[i]) ;
	const float b =   0.5f * (c1 - MATRIX_TRACE(c2));
	const float c = - det(mls_fluid[i]);

	const float p = (3.f * b - a*a) / 3.f;
	const float q = (2.f*a*a*a-9.f*a*b+27.f*c)/27.f;
	const float delta = pow( q/2.f , 2.f) + pow( p/3.f , 3.f);
	const float phi = acos( (-q/2.f) / sqrt(-pow(p/3.f , 3.f) ) );
	const float pi = 3.14159f;

	if(delta == 0.00f){
		if((p == 0.f)  && (q == 0.f)){
			const float x1 = - a/3.f;
			const float x2 = - a/3.f;
			const float x3 = - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda[i] = fmin(inter, 1.f/x3);
		}
		else {
 			const float x1 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x2 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x3 = - (4.f*p*p)/(9.f*q) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda[i] = fmin(inter, 1.f/x3); 
		}
	}
	else if (delta > 0.00f){
		const float xi = -q/2.f - sqrt(delta);
		const float yi = -q/2.f + sqrt(delta);

		const float x1 = copysign(pow( fabs(yi) , 1.f/3.f ), yi) + copysign(pow( fabs(xi) , 1.f/3.f ), xi) - a/3.f;

	lambda[i] = 1.f/x1;
	}
	else {
		const float x1 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*0*pi)/3.f) - a/3.f;
		const float x2 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*1*pi)/3.f) - a/3.f;
		const float x3 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*2*pi)/3.f) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda[i] = fmin(inter, 1.f/x3);

	}

    #endif
}

__kernel void first(const __global int* imove,          
					__global unsigned int* frees,
                   const __global float* lambda_bi,
					const __global vec* dshepard,
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
	if (imove[i] != 1)
		return;

	if (lambda_bi[i] < 0.2f){
		frees[i] = 1;
		return;
	}
}

__kernel void normal(const __global int* imove,
                    const __global matrix* mls_fluid,
                    const __global matrix* renorm,
                    const __global float* lambda,
                    const __global float* lambda_bi,
					const __global float* m,
					const __global float* rho,
					const __global vec* r,
					const __global vec* dshepard,
                    const __global unsigned int* counter_mls_in,
					__global vec* inormal,
					__global int* na,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells)
{
    const uint i = get_global_id(0);
//    const uint it = get_local_id(0);
    if(i >= N)
        return;
	if(imove[i] != 1){
		inormal[i].XYZ = VEC_ZERO.XYZ;
		return;
	}

	const vec_xyz r_i = r[i].XYZ;
/*
	// Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _INORMAL_ inormal[i]
    #else
        #define _INORMAL_ inormal_l[it]
        __local vec_xyz inormal_l[LOCAL_MEM_SIZE];
		_INORMAL_ = VEC_ZERO;
    #endif

*/	
/*
	const vec nnormal = - MATRIX_DOT(renorm[i], dshepard[i]);
		if(length(nnormal) > 0.01f){
			na[i] = 1;
		}
		else{
			na[i] = 0;
		}

	if(na[i] == 1){
		inormal[i] = normalize(nnormal);
	}
	else{
		inormal[i] = VEC_ZERO;
		}
*/
	inormal[i].XYZ = VEC_ZERO.XYZ;


	if(counter_mls_in[i] > 5){
		inormal[i] = - MATRIX_DOT(renorm[i], dshepard[i]);
	}

	inormal[i] = fast_normalize(inormal[i]); 

/*
	if(lambda_bi[i] > 0.75f){
		inormal[i].XYZ = VEC_ZERO.XYZ;
		return;
	}


	BEGIN_LOOP_OVER_NEIGHS(){
       if((imove[j] != 1)){
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
			_INORMAL_ -= (lambda[j] - lambda[i]) * MATRIX_DOT(renorm[i], dshepard[i]) ;
		
        }
    }END_LOOP_OVER_NEIGHS()


	const vec_xyz norm = normalize(_INORMAL_.XYZ) ; 

    #ifdef LOCAL_MEM_SIZE
	    inormal[i].XYZ = norm;
    #endif

*/
}


__kernel void second(const __global int* imove,
                    const __global vec* r,
                    const __global float* m,
                    const __global vec* inormal,
					const __global vec* normal,
					const __global vec* tangent,
					const __global float* lambda,
					const __global float* lambda_bi,
					const __global float* shepard,
					const __global vec* dshepard,
					float h,
					float dr,
					__global unsigned int* frees,
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

    if(frees[i] == 1)
        return;

    const vec_xyz r_i = r[i].XYZ;
	const vec_xyz n_i = inormal[i].XYZ;

	const vec_xyz dt_i = h * n_i;
	const vec_xyz t_i = r_i + dt_i;


    #ifndef HAVE_3D
    vec_xyz tau_i = VEC_ZERO;
	tau_i.x = -n_i.y;
	tau_i.y = n_i.x;
	#else
    vec_xyz tau_i = cross(r_i, n_i);
	#endif
 
    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _FREES_ frees[i]
    #else
        #define _FREES_ frees_l[it]
        __local unsigned int frees_l[LOCAL_MEM_SIZE];
    #endif

    _FREES_ = 1;

	if(lambda_bi[i] > 0.75f){
	frees[i] = 0;
	return;
	}

    BEGIN_LOOP_OVER_NEIGHS(){

        const vec_xyz r_ij = r[j].XYZ - r_i;
	const vec_xyz r_jT = r[j].XYZ - t_i;
        const float q = length(r_ij) / H;

        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
	if(imove[j] != 1){
	    j++;
	    continue;
	}

    #ifndef HAVE_3D
        {
			if((length(r_ij) >= sqrt(2.0) * h) && (length(r_jT) < h)){
		//	if(length(r_jT) < 0.99 * h){ 
            	 _FREES_ = 0;
			}
			if((length(r_ij) < sqrt(2.0) * h) && (fabs(dot(n_i, r_jT)) + fabs(dot(tau_i, r_jT)) < h) ){
		//	if(fabs(dot(n_i, r_jT)) + fabs(dot(tau_i, r_jT)) < 0.99 * h){
				 _FREES_ = 0;
			}
        }
    #else
        {
			if((length(r_ij) >= sqrt(2.f) * h) && (length(r_jT) < h)){ 
		//	if(length(r_jT) < 0.99 * h){ 
            	 _FREES_ = 0;
			}
			if((length(r_ij) < sqrt(2.f) * h) && (acos( dot(inormal[i].XYZ, r_ij) / length(r_ij) ) < 0.25f * 3.14159265359f )){
		//	if((acos( dot(inormal[i].XYZ, r_ij) / length(r_ij) ) < 0.4f * 3.14159265359f )){
				 _FREES_ = 0;
			}


        }
    #endif
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        frees[i] = _FREES_;
    #endif
}

__kernel void region(const __global int* imove,
                    const __global vec* r,
                    const __global float* m,
                    const __global vec* inormal,
					const __global float* lambda,
					const __global unsigned int* frees,
					__global unsigned int* region,
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
        #define _REGION_ region[i]
    #else
        #define _REGION_ region_l[it]
        __local unsigned int region_l[LOCAL_MEM_SIZE];
    #endif
    _REGION_ = 0;


    BEGIN_LOOP_OVER_NEIGHS(){
        if((imove[j] != 1)){
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
			if( (frees[j] == 1) ){
			_REGION_ = 1;
			break;
			}
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        region[i] = _REGION_;
    #endif
}