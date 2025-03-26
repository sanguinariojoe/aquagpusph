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
__kernel void lambdabi(const __global int* imove,
                   const __global matrix* renorm,
                   __global float* lambda_bi,
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
	if (imove[i] != 1)
		return;
    #ifndef HAVE_3D

	float a = renorm[i].s0;
	float b = renorm[i].s1;
	float c = renorm[i].s2;
	float d = renorm[i].s3;

	const float lambdabi1 = 1.f / ( 0.5f * (a + d) + 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) );
	const float lambdabi2 = 1.f / ( 0.5f * (a + d) - 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) );
	
    lambda_bi[i] = fmin(lambdabi1,lambdabi2);

	#else

	const float c1 = MATRIX_TRACE(renorm[i]) * MATRIX_TRACE(renorm[i]);
	const float16 c2 = MATRIX_MUL(renorm[i], renorm[i]);

	const float a = - MATRIX_TRACE(renorm[i]) ;
	const float b =   0.5f * (c1 - MATRIX_TRACE(c2));
	const float c = - det(renorm[i]);

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

	lambda_bi[i] = min(inter, 1.f/x3);
		}
		else {
 			const float x1 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x2 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x3 = - (4.f*p*p)/(9.f*q) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda_bi[i] = min(inter, 1.f/x3); 
		}
	}
	else if (delta > 0.00f){
		const float xi = -q/2.f - sqrt(delta);
		const float yi = -q/2.f + sqrt(delta);

		const float x1 = copysign(pow( fabs(yi) , 1.f/3.f ), yi) + copysign(pow( fabs(xi) , 1.f/3.f ), xi) - a/3.f;

	lambda_bi[i] = 1.f/x1;
	}
	else {
		const float x1 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*0*pi)/3.f) - a/3.f;
		const float x2 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*1*pi)/3.f) - a/3.f;
		const float x3 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*2*pi)/3.f) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda_bi[i] = min(inter, 1.f/x3);

	}

    #endif
}
