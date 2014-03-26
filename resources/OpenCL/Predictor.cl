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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** Quasi-second order time integration predictor stage.
 * @param imove Fix particles flag.
 * @param ifluid Fluid identifiers.
 * @param pos Positions.
 * @param v Velocities.
 * @param f Accelerations.
 * @param dens Densities.
 * @param mass Masses.
 * @param drdt Density varaition rates.
 * @param posin Positions from the corrector stage.
 * @param vin Velocities from the corrector stage.
 * @param fin Accelerations from the corrector stage.
 * @param densin Densities from the corrector stage.
 * @param massin Masses from the corrector stage.
 * @param drdtin Density varaition rates from the corrector stage.
 * @param press Pressure from the corrector stage.
 * @param refd Reference density.
 * @param gamma Gamma.
 * @param N Number of particles.
 * @param t Simulation time.
 * @param dt Time step.
 * @param cs Sound speed.
 * @param g Gravity force.
 * @param minDens Minimum tolerated density value.
 * @param maxDens Maximum tolerated density value.
 */
__kernel void Predictor(__global int* imove, __global int* ifluid,
                        __global vec* pos, __global vec* v, __global vec* f,
                        __global float* dens, __global float* mass,
                        __global float* drdt, __global vec* posin,
                        __global vec* vin, __global vec* fin,
                        __global float* densin, __global float* massin,
                        __global float* drdtin, __global float* press, 
                        __constant float* refd, __constant float* gamma,
                        unsigned int N, float t, float dt, float cs, vec g,
                        float minDens, float maxDens)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Stabilization time
	if(t < 0.f){
		if(imove[i] > 0){
			vec dr = 0.5f * dt * dt * (fin[i] + g);
			posin[i] += dr;
		}
		dt = 0.f;
	}
	// Time step modified by imove flag
	float DT = dt;
	if(imove[i] <= 0)
		DT = 0.f;

	// Predictor step for the fluid and walls
	v[i] = vin[i] + DT * (fin[i] + g);
	// mass[i]    = massin[i]*(1.f + dt*drdtin[i]/densin[i]);
	#if __BOUNDARY__ == 1
		// Continuity equation must be solved for the fixed particles too
		dens[i] = densin[i] + dt * drdtin[i];
	#else
		dens[i] = densin[i] + DT * drdtin[i];
	#endif
	if(dens[i] < minDens) dens[i] = minDens;
	if(dens[i] > maxDens) dens[i] = maxDens;
	pos[i] = posin[i] + DT * vin[i] + 0.5f * DT * DT * (fin[i] + g);
	// Batchelor 1967
	{
		const float ddenf = dens[i] / refd[ifluid[i]];
		const float prb = cs * cs * refd[ifluid[i]] / gamma[ifluid[i]];
		press[i] = prb * (pow(ddenf, gamma[ifluid[i]]) - 1.f);
	}
	// Output variables reinitialization
	f[i] = VEC_ZERO;
	drdt[i] = 0.f;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
