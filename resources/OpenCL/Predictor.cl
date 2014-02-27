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

#ifndef M_PI
	#define M_PI 3,14159265359
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

/** Quasi-second order time integration predictor stage.
 * @param imove Fix particles flag.
 * @param ifluid Fluid identifier.
 * @param pos Position of particles.
 * @param v Velocity of particles.
 * @param f Forces over particles.
 * @param dens Density of particles.
 * @param mass Mass of particles.
 * @param drdt Density evolution of particles.
 * @param posin Position of particles.
 * @param vin Velocity of particles.
 * @param fin Forces over particles.
 * @param densin Density of particles.
 * @param massin Mass of particles.
 * @param drdtin Density evolution of particles.
 * @param press Pressure.
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
__kernel void Predictor(_g int* imove, _g int* ifluid, _g vec* pos, _g vec* v,
                        _g vec* f, _g float* dens, _g float* mass,
                        _g float* drdt, _g vec* posin,
                        _g vec* vin, _g vec* fin, _g float* densin,
                        _g float* massin, _g float* drdtin,
                        _g float* press, _c float* refd, _c float* gamma,
                        unsigned int N, float t, float dt, float cs, vec g,
                        float minDens, float maxDens)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	float ddenf, prb;
	// Stabilization time
	if(t < 0.f){
		if(imove[i] > 0){
			vec dr    = 0.5f*dt*dt*(fin[i] + g);
			posin[i] += dr;
		}
		dt = 0.f;
	}
	// Time step modified by imove flag
	float DT = dt;
	if(imove[i] <= 0)
		DT = 0.f;

	// Predictor step for the fluid and walls
	v[i]       = vin[i] + DT*( fin[i] + g );
	// mass[i]    = massin[i]*(1.f + dt*drdtin[i]/densin[i]);
	#if __BOUNDARY__ == 1
		// Continuity equation must be solved for fixed particles too
		dens[i]    = densin[i] + dt*drdtin[i];
	#else
		dens[i]    = densin[i] + DT*drdtin[i];
	#endif
	if(dens[i] < minDens) dens[i] = minDens;
	if(dens[i] > maxDens) dens[i] = maxDens;
	pos[i]     = posin[i] + DT*vin[i] + 0.5f*DT*DT*(fin[i] + g);
	// Batchelor 1967
	{
		ddenf    = dens[i]/refd[ifluid[i]];
		// Pressure term
		prb      = cs*cs*refd[ifluid[i]]/gamma[ifluid[i]];
		press[i] = prb*(pow(ddenf,gamma[ifluid[i]]) - 1.f);
	}
	// Reinitializating of output variables
	f[i]     = VEC_ZERO;
	drdt[i]  = 0.f;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
