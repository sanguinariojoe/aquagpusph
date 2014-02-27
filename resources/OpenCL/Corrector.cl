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
	#define M_PI 3.14159265359
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

/** Convinient velocity clamping correction, that ensures that
 * a particle never can move more than 0.1h in one time step.
 * Requires that minimum time step was set greather than 0
 * seconds, and of course, velocity clamping was set true.
 * @param imove Fix particles flag. Fixed particles must move
 * as required by solid motions, so their velocity can't be clamped.
 * @param v Velocity of particles.
 * @param f Forces over particles.
 * @param hp Kernel height of particles.
 * @param fin Forces over particles.
 * @param mindt Minimum time step.
 * @param grav Gravity acceleration.
 * @param N Number of particles.
 */
__kernel void ClampVel(_g int* imove, _g vec* v, _g vec* f, _g vec* fin,
                        float mindt, vec grav, unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;
	if(imove[i] <= 0)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Compute the distance moved if dt is equal to mindt
	vec r = mindt*v[i] + mindt*mindt*(f[i] + 0.5f*(grav - fin[i]));
	// Ensure that the particle is not moving too much
	float dr2   = dot(r,r);
	float maxdr = 0.1f*h;
	if( dr2 > maxdr*maxdr){
		// Correct the displacement vector in order to fit to 0.1h
		r   *= maxdr / sqrt(dr2);
		// Get the acceleration that achieves this displacement
		f[i] = (r - mindt*v[i])/(mindt*mindt) + 0.5f*(fin[i] - grav);
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

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
 * @param sigma Viscosity time step term.
 * @param N Number of particles.
 * @param t Simulation time.
 * @param dt Time step.
 */
__kernel void Corrector(_g int* imove, _g vec* pos, _g vec* v, _g vec* f,
                        _g float* dens, _g float* mass, _g float* drdt,
                        _g vec* posin, _g vec* vin, _g vec* fin,
                        _g float* densin, _g float* massin, _g float* drdtin,
                        _g float* sigma,
                        unsigned int N, float t, float dt)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Stabilization time
	if(t < 0.f)
		dt = 0.f;
	float DT, HDT;
	DT = dt;
	if(imove[i] <= 0)
	        DT = 0.f;
	HDT = 0.5f*DT;

	// Corrector step for the fluid
	v[i]      = v[i] + HDT*(f[i] - fin[i]);
	// mass[i]   = mass[i] * (1.f +  hdt*(drdt[i] - drdtin[i]) / dens[i]);
	#if __BOUNDARY__ == 1
		// Continuity equation must be solved for fixed particles
		dens[i]   = dens[i] + 0.5f*dt*(drdt[i] - drdtin[i]);
	#else
		dens[i]   = dens[i] + HDT*(drdt[i] - drdtin[i]);
	#endif
	/* Calculate initial positions, mass and thermal energy
	 * for the next step.
	 */
	posin[i]  = pos[i];
	vin[i]    = v[i];
	massin[i] = mass[i];
	densin[i] = dens[i];
	fin[i]    = f[i];
	drdtin[i] = drdt[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
