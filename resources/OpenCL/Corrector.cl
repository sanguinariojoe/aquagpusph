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
 * @brief Leap-frog time integration scheme corrector stage.
 * (See Aqua::CalcServer::Corrector for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Velocity clamping correction.
 *
 * This tool asserts that a particle never can move more than a distance
 * \f$ 0.1 h \f$ in one time step.
 * Requires that a minimum time step was set greater than 0 seconds, and the
 * velocity clamping has been activated.
 * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters::dt_min
 * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters::velocity_clamp
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param v Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dvdt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param dvdtin Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param mindt Minimum time step \f$ \Delta t \f$.
 * @param grav Gravity acceleration \f$ \mathbf{g} \f$.
 * @param N Number of particles.
 */
__kernel void ClampVel(__global int* imove,
                       __global vec* v,
                       __global vec* dvdt,
                       __global vec* dvdtin,
                       float mindt,
                       vec grav,
                       unsigned int N)
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
    const float dt = mindt;
	vec r = dt * v[i] + dt * dt * (dvdt[i] + 0.5f * (grav - dvdtin[i]));
	// Ensure that the particle is not moving too far
	const float dr2 = dot(r, r);
	const float maxdr = 0.1f * h;
	if( dr2 > maxdr * maxdr){
		// Correct the displacement vector in order to fit to 0.1h
		r *= maxdr / sqrt(dr2);
		// Get the acceleration that achieves this displacement
		dvdt[i] = (r - dt * v[i]) / (dt * dt) + 0.5f * (dvdtin[i] - grav);
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

/** @brief Leap-frog time integration scheme corrector stage.
 *
 * Time integration is based in the following quasi-second order
 * Predictor-Corrector integration scheme:
 *   - \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t \left(
        \mathbf{g} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
     \right)
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n + 1/2} -
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n - 1/2}
     \right)
     \f$
 *   - \f$ \mathbf{r}_{n+1} = \mathbf{r}_{n} + \Delta t \, \mathbf{u}_{n}
     + \frac{\Delta t^2}{2} \left(
        \mathbf{g} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
     \right)
     \f$
 *   - \f$ \rho_{n+1} = \rho_{n} + \Delta t
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n+1/2}
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n + 1/2} -
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n - 1/2}
     \right)
     \f$
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param pos Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param v Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dvdt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param dens Density \f$ \rho_{n+1/2} \f$.
 * @param mass Mass \f$ m_{n+1/2} \f$.
 * @param drdt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param drdt_F Density rate of change restricted to the diffusive term
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2, F} \f$.
 * @param posin Position \f$ \mathbf{r}_{n} \f$.
 * @param vin Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dvdtin Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param densin Density \f$ \rho_{n} \f$.
 * @param massin Mass \f$ m_{n} \f$.
 * @param drdtin Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param t Simulation time \f$ t \f$.
 * @param dt Time step \f$ \Delta t \f$.
 * @see Predictor.cl
 * @see Aqua::CalcServer::Predictor
 * @see Aqua::CalcServer::Corrector
 */
__kernel void Corrector(__global int* imove,
                        __global vec* pos,
                        __global vec* v,
                        __global vec* dvdt,
                        __global float* dens,
                        __global float* mass,
                        __global float* drdt,
                        __global float* drdt_F,
                        __global vec* posin,
                        __global vec* vin,
                        __global vec* dvdtin,
                        __global float* densin,
                        __global float* massin,
                        __global float* drdtin,
                        unsigned int N,
                        float t,
                        float dt)
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
	v[i] = v[i] + HDT * (dvdt[i] - dvdtin[i]);
	// mass[i]   = mass[i] * (1.f +  hdt*(drdt[i] - drdtin[i]) / dens[i]);
	#if __BOUNDARY__ == 1
		// Continuity equation must be solved for fixed particles
		dens[i] = dens[i] + 0.5f * dt * (drdt[i] + drdt_F[i] - drdtin[i]);
	#else
		dens[i] = dens[i] + HDT * (drdt[i] + drdt_F[i] - drdtin[i]);
	#endif
	/* Calculate initial positions, mass and thermal energy
	 * for the next step.
	 */
	posin[i] = pos[i];
	vin[i] = v[i];
	massin[i] = mass[i];
	densin[i] = dens[i];
	dvdtin[i] = dvdt[i];
	drdtin[i] = drdt[i] + drdt_F[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
