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
 * @brief Leap-frog time integration scheme predictor stage.
 * (See Aqua::CalcServer::Predictor for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Leap-frog time integration scheme predictor stage.
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
 * @param ifluid Fluid index.
 * @param pos Position \f$ \mathbf{r}_{n+1} \f$.
 * @param v Velocity \f$ \mathbf{u}_{n+1} \f$.
 * @param dvdt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param dens Density \f$ \rho_{n+1} \f$.
 * @param mass Mass \f$ m_{n+1} \f$.
 * @param drdt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param posin Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param vin Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dvdtin Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param densin Density \f$ \rho_{n+1/2} \f$.
 * @param massin Mass \f$ m_{n+1/2} \f$.
 * @param drdtin Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param press Pressure \f$ p \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param gamma Eq. of state exponent \f$ \gamma \f$.
 * @param N Number of particles.
 * @param t Simulation time \f$ t \f$.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param grav Gravity acceleration \f$ \mathbf{g} \f$.
 * @param minDens Minimum tolerated density value \f$ \rho_{min} \f$.
 * @param maxDens Maximum tolerated density value \f$ \rho_{max} \f$.
 * @see Corrector.cl
 * @see Aqua::CalcServer::Predictor
 * @see Aqua::CalcServer::Corrector
 */
__kernel void Predictor(__global int* imove,
                        __global int* ifluid,
                        __global vec* pos,
                        __global vec* v,
                        __global vec* dvdt,
                        __global float* dens,
                        __global float* mass,
                        __global float* drdt,
                        __global vec* posin,
                        __global vec* vin,
                        __global vec* dvdtin,
                        __global float* densin,
                        __global float* massin,
                        __global float* drdtin,
                        __global float* press, 
                        __constant float* refd,
                        __constant float* gamma,
                        unsigned int N,
                        float t,
                        float dt,
                        float cs,
                        vec grav,
                        float minDens,
                        float maxDens)
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
			vec dr = 0.5f * dt * dt * (dvdtin[i] + grav);
			posin[i] += dr;
		}
		dt = 0.f;
	}
	// Time step modified by imove flag
	float DT = dt;
	if(imove[i] <= 0)
		DT = 0.f;

	// Predictor step for the fluid and walls
	v[i] = vin[i] + DT * (dvdtin[i] + grav);
	// mass[i]    = massin[i]*(1.f + dt*drdtin[i]/densin[i]);
	#if __BOUNDARY__ == 1
		// Continuity equation must be solved for the fixed particles too
		dens[i] = densin[i] + dt * drdtin[i];
	#else
		dens[i] = densin[i] + DT * drdtin[i];
	#endif
	if(dens[i] < minDens) dens[i] = minDens;
	if(dens[i] > maxDens) dens[i] = maxDens;
	pos[i] = posin[i] + DT * vin[i] + 0.5f * DT * DT * (dvdtin[i] + grav);
	// Batchelor 1967
	{
		const float ddenf = dens[i] / refd[ifluid[i]];
		const float prb = cs * cs * refd[ifluid[i]] / gamma[ifluid[i]];
		press[i] = prb * (pow(ddenf, gamma[ifluid[i]]) - 1.f);
	}
	// Output variables reinitialization
	dvdt[i] = VEC_ZERO;
	drdt[i] = 0.f;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
