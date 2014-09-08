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
 * @brief Quaternion based motions.
 * (See Aqua::CalcServer::Movement::Quaternion for details)
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Move the boundary particles/elements according to the local
 * reference system (Quaternion).
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ifluid Fluid index.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param dens Density \f$ \rho \f$.
 * @param mass Mass \f$ m \f$.
 * @param relPos Local position of each particle \f$ \mathbf{r}_{local} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param CoR Current time step quaternion center.
 * @param X Current time step quaternion x axis vector.
 * @param Y Current time step quaternion y axis vector.
 * @param Z Current time step quaternion z axis vector.
 * @param cor Previous time step quaternion center.
 * @param x Previous time step quaternion x axis vector.
 * @param y Previous time step quaternion y axis vector.
 * @param z Previous time step quaternion z axis vector.
 */
__kernel void Movement( __global int* imove, __global int* ifluid, __global vec* pos,
                        __global vec* normal, __global vec* v, __global float* dens,
                        __global float* mass, __global vec* relPos,
                        __global vec* relNormal, unsigned int N, float dt, vec CoR,
                        vec X, vec Y, vec Z, vec cor, vec x, vec y, vec z)
{
	// find position in global arrays
	int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	if(imove[i]>0)		// Then is a fluid particles
		return;

	// Position at actual time step computation
	vec newX, newY, newPos;
	vec oldX, oldY, oldPos;
	vec nx, ny;
	#ifdef HAVE_3D
		vec oldZ, newZ;
		vec nz;
	#endif
	newX = relPos[i].x * X;
	newY = relPos[i].y * Y;
	#ifndef HAVE_3D
		newPos = CoR + newX + newY;
	#else
		newZ = relPos[i].z * Z;
		newPos = CoR + newX + newY + newZ;
	#endif
	// Position at previous time step computation
	oldX = relPos[i].x * x;
	oldY = relPos[i].y * y;
	#ifndef HAVE_3D
		oldPos = cor + oldX + oldY;
	#else
		oldZ = relPos[i].z * z;
		oldPos = cor + oldX + oldY + oldZ;
	#endif
	// Set output (Sensor must preserve velocity field as output value)
	pos[i]    = newPos;
	if(dt <= 0.f){
		v[i]   = VEC_ZERO;
	}
	else if(imove[i]){
		v[i]   = (newPos - oldPos)/dt;
	}
	// Normal computation
	nx = relNormal[i].x * x;
	ny = relNormal[i].y * y;
	#ifndef HAVE_3D
		normal[i] = nx + ny;
	#else
		nz = relNormal[i].z * z;
		normal[i] = nx + ny + nz;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

