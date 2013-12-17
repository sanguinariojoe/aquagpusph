/*
 *  This file is part of AQUA-gpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUA-gpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUA-gpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
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

/** Method called outside to fill with particles the fluid.
 * @param imove (output) Fix particles flag.
 * @param ifluid (output) Fluid identifier.
 * @param pos (output) Position of particles.
 * @param normal (output) Normal of particles (used for boundary particles/vertexes).
 * @param v (output) Velocity of particles.
 * @param dens (output) dens[i] of particles.
 * @param mass (output) Mass of particles.
 * @param hp (output) Kernel height of particles.
 * @param relPos Relative position of particles to quaternion.
 * @param relNormal Particle normal at t=0 (used for boundary particles/vertexes).
 * @param N Number of particles.
 * @param dt Time step.
 * @param CoR Center of rotation of quaternion.
 * @param X x axis vector of quaternion.
 * @param Y y axis vector of quaternion.
 * @param Z z axis vector of quaternion.
 * @param cor Center of rotation of quaternion at previous time step.
 * @param x x axis vector of quaternion at previous time step.
 * @param y y axis vector of quaternion at previous time step.
 * @param z z axis vector of quaternion at previous time step.
 */
__kernel void Movement( _g int* imove, _g int* ifluid, _g vec* pos,
                        _g vec* normal, _g vec* v, _g float* dens,
                        _g float* mass, _g float* hp, _g vec* relPos,
                        _g vec* relNormal, unsigned int N, float dt, vec CoR,
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

