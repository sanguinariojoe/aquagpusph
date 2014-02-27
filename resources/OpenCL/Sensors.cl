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
#ifndef iM_PI
	#define iM_PI 0.318309886f
#endif

#ifndef uint
	#define uint unsigned int
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

#ifdef _l
	#error '_l' is already defined.
#endif
#define _l __local

/** Set the sensors interpolated values.
 * @param f Pressure components, it means the interpolated value and
 * the hydrostatic compensation.
 * @param drdt Interpolated density value.
 * @param press Pressure.
 * @param pressin Pressure backup.
 * @param dens Density.
 * @param densin Density backup.
 * @param refd Density of reference.
 * @param ifluid Fluid identifier.
 * @param gamma Gamma EOS exponent.
 * @param shepard Shepard term (0th correction).
 * @param cs Speed of sound.
 * @param i0 First particle which is a sensor.
 * @param N Number of sensors.
 */
__kernel void Sensors( _g vec* f, _g float* drdt, _g float* press,
                       _g float* pressin, _g float* dens, _g float* densin,
                       _g float* refd, _g uint* ifluid, _c float* gamma,
                       _g float* shepard, float cs,
                       uint i0, uint N )
{
	uint i = get_global_id(0);
	if(i >= N)
		return;
	i += i0;	

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	float iShepard = shepard[i];
	if(iShepard < 0.01f){
		// It will be considered that there are not enough
		// particles to interpolate
		iShepard   = 1.f;
		f[i]     = VEC_ZERO;
		drdt[i]  = 0.f;
	}
	dens[i]   = drdt[i]/iShepard;
	densin[i] = dens[i];
	// Batchelor 1967
	float rdenf = refd[ifluid[i]];
	float gammf = gamma[ifluid[i]];
	float ddenf = dens[i]/rdenf;
	float prb   = cs*cs*rdenf/gammf;
	// Based on the density value
	// press[i]    = prb*(pow(ddenf,gammf) - 1.f) + f[i].y/iShepard;
	// All interpolated
	// press[i]    = (f[i].x + f[i].y)/iShepard;
	// Maximum value (Batchelor vs interpolation)
	press[i]    = max(prb*(pow(ddenf,gammf)-1.f), f[i].x/iShepard) + f[i].y/iShepard;
	pressin[i]  = press[i];
	// Restore zero values
	f[i]     = VEC_ZERO;
	drdt[i]  = 0.f;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

