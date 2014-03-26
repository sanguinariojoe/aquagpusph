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

/** Permutate the data arrays from the unsorted space to the sorted one, and
 * viceversa.
 * @param ifluid permuted output fluid identifier.
 * @param ifluid_in Unpermuted fluid identifier.
 * @param imove permuted output moving flag.
 * @param imove_in Unpermuted moving flag.
 * @param pos permuted output position.
 * @param pos_in Unpermuted position.
 * @param normal permuted output position.
 * @param normal_in Unpermuted position.
 * @param v permuted output velocity.
 * @param v_in Unpermuted velocity.
 * @param f permuted output acceleration.
 * @param f_in Unpermuted acceleration.
 * @param dens permuted output density.
 * @param dens_in Unpermuted density.
 * @param press permuted output pressure.
 * @param press_in Unpermuted pressure.
 * @param mass permuted output mass.
 * @param mass_in Unpermuted mass.
 * @param drdt permuted output density variation rate.
 * @param drdt_in Unpermuted density variation rate.
 * @param drdt_F permuted output density variation rate (due to Delta-SPH).
 * @param drdt_F_in Unpermuted density variation rate (due to Delta-SPH).
 * @param shepard permuted output Shepard factor.
 * @param shepard_in Unpermuted Shepard factor.
 * @param permut Permutations array.
 * @param N Number of particles.
 */
__kernel void Permutate(_g int* ifluid, _g int* ifluid_in,
                        _g int* imove, _g int* imove_in,
                        _g vec* pos, _g vec* pos_in,
                        _g vec* normal, _g vec* normal_in,
                        _g vec* v, _g vec* v_in,
                        _g vec* f, _g vec* f_in,
                        _g float* dens, _g float* dens_in,
                        _g float* press, _g float* press_in,
                        _g float* mass, _g float* mass_in,
                        _g float* drdt, _g float* drdt_in,
                        _g float* drdt_F, _g float* drdt_F_in,
                        _g float* shepard, _g float* shepard_in,
                        _g uint *permut,
                        uint N)
{
	uint i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// We assume i in the unsorted space, and i_out in the sorted one 
	const uint i_out = permut[i];

	ifluid[i_out] = ifluid_in[i];
	imove[i_out] = imove_in[i];
	pos[i_out] = pos_in[i];
	normal[i_out] = normal_in[i];
	v[i_out] = v_in[i];
	f[i_out] = f_in[i];
	dens[i_out] = dens_in[i];
	press[i_out] = press_in[i];
	mass[i_out] = mass_in[i];
	drdt[i_out] = drdt_in[i];
	drdt_F[i_out] = drdt_F_in[i];
	shepard[i_out] = shepard_in[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
