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
 * @brief Boundary element - Fluid particle interaction.
 * (See ElasticBounce.cl for details)
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for thye particles interaction, which imply more registries
 * consumption.
 */

if(imove[j] >= 0){
	j++;
	continue;
}

// ------------------------------------------------------------------
// face properties
// ------------------------------------------------------------------
vec n_j = normal[j];
const vec v_j = v[j];
// ------------------------------------------------------------------
// Vertex relation
// ------------------------------------------------------------------
const vec r  = pos_i - pos[j];
float r0 = dot(r, n_j);
const vec rt = r - r0 * n_j;
if(dot(rt, rt) >= r_element * r_element){
    // The particle is passing too far from the wall element
	j++;
	continue;
}

// Test for swap normal (that must be internal oriented)
if(r0 < 0.f){
	n_j = -n_j;
	r0 = -r0;
}
// ------------------------------------------------------------------
// Movement data
// ------------------------------------------------------------------
const float v_n = dot(v_i - v_j, n_j);
const float dvdt_n = dot(dvdt_i, n_j);
const float g_n = dot(grav, n_j);
const float dist = dt * v_n + 0.5f * dt * dt * (dvdt_n + g_n);
// ------------------------------------------------------------------
// Since normal has been internally oriented, if dist < 0, the
// particle is running against the wall, and then two cases can be
// discriminated:
//   - The particle is placed in the effect zone of the wall.
//   - The partilce is placed outside the effect zone, but will enter inside it.
// ------------------------------------------------------------------
if((dist < 0.f) && (r0 + dist <= __MIN_BOUND_DIST__ * h)){
	// ------------------------------------------------------------------
	// Reflect particle velocity (using elastic factor)
	// ------------------------------------------------------------------
	// As first approach, particle can be just fliped
	// v[i] = v_i - (1.f + __ELASTIC_FACTOR__) * (
        // v_n + 0.5f * dt * (dvdt_n + g_n)) * n_j;

	// A second approach is setting an acceleration equal to the gravity
	// Just trying to don't perturbate the moments meassurement, fliping
	// later the velocity
	dvdt[i] = dvdt_i - (dvdt_n + g_n) * n_j;
	v[i] = v_i - (1.f + __ELASTIC_FACTOR__) * v_n * n_j;

	// Modify the value for the next walls test.
	v_i = v[i];
	dvdt_i = dvdt[i];
}
