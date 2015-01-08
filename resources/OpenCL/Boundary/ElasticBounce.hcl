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

// ------------------------------------------------------------------
// Movement data
// ------------------------------------------------------------------
const float v_n = dot(v_i - v[j].XYZ, n_j);
const float dvdt_n = dot(dvdt_i, n_j);
const float dist = dt * v_n + 0.5f * dt * dt * dvdt_n;
if(dist < 0.f){
    // The particle is already running away from the boundary
    j++;
    continue;
}

// ------------------------------------------------------------------
// The particle should be corrected if:
//   - It is already placed in the effect zone.
//   - It is entering inside the effect zone.
// ------------------------------------------------------------------
if(r0 - dist <= __MIN_BOUND_DIST__ * dr){
    // ------------------------------------------------------------------
    // Reflect particle velocity (using elastic factor)
    // ------------------------------------------------------------------
    // As first approach, particle can be just fliped
    // v[i] = v_i - (1.f + __ELASTIC_FACTOR__) * (
        // v_n + 0.5f * dt * (dvdt_n + g_n)) * n_j;

    // A second approach is setting an acceleration equal to the gravity
    // Just trying to don't perturbate the moments meassurement, fliping
    // later the velocity
    dvdt[i].XYZ = dvdt_i - dvdt_n * n_j;
    v[i].XYZ = v_i - (1.f + __ELASTIC_FACTOR__) * v_n * n_j;

    // Modify the value for the next walls test.
    v_i = v[i].XYZ;
    dvdt_i = dvdt[i].XYZ;
}
