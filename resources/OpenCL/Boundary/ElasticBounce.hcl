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
const float u_n = dot(u_i - u[j].XYZ, n_j);
const float dudt_n = dot(dudt_i, n_j);
const float dist = dt * u_n + 0.5f * dt * dt * dudt_n;
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
    // u[i] = u_i - (1.f + __ELASTIC_FACTOR__) * (
        // u_n + 0.5f * dt * (dudt_n + g_n)) * n_j;

    // A second approach is setting an acceleration equal to the gravity
    // Just trying to don't perturbate the moments meassurement, fliping
    // later the velocity
    dudt[i].XYZ = dudt_i - dudt_n * n_j;
    u[i].XYZ = u_i - (1.f + __ELASTIC_FACTOR__) * u_n * n_j;

    // Modify the value for the next walls test.
    u_i = u[i].XYZ;
    dudt_i = dudt[i].XYZ;
}
