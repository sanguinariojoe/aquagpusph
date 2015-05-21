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
 * @brief Boundary element - Fluid particle interaction (friction force).
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for thye particles interaction, which imply more registries
 * consumption.
 *
 * @see https://hal.archives-ouvertes.fr/hal-00691603/document (Martin Ferrand,
 * Dominique Laurence, Benedict Rogers, Damien Violeau, Christophe Kassiotis.
 * Unified semi-analytical wall boundary conditions for inviscid, laminar or
 * turbulent flows in the meshless SPH method. International Journal for 
 * Numerical Methods in Fluids, Wiley-Blackwell, 2013, 71 (476-472).
 */

#if __LAP_FORMULATION__ != __LAP_MORRIS__ && \
    __LAP_FORMULATION__ != __LAP_MONAGHAN__
    #error Unknown Laplacian formulation: __LAP_FORMULATION__
#endif

const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
const float area_j = m[j];

{
    const float w_ij = 2.f * kernelW(q) * CONW * area_j;

    const vec_xyz du = u[j].XYZ - u_i;
    const vec_xyz du_t = du - dot(du, n_j) * n_j;

    #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
        const float r2 = (q * q + 0.01f) * H * H;
        _LAPU_ += w_ij * dot(du, r_ij) / (r2 * rho_i) * n_j;
    #endif
    #if __LAP_FORMULATION__ == __LAP_MORRIS__ || \
          __LAP_FORMULATION__ == __LAP_MONAGHAN__
        const float dr_n = max(fabs(dot(r_ij, n_j)), dr);
        _LAPU_ += w_ij / (rho_i * dr_n) * du_t;
    #endif
}
