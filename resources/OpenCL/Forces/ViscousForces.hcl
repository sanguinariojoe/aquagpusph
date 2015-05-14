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
 * @brief Fluid particle - Boundary element interaction (friction force).
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for thye particles interaction, which imply more registries
 * consumption.
 *
 * @see NoSlipBI.cl.
 */

const float rho_j = rho[j];
const float m_j = m[j];

{
    const float dr_n = max(fabs(dot(r_ij, n_i)), dr);
    const vec_xyz du = u[j].XYZ - u_i;
    const vec_xyz du_t = du - dot(du, n_i) * n_i;

    const float w_ij = kernelW(q) * CONW * area_i;

    _F_ += 2.f * m_j * w_ij / (rho_j * dr_n) * du_t;
}
