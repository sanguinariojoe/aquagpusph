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

const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
const float rho_j = rho[j];
if(rho_j <= 0.01f * refd_i){
    j++;
    continue;
}

{
    const float area_j = m[j];
    const float p_j = p[j];
    const vec_xyz du = u[j].XYZ - u_i;
    const float w_ij = kernelW(q) * CONW * area_j;

    _GRADP_ += (p_i + p_j) / rho_i * w_ij * n_j;
    _DIVU_ += rho_i * dot(du, n_j) * w_ij;
}
