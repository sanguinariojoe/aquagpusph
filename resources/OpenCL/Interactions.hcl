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

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

const float rho_j = rho[j];
const float p_j = p[j];
const float f_ij = kernelF(q) * CONF * m[j];
const float udr = dot(u[j].XYZ - u_i, r_ij);

const float prfac = (p_i + p_j) / (rho_i * rho_j);
_GRADP_ += prfac * f_ij * r_ij;

vec_xyz lapufac = VEC_ZERO.XYZ;
#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    const float r2 = (q * q + 0.01f) * H * H;
    lapufac = __CLEARY__ * udr / (r2 * rho_i * rho_j) * r_ij;
#elif __LAP_FORMULATION__ == __LAP_MORRIS__
    lapufac = 2.f / (rho_i * rho_j) * (u[j].XYZ - u_i);
#else
    #error Unknown Laplacian formulation: __LAP_FORMULATION__
#endif
_LAPU_ += f_ij * lapufac;

_DIVU_ += udr * f_ij * rho_i / rho_j;

const float drfac = (p_j - p_i) - refd_i * dot(g.XYZ, r_ij);
_LAPP_ += drfac * f_ij / rho_j;
