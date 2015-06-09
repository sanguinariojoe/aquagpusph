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
 * @brief Particles pair viscous term energy computation
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for the particles interaction, which may imply more registers
 * consumption.
 */

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 4.f
    #else
        #define __CLEARY__ 5.f
    #endif
#endif

const float rho_j = rho[j];
const float m_j = m[j];
const float f_ij = kernelF(q) * CONF * m_j;

float lapufac = 0.f;
if(move_j > 0){
#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    const float udr = dot(u[j].XYZ - u_i, r_ij);
    const float r2 = (q * q + 0.01f) * H * H;
    lapufac = __CLEARY__ * udr * udr / (r2 * rho_j * rho_i);
#elif __LAP_FORMULATION__ == __LAP_MORRIS__
    const vec_xyz du = u[j].XYZ - u_i;
    lapufac = dot(du, du) / (rho_i * rho_j);
#else
    #error Unknown Laplacian formulation: __LAP_FORMULATION__
#endif
}
_DWDT_ -= f_ij * lapufac;
