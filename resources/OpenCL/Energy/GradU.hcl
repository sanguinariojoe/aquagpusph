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

const vec_xyz du = u[j].XYZ - u_i;
const vec_xyz f_ij = kernelF(q) * CONF * m[j] / rho[j] * r_ij;

{
    _GRADUX_ +=  du.x * f_ij;
    _GRADUX_ +=  du.x * f_ij;
    #ifdef HAVE_3D
        _GRADUz_ +=  du.z * f_ij;
    #endif
}
