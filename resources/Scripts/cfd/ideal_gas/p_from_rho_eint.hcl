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

#ifndef _PERFECT_GAS_P_INCLUDED_
#define _PERFECT_GAS_P_INCLUDED_

/// @brief 
/// @param gamma polytropic coefficient
/// @param p pressure
/// @param eint internal energy
/// @param rho density
/// @return speed of sound

float p_from_rho_eint(float gamma, float rho, float eint)
{
    return (gamma - 1.0f) * rho * eint;
}

float eint_from_rho_p(float gamma, float rho, float p)
{
    return p /((gamma - 1.0f) * rho);
}

#endif //_PERFECT_GAS_P_INCLUDED_