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

#ifndef _ARRHENIUS_DETONATION_H_INCLUDED_
#define _ARRHENIUS_DETONATION_H_INCLUDED_

#ifndef SPECIES_HEADER
#error "working with species requires to load a backend module"
#endif
#include SPECIES_HEADER

#include "resources/Scripts/cfd/reactive/reaction_generic.hcl"

#define E_ch 7000.0f
#define K_ch 1.0e7f

/** @brief Compute a measure of the compression.
 *
 * @param zeta progress of reaction.
 * @param T temperature.
 *
 */
inline float
arrhenius_detonation(float zeta, float T)
{

	return K_ch * (1.0f - zeta) * exp(-E_ch / T);
}

/** @brief Compute rate of progress of the reaction. Inside of each particle.
 *
 * @param z tracer.
 * @param y_0 same component as tracer bt burning.
 * @param T temperature.
 * @param MMix molar mass of the mixture.
 *  */
#define LOW_CONC 4.0e-2f
inline float
zeta_dot_calc_arrhenius(float z, float T, float y_0, float MMix)
{
	const species_t Mis = MIS;
	float zx = MMix / Mis.SPECIES_COMPONENT0 * z;

	float zeta;

	// check that there is mroe than 4%
	if (zx > LOW_CONC) {
		zeta = give_zeta(z, y_0);
		return arrhenius_detonation(zeta, T);
	} else {
		return 0.0f;
	}
}

/** @brief Compute rate of progress of the reaction for each mass fraction for
 * each particle.
 *
 * @param z tracer.
 * @param T temperature.
 * @param ys chemical species.
 * @param w_rhos sink sources of the mass fractions.
 * @param deintdt change of the internal energy.
 * @param zeta_dot rate of change of the dimensionless progress of the reaction.
 *  */
void
w_rhos_arrhenius_det(float z,
                     float T,
                     species_t ys,
                     __global species_t* w_rhos,
                     __global float* deintdt,
                     __global float* zeta_dot)
{

	float MMix = molar_mass_mixture(ys);

	*zeta_dot = zeta_dot_calc_arrhenius(z, T, ys.SPECIES_COMPONENT0, MMix);

	w_from_zeta_dot(w_rhos, deintdt, zeta_dot, MMix);

	return;
}

#endif // _ARRHENIUS_DETONATION_H_INCLUDED_