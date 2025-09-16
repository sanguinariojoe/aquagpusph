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
// #include "resources/Scripts/cfd/species/species_auxiliary.hcl"
#include "resources/Scripts/cfd/reaction/reaction_generic.hcl"

__constant float S_L = 3.0f;
__constant float sigma = 5.0f;

/** @brief Compute factor for the Zimont's model.
 *
 * @param mod_grad_zeta module of the gradient of z
 *  */
inline float
zimont_deflagration(float mod_grad_zeta)
{

	return sigma * S_L * mod_grad_zeta;
}

/** @brief Compute rate of progress of the reaction. Inside of each particle.
 *
 * @param z tracer.
 * @param MMix molar mass of the mixture.
 * @param mod_grad_zeta module of the gradient of z
 *  */
float
zeta_dot_calc_zimont(float z, float MMix, float mod_grad_zeta)
{

	float zx = MMix / Mis[0] * z;

	// check that there is mroe than 4%
	if (zx > 4.e-2f) {
		return zimont_deflagration(mod_grad_zeta)
	} else {
		//        zeta = 1.0;
		return 0.0f;
	}
}

/** @brief Compute rate of progress of the reaction for each mass fraction for each particle.
 *
 * @param z tracer.
 * @param T temperature.
 * @param ys chemical species.
 * @param w_rhos sink sources of the mass fractions.
 * @param deintdt change of the internal energy.
 * @param zeta_dot rate of change of the dimensionless progress of the reaction.
 *  */
void
w_rhos_arrhenius_def(float z,
                     float T,
                     species_t ys,
                     __flobal species_t* w_rhos,
                     __global float* deintdt,
                     __global float* zeta_dot,
                     vec* grad_zeta)
{

	float MMix = molar_mass_mixture(ys);

	float mod_grad_zeta = sqrt(dot(grad_zeta, grad_zeta));

	*zeta_dot = zeta_dot_calc_zimont(z, MMix, mod_grad_zeta);

	w_from_zeta_dot(w_rhos, deintdt, zeta_dot, MMix);

	return;
}

#endif // _ARRHENIUS_DETONATION_H_INCLUDED_