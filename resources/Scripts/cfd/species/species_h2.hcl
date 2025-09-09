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

#ifndef _SPECIES_H2_INCLUDED_
#define _SPECIES_H2_INCLUDED_

#define PASTE(x, y) x##_##y
#define EVALUATE_AND_PASTE(x, y) PASTE(x, y)
#define reduce_sum EVALUATE_AND_PASTE(reduce_sum, species_t)

#define R_gas 8.31f
#define CPS (species_t)(14200.f, 918.f, 1040.f, 2050.f)
#define MIS (species_t)(0.002f, 0.032f, 0.028f, 0.018f)
#define NUS (species_t)(-1.0f, -0.5f, 0.0f, 1.0f)
#define HPLUS_MASS (species_t)(0.0f, 0.0f, 0.0f, -285.83e3f / 0.018f)

inline float
molar_mass_mixture(const species_t y)
{
	const species_t Mis = MIS;
	return reduce_sum(y * MIS);
}

inline float
calc_cp_mix(species_t y)
{
	const species_t cps = CPS;
	return reduce_sum(y * cps);
}

void
X_from_Y(const __global species_t* ys, __global species_t* xs)
{
	float MMix;
	const species_t Mis = MIS;

	MMix = molar_mass_mixture(*ys);

	*xs = MMix / Mis * *ys;

	return;
}

void
calc_gamma_cv(const __global species_t* ys,
              __global float* gamma,
              __global float* cv)
{

	float MMix, R_mix, cp_local, cv_local;

	MMix = molar_mass_mixture(*ys);
	cp_local = calc_cp_mix(*ys);

	R_mix = R_gas / MMix;

	cv_local = cp_local - R_mix;

	*gamma = cp_local / cv_local;
	*cv = cv_local;

	return;
}

void
calc_gamma_cp_cv(const __global species_t* ys,
                 __global float* gamma,
                 __global float* cv,
                 __global float* cp)
{
	float MMix, R_mix, cp_local, cv_local;

	MMix = molar_mass_mixture(*ys);
	cp_local = calc_cp_mix(*ys);

	R_mix = R_gas / MMix;

	cv_local = cp_local - R_mix;

	*gamma = cp_local / cv_local;
	*cv = cv_local;
	*cp = cp_local;

	return;
}

#endif // _SPECIES_H2_INCLUDED_
