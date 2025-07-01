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

/** @defgroup ideal_gas Preset to model ideal gases
 *
 * @brief A preset to model ideal gases within @ref cfd preset
 *
 * @{
 */

/** @file
 * @brief Equation Of State (EOS) for ideal gases.
 *
 * This is an replacement for resources/Scripts/basic/EOS.cl
 */

#ifndef EXCLUDED_PARTICLE
#define EXCLUDED_PARTICLE(index) (imove[index] <= 0) && (imove[index] != -1)
#endif

#include "resources/Scripts/types/types.h"
//#include "resources/Scripts/cfd/species/species_auxiliary.hcl"
#include SPECIES_HEADER
//#define TEST_HEADER_MACRO "species_h2.hcl"
//#include TEST_HEADER_MACRO
/** @brief Ideal gas Equation Of State (EOS) computation
 *
 * The equation of state relates the pressure, density and internal energy
 * fields,
 * \f$ p = \rho \left( \gamma - 1 \right) e \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param eint Internal energy \f$ e_{n+1/2} \f$.
 * @param p Pressure \f$ p_{n+1/2} \f$.
 * @param gamma Heat capacity ratio \f$ \gamma \f$.
 * @param N Number of particles.
 * @param xs molar fraction
 * @param ys mass fraction
 * @param p pressure
 * @param T temperature
 * @param gamma polytropic coefficient
 * @param cp heat at constant pressure
 * @param cv heat at constant volume
 * @param nu kinematic viscosity
 * @param xi thermal diffisivity
 * @param mu dynamic viscosity
 * @param lambda thermal conductivity
 */

__kernel void
entry(const __global unsigned int* iset,
      const __global int* imove,
      const __global float* rho,
      const __global float* eint,
	  const __global vec16* ys,
/*      const __global float* y_H2,
      const __global float* y_O2,
      const __global float* y_N2,
      const __global float* y_H2O,*/
	  const __global vec16* xs,
/*      __global float* x_H2,
      __global float* x_O2,
      __global float* x_N2,
      __global float* x_H2O,*/
      __global float* p,
      __global float* T,
      __global float* gamma,
      __global float* cp,
      __global float* cv,
      const __global float* nu,
      const __global float* xi,
      __global float* mu,
      __global float* lambda,
      usize N)
{
	usize i = get_global_id(0);
	if (i >= N)
		return;
	if (EXCLUDED_PARTICLE(i))
		return;

	//__global float cp[]={0.0f,};

	calc_gamma_cp_cv(
	    //y_H2[i], y_O2[i], y_N2[i], y_H2O[i], 
		ys+i,
		gamma + i, cv + i, cp + i);
	X_from_Y(/*y_H2[i],
	         y_O2[i],
	         y_N2[i],
	         y_H2O[i],*/
			 ys+i,
	         /*x_H2 + i,
	         x_O2 + i,
	         x_N2 + i,
	         x_H2O + i*/
			 xs+i);

	p[i] = (gamma[i] - 1.0f) * rho[i] * eint[i];
	T[i] = eint[i] / cv[i];
	// printf("cv = %f\n", cv[i]);
	// printf("T = %f\n", T[i]);
	lambda[i] = cp[i] * rho[i] * xi[i] * sqrt(T[i] / 298.0f);

	// This line is for debug!!!!
	lambda[i] = 1012.0f * xi[i];

	mu[i] = rho[i] * nu[i] * sqrt(T[i] / 298.0f);
}

/*
 * @}
 */
