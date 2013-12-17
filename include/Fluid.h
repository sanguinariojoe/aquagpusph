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

#ifndef FLUID_H_INCLUDED
#define FLUID_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------------
// Include OpenCL libraries
// ----------------------------------------------------------------------------
#include <CL/cl.h>

// ----------------------------------------------------------------------------
// Include Singleton abstract class
// ----------------------------------------------------------------------------
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** @class Fluid Fluid.h Fluid.h
 * @brief Host - Server transfer layer for fluid entity. This transfer layer
 * builds the fluid, and send to Server, when server end a round of
 * iterations, this class recovers all info, and send it to the output manager.
 */
class Fluid : public Aqua::Singleton<Aqua::InputOutput::Fluid>
{
public:
	/** Constructor.
	 */
	Fluid();

	/** Destructor.
	 */
	~Fluid();

	/** Retrieves all data from server.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool retrieveData();

	/** Get the number of particles
	 * @return Number of particles.
	 */
	unsigned int n(){return nParticle;}
	/** Get the number of particles
	 * @return Number of particles.
	 */
	unsigned int nParticles(){return n();}

	/** Get the number of fluids
	 * @return Number of fuid species.
	 */
	unsigned int nFluids(){return nFluid;}

	/** Fix particles & fluid flag. \n
	 * imove>0 for every fluid. \n
	 * imove=0 for sensors. \n
	 * imove<0 for boundary particles/vertex.
	 * @note different imove number imply different fluid species.
	 */
	int *imove;
	/// Fluid of any particle
	int *ifluid;
	/// Position of all particles.
	vec *pos;
	/// Normal at each particle. (Used for boundary particles/vertexes)
	vec *normal;
	/// Velocity of all particles.
	vec *v;
	/// Density of all particles.
	float *dens;
	/// Density variation of all particles.
	float *drdt;
	/// hp of all particles.
	float *hp;
	/// Acceleration of all particles.
	vec *f;
	/// Mass of all particles.
	float *mass;
	/// Pressure of all particles.
	float *press;
	/// Shepard term.
	float *shepard;
	/// Shepard term gradient.
	vec *gradShepard;
	/// grad(p) output
	vec *gradP;
private:
	/// Number of particles
	unsigned int nParticle;
	/// Number of fluids
	unsigned int nFluid;
};

}}  // namespace

#endif // FLUID_H_INCLUDED
