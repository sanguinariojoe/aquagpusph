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

#ifndef ASCII_H_INCLUDED
#define ASCII_H_INCLUDED

#include <sphPrerequisites.h>
#include <Fluid.h>

namespace Aqua{ namespace InputOutput{ namespace Input{

/** Reads an ASCII data file to setup the particles.
 * @param path Input file path.
 * @param ifluid Index of the fluid.
 * @param i0 Starting index of the particle.
 * @param n number of particles to read.
 * @param refd Reference density.
 * @param h Reference kernel height.
 * @param F Fluid host instance.
 * @return false if all gone right, true otherwise.
 */
bool loadASCII(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F);

}}} // namespaces

#endif // ASCII_H_INCLUDED
