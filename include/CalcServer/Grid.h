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

#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include Reduction tool
// ----------------------------------------------------------------------------
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class Grid Grid.h CalcServer/Grid.h
 * @brief Grid allocation stage. All particle interaction are limited
 * in distance with a know parameter, so the space can be discretized
 * in cells using this charasteristic length, and guarantee that a
 * particle only can interact with other particles included into neighbour
 * cells. Grid stage builds the discretized space grid.
 */
class Grid : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Grid();

	/** Destructor.
	 */
	~Grid();

	/** Performs space discretization.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
    /// Minimum position reduction tool
    Reduction *maximum;
    /// Maximum position reduction tool
    Reduction *minimum;
};

}}  // namespace

#endif // GRID_H_INCLUDED
