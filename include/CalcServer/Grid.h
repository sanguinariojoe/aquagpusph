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
 * @brief Link-list grid configuration.
 * (See Aqua::CalcServer::Grid for details)
 */

#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class Grid Grid.h CalcServer/Grid.h
 * @brief Link-list grid configuration tool.
 *
 * This tool will compute the global bounds box (including boundary elements
 * and sensors), setting up the lattice mesh where the particles will be
 * allocated.
 *
 * In fact this tool is consisting in a pair of reductions to compute the bounds
 * box.
 * @see Aqua::CalcSercer::Reduction
 * @see Aqua::CalcServer::LinkList
 */
class Grid : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Grid();

    /// Destructor.
    ~Grid();

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /// Minimum position reduction tool
    Reduction *_maximum;
    /// Maximum position reduction tool
    Reduction *_minimum;
};

}}  // namespace

#endif // GRID_H_INCLUDED
